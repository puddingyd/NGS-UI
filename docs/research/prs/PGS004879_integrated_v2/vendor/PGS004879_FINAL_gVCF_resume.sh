#!/usr/bin/env bash
# PGS004879 FINAL: sample-independent single-sample GRCh38 gVCF pipeline.
# Research validation for one known sample does NOT establish performance in other samples.
set -Eeuo pipefail

usage() {
cat <<'HELP'
PGS004879 FINAL | single-sample GRCh38 gVCF -> variant QC -> strand verification -> PGEN -> PGS

Required:
  --gvcf FILE             GRCh38 gVCF, bgzip indexed (one sample)
  --out DIR               New output directory, never overwrite
Optional:
  --work DIR              Tool/data root (default: /home/n102968/2026_INTERVENE_PRS/26T00970)
  --scorefile FILE        Harmonized PGS004879_hmPOS_GRCh38.txt.gz
  --ref FILE              GRCh38 FASTA and .fai
  --dbsnp FILE            NCBI dbSNP157 GRCh38 VCF and index
  --plink-img FILE        pgsc_calc PLINK2 2.00a5.10 Apptainer image
  --container-config FILE  Existing local container overrides for pgsc_calc v2.3.0
  --sampleset NAME        Sampleset label (default: case<SAMPLE>)
  --min-dp INT            Default: 10
  --min-gq INT            Default: 20
  --min-overlap FRACTION  Default: 0.75
  --palindromic-mode assume-direct|exclude   Default: assume-direct (historic compatibility; WARN)
  --check-only            Read header, verify inputs, no output directory created
  --discover-only         Build all score-locus audit + inputs, but no PLINK2/Nextflow
  --expect-sum NUMBER     Optional external benchmark, checked only after scoring
  --expect-denom INT      Optional external benchmark, checked only after scoring
  --threads INT           Default: 2
  --memory-mb INT         Default: 8192
  -h, --help              Print this message

Scope: PGS004879, GRCh38, autosomal biallelic SNP scores, single-sample gVCF.
No historical 26T00970 VCF/audit/score is used to determine this case's genotype.
All exclusions, uncertain strand and score outputs are reported per sample.
HELP
}
WORK=/home/n102968/2026_INTERVENE_PRS/26T00970
GVCF='' OUT='' SCORE='' REF='' DB='' PLINK_IMG='' CONFIG=''
SAMPLESET='' EXPECT_SUM='' EXPECT_DENOM=''
MIN_DP=10 MIN_GQ=20 MIN_OVERLAP=0.75 PALINDROMIC_MODE=assume-direct
THREADS=2 MEMORY_MB=8192 CHECK_ONLY=0 DISCOVER_ONLY=0
while (($#));do
  case "$1" in
    --work) WORK="${2:?}";shift 2;;
    --gvcf) GVCF="${2:?}";shift 2;;
    --out) OUT="${2:?}";shift 2;;
    --scorefile) SCORE="${2:?}";shift 2;;
    --ref) REF="${2:?}";shift 2;;
    --dbsnp) DB="${2:?}";shift 2;;
    --plink-img) PLINK_IMG="${2:?}";shift 2;;
    --container-config) CONFIG="${2:?}";shift 2;;
    --sampleset) SAMPLESET="${2:?}";shift 2;;
    --min-dp) MIN_DP="${2:?}";shift 2;;
    --min-gq) MIN_GQ="${2:?}";shift 2;;
    --min-overlap) MIN_OVERLAP="${2:?}";shift 2;;
    --palindromic-mode) PALINDROMIC_MODE="${2:?}";shift 2;;
    --threads) THREADS="${2:?}";shift 2;;
    --memory-mb) MEMORY_MB="${2:?}";shift 2;;
    --expect-sum) EXPECT_SUM="${2:?}";shift 2;;
    --expect-denom) EXPECT_DENOM="${2:?}";shift 2;;
    --check-only) CHECK_ONLY=1;shift;;
    --discover-only) DISCOVER_ONLY=1;shift;;
    --help|-h) usage;exit 0;;
    *) echo "Unknown option: $1" >&2; usage >&2;exit 2;;
  esac
done
[[ -n "$GVCF" && -n "$OUT" ]] || { usage >&2; exit 2; }
[[ "$PALINDROMIC_MODE" == assume-direct || "$PALINDROMIC_MODE" == exclude ]] || { echo 'Invalid --palindromic-mode' >&2; exit 2; }
[[ "$MIN_DP" =~ ^[0-9]+$ && "$MIN_GQ" =~ ^[0-9]+$ && "$THREADS" =~ ^[0-9]+$ && "$MEMORY_MB" =~ ^[0-9]+$ ]] || { echo 'Invalid numeric option' >&2;exit 2; }
WORK="$(realpath -m "$WORK")"
GVCF="$(realpath -m "$GVCF")"
OUT="$(realpath -m "$OUT")"
SCORE="$(realpath -m "${SCORE:-$WORK/scorefiles/PGS004879_hmPOS_GRCh38.txt.gz}")"
REF="$(realpath -m "${REF:-/home/pipeline/reference/hg38/Homo_sapiens_assembly38.fasta}")"
DB="$(realpath -m "${DB:-$WORK/reference/dbsnp157/GCF_000001405.40.gz}")"
PLINK_IMG="$(realpath -m "${PLINK_IMG:-$WORK/pgsc_calc/apptainer_cache/ghcr.io-pgscatalog-plink2-2.00a5.10.img}")"
CONFIG="$(realpath -m "${CONFIG:-$WORK/corrected_v1_chrwise/use_local_containers.config}")"
[[ "$OUT" != "$WORK" ]] || { echo 'Refusing output==work root' >&2;exit 2; }
for input in "$GVCF" "$SCORE" "$REF" "$DB" "$PLINK_IMG" "$CONFIG";do
  [[ "$OUT" != "$input" && "$OUT" != "$input"/* ]] || { echo "Unsafe output overlaps input: $input" >&2;exit 2; }
done
for c in python3 bcftools samtools apptainer nextflow realpath gzip sha256sum awk;do
  command -v "$c" >/dev/null || { echo "Missing tool: $c" >&2;exit 2; }
done
for p in "$GVCF" "$SCORE" "$REF" "$DB" "$PLINK_IMG" "$CONFIG";do
  [[ -s "$p" ]] || { echo "Missing input: $p" >&2;exit 2; }
done
[[ -s "${REF}.fai" ]] || { echo 'Missing FASTA index' >&2;exit 2; }
[[ -s "${DB}.tbi" || -s "${DB}.csi" ]] || { echo 'Missing dbSNP index' >&2;exit 2; }
[[ -s "${GVCF}.tbi" || -s "${GVCF}.csi" ]] || { echo 'Missing gVCF index' >&2;exit 2; }
SAMPLE="$(bcftools query -l "$GVCF")"
[[ $(printf '%s\n' "$SAMPLE" | sed '/^$/d' | wc -l) -eq 1 ]] || { echo 'Requires exactly one gVCF sample' >&2;exit 2; }
[[ "$SAMPLE" =~ ^[a-zA-Z0-9_.-]+$ ]] || { echo 'Unsafe sample identifier in gVCF' >&2;exit 2; }
SAMPLESET="${SAMPLESET:-case${SAMPLE}}"
[[ "$SAMPLESET" =~ ^[a-zA-Z0-9_.-]+$ ]] || { echo 'Unsafe sampleset label' >&2;exit 2; }
PGS_ID=PGS004879
export WORK GVCF OUT SCORE REF DB PLINK_IMG CONFIG SAMPLE SAMPLESET MIN_DP MIN_GQ MIN_OVERLAP PALINDROMIC_MODE EXPECT_SUM EXPECT_DENOM THREADS MEMORY_MB PGS_ID
python3 - <<'PYOPT'
import os
from decimal import Decimal
v=Decimal(os.environ['MIN_OVERLAP'])
if not Decimal(0)<v<=1:raise ValueError('min-overlap must be (0,1]')
for k in ('THREADS','MEMORY_MB'):
 if int(os.environ[k])<=0:raise ValueError(f'{k} must be positive')
PYOPT
if (( CHECK_ONLY ));then
  echo "CHECK_ONLY PASS: SAMPLE=$SAMPLE GVCF=$GVCF SCORE=$SCORE DBSNP=$DB OUT=$OUT"
  exit 0
fi
[[ ! -e "$OUT" ]] || { echo "Output already exists (refuse overwrite): $OUT" >&2; exit 2; }
mkdir -p "$OUT"/{audit,prep,scorefiles,discovery,scripts,logs,roundtrip_22chr/{input,pgen,export,logs},results}
exec > >(tee "$OUT/logs/run.log") 2>&1
trap 'echo "FAILED at line $LINENO; inspect $OUT/logs/run.log" >&2' ERR
echo "START $(date -Is) | sample=$SAMPLE | sampleset=$SAMPLESET | output=$OUT"
echo "Tool versions:";bcftools --version | head -2;samtools --version | head -2
apptainer exec --bind /home:/home "$PLINK_IMG" plink2 --version
nextflow -version | head -8
sha256sum "$GVCF" "$SCORE" > "$OUT/audit/source_sha256.txt"
stat -c '%n %s bytes %y' "$REF" "$DB" "$PLINK_IMG" "$CONFIG" > "$OUT/audit/references_file_manifest.txt"
cp "$CONFIG" "$OUT/local_containers.config"
cat > "$OUT/scripts/gvcf_to_pgs.py" <<'PY_GVCF'
#!/usr/bin/env python3
"""Rebuild PGS004879 score-directed loci from gVCF, never copying historical GTs."""
import csv, gzip, json, os, re, subprocess, sys
from collections import Counter, defaultdict
from decimal import Decimal
from pathlib import Path
from bisect import bisect_left

ROOT=Path(os.environ['OUT']); SRC=ROOT/'discovery'; AUD=ROOT/'audit'
GVCF=Path(os.environ['GVCF']); SCORE=Path(os.environ['SCORE'])
REF=Path(os.environ['REF']); DB=Path(os.environ['DB'])
SAMPLE=os.environ['SAMPLE']; PGS_ID=os.environ['PGS_ID']
MIN_DP=int(os.environ['MIN_DP']); MIN_GQ=int(os.environ['MIN_GQ'])
PALINDROMIC_MODE=os.environ['PALINDROMIC_MODE']
SCORE_OUT=ROOT/'scorefiles/PGS004879_hmPOS_GRCh38.txt.gz'
UNSORTED=ROOT/'prep/score_ready.unsorted.vcf'
AUDIT=SRC/'strand_flip_loci.tsv'
SUMMARY=AUD/'upstream_reconstruction_summary.json'
DECISIONS=AUD/'all_score_locus_decisions.tsv'
DIFF=AUD/'unmatched_loci.tsv'
SCORE_KEYS={'hm_chr','hm_pos','effect_allele','other_allele','effect_weight'}
VALID=set('ACGT'); COMP=str.maketrans('ACGT','TGCA')

def fail(msg): raise RuntimeError(msg)
def key(c,p): return (str(c).removeprefix('chr'),int(p))
def comp(s): return s.translate(COMP)
def ids(s): return {x for x in re.split(r'[;,\s]+',s or '') if re.fullmatch(r'rs\d+',x)}
def gt_dose(gt):
    a=re.split(r'[/|]',gt)
    if len(a)!=2 or any(v not in ('0','1') for v in a): return None
    return a.count('1')
def parse_info_end(s,pos):
    m=re.search(r'(?:^|;)END=(\d+)(?:;|$)',s)
    return int(m.group(1)) if m else pos

def score_targets():
    targets={}; raw=[]; skipped=[]
    with gzip.open(SCORE,'rt') as f:
        for meta in f:
            if meta.startswith('#'):
                if re.match(r'^#pgs_id[=:]\s*PGS[0-9]+',meta,flags=re.I):
                    pid=re.search(r'PGS[0-9]+',meta,re.I).group(0).upper()
                    if pid!=PGS_ID:fail(f'Incorrect scoring file: {pid}, expected {PGS_ID}')
                continue
            break
    with gzip.open(SCORE,'rt') as f:
        reader=csv.DictReader((s for s in f if not s.startswith('#')),delimiter='\t')
        if not SCORE_KEYS <= set(reader.fieldnames or []): fail('Required score columns not present')
        for i,r in enumerate(reader,1):
            raw.append(r)
            chrom=str(r['hm_chr']).removeprefix('chr')
            pos=r['hm_pos']
            if not chrom.isdecimal() or not 1<=int(chrom)<=22 or not pos.isdecimal():
                skipped.append((i,'NO_AUTOSOMAL_COORDINATE'));continue
            k=(chrom,int(pos))
            if k in targets: fail(f'Duplicate score coordinate: {k}')
            targets[k]=r
    if not raw: fail('Empty scoring file')
    if len(targets)==0: fail('No autosomal SNPs in scoring file')
    return targets,skipped

def fasta_at_targets(targets):
    p=SRC/'score_target_fasta_regions.txt'
    with p.open('w') as f:
        for c,pos in sorted(targets,key=lambda k:(int(k[0]),k[1])):
            f.write(f'chr{c}:{pos}-{pos}\n')
    # Stream output, avoid holding 610K subprocess results as one giant string.
    bases={}; current=None
    proc=subprocess.Popen(['samtools','faidx','-r',str(p),str(REF)],stdout=subprocess.PIPE,text=True)
    for line in proc.stdout:
        line=line.strip()
        if line.startswith('>'):
            m=re.match(r'^>chr(\d+):(\d+)-(\d+)',line)
            if not m: fail(f'Unexpected FASTA header: {line}')
            current=(m.group(1),int(m.group(2)))
        elif line:
            if current is None or len(line)!=1: fail('Unexpected FASTA sequence format')
            bases[current]=line.upper();current=None
    if proc.wait()!=0: fail('samtools faidx failed')
    if set(bases)!=set(targets): fail('FASTA did not return all score loci')
    return bases

def get_contigs():
    h=subprocess.check_output(['bcftools','view','-h',str(GVCF)],text=True)
    contigs=set(re.findall(r'^##contig=<ID=([^,>]+)',h,re.M))
    mapping={}
    for i in range(1,23):
        mapping[str(i)]=f'chr{i}' if f'chr{i}' in contigs else str(i) if str(i) in contigs else None
        if mapping[str(i)] is None: fail(f'gVCF missing chr{i} contig')
    return mapping

def source_record(parts):
    if len(parts)!=10: fail('gVCF must have exactly one sample (10 columns)')
    c,p,rs,ref,alt,qual,flt,info,fmt,sample=parts
    pos=int(p); end=parse_info_end(info,pos)
    fmtkeys=fmt.split(':'); vals=sample.split(':'); samplefields=dict(zip(fmtkeys,vals))
    real=[a.upper() for a in alt.split(',') if a not in ('.','<NON_REF>','<*>')]
    is_block=(len(real)==0 and (end>pos or alt in ('<NON_REF>','<*>','.')))
    if not is_block and end<pos+len(ref)-1: end=pos+len(ref)-1
    return {'c':c,'pos':pos,'end':end,'ref':ref.upper(),'alt':alt.upper(), 'real':real,
            'filter':flt,'gt':samplefields.get('GT','.'),'dp':samplefields.get('DP','.'),
            'min_dp':samplefields.get('MIN_DP','.'),'gq':samplefields.get('GQ','.'),'block':is_block}

def observe_at(k,source,fasta,score):
    """Project the target PGS SNV from a (possibly multiallelic) gVCF record.

    Non-target alleles are never silently collapsed into REF.  For explicit
    variants, only source GT alleles {0, target_index} are allowed; the resulting
    biallelic VCF GT is renumbered into {0,1} without altering dosage or phase.
    """
    c,pos=k; ea=score['effect_allele'].upper();oa=score['other_allele'].upper()
    r={'CHROM':f'chr{c}','POS':pos,'RSID':score.get('rsID','.'),
       'HM_RSID':score.get('hm_rsID','.'),'EFFECT_ALLELE':ea,'OTHER_ALLELE':oa,
       'WEIGHT':score['effect_weight'],'FASTA_REF':fasta[k],
       'COMP_EFFECT':comp(ea),'COMP_OTHER':comp(oa),'SOURCE_GT':'',
       'SOURCE_FILTER':'','SOURCE_POS':'','SOURCE_END':'','SOURCE_REF':'',
       'SOURCE_ALT':'','SOURCE_DP':'','SOURCE_MIN_DP':'','SOURCE_GQ':'',
       'FLIP_TRIAGE':'','AUDIT_FLAGS':'','PRELIMINARY_STATUS':'',
       'RECALC_EFFECT_DOSAGE':'','RECALC_CONTRIBUTION':'',
       'ORIENTATION':'','STATUS':'','VCF_REF':'','VCF_ALT':'','VCF_GT':'',
       'ALT_DOSAGE':'','SOURCE_TARGET_ALT_INDEX':'','GVCF_TYPE':'','STRAND_AMBIGUITY':''}
    def reject(status): r['STATUS']=status;return r
    if len(ea)!=1 or len(oa)!=1 or ea not in VALID or oa not in VALID or ea==oa:
        return reject('UNSUPPORTED_SCORE_ALLELES')
    if fasta[k] not in VALID: return reject('NON_ACGT_FASTA_BASE')
    if source is None: return reject('NO_GVCF_COVERAGE')
    r.update(SOURCE_GT=source['gt'],SOURCE_FILTER=source['filter'],SOURCE_POS=source['pos'],
             SOURCE_END=source['end'],SOURCE_REF=source['ref'],SOURCE_ALT=source['alt'],
             SOURCE_DP=source['dp'],SOURCE_MIN_DP=source['min_dp'],SOURCE_GQ=source['gq'],
             GVCF_TYPE='BLOCK' if source['block'] else 'VARIANT')
    ref=fasta[k];r['VCF_REF']=ref
    direct_pair={ea,oa};flip_pair={comp(ea),comp(oa)}
    if source['block']:
        if source['gt'].replace('|','/')!='0/0': return reject('NON_HOMREF_REFERENCE_BLOCK')
        if ref in direct_pair: orient='DIRECT';forward=(ea,oa)
        elif ref in flip_pair: orient='FLIP';forward=(comp(ea),comp(oa))
        else:return reject('ALLELE_MISMATCH')
        projected_gt=source['gt'];alt_dosage=0
        target_index=0
    else:
        if source['pos']!=pos or source['ref']!=ref or len(ref)!=1:
            return reject('COMPLEX_MULTIALLELIC_OR_REF_MISMATCH')
        alt_list=[x.upper() for x in source['alt'].split(',')]
        direct_alt=next(iter(direct_pair-{ref}),None) if ref in direct_pair else None
        flip_alt=next(iter(flip_pair-{ref}),None) if ref in flip_pair else None
        if direct_alt is not None and direct_alt in alt_list:
            orient='DIRECT';forward=(ea,oa);target_alt=direct_alt
        elif flip_alt is not None and flip_alt in alt_list:
            orient='FLIP';forward=(comp(ea),comp(oa));target_alt=flip_alt
        else:
            return reject('ALLELE_MISMATCH')
        if target_alt not in VALID or alt_list.count(target_alt)!=1:
            return reject('UNSUPPORTED_TARGET_ALT')
        target_index=alt_list.index(target_alt)+1
        original_gt=source['gt']
        gtparts=re.split(r'[/|]',original_gt)
        if len(gtparts)!=2 or any(not x.isdigit() for x in gtparts):
            return reject('NON_DIPLOID_OR_MISSING_GT')
        nums=[int(x) for x in gtparts]
        if any(x not in (0,target_index) for x in nums):
            return reject('NON_TARGET_ALT_IN_GT')
        delimiter='|' if '|' in original_gt else '/'
        projected_gt=delimiter.join(str(int(x==target_index)) for x in nums)
        alt_dosage=sum(x==target_index for x in nums)
    r['ORIENTATION']=orient
    r['VCF_ALT']=forward[0] if ref==forward[1] else forward[1]
    r['VCF_GT']=projected_gt
    r['SOURCE_TARGET_ALT_INDEX']=target_index
    if not r['RSID'].startswith('rs') and r['HM_RSID'].startswith('rs'):
        r['RSID']=r['HM_RSID']
    if orient=='FLIP':
        r['FLIP_TRIAGE']='ALLELE_COMPATIBLE_BLOCK' if source['block'] else 'ALLELE_COMPATIBLE_VARIANT'
    if {ea,oa}=={comp(ea),comp(oa)}:
        r['STRAND_AMBIGUITY']='PALINDROMIC_ORIENTATION_NOT_RESOLVED'
        if PALINDROMIC_MODE=='exclude':return reject('PALINDROMIC_STRAND_AMBIGUOUS')
    if source['filter']!='PASS':return reject('NON_PASS_FILTER')
    dp=(source['min_dp'] if source['block'] and source['min_dp'] not in ('','.','NA') else source['dp'])
    if dp in ('','.','NA') or source['gq'] in ('','.','NA'):
        return reject('MISSING_DP_OR_GQ')
    try: dpn=int(dp);gqn=int(source['gq'])
    except ValueError:return reject('MISSING_DP_OR_GQ')
    if dpn<MIN_DP:return reject('LOW_DP')
    if gqn<MIN_GQ:return reject('LOW_GQ')
    # Effect allele dosage uses a target-specific biallelic projected GT.
    dose=alt_dosage if forward[0]==r['VCF_ALT'] else 2-alt_dosage
    r['ALT_DOSAGE']=alt_dosage;r['RECALC_EFFECT_DOSAGE']=dose
    try: con=Decimal(score['effect_weight'])*Decimal(dose)
    except Exception:return reject('INVALID_EFFECT_WEIGHT')
    r['RECALC_CONTRIBUTION']=str(con)
    r['STATUS']='CANDIDATE_'+orient
    return r

def choose_observation(k,active,bases,score):
    """Prefer exact target SNVs over overlapping indels/blocks.

    Finalization is delayed until the stream advances BEYOND the target POS,
    so same-position records cannot be accidentally skipped.  A preceding
    REF-spanning indel cannot mask an exact-position SNV.
    """
    c,pos=k;ref=bases[k]
    ea=score['effect_allele'].upper();oa=score['other_allele'].upper()
    allele_pairs=({ea,oa},{comp(ea),comp(oa)})
    applicable=[s for s in active if s['pos']<=pos<=s['end']]
    if not applicable:return observe_at(k,None,bases,score)
    def relevance(s):
        exact=s['pos']==pos
        if not s['block'] and exact and s['ref']==ref:
            alts=s['alt'].split(',')
            compatible=any(ref in pair and any(a in pair and a!=ref for a in alts)
                           for pair in allele_pairs)
            return 50 if compatible else 40
        if exact and not s['block']:return 35
        if s['block'] and s['gt'].replace('|','/')=='0/0':return 30
        if exact:return 20
        return 10
    best=max(relevance(s) for s in applicable)
    chosen=[s for s in applicable if relevance(s)==best]
    judged=[observe_at(k,s,bases,score) for s in chosen]
    # Distinct genotype claims at identical priority cannot be resolved silently.
    good=[r for r in judged if r['STATUS'] in ('CANDIDATE_DIRECT','CANDIDATE_FLIP')]
    if len(good)>1:
        signature={(r['VCF_REF'],r['VCF_ALT'],r['ALT_DOSAGE'],r['ORIENTATION']) for r in good}
        if len(signature)>1:
            fail(f'Conflicting gVCF records at chr{c}:{pos}: {signature}')
    if good:return good[0]
    # Preserve the most informative failed observation, not an unrelated block.
    for s in ('LOW_GQ','LOW_DP','NON_PASS_FILTER','NON_TARGET_ALT_IN_GT'):
        for r in judged:
            if r['STATUS']==s:return r
    return judged[0]

def find_observations(targets,bases):
    contig=get_contigs();chromosome=defaultdict(list)
    for k in targets:chromosome[k[0]].append(k)
    out={}
    for chrom in sorted(chromosome,key=int):
        ordered=sorted(chromosome[chrom],key=lambda k:k[1]);i=0
        log=ROOT/f'logs/gvcf_chr{chrom}.stderr.log'
        active=[]
        with log.open('w') as err:
            proc=subprocess.Popen(['bcftools','view','-H','-r',contig[chrom],str(GVCF)],
                                  stdout=subprocess.PIPE,stderr=err,text=True)
            last_pos=0
            for line in proc.stdout:
                if not line.strip():continue
                src=source_record(line.rstrip('\n').split('\t'))
                if src['c']!=contig[chrom]:fail('Wrong chromosome returned by bcftools')
                if src['pos']<last_pos:fail(f'gVCF is not position sorted on chr{chrom}')
                last_pos=src['pos']
                # Process targets strictly BELOW current POS. Same-position
                # records are all considered before processing that target.
                while i<len(ordered) and ordered[i][1]<src['pos']:
                    k=ordered[i]
                    out[k]=choose_observation(k,active,bases,targets[k])
                    i+=1
                active=[a for a in active if a['end']>=src['pos']]
                active.append(src)
            if proc.wait()!=0:fail(f'bcftools gVCF streaming failed: chr{chrom}; see {log}')
        while i<len(ordered):
            k=ordered[i]
            out[k]=choose_observation(k,active,bases,targets[k])
            i+=1
        print(f'chr{chrom}: inspected {len(ordered)} score positions',flush=True)
    if set(out)!=set(targets):fail('Score observations incomplete')
    return out

def check_flips(outs):
    # Query all orientation-compatible FLIPs, including low-QC site, in dbSNP157.
    flips=[r for r in outs.values() if r['ORIENTATION']=='FLIP']
    header=subprocess.check_output(['bcftools','view','-h',str(DB)],text=True)
    accessions={}
    for acc in re.findall(r'^##contig=<ID=(NC_\d{6}\.\d+)',header,re.M):
        idx=int(acc[3:9]);
        if 1<=idx<=22:accessions[str(idx)]=acc
    if len(accessions)!=22:fail('dbSNP157 is not compatible with expected GRCh38 RefSeq contigs')
    regions=AUD/'flip_dbsnp_regions.tsv'
    with regions.open('w') as f:
        for r in sorted(flips,key=lambda r:(int(key(r['CHROM'],r['POS'])[0]),int(r['POS']))):
            c,p=key(r['CHROM'],r['POS']);f.write(f'{accessions[c]}\t{p}\n')
    dbs=defaultdict(list)
    if flips:
        query=subprocess.check_output(['bcftools','query','--regions-overlap','pos','-R',str(regions),
                                       '-f','%CHROM\t%POS\t%ID\t%REF\t%ALT\n',str(DB)],text=True)
        for line in query.splitlines():
            c,p,i,ref,alt=line.split('\t')
            dbs[(str(int(c[3:9])),int(p))].append((ids(i),ref.upper(),alt.upper().split(',')))
    for r in flips:
        if r['STATUS']!='CANDIDATE_FLIP':continue
        k=key(r['CHROM'],r['POS']); ref=r['FASTA_REF'];pair={r['COMP_EFFECT'],r['COMP_OTHER']}
        names=ids(r['RSID'])|ids(r['HM_RSID'])
        matched=any(names & i and dref==ref and any({ref,a}==pair for a in alts) for i,dref,alts in dbs.get(k,[]))
        if not matched:r['STATUS']='FLIP_DBSNP_IDENTITY_UNCONFIRMED'
        else:r['STATUS']='FLIP_CONFIRMED'
    return flips



def emit_and_build(targets, observations, skipped, flips):
    names=list(next(iter(observations.values())).keys())
    with DECISIONS.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=names,delimiter='\t');w.writeheader()
        for k in sorted(observations,key=lambda t:(int(t[0]),t[1])):
            w.writerow(observations[k])

    counts=Counter(r['STATUS'] for r in observations.values())
    direct=[r for r in observations.values() if r['STATUS']=='CANDIDATE_DIRECT']
    good_flips=[r for r in observations.values() if r['STATUS']=='FLIP_CONFIRMED']
    accepted=direct+good_flips
    selected={key(r['CHROM'],r['POS']):r for r in accepted}
    if len(selected)!=len(accepted):fail('Duplicate accepted score loci')
    expected=sum((Decimal(r['RECALC_CONTRIBUTION']) for r in accepted),Decimal(0))
    denom=2*len(accepted)
    denom_score=len(targets)+len(skipped)
    overlap=len(accepted)/denom_score

    with DIFF.open('w',newline='') as f:
        w=csv.writer(f,delimiter='\t')
        w.writerow(['CHROM','POS','RSID','STATUS','GVCF_TYPE','SOURCE_GT','SOURCE_DP','SOURCE_GQ','SOURCE_FILTER','STRAND_AMBIGUITY'])
        for k in sorted(observations,key=lambda t:(int(t[0]),t[1])):
            r=observations[k]
            if k not in selected:
                w.writerow([r[c] for c in ['CHROM','POS','RSID','STATUS','GVCF_TYPE','SOURCE_GT','SOURCE_DP','SOURCE_GQ','SOURCE_FILTER','STRAND_AMBIGUITY']])
        for idx, why in skipped:
            w.writerow(['','','',f'SCORE_ROW_{idx}:{why}','','','','','',''])

    flags=Counter(r['STRAND_AMBIGUITY'] for r in accepted if r['STRAND_AMBIGUITY'])
    summary={
        'sample':SAMPLE,'score_id':PGS_ID,'score_rows':denom_score,
        'target_coordinate_loci':len(targets),'unsupported_score_rows':len(skipped),
        'direct':len(direct),'flip_confirmed':len(good_flips),
        'flip_candidates_total':len(flips),
        'matched':len(accepted),'unmatched':denom_score-len(accepted),
        'match_fraction':overlap,'predicted_sum':str(expected),
        'predicted_denom':denom,'predicted_avg':str(expected/Decimal(denom)) if denom else None,
        'orientation_counts':dict(sorted(counts.items())),
        'strand_ambiguity_warnings':dict(flags),
        'min_dp':MIN_DP,'min_gq':MIN_GQ,'palindromic_mode':PALINDROMIC_MODE,
        'flip_requires_dbsnp157_rsid_position_alleles':True,
        'historical_sample_regression':False}
    SUMMARY.write_text(json.dumps(summary,indent=2,ensure_ascii=False)+'\n')
    print('UPSTREAM SUMMARY:',json.dumps(summary,ensure_ascii=False),flush=True)
    if len(accepted)==0:fail('No qualified variants; cannot score')
    if overlap<float(os.environ['MIN_OVERLAP']):
        fail(f'Match rate {overlap:.5%} below minimum {float(os.environ["MIN_OVERLAP"]):.1%}')

    # Scoring allele changes are determined anew for THIS gVCF; no historical
    # list or pre-known 330 sites. Only FLIP_CONFIRMED loci are complemented.
    flipped={key(r['CHROM'],r['POS']):r for r in good_flips}
    changed=set(); seen_rows=0
    with gzip.open(SCORE,'rt') as src,gzip.open(SCORE_OUT,'wt') as dst:
        cols=None
        for line in src:
            if line.startswith('#'):
                dst.write(line);continue
            if cols is None:
                cols=line.rstrip('\n').split('\t');dst.write(line);continue
            x=line.rstrip('\n').split('\t')
            if len(x)!=len(cols):fail('Malformed scoring file')
            r=dict(zip(cols,x));seen_rows+=1
            p=r.get('hm_pos',''); c=r.get('hm_chr','')
            k=key(c,p) if c.removeprefix('chr').isdigit() and p.isdigit() else None
            if k in flipped:
                srcrow=flipped[k]
                if k in changed:fail('Duplicate flip score row')
                if r['effect_allele'].upper()!=srcrow['EFFECT_ALLELE'] or r['other_allele'].upper()!=srcrow['OTHER_ALLELE']:
                    fail(f'Original score effect/other mismatch: {k}')
                x[cols.index('effect_allele')]=srcrow['COMP_EFFECT']
                x[cols.index('other_allele')]=srcrow['COMP_OTHER']
                changed.add(k)
            dst.write('\t'.join(x)+'\n')
    if seen_rows!=denom_score or set(flipped)!=changed:
        fail('Scoring file row/correction count inconsistency')

    with UNSORTED.open('w') as dst:
        dst.write('##fileformat=VCFv4.2\n')
        for i in range(1,23):dst.write(f'##contig=<ID=chr{i}>\n')
        dst.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        dst.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth or minimum block depth">\n')
        dst.write('##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">\n')
        dst.write(f'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{SAMPLE}\n')
        for k in sorted(selected,key=lambda t:(int(t[0]),t[1])):
            r=selected[k]
            dp=r['SOURCE_MIN_DP'] if r['GVCF_TYPE']=='BLOCK' and r['SOURCE_MIN_DP'] not in ('','.','NA') else r['SOURCE_DP']
            ident=r['RSID'] if str(r['RSID']).startswith('rs') else '.'
            dst.write('\t'.join(map(str,[r['CHROM'],r['POS'],ident,r['VCF_REF'],r['VCF_ALT'],'.','PASS','.',
                          'GT:DP:GQ',f'{r["VCF_GT"]}:{int(dp)}:{int(r["SOURCE_GQ"])}']))+'\n')
    print(f'INPUT BUILD PASS: {len(accepted)} variants; {len(flipped)} dbSNP-confirmed FLIPs; {denom_score-len(accepted)} unmatched')
    if flags:
        print('CAUTION: palindromic strand orientation not independently resolved; see strand_ambiguity_warnings in summary')
    return summary

if __name__=='__main__':
    targets,skipped=score_targets()
    print('Score rows:',len(targets)+len(skipped),'unique target loci:',len(targets),flush=True)
    bases=fasta_at_targets(targets)
    observations=find_observations(targets,bases)
    flips=check_flips(observations)
    emit_and_build(targets,observations,skipped,flips)

PY_GVCF

python3 "$OUT/scripts/gvcf_to_pgs.py"
if (( DISCOVER_ONLY ));then
  echo "DISCOVERY COMPLETE: $OUT/audit/upstream_reconstruction_summary.json"
  exit 0
fi

printf '\n### Verify VCF and genotype round-trip\n'
VCF="$OUT/prep/score_ready.vcf.gz"
bcftools sort -Oz -o "$VCF" "$OUT/prep/score_ready.unsorted.vcf"
bcftools index -f -t "$VCF"

cat > "$OUT/scripts/verify_and_reconcile.py" <<'PY_VERIFY'
#!/usr/bin/env python3
import csv,gzip,json,os,re,sys
from pathlib import Path
from decimal import Decimal
P=Path(os.environ['OUT']); SAMPLE=os.environ['SAMPLE'];SAMPLESET=os.environ['SAMPLESET']
SCORE=P/'scorefiles/PGS004879_hmPOS_GRCh38.txt.gz'
VCF=P/'prep/score_ready.vcf.gz'
INFO=json.loads((P/'audit/upstream_reconstruction_summary.json').read_text())
TOL=Decimal('1e-12')
def bad(s):raise RuntimeError(s)
def key(c,pos):return(str(c).removeprefix('chr'),int(pos))
def dosage(fmt,s):
    pairs=dict(zip(fmt.split(':'),s.split(':')))
    gt=pairs.get('GT','')
    a=re.split(r'[/|]',gt)
    if len(a)!=2 or any(x not in ('0','1') for x in a):bad(f'Invalid projected GT {gt}')
    return a.count('1')
def vcf(path):
    with gzip.open(path,'rt') as f:
        for line in f:
            if not line.startswith('#'):
                x=line.rstrip('\r\n').split('\t')
                if len(x)!=10:bad('VCF has wrong number of sample columns')
                yield x

def input_qc():
    score={}
    with gzip.open(SCORE,'rt') as f:
        rd=csv.DictReader((x for x in f if not x.startswith('#')),delimiter='\t')
        for r in rd:
            c=r.get('hm_chr','').removeprefix('chr');p=r.get('hm_pos','')
            if not (c.isdigit() and p.isdigit() and 1<=int(c)<=22):continue
            k=key(c,p)
            if k in score:bad(f'Duplicate score locus {k}')
            score[k]=r
    tally=Decimal(0);found=set();records={}
    for x in vcf(VCF):
        k=key(x[0],x[1]);ref=x[3];alt=x[4]
        if k in found:bad(f'Duplicate genotype at {k}')
        if k not in score:bad(f'VCF locus missing score: {k}')
        found.add(k)
        sc=score[k];ea=sc['effect_allele'];oa=sc['other_allele']
        if {ref,alt}!={ea,oa} or len(ref)!=1 or len(alt)!=1:bad(f'Allele mismatch at {k}')
        d=dosage(x[8],x[9]); eff=d if ea==alt else 2-d
        tally+=Decimal(sc['effect_weight'])*Decimal(eff)
        records[k]=(ref,alt,d)
    if len(found)!=INFO['matched']:bad('Input record count differs from discovery')
    if abs(tally-Decimal(INFO['predicted_sum']))>TOL:bad('Recalculated additive SUM differs from discovery')
    if 2*len(found)!=INFO['predicted_denom']:bad('Expected denominator mismatch')
    print(f'INPUT QC PASS: {len(found)} exact loci; independently recomputed SUM {tally}')
    return records

def roundtrip():
    orig=input_qc();seenpv=set();seenvcf=set();logs=[]
    for n in range(1,23):
        pref=P/f'roundtrip_22chr/pgen/chr{n}'
        if not (pref.with_suffix('.pgen').exists()):
            if any(k[0]==str(n) for k in orig):bad(f'Missing chr{n} PGEN')
            continue
        n1=n2=0
        with open(str(pref)+'.pvar') as f:
            for line in f:
                if line.startswith('#'):continue
                x=line.rstrip('\n').split('\t');k=key(x[0],x[1])
                if k in seenpv or k not in orig or (x[3],x[4])!=orig[k][:2]:bad(f'PVAR allele error {k}')
                seenpv.add(k);n1+=1
        for x in vcf(P/f'roundtrip_22chr/export/chr{n}.vcf.gz'):
            k=key(x[0],x[1])
            if k in seenvcf or k not in orig or (x[3],x[4],dosage(x[8],x[9]))!=orig[k]:
                bad(f'Roundtrip mismatch at {k}')
            seenvcf.add(k);n2+=1
        if n1!=n2:bad(f'chr{n} PVAR/export count mismatch')
        logs.append((f'chr{n}',n1,n2,'PASS'))
    if seenpv!=set(orig) or seenvcf!=set(orig):bad('Missing chromosome roundtrip loci')
    with open(P/'audit/roundtrip_22chr.tsv','w') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(['CHROM','PVAR','EXPORTED','STATUS']);w.writerows(logs)
    print(f'ROUNDTRIP PASS: {len(orig)}/{len(orig)}; {len(logs)} chromosome(s); zero errors')

def reconcile():
    path=P/f'results/{SAMPLESET}/score/aggregated_scores.txt.gz'
    if not path.exists():bad(f'Missing pgsc_calc aggregated scores: {path}')
    with gzip.open(path,'rt') as f: rows=list(csv.DictReader(f,delimiter='\t'))
    if len(rows)!=1:bad(f'Expected one sample and one score; got {len(rows)} rows')
    row=rows[0]
    if row['IID']!=SAMPLE or row['FID']!=SAMPLE:bad('Sample ID changed during pgsc_calc')
    got=Decimal(row['SUM']); predicted=Decimal(INFO['predicted_sum'])
    if abs(got-predicted)>TOL:bad(f'SUM mismatch: observed={got} predicted={predicted}')
    den=int(row['DENOM'])
    if den!=INFO['predicted_denom']:bad(f'DENOM mismatch {den} != {INFO["predicted_denom"]}')
    if abs(Decimal(row['AVG'])-got/Decimal(den))>Decimal('1e-16'):bad('AVG mismatch')
    expected_sum=os.environ.get('EXPECT_SUM')
    expected_denom=os.environ.get('EXPECT_DENOM')
    if expected_sum and abs(got-Decimal(expected_sum))>TOL:bad('Optional external expected SUM mismatch')
    if expected_denom and den!=int(expected_denom):bad('Optional external expected DENOM mismatch')
    html=P/f'results/{SAMPLESET}/score/report.html'
    if not html.is_file() or html.stat().st_size==0:bad('Missing HTML report')
    result={**INFO,'status':'PASS','observed_sum':row['SUM'],
            'observed_denom':den,'observed_avg':row['AVG'],'iid':row['IID'],
            'sum_discrepancy':str(abs(got-predicted)),'roundtrip_loci':INFO['matched'],
            'score_report':str(html),'algorithm_validation':'per_case_self_consistency',
            'independent_other_sample_validation':'NOT_YET_ESTABLISHED'}
    (P/'audit/final_reconciliation.json').write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n')
    print('FINAL RECONCILIATION PASS:',json.dumps({'sum':row['SUM'],'denom':den,'matched':INFO['matched'],'unmatched':INFO['unmatched'],'sum_error':str(abs(got-predicted))}))

if __name__=='__main__':
    {'input':input_qc,'roundtrip':roundtrip,'reconcile':reconcile}[sys.argv[1]]()
PY_VERIFY

python3 "$OUT/scripts/verify_and_reconcile.py" input
RUN="$OUT/roundtrip_22chr"
for N in $(seq 1 22);do
  CHR="chr${N}"; CV="$RUN/input/$CHR.vcf.gz";P="$RUN/pgen/$CHR";EXP="$RUN/export/$CHR"
  bcftools view -r "$CHR" -Oz -o "$CV" "$VCF"
  CNT="$(bcftools view -H "$CV" | wc -l)"
  if [[ "$CNT" -eq 0 ]];then echo "Skip empty $CHR";rm -f "$CV";continue;fi
  echo "Processing $CHR ($CNT variants)"
  if ! apptainer exec --bind /home:/home "$PLINK_IMG" plink2 \
     --vcf "$CV" --make-pgen --threads "$THREADS" --memory "$MEMORY_MB" \
     --out "$P" > "$RUN/logs/$CHR.import.log" 2>&1;then
     tail -70 "$RUN/logs/$CHR.import.log";exit 1
  fi
  if ! apptainer exec --bind /home:/home "$PLINK_IMG" plink2 \
     --pfile "$P" --export vcf bgz id-paste=iid --threads "$THREADS" --memory "$MEMORY_MB" \
     --out "$EXP" > "$RUN/logs/$CHR.export.log" 2>&1;then
     tail -70 "$RUN/logs/$CHR.export.log";exit 1
  fi
done
python3 "$OUT/scripts/verify_and_reconcile.py" roundtrip

SHEET="$OUT/samplesheet.csv"
printf 'sampleset,path_prefix,chrom,format\n' > "$SHEET"
for N in $(seq 1 22);do
  if [[ -s "$RUN/pgen/chr${N}.pgen" ]];then
    printf '%s,%s,%s,pfile\n' "$SAMPLESET" "$RUN/pgen/chr${N}" "$N" >> "$SHEET"
  fi
done
[[ $(wc -l < "$SHEET") -ge 2 ]] || { echo 'No chromosomes in samplesheet' >&2;exit 2; }
echo 'Starting pgsc_calc v2.3.0'
(
  cd "$OUT"
  nextflow run pgscatalog/pgsc_calc -r v2.3.0 -profile singularity -resume \
     -c "$OUT/local_containers.config" -w "$OUT/nextflow_work" \
     --input "$SHEET" --target_build GRCh38 \
     --scorefile "$OUT/scorefiles/PGS004879_hmPOS_GRCh38.txt.gz" \
     --outdir "$OUT/results" --skip_ancestry \
     > "$OUT/logs/pgsc_calc.log" 2>&1
) || { echo 'pgsc_calc failed:';tail -90 "$OUT/logs/pgsc_calc.log";exit 1; }
python3 "$OUT/scripts/verify_and_reconcile.py" reconcile
sha256sum "$VCF" "$OUT/scorefiles/PGS004879_hmPOS_GRCh38.txt.gz" > "$OUT/audit/output_sha256.txt"
echo "FINAL PASS: $OUT/audit/final_reconciliation.json"
