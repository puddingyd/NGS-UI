"""Attempt-scoped runtime evidence; unknown versions are never inferred."""
import datetime
import os
from pathlib import Path
import re
import shutil
import subprocess

from .common import sha, plink


def probe(argv):
    executable = shutil.which(str(argv[0]))
    if not executable:
        return dict(status='unavailable', argv=list(map(str,argv)))
    try:
        p = subprocess.run([executable, *map(str,argv[1:])], capture_output=True, timeout=30)
        text = (p.stdout+p.stderr).decode('utf-8',errors='replace')[:12000]
        return dict(status='observed' if p.returncode == 0 else 'failed', executable=executable,
                    exit_code=p.returncode, output=text, argv=list(map(str,argv)))
    except (OSError,subprocess.TimeoutExpired) as exc:
        return dict(status='failed', executable=executable, error=type(exc).__name__)


def collect(cfg, gvcf, stage, root):
    info = dict(schema='PRS_RUNTIME_PROVENANCE_V1', stage=stage,
        inspected_at=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        meaning='Tools inspected for this attempt; not retroactive evidence of an older execution',
        caller=dict(status='unknown', header=[]), tools={})
    info['code_sha256'] = {p.name:sha(p) for p in Path(__file__).parent.glob('*.py')}
    info['repository_commit'] = probe(['git','-C',str(Path(__file__).resolve().parents[3]),'rev-parse','HEAD'])
    for name,args in [('bcftools',['--version']),('samtools',['--version']),
                      ('nextflow',['-version']),('apptainer',['--version']),('singularity',['--version'])]:
        info['tools'][name] = probe([name,*args])
    if gvcf and shutil.which('bcftools'):
        try:
            p = subprocess.run(['bcftools','view','-h',str(gvcf)],capture_output=True,timeout=30)
            if p.returncode == 0:
                # Do not collect sample IDs or command lines containing identifiers.
                lines = [s[:2000] for s in p.stdout.decode(errors='replace').splitlines()
                         if re.match(r'##(?:source|reference|DRAGENVersion|GATKVersion|DeepVariantVersion)=',s,re.I)]
                info['caller'] = dict(status='observed' if lines else 'unknown',header=lines,
                    version_status='observed' if any(re.search(r'\d+\.\d+',s) for s in lines if 'reference=' not in s.lower()) else 'unknown')
        except (OSError,subprocess.TimeoutExpired):
            pass
    info['resources'] = {}
    for name in ('plink_img','ancestry_nextflow_config'):
        p = Path(cfg.get(name,''))
        info['resources'][name] = (dict(status='observed',path=str(p),sha256=sha(p))
            if p.is_file() else dict(status='unavailable'))
    info['tools']['plink_container'] = probe(plink(cfg,'--version')) if cfg.get('plink_img') else dict(status='unavailable')
    assets = Path(os.environ.get('NXF_ASSETS',str(Path(os.environ.get('NXF_HOME',str(Path.home()/'.nextflow')))/'assets')))
    checkout = assets/'pgscatalog/pgsc_calc'
    info['nextflow_pipeline'] = dict(requested='v2.3.0', resolved_commit=probe(['git','-C',str(checkout),'rev-parse','HEAD'])
        if (checkout/'.git').exists() else dict(status='unavailable'))
    return info


def record_launch(info, root):
    log = Path(root)/'ancestry/.nextflow.log'
    if log.is_file():
        text=log.read_text(errors='replace')
        launches = re.findall(r'Launching[^\n]+revision:\s*([^\n]+)',text)
        info['nextflow_pipeline']['launch_records'] = launches[-1:]
        images=sorted(set(re.findall(r'''(/[^\s"'<>]+\.(?:sif|img))\b''',text)))
        info['nextflow_pipeline']['logged_containers']=[dict(path=p,sha256=sha(p)) for p in images if Path(p).is_file()]
        info['nextflow_pipeline']['container_scope']='Only image paths recorded in this launch log; absent entries do not prove an exhaustive container inventory'
