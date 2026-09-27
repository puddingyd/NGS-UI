# In-house allele frequency (in-house AF)

Build a cohort allele-frequency database from our own NovaSeq + DRAGEN WGS
gVCFs, and annotate the pipeline SNV TSV with `INHOUSE_AF` the same way
we annotate `GNOMAD_G_AF`. Goal: an in-house AF (incl. **rare** variants) for
annotation, and a filter for locally common polymorphisms that gnomAD
under-represents.

> **Scope:** SNV/indel only. SV in-house frequency (breakpoint /
> reciprocal-overlap matching) is a later phase. Mito is separate.

**目前狀態（2026-09）**：cohort **1397 隻**、**57,171,081** 個站點、AN = **2794**。

---

# 批次更新 SOP（新樣本下機後照這個做）

> 這是**操作手冊**。下面「Cohort / Files / Phase …」是設計與歷史紀錄，平常不用看。

## 機器角色

| 機器 | 位置 | 負責 |
|---|---|---|
| **DGX**（`n102968@dgx2`） | DB 在 `/raid/DGM/n102968/inhouse_af`、script 在 `~/dgx_stage/inhouse_af` | 建資料庫（ingest / accumulate / publish）。有 1.5 TB RAM，無外網 |
| **DGM**（`n102968@server`） | git checkout `~/NGS_UI/inhouse-af` | 從 GitHub 拉 code，再中繼給 DGX |
| **NGS-UI 主機** | `~/NGS_UI/NGS-UI`、DB 裝在 `~/NGS_UI/biotools/inhouse_af/` | 部署 DB、backfill 既有樣本 |

DGM 與 DGX 之間**沒有 ssh**，靠共用的 datalake 當中繼。

## Step 0｜環境變數（每次開新 shell 都要設）

**DGX：**
```bash
export PATH=$HOME/bin:$PATH          # bcftools/bgzip/tabix 的 apptainer wrapper
DB=/raid/DGM/n102968/inhouse_af
SCR=~/dgx_stage/inhouse_af
for c in /datalake_Intermediate/pipeline/reference/hg38/Homo_sapiens_assembly38.fasta \
         /datalake_Intermediate/datalake_Intermediate/pipeline/reference/hg38/Homo_sapiens_assembly38.fasta; do
  [ -f "$c" ] && REF="$c" && break
done; echo "REF=$REF"
```
- ✅ 正常：印出 `REF=/datalake_Intermediate/...fasta`
- ❌ `REF=` 空的 → datalake 掛載點又變了，用 `find /datalake* -name 'Homo_sapiens_assembly38.fasta' 2>/dev/null | head` 找

## Step 1｜同步最新 script 到 DGX

**DGM：**
```bash
cd ~/NGS_UI/inhouse-af
git status --short | head                       # 先確認沒有未提交的修改
git checkout claude/plan-ngs-ui-RQW8J
git pull --ff-only origin claude/plan-ngs-ui-RQW8J
scripts/inhouse_af/sync_to_dgx.sh --repo ~/NGS_UI/inhouse-af \
  --via /home/datalake_Intermediate/n102968/_sync_inhouse_af
```
- ✅ 會印 `branch : claude/plan-ngs-ui-RQW8J`、`HEAD -> <commit>`，最後印一行 DGX 要跑的 `cp`
- ❌ `Not possible to fast-forward, aborting` → checkout 停在舊分支。用上面的 `git checkout` 切過去再 pull。**不要用 `--no-pull` 硬跳過**，那會把舊 script 送出去（踩過一次）

**DGX**（貼上它印出來的那行）：
```bash
mkdir -p ~/dgx_stage/inhouse_af && cp -f /datalake_Intermediate/n102968/_sync_inhouse_af/*.{py,sh,yml,txt} ~/dgx_stage/inhouse_af/
python3 "$SCR/accumulate.py" --selftest
python3 "$SCR/annotate_inhouse_af.py" --selftest
```
- ✅ 兩個都印 `selftest OK`
- ❌ `error: the following arguments are required: --db-dir` → 拿到的是舊版，回 DGM 確認分支

## Step 2｜掃描有哪些新檢體（唯讀，不動 DB）

**DGX：**
```bash
find /datalake_Raw/datalake_Raw/Novaseq -maxdepth 5 -name '*.hard-filtered.gvcf.gz' \
  -printf '%s\t%p\n' > "$DB/gvcf_sizes.new.txt"

python3 "$SCR/select_cohort.py" --sizes "$DB/gvcf_sizes.new.txt" \
  --exclude-range 'VAL-:37-54' \
  --out-manifest "$DB/cohort_manifest.new.tsv" --out-list "$DB/cohort_gvcfs.new.txt"

sqlite3 "$DB/counts.sqlite" "SELECT sample_id FROM samples;" | sort > "$DB/.ingested.txt"
awk -F'\t' 'NR>1 && $4=="include"{print $1}' "$DB/cohort_manifest.new.tsv" | sort > "$DB/.included.txt"
comm -23 "$DB/.included.txt" "$DB/.ingested.txt" > "$DB/new_samples.txt"

echo "新 cohort: $(wc -l < "$DB/.included.txt")  DB 現有: $(wc -l < "$DB/.ingested.txt")  要新增: $(wc -l < "$DB/new_samples.txt")"
echo "--- DB 有但 cohort 沒有（必須是空的）---"; comm -13 "$DB/.included.txt" "$DB/.ingested.txt"
```
- ✅ `EXCLUDE_standard_ref 18`、`EXCLUDE_broken_empty 1`（排除規則有作用）；最後一行**空的**
- ⚠️ 最後一行**不是空的** → 有樣本從 datalake 消失或改名。**先停下來釐清**：AN track 是累加的，沒辦法用增量移除樣本
- ⚠️ 新增清單裡有不該收的對照品 → 加進 `--exclude-range` 或 `--exclude-id-file` 再跑一次

## Step 3｜ingest 新檢體

**DGX：**
```bash
awk -F'\t' 'NR==FNR{new[$1];next} FNR>1 && $4=="include" && ($1 in new){
  if (match($5, /\/(vcf\.gz|other)\//)) print substr($5,1,RSTART-1) "/other/" $1;
}' "$DB/new_samples.txt" "$DB/cohort_manifest.new.tsv" > "$DB/new_other_dirs.txt"

# 驗證每個目錄都有 gVCF 與 ploidy CSV（遞迴找，因新 run 用 germline_seq/ 巢狀）
miss_g=0; miss_p=0
while read -r d; do
  [ -n "$(find "$d" -name '*.hard-filtered.gvcf.gz' -print -quit 2>/dev/null)" ] || miss_g=$((miss_g+1))
  [ -n "$(find "$d" -name '*ploidy_estimation_metrics.csv' -print -quit 2>/dev/null)" ] || miss_p=$((miss_p+1))
done < "$DB/new_other_dirs.txt"
echo "dirs=$(wc -l < "$DB/new_other_dirs.txt")  missing gVCF=$miss_g  ploidy CSV=$miss_p"

nohup "$SCR/ingest_batch.sh" --dirs-file "$DB/new_other_dirs.txt" \
  --ref "$REF" --out-dir "$DB/per_sample" --jobs 16 > "$DB/ingest.log" 2>&1 &
```
- ✅ 兩個 missing 都是 **0**；log 持續出現 `[ok] <id>`
- ❌ `missing ploidy CSV` > 0 → 性別判不出來會被當 ambiguous（X/Y 不計入），先查那些目錄的結構
- 進度：`echo "$(ls "$DB/per_sample"/*/qc.json | wc -l) / <新 cohort 總數>"`；`grep -c '\[FAIL\]' "$DB/ingest.log"` 要維持 0
- 可中斷重跑（已有 `qc.json` 的會跳過）
- ⏱ 720 隻約數小時

ingest 完做一次性別檢查：
```bash
for s in "$DB"/per_sample/*/qc.json; do grep -o '"sex_class": *"[a-z]*"' "$s"; done | sort | uniq -c
```
- ✅ female/male 各佔多數，`ambiguous` 只有零星幾隻且散落在不同 run
- ⚠️ `ambiguous` 集中成一整批（幾十上百隻）→ 那批的 ploidy CSV 沒被讀到，**先別 accumulate**

## Step 4｜accumulate + publish

**DGX：**
```bash
nohup bash -c "
  '$SCR/accumulate.py' --db-dir '$DB' --sort-tmp '$DB/.sorttmp' --jobs 16 &&
  '$SCR/publish_af.py' --db-dir '$DB' --ref '$REF'
" > "$DB/rebuild.log" 2>&1 &
tail -f "$DB/rebuild.log"
```
依序會看到：
```
[accumulate] adding N sample(s) (have M)
[accumulate]   + <sample>  (…… variant rows)      ← 每隻一行
[accumulate] rebuilding AN track (old + N); sort tmp=…
[accumulate]   phase 1/3 demux: N inputs, 16 workers (…)
[accumulate]   phase 2/3 reduce: 25 chromosomes, 16 in parallel
[accumulate]   phase 3/3 concat + bgzip
[accumulate] AN track -> …/an_track.bg.gz  (NNNNs)
[accumulate] cohort=<總數> samples, <變異數> distinct variants
[publish] <站點數> sites written, <少量> dropped (AN=0/AC=0)
```
- ⏱ 增量（`old + N`）比全量快很多；全量 1397 隻約 8 小時
- **phase 1 沒有輸出不代表當掉**：進度看 `du -sh "$DB"/.sorttmp/an_track.*`；讀到第幾份 BED 用
  ```bash
  WD=$(ls -d "$DB"/.sorttmp/an_track.* | head -1)
  for p in $(pgrep -x gzip); do tr '\0' '\n' < /proc/$p/cmdline 2>/dev/null | grep -o 'per_sample/[^/]*'; done | head
  ```
- ❌ **中途失敗（OOM／磁碟滿／被 kill）絕對不要直接重跑** `accumulate.py`：`counts.sqlite` 已經寫入新樣本、`an_track` 還是舊的，直接重跑會印 `AN track unchanged` 然後**所有 AF 靜靜地偏高**。正確復原：
  ```bash
  "$SCR/accumulate.py" --db-dir "$DB" --sort-tmp "$DB/.sorttmp" --rebuild-an-track --jobs 16
  ```
  （`--rebuild-an-track` 會從全部樣本重建，較久但正確）

## Step 5｜驗收（重要，2 分鐘）

**DGX：**
```bash
ls -lah "$DB/an_track.bg.gz" "$DB/inhouse_af.hg38.vcf.gz"
tabix "$DB/inhouse_af.hg38.vcf.gz" chr1:1000000-1100000 | head -3
tabix "$DB/inhouse_af.hg38.vcf.gz" chr1:45330228-45330228   # MUTYH 常見 indel
```
**判準**：`INHOUSE_AN` 要接近 **2 × 樣本數**；而且拿一個**常見變異**對照上一版——
> **AF 應該幾乎不變，AC 和 AN 各自等比例增加。**

這是最有力的檢查：常見變異的頻率不該因為樣本數變多而改變。若 AF 大幅偏移，代表 AN track 與 counts 不同步（見 Step 4 的復原）。

## Step 6｜部署到 NGS-UI 主機

**DGX：**
```bash
mkdir -p /datalake_Intermediate/n102968/_deploy_inhouse_af
cp -f "$DB/inhouse_af.hg38.vcf.gz" "$DB/inhouse_af.hg38.vcf.gz.tbi" \
      /datalake_Intermediate/n102968/_deploy_inhouse_af/
```
**NGS-UI 主機：**
```bash
cd ~/NGS_UI/NGS-UI && git pull --ff-only origin claude/plan-ngs-ui-RQW8J
scripts/inhouse_af/deploy_inhouse_af_db.sh \
  /home/datalake_Intermediate/n102968/_deploy_inhouse_af/inhouse_af.hg38.vcf.gz
```
- ✅ 會驗 bgzip、沿用或重建 `.tbi`、探測一個站點、印站點數，最後**原子換檔**
- ✅ 裝到預設路徑就**不用設 `NGS_UI_INHOUSE_AF_DB`**
- ❌ `compression ... is not BGZF` → 來源檔壞了，回 DGX 重新 publish

## Step 7｜backfill 既有樣本

**cohort 一變，所有既有樣本的 AF 都過期了，一定要重跑。**

**NGS-UI 主機：**
```bash
python3 scripts/annotate_inhouse_af.py --selftest      # 先確認 code 是新的
python3 scripts/backfill_inhouse_af.py --selftest      # overlay 直接更新 == 參考流程

python3 scripts/backfill_inhouse_af.py --dry-run | head     # 看解析出的路徑對不對
time python3 scripts/backfill_inhouse_af.py <某一隻 SID>    # 先測一隻並計時

nohup python3 scripts/backfill_inhouse_af.py --continue-on-error > ~/NGS_UI/backfill.log 2>&1 &
tail -f ~/NGS_UI/backfill.log
```
腳本會**依「解析出來的那個檔案是什麼」**自動選模式（不是看 layout 標記）：

**Overlay 模式**（解析到的是 pipeline 的 `03_acmg/*.snv_indel.acmg.tsv`，不論在 unified 還是舊 pipeline root）
——`03_acmg` 是**唯讀的原始真相**，匯出報告讀它，**絕不寫入**：
1. 讀 `03_acmg` TSV（唯讀）建 allele key → 串流掃一次 in-house DB 找命中
2. **直接在** `08_postprocessing/<sid>.snv_annotations.sqlite` 裡、同一個 SQLite transaction：
   先把所有舊的 `INHOUSE_*` 拿掉（DB 更新後已過期；payload 變空的列整列刪），
   再把新值 **merge 進每列既有的 payload**——GeneBe / SpliceAI / MANE / LitVar2 等其他欄位完全不動
3. 更新 overlay meta 的欄位清單；raw 沒變，所以 source signature 仍有效（`is_current()` 仍成立）
4. 重建 review TSV
5. **gene index 不用重建**（從未變動的 raw 建的，byte offset 沒動）

> 舊做法是「raw + overlay 還原成完整 working TSV → 註解 → 用 `build_overlay()` 重新 diff」，結果相同，
> 但每隻 WGS 要在 NFS 上寫出／讀回好幾倍 TSV 大小的資料，一隻 15 分鐘以上且中途完全沒輸出。
> 現在不產生 working TSV。`--selftest` 會拿舊做法當參考，驗證直接更新的 overlay 逐列相同、
> raw md5 不變、其他註解保留、過期的 `INHOUSE_*` 被清掉、重跑結果相同（idempotent）。

**In-place 模式**（解析到的是舊 UI 副本 `<NGS_UI_HOME>/tertiary_output/<sample>/snv_indel.annotated.tsv`）
——那份本來就是註解後的副本，所以原地改寫，而且**必須一起重建 gene index**。

> ⚠️ 判斷依據是**檔案本身**，不是 `uses_unified_layout()`。有些樣本標記為 legacy，但 `snv_raw_tsv()` 仍會回傳舊 pipeline root 下的 `03_acmg` 檔——那還是原始真相，不能改。腳本另外有硬性防護：`do_inplace()` 一旦看到路徑含 `03_acmg` 就直接報錯中止。

`--dry-run` 會逐隻印出模式與路徑，`raw` 那行標 `(read-only)` 或 `(REWRITTEN IN PLACE)`。**全跑前先確認沒有任何 `03_acmg` 的路徑被標成 REWRITTEN IN PLACE。**

每隻樣本會印帶累計秒數的進度行，正常的 overlay 模式長這樣：
```
  • 26T00028-dragen  [overlay]
      [    0.0s] scan raw TSV (4.21 GB)
      [   55.3s] 5,312,004 rows, 5,401,877 allele keys; join in-house DB
      [  140.2s] 5,301,550 rows matched; update overlay
      [  260.8s] overlay updated (5,301,550 rows with INHOUSE_AF); rebuild review TSV
      [  330.1s] done
```
（數字只是示意；秒數是從腳本啟動起算的累計值。）

- ✅ matched / rows 比例約 **99.9%**（DRAGEN 樣本）；in-place 模式另印 `[inhouse-af] N variants, M matched in-house AF DB`
- ⚠️ `-nckuh` 樣本配對率較低（約 87%）是**正常的**：DB 用 DRAGEN gVCF 建的，in-house pipeline 的 variant caller 表示法與變異集合不同
- ❌ `overlay is stale (raw TSV changed after post-processing)` → `03_acmg` 在上次 post-processing 之後被重跑過，UI 本來就會忽略這份 overlay；先重跑該樣本的三級 post-processing，再 backfill。overlay 不會被動到
- ❌ `no SNV TSV (...), skip` → 該樣本沒有 `03_acmg/*.snv_indel.acmg.tsv`，通常是三級分析沒跑完
- ❌ `in-house AF DB not found` → Step 6 沒做或路徑不對
- 預設**一失敗就停**；要跳過壞樣本繼續用 `--continue-on-error`
- 中途被砍（Ctrl-C / kill）是安全的：overlay 更新是單一 transaction，沒 commit 就整個 rollback，
  不會留下一半新一半舊的 AF；重跑即可
- ⏱ 先用單隻 `time` 的結果乘以樣本數估總時間（本機 1M 列合成資料約 35 秒；WGS 主要花在讀 raw TSV 兩次與掃 DB 一次）

## Step 8｜重啟與確認

```bash
# 重啟 NGS-UI 服務（載入新的 adapter / config）
```
瀏覽器 **Ctrl-Shift-R** 強制重新整理，開一個個案搜常見基因：
- ✅ SNV 卡片的 `AF_nckuh` 顯示成 `0.44202 (1235/2794)`，括號裡的 AN = 2 × 樣本數
- ✅ 不在 DB 的變異顯示乾淨的 `—`
- ❌ 顯示 `— (0/0)` → adapter 是舊版（缺 `_first_num`），確認主機 checkout 有 pull 到且服務已重啟
- ❌ 數字還是舊的 → 那隻樣本還沒 backfill，或瀏覽器快取沒清

---

## Cohort

613 WGS samples (one DRAGEN gVCF each), selected from the raw datalake by
`select_cohort.py`:

- de-duplicated to one path per sample (same sample is staged under both
  `.../<run>/other/<sample>/` and `.../<run>/vcf.gz/`; we prefer the
  `vcf.gz` copy),
- excluded 1 broken/empty gVCF (`< 50 MB`),
- excluded 18 standard-reference controls (`VAL-37`..`VAL-54`).

No MRN/family de-duplication yet — a patient sequenced twice under different
ids, or related individuals, still count more than once. Acceptable for now
(small effect on *common* variants; revisit before using in-house AF as
population-frequency evidence). The manifest carries `sample_id` + `run` so
MRN/family columns can be bolted on later.

## Files

| file | what |
|---|---|
| `select_cohort.py` | gVCF path list → deduplicated, filtered cohort (`--out-list`) + audit manifest. Stdlib only. |
| `build_inhouse_af.sh` | cohort gVCFs → joint genotyping (GLnexus) → `inhouse_af.hg38.vcf.gz` (sites-only, `INHOUSE_{AC,AN,AF,NHOM}`). |
| `update_inhouse_af.sh` | per-batch SOP: append new samples (dedup), re-genotype the full cohort, atomic-swap the AF DB. Maintains `cohort_manifest.tsv` + `updates.log`. |
| `ingest_sample.py` | **incremental Phase A** — per-sample tool: DRAGEN ploidy → sex, stream gVCF → ploidy-weighted callable BED (DP≥10) + AC/hom/het/hemi + chrM carrier. Writes `per_sample/{id}/`. See `DESIGN_incremental.md`. |
| `accumulate.py` | **incremental Phase B (1/2)** — fold per-sample output into `counts.sqlite` (UPSERT) + `an_track.bg.gz` (event/delta rebuild). Idempotent dedup via the `samples` table. |
| `publish_af.py` | **incremental Phase B (2/2)** — merge-join counts × AN track → `inhouse_af.hg38.vcf.gz` (the `INHOUSE_*` sites VCF). |
| `DESIGN_incremental.md` | full design of the incremental (per-sample) production path + reconciliation. |

The cohort list / manifest / AF DB all contain patient sample ids and
datalake paths — **keep them out of git** (put them under
`$NGS_UI_HOME/biotools/inhouse_af/`, like the gnomAD/GeneBe DBs). Only these
scripts are committed.

## Phase 0 runbook (validate the toolchain first)

Run on the cluster head node where the datalake + reference live.

### 1. Select the cohort

```bash
# the size list is `du -h .../*.hard-filtered.gvcf.gz` over the run dirs
scripts/inhouse_af/select_cohort.py \
  --sizes hard_filtered_gvcf_sizes.txt \
  --exclude-range 'VAL-:37-54' \
  --out-manifest $NGS_UI_HOME/biotools/inhouse_af/cohort_manifest.tsv \
  --out-list     $NGS_UI_HOME/biotools/inhouse_af/cohort_gvcfs.txt
# → 613 included across 10 runs
```

### 2. Smoke test on ONE run (~64 samples) before the full cohort

```bash
# GLnexus — prefer the static binary (single file, zero deps, air-gap friendly):
#   download glnexus_cli from https://github.com/dnanexus-rnd/GLnexus/releases
export GLNEXUS_BIN=$NGS_UI_HOME/biotools/glnexus/glnexus_cli   # OR GLNEXUS_SIF=...
# bcftools: host binary on PATH is fine; or BCFTOOLS_SIF=<bcftools.sif> on DGM.

# restrict to one run by grepping the list, OR use --max-samples for a quick run
grep 20251118_LH00873_0004 $NGS_UI_HOME/biotools/inhouse_af/cohort_gvcfs.txt \
  > /tmp/run1_gvcfs.txt

scripts/inhouse_af/build_inhouse_af.sh \
  --list /tmp/run1_gvcfs.txt \
  --out  /tmp/inhouse_af.run1.vcf.gz \
  --threads 16 --mem-gbytes 96
```

The reference defaults to
`/home/datalake_Intermediate/pipeline/reference/hg38/Homo_sapiens_assembly38.fasta`
(DGM). Override with `--ref` when it lives elsewhere.

**What to check (this is the Phase 0 validation):**

1. **Config: use `glnexus_config_dragen.yml` (the default).** The built-in
   `gatk` preset sets `revise_genotypes: true` and aborts on DRAGEN gVCFs with
   *"couldn't find genotype likelihoods (NotFound)"* — DRAGEN records don't
   always carry PL. Our config keeps the `gatk` unifier/quality thresholds but
   sets `revise_genotypes: false`, so GLnexus trusts DRAGEN's GT calls. (This is
   what Phase 0 resolved; `gatk` fails at the genotyping step even though
   discovery succeeds.)
2. **AN denominator is right** — the verify block prints the AN distribution;
   at common sites it should peak near `2 × N_samples` (≈128 for one run).
   If AN is way below that, the gVCF reference blocks aren't being read
   (wrong config / gVCFs lack `<NON_REF>` blocks).
3. **Known common variant sanity** — pick a few well-known common SNPs and
   confirm `INHOUSE_AF` is in the right ballpark vs gnomAD EAS.

### 3. Demonstrate accumulation (the "every batch" requirement)

```bash
# run1 alone, then run1+run2 — AC should grow, AN ≈ 2×N, AF stays sane
grep -E '20251118_LH00873_0004|20251121_LH00873_0005' \
  $NGS_UI_HOME/biotools/inhouse_af/cohort_gvcfs.txt > /tmp/run12_gvcfs.txt
scripts/inhouse_af/build_inhouse_af.sh --list /tmp/run12_gvcfs.txt \
  --out /tmp/inhouse_af.run12.vcf.gz --threads 16 --mem-gbytes 96
```

### 4. Full cohort

**Restrict to primary contigs.** DRAGEN gVCFs carry ~3340 decoy/alt/HLA
contigs; variants there are useless for annotating `snv_indel.annotated.tsv`
(primary assembly only) and processing them is slow. Build a one-time BED of
chr1-22, X, Y, M from the reference `.fai` and pass it with `--bed`:

```bash
REF=/home/datalake_Intermediate/pipeline/reference/hg38/Homo_sapiens_assembly38.fasta
awk 'BEGIN{OFS="\t"} $1 ~ /^chr([0-9]+|X|Y|M)$/ {print $1,0,$2}' "$REF.fai" \
  > $NGS_UI_HOME/biotools/inhouse_af/primary_contigs.bed

scripts/inhouse_af/build_inhouse_af.sh \
  --list $NGS_UI_HOME/biotools/inhouse_af/cohort_gvcfs.txt \
  --out  $NGS_UI_HOME/biotools/inhouse_af/inhouse_af.hg38.vcf.gz \
  --bed  $NGS_UI_HOME/biotools/inhouse_af/primary_contigs.bed \
  --threads 32 --mem-gbytes 192
```

**Timing note.** GLnexus bulk-loads every gVCF into a scratch DB before
genotyping, and that load is data-volume bound (~2.5 h / 64 WGS observed →
roughly a day for 600+). `--bed` speeds the discover/genotype tail but not the
bulk load. If the full-rebuild-per-batch time becomes painful at scale, that is
the signal to switch to the DRAGEN iterative gVCF Genotyper (see below).

## Running on DGM now, air-gapped DGX later

Both machines can read the reference and the gVCFs; the only difference is
that **on the DGX the datalake paths drop the leading `/home`** (reference at
`/datalake_Intermediate/.../hg38/...`, gVCFs at `/datalake_Raw/Novaseq/...`).

Every external tool is pluggable, so the same script runs in both places:

| | DGM (now) | air-gapped DGX (later) |
|---|---|---|
| GLnexus | `GLNEXUS_BIN` (static) or `GLNEXUS_SIF` | `GLNEXUS_BIN` static binary (no network) |
| bcftools | host `bcftools` or `BCFTOOLS_SIF` | host `bcftools` binary |
| reference | `--ref /home/datalake_Intermediate/.../hg38/Homo_sapiens_assembly38.fasta` | `--ref /datalake_Intermediate/.../hg38/Homo_sapiens_assembly38.fasta` |
| cohort list | absolute `/home/...` paths | reuse same list + `--strip-path-prefix /home` |

To move to the DGX, stage offline once: the `glnexus_cli` binary, a `bcftools`
binary, and (already present) the reference + gVCFs. GLnexus is **CPU-only**
(the DGX GPUs are irrelevant) — it just needs many cores, RAM, and scratch
disk. Example DGX invocation:

```bash
GLNEXUS_BIN=/opt/glnexus/glnexus_cli BCFTOOLS_BIN=bcftools \
scripts/inhouse_af/build_inhouse_af.sh \
  --list cohort_gvcfs.txt --strip-path-prefix /home \
  --ref /datalake_Intermediate/pipeline/reference/hg38/Homo_sapiens_assembly38.fasta \
  --out /path/inhouse_af.hg38.vcf.gz --threads 32 --mem-gbytes 192
```

When a `*_SIF` is set the script apptainer-execs it (binding `APPTAINER_BIND`,
default `/home`); otherwise it runs the host binary with no container at all.

## Incremental update design (each new run)

Two options; Phase 0 uses **A** as the always-correct baseline.

**A. Full re-genotype every batch (baseline, what these scripts do).**
Append the new run's gVCFs to `cohort_gvcfs.txt`, re-run `build_inhouse_af.sh`
over the whole cohort, atomic-swap the output. Simple and always correct.
GLnexus scales to thousands of WGS, so for ~600→growing this is acceptable
(hours, run off-peak). The pitfall it avoids: you cannot just genotype the new
batch alone and add counts, because a site that is variant only in old batches
has no record in the new batch's joint VCF — its AN contribution from the new
(hom-ref) samples would be lost.

**B. DRAGEN iterative gVCF Genotyper (true incremental, later).**
DRAGEN's native iterative genotyper maintains a cohort "census" and folds in
new gVCFs without re-processing the old ones, back-filling AN at known sites
from each gVCF's reference blocks. This is the right long-term engine if
re-genotyping time becomes a problem, but it needs DRAGEN hardware/license
time and its own validation. Keep the same `select_cohort.py` front-end and
the same `INHOUSE_*` sites-VCF output contract so the UI side doesn't change.

Either way, maintain `cohort_manifest.tsv` as the audit record of which
samples (and runs) are in the current DB, and **dedup on append** so a re-run
batch is never counted twice.

## Phase 1: each new batch (after Phase 0 validates the config)

Once `--config` is confirmed, you never call `build_inhouse_af.sh` by hand
again — `update_inhouse_af.sh` is the single per-batch entry point:

```bash
COHORT=$NGS_UI_HOME/biotools/inhouse_af

# new run lands → list its gVCFs (exact bytes, tab-separated)
find /home/datalake_Raw/Novaseq/<NEW_RUN> -name '*.hard-filtered.gvcf.gz' \
  -printf '%s\t%p\n' > /tmp/new_batch_sizes.txt

scripts/inhouse_af/update_inhouse_af.sh \
  --cohort-dir "$COHORT" \
  --new-sizes  /tmp/new_batch_sizes.txt \
  --exclude-range 'VAL-:37-54' \
  --threads 32 --mem-gbytes 192
# → appends only unseen samples, re-genotypes the whole cohort, atomic-swaps
#   $COHORT/inhouse_af.hg38.vcf.gz, logs to $COHORT/updates.log
```

Re-running the same batch is safe (dedup by `sample_id` → "added=0, no
rebuild"). Use `--manifest-only` to update bookkeeping without rebuilding, and
`--strip-path-prefix /home` on the DGX.

## Incremental Phase A — per-sample ingest (validate before wiring up Phase B)

`ingest_sample.py` turns one DRAGEN sample into its per-sample contributions
(`callable.weighted.bed.gz` + `counts.tsv` + `qc.json`). Validate it on a few
samples before building the accumulator/publish (Phase B).

```bash
# pure-logic unit checks (no I/O / no bcftools)
scripts/inhouse_af/ingest_sample.py --selftest

# one sample (point at its DRAGEN other/{id}/ dir; auto-finds gVCF + ploidy CSV)
scripts/inhouse_af/ingest_sample.py \
  --sample-dir /home/datalake_Raw/Novaseq/<run>/other/25G00042 \
  --ref /home/datalake_Intermediate/pipeline/reference/hg38/Homo_sapiens_assembly38.fasta \
  --out-dir $NGS_UI_HOME/biotools/inhouse_af/per_sample
# bcftools: host binary, or export BCFTOOLS_SIF=...

cat $NGS_UI_HOME/biotools/inhouse_af/per_sample/25G00042/qc.json
```

**What to sanity-check on a handful of samples (mix of XX / XY):**
- `qc.json` `sex_class` matches the DRAGEN karyotype; `callable_bp` for chrX is
  ~2× higher in XX than XY, chrY ~0 in XX; chrM small.
- `counts.tsv`: XY samples produce `hemi` rows on chrX-nonPAR/chrY; XX produce
  none; chrM rows are `mt_hom`/`mt_het` with an `af`.
- BED weights: 2 on autosomes/PAR, 1 on male chrX-nonPAR & chrY, 1 on chrM,
  nothing on chrY for XX.

## Incremental Phase B — accumulate + publish

After ingesting samples (Phase A → `per_sample/{id}/`), fold them in and render
the sites VCF:

```bash
DB=$NGS_UI_HOME/biotools/inhouse_af

# accumulate every per_sample/ dir not yet ingested (idempotent; safe to re-run)
scripts/inhouse_af/accumulate.py --db-dir "$DB"
#   -> counts.sqlite (UPSERT) + an_track.bg.gz (event/delta rebuild) + samples table

# render the INHOUSE_* sites VCF
scripts/inhouse_af/publish_af.py --db-dir "$DB" \
  --ref /home/datalake_Intermediate/pipeline/reference/hg38/Homo_sapiens_assembly38.fasta
#   -> $DB/inhouse_af.hg38.vcf.gz  (bgzip+tabix when available)
```

- `accumulate.py` is incremental: each call only adds samples not already in the
  `samples` table, and the AN track is rebuilt from `old track + new BEDs` (no
  gVCF re-processing). Re-running a batch is a no-op.
- `publish_af.py` merge-joins in one streaming pass; AN=0 / AC=0 sites are
  dropped. Output is the same contract as the GLnexus path.

**Validation against GLnexus (ground truth):** ingest + accumulate + publish the
same 8 samples, then compare `INHOUSE_AF` to the 8-sample GLnexus VCF (scatter /
correlation). Differences are expected only at borderline calls / complex indels;
AN should agree closely since both see the same coverage.

The whole Phase A→B chain was end-to-end verified on a synthetic 3-sample cohort
(hand-checked AC/AN/AF, incremental add, and dedup).

## Phase C — validation on the DGX (incremental vs GLnexus)

The production target is the **air-gapped DGX** (1.5 TB RAM, no OOM headroom
worries; no network). Everything — the GLnexus baseline *and* the incremental
path — runs there. Two deployment details make the same scripts work unchanged:

**1. apptainer wrappers for the htslib tools.** The DGX has no `bcftools` /
`bgzip` / `tabix` on PATH, but it does have `apptainer` and a samtools/bcftools
`.sif`. Drop thin wrappers in `~/bin` so every script that calls `bcftools`
(via `BCFTOOLS_BIN=bcftools`) just works:

```bash
mkdir -p ~/bin
cat > ~/bin/bcftools <<'EOF'
#!/usr/bin/env bash
exec apptainer exec --bind /datalake_Raw,/datalake_Intermediate,/home,/raid \
  /path/to/bcftools.sif bcftools "$@"
EOF
# identical wrappers for bgzip and tabix (same sif, last arg = tool name)
chmod +x ~/bin/{bcftools,bgzip,tabix}
export PATH=$HOME/bin:$PATH
```

The wrappers must `--bind` every filesystem the tool touches (reference,
gVCFs, and the scratch/output dir) or you get "No such file" inside the
container.

**2. Paths drop the leading `/home`, and scratch must be local.**

| | value on the DGX |
|---|---|
| reference | `--ref /datalake_Intermediate/pipeline/reference/hg38/Homo_sapiens_assembly38.fasta` |
| cohort list | regenerate with DGX paths, or reuse the DGM list + `--strip-path-prefix /home` |
| GLnexus scratch (`--scratch`) | a **fast local** disk, e.g. `/raid/DGM/<user>/glnexus_scratch` — **never** NFS (`/datalake*`). GLnexus's RocksDB bulk-load thrashes random I/O and NFS makes it pathologically slow. |
| accumulate sort tmp (`--sort-tmp`) | likewise a big local disk, not `/tmp` (the event stream is tens of GB). |

A 64-sample run reproduced **bit-for-bit identical** output on the DGX vs DGM
(same 19,009,919 GLnexus variant records), confirming the port is clean.

### Validation result (64-sample run)

Pearson correlation of `INHOUSE_AF` between the **incremental** path and the
**GLnexus** ground truth, at shared normalized sites:

| stratum | sites | Pearson r |
|---|---|---|
| SNP, AN ≥ 120 | — | **0.9998** |
| indel, AN ≥ 120 | — | **0.9962** |
| all shared sites | 13.1 M | ~0.5 |

The headline ~0.5 over *all* sites is **not a bug** — it is entirely the
low-AN tail (~1.87 M sites, ~14%) where GLnexus emits a no-call for samples it
is unsure about, so its AN/AF there is based on fewer samples than ours. Our
incremental method counts every sample whose coverage is callable (DP≥10), so
it is *more complete* for in-house AF, not less accurate. On well-called sites
(AN ≥ 120, i.e. ≥60/64 samples) the two methods agree to **r ≈ 0.9998 (SNP) /
0.9962 (indel)**. See `DESIGN_incremental.md` §6a.

**Re-confirmed at full cohort (677 incremental vs 675 GLnexus):** well-called
sites (AN ≥ 1300 of 1350) agree to **r = 0.9999 (SNP) / 0.9960 (indel)**, and
the incremental DB is a **superset** of GLnexus (99.4 % of GLnexus sites present,
plus ~18 M more rare/low-quality sites). The low overall r (~0.15) is again the
rare/low-AN tail, magnified at N=677 because singletons dominate the site count.
See `DESIGN_incremental.md` §6b.

> **Indel AN double-count — fixed during Phase C.** An early run showed
> r ≈ 0.5 even at high AN. Cause: a variant's callable interval was the full
> REF span `[pos-1, pos-1+len(ref))`, which for multi-base REF / deletions
> overlapped the *same sample's* adjacent reference block and double-counted
> AN; a compounding off-by-one in the publish AN lookup landed on the
> double-counted base. Fixed by (a) recording a **single anchor base**
> `[pos-1, pos)` per variant in `ingest_sample.py`, and (b) looking AN up at
> `pos-1` in `publish_af.py`. Both fixes are in the committed scripts.

## Quarterly GLnexus audit (ground-truth baseline) — SOP

The incremental path is production; **GLnexus is only re-run periodically** as an
independent ground-truth to re-confirm the incremental AF. GLnexus is
all-or-nothing: **one malformed gVCF aborts the whole multi-hour load.** Two
guards make that survivable:

1. **`known_bad_gvcfs.txt`** — sample ids GLnexus can't load (with the reason).
   `make_glnexus_list.sh` drops them from the GLnexus list. These samples stay
   in the *incremental* cohort (its parser is tolerant) — only the GLnexus
   baseline excludes them.
2. **`preflight_glnexus.sh`** — a fast, parallel pre-scan that finds *new*
   malformed gVCFs **before** the load (not 11 h in), and appends them to
   `known_bad_gvcfs.txt`.

So you do **not** discover bad samples by losing a full run — you scan first:

```bash
DB=/raid/DGM/n102968/inhouse_af

# 1. pre-flight the whole cohort in parallel (~1–2 h; NFS-read-bound, not a full load)
scripts/inhouse_af/preflight_glnexus.sh \
  --list "$DB/cohort_gvcfs.txt" --out-prefix "$DB/glnexus_preflight" --jobs 32
#   -> $DB/glnexus_preflight.clean.txt  + appends any new bad ids to known_bad_gvcfs.txt

# 2. build the GLnexus-safe list (drops known-bad ids)
scripts/inhouse_af/make_glnexus_list.sh \
  --in "$DB/cohort_gvcfs.txt" --out "$DB/cohort_gvcfs.glnexus.txt"

# 3. joint-genotype (scratch on LOCAL /raid, keep the cohort BCF as insurance)
scripts/inhouse_af/build_inhouse_af.sh \
  --list "$DB/cohort_gvcfs.glnexus.txt" \
  --ref "$REF" --bed "$DB/primary_contigs.bed" \
  --scratch "$DB/glnexus_scratch" --keep-cohort \
  --out "$DB/inhouse_af.glnexus.vcf.gz" --threads 32 --mem-gbytes 256

# 4. compare incremental vs GLnexus (Phase C)
scripts/inhouse_af/compare_inhouse_af.py \
  --incremental "$DB/inhouse_af.hg38.vcf.gz" \
  --glnexus     "$DB/inhouse_af.glnexus.vcf.gz"
```

`preflight_glnexus.sh` is stdlib-only (fast C decompression via `bgzip -dc` piped
into `validate_gvcf_glnexus.py`, which flags the "wrong # of GT entries" defect
GLnexus rejects). `compare_inhouse_af.py` merge-joins the two sites VCFs and
prints Pearson r stratified by SNP/indel and AN floor — the low-AN tail (where
GLnexus no-calls low-confidence samples) is separated from the well-called sites
the method is judged on.

**Adding a batch does NOT need any of this.** Routine per-batch updates run the
incremental path only (`ingest_sample.py` → `accumulate.py` → `publish_af.py`),
which never re-genotypes the cohort and never trips on a malformed gVCF. The
GLnexus audit is a periodic cross-check, not part of the update loop.

## Integrating into NGS-UI (Phase 2 — DONE)

Implemented on the main line (mirrors the GeneBe/GIAB annotate pattern — the
sites VCF is joined into the TSV as columns, NOT `bcftools annotate`, because
`snv_indel.annotated.tsv` is a TSV):

1. **`config.py`** — `INHOUSE_AF_DB` (env `NGS_UI_INHOUSE_AF_DB`, default
   `$NGS_UI_HOME/biotools/inhouse_af/inhouse_af.hg38.vcf.gz`; missing = off).
2. **`scripts/annotate_inhouse_af.py`** — single streaming pass over the sorted
   sites VCF (not per-variant tabix seeks), fills `INHOUSE_AC/AN/AF` on matching
   `(chrom,pos,ref,alt)`, atomic replace, no-op when the DB is absent.
3. **`run_stopgaps.sh`** — runs `annotate_inhouse_af.py` after `giab-strata` and
   before `review-tsv`/`gene-index` (so the columns reach the UI and the gene
   index is rebuilt after the whole-TSV rewrite). `--skip-inhouse-af` /
   `--inhouse-af-db` flags.
4. **`adapters/snv_tsv.py`** — surfaces `inhouse_af` / `inhouse_ac` / `inhouse_an`.
5. **`services/inhouse_af_mito.py` + `adapters/mito_tsv.py`** — at backend
   startup, query only the indexed chrM slice and cache
   `INHOUSE_AC/AN/AF/NHOM/HET_MT`. For chrM, AC is the number of carriers, AN
   is callable mitochondrial genomes, and AF is carrier frequency; the source
   Mito TSV is not rewritten.
6. **Front-end** — SNV card shows an **`AF_nckuh`** row `AF (AC/AN)` (e.g.
   `0.045 (61/1352)`) under `AF`/`AF_eas`; `1000G EAS` moved into *More*. No
   display filter (deliberately). Mito cards show
   `AF (carrier AC/callable mito AN; hom/het carriers)`. "—" when the variant
   isn't in the DB.
7. **`scripts/backfill_inhouse_af.py`** — for existing samples / after a batch
   refresh. Pipeline `03_acmg` sources (read-only): update `INHOUSE_*` directly
   in the `08_postprocessing` overlay SQLite (one transaction, other fields
   kept) → rebuild review TSV; gene index untouched. Legacy UI copies:
   annotate in place → rebuild review TSV → **rebuild gene index** (offsets shift).
8. **`scripts/inhouse_af/deploy_inhouse_af_db.sh`** — atomic-install the sites VCF
   onto the NGS-UI host under `biotools/inhouse_af/`.

Case list / DOCX intentionally NOT wired (main-screen display only). See
`PHASE2_PLAN.md` for the decisions.

### Verify chrM records and format on the NGS-UI server

```bash
ngs_af_db="${NGS_UI_INHOUSE_AF_DB:-$HOME/NGS_UI/biotools/inhouse_af/inhouse_af.hg38.vcf.gz}"

ls -lh "$ngs_af_db" "$ngs_af_db.tbi"

# The matching contig row should have a non-zero record count.
bcftools index -s "$ngs_af_db" \
  | awk '$1 ~ /^(chrM|MT|M|chrMT)$/ {print}'

# Confirm the published field definitions.
bcftools view -h "$ngs_af_db" \
  | grep -E '^##(source|contig=<ID=(chrM|MT|M|chrMT)|INFO=<ID=INHOUSE_)'

# Inspect five real records. Expected columns:
# CHROM POS REF ALT carrier_AC callable_mito_AN carrier_AF homoplasmic heteroplasmic
bcftools query -r chrM \
  -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/INHOUSE_AC\t%INFO/INHOUSE_AN\t%INFO/INHOUSE_AF\t%INFO/INHOUSE_NHOM\t%INFO/INHOUSE_HET_MT\n' \
  "$ngs_af_db" | head -n 5

# Every current chrM row should satisfy AC = NHOM + HET_MT.
bcftools query -r chrM \
  -f '%INFO/INHOUSE_AC\t%INFO/INHOUSE_NHOM\t%INFO/INHOUSE_HET_MT\n' \
  "$ngs_af_db" \
  | awk '$1 != $2 + $3 {bad++} END {print "bad_rows=" (bad+0)}'
```

The incremental publisher writes `chrM` and defines:

- `INHOUSE_AC` (`Number=A`): samples carrying the ALT.
- `INHOUSE_AN` (`Number=1`): callable mitochondrial genomes, one per sample at
  the site.
- `INHOUSE_AF` (`Number=A`): `AC / AN`, i.e. carrier frequency.
- `INHOUSE_NHOM` (`Number=A`): homoplasmic carriers (`FORMAT/AF >= 0.95` at
  ingest).
- `INHOUSE_HET_MT` (`Number=A`): heteroplasmic carriers.

## Caveats to keep in mind

- **Disease-referral cohort, not population controls.** Use in-house AF to
  flag recurrent artifacts / locally common variants; do **not** use it as
  ACMG BA1/BS1 population-frequency evidence (that stays gnomAD).
- **N=677** → minimum resolvable allele frequency ≈ 1/1354 ≈ 0.07%; show AC/AN.
- **No relatedness/MRN dedup yet** (see Cohort note).
