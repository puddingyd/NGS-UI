#!/usr/bin/env bash
# =========================================================
# deploy_inhouse_af_db.sh — 把 in-house AF sites VCF 原子安裝到部署機
# =========================================================
# 在「跑 NGS-UI 的那台機器」上執行。把增量版 publish 出來的
# inhouse_af.hg38.vcf.gz(+.tbi) 驗證後原子搬進 biotools/inhouse_af/，
# 讓 config.py 的 NGS_UI_INHOUSE_AF_DB 找得到（缺檔會靜默 disable）。
#
#   1. 驗證來源 bgzip 完整性（bgzip -t）
#   2. 確保 .tbi：來源有且比 .gz 新就沿用；否則 tabix -p vcf 重建
#   3. 探測一個已知位點（sanity）
#   4. 原子安裝到 DEST（先 .gz 再 .tbi，同檔系統 mv = rename）
#
# 來源檔在 DGX（例 /raid/DGM/n102968/inhouse_af/inhouse_af.hg38.vcf.gz）；
# 若部署機與 DGX 共用 datalake，先把 .gz(+.tbi) 複製到共享路徑，再用它當 SRC。
#
# 用法（在部署機）：
#   scripts/inhouse_af/deploy_inhouse_af_db.sh <SRC.vcf.gz> [--dest DEST]
#
#   SRC   來源 inhouse_af.hg38.vcf.gz（會一併找同名 .tbi）
#   DEST  預設 $NGS_UI_INHOUSE_AF_DB，否則
#         $NGS_UI_HOME/biotools/inhouse_af/inhouse_af.hg38.vcf.gz，
#         否則 $HOME/NGS_UI/biotools/inhouse_af/inhouse_af.hg38.vcf.gz
# =========================================================
set -euo pipefail

NGS_UI_HOME="${NGS_UI_HOME:-$HOME/NGS_UI}"
DEST="${NGS_UI_INHOUSE_AF_DB:-$NGS_UI_HOME/biotools/inhouse_af/inhouse_af.hg38.vcf.gz}"
SRC=""
PROBE_REGION="chr1:10146-10146"   # 早期常見位點，用來確認索引可查

while [ $# -gt 0 ]; do
  case "$1" in
    --dest)    DEST="$2"; shift 2;;
    --probe)   PROBE_REGION="$2"; shift 2;;
    -h|--help) sed -n '2,36p' "$0"; exit 0;;
    -*)        echo "unknown arg: $1" >&2; exit 2;;
    *)         SRC="$1"; shift;;
  esac
done

[ -n "$SRC" ] || { echo "ERROR: 需要 SRC（來源 inhouse_af.hg38.vcf.gz）" >&2; exit 2; }
[ -s "$SRC" ] || { echo "ERROR: 來源不存在／為空：$SRC" >&2; exit 1; }
for bin in bgzip tabix; do
  command -v "$bin" >/dev/null 2>&1 || { echo "ERROR: PATH 上需要 '$bin'" >&2; exit 1; }
done

echo "in-house AF DB deploy"
echo "  source: $SRC"
echo "  dest:   $DEST"

echo "[1/4] 驗證來源 bgzip 完整性…"
bgzip -t "$SRC"

DEST_DIR="$(dirname "$DEST")"
mkdir -p "$DEST_DIR"
# 在目的資料夾裡建暫存，最後 mv 才是同檔系統的原子 rename
TMP_GZ="$(mktemp "$DEST_DIR/.inhouse_af_deploy.XXXXXX").vcf.gz"
cleanup() { rm -f "$TMP_GZ" "$TMP_GZ.tbi"; }
trap cleanup EXIT

cp -f "$SRC" "$TMP_GZ"

echo "[2/4] 準備 .tbi 索引…"
if [ -f "$SRC.tbi" ] && [ "$SRC.tbi" -nt "$SRC" ]; then
  echo "      沿用來源 .tbi（比 .gz 新）"
  cp -f "$SRC.tbi" "$TMP_GZ.tbi"
else
  echo "      來源 .tbi 缺失或過舊 → tabix -p vcf 重建"
  tabix -p vcf "$TMP_GZ"
fi

echo "[3/4] sanity 探測：$PROBE_REGION"
if ! tabix "$TMP_GZ" "$PROBE_REGION" | head -1; then
  echo "      WARNING: 探測沒回傳（該位點可能本來就沒有，非致命）" >&2
fi
echo -n "      站點數（bcftools index -n，若有 bcftools）："
if command -v bcftools >/dev/null 2>&1; then bcftools index -n "$TMP_GZ" 2>/dev/null || echo "n/a"; else echo "略過（無 bcftools）"; fi

echo "[4/4] 原子安裝…"
mv -f "$TMP_GZ"     "$DEST"
mv -f "$TMP_GZ.tbi" "$DEST.tbi"
trap - EXIT

echo "done → $DEST (+ .tbi)"
ls -la "$DEST" "$DEST.tbi"
echo
echo "接著在 NGS-UI 設定 NGS_UI_INHOUSE_AF_DB=$DEST（若 DEST 已是預設路徑則免）。"
echo "注意：畫面出現 AF_nckuh 還需要 Phase 2 的 annotate/adapter/前端變更（見 PHASE2_PLAN.md）。"
