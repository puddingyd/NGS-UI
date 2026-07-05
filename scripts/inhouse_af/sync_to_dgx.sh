#!/usr/bin/env bash
# =========================================================
# sync_to_dgx.sh — 在 DGM 上執行：pull 分支 → 把 inhouse_af scripts 送到 DGX
# =========================================================
# 用途：DGX 沒有網路、不能 git pull，所以在有網路的 DGM 上更新 repo，
# 再用 rsync/scp 把 scripts/inhouse_af/ 的「程式碼」複製到 DGX。
# 只送程式碼（*.py *.sh *.yml *.md *.txt），不送 __pycache__，也不碰任何
# 病患資料 / DB（那些本來就不在這個資料夾）。
#
# 用法（在 DGM 上）：
#   scripts/inhouse_af/sync_to_dgx.sh                 # 用下面的預設
#   scripts/inhouse_af/sync_to_dgx.sh --no-pull       # 跳過 git，只送檔
#   scripts/inhouse_af/sync_to_dgx.sh \
#     --repo ~/NGS_UI/NGS-UI --branch claude/pensive-johnson-d68tmi \
#     --dgx n102968@dgx2 --dest '~/dgx_stage/inhouse_af'
#
# 可用環境變數覆寫：REPO / BRANCH / DGX / DEST
# =========================================================
set -euo pipefail

# ---- 預設值（依目前環境；用旗標或環境變數覆寫）--------------------------
REPO="${REPO:-$HOME/NGS_UI/NGS-UI}"
BRANCH="${BRANCH:-claude/pensive-johnson-d68tmi}"
DGX="${DGX:-n102968@dgx2}"
DEST="${DEST:-~/dgx_stage/inhouse_af}"     # DGX 上的目的資料夾（注意是 flat，不含 scripts/）
VIA="${VIA:-}"                             # 若設了：走「共享 datalake 中繼」而非 ssh/scp
DGX_VIEW="${DGX_VIEW:-}"                   # --via 在 DGX 上看到的路徑（預設把開頭 /home 去掉）
DO_PULL=1

while [ $# -gt 0 ]; do
  case "$1" in
    --repo)     REPO="$2"; shift 2;;
    --branch)   BRANCH="$2"; shift 2;;
    --dgx)      DGX="$2"; shift 2;;
    --dest)     DEST="$2"; shift 2;;
    --via)      VIA="$2"; shift 2;;
    --dgx-view) DGX_VIEW="$2"; shift 2;;
    --no-pull)  DO_PULL=0; shift;;
    -h|--help)  sed -n '2,28p' "$0"; exit 0;;
    *) echo "unknown arg: $1" >&2; exit 2;;
  esac
done

SRC="$REPO/scripts/inhouse_af"
[ -d "$SRC" ] || { echo "ERROR: 找不到 $SRC （--repo 指到 NGS-UI checkout 根目錄）" >&2; exit 2; }

echo "== sync_to_dgx =="
echo "  repo   : $REPO"
echo "  branch : $BRANCH"
echo "  src    : $SRC"
echo "  dgx    : $DGX"
echo "  dest   : $DEST"
echo

# ---- 1. 更新 repo 到指定分支（fast-forward；不會動到你未提交的東西）-----
if [ "$DO_PULL" -eq 1 ]; then
  echo "[1/3] git fetch + checkout + pull --ff-only ..."
  git -C "$REPO" fetch origin "$BRANCH"
  git -C "$REPO" checkout "$BRANCH"
  # 只允許快轉；若本地有分歧會停下來讓你自己處理，不會硬合併
  git -C "$REPO" pull --ff-only origin "$BRANCH"
  echo "      HEAD -> $(git -C "$REPO" rev-parse --short HEAD)  $(git -C "$REPO" log -1 --pretty=%s)"
else
  echo "[1/3] 略過 git pull（--no-pull）"
fi
echo

# ---- 模式 A：共享 datalake 中繼（--via），不需 ssh/主機名 ---------------
if [ -n "$VIA" ]; then
  echo "[2/2] 複製到共享中繼資料夾（DGX 也掛載得到）..."
  mkdir -p "$VIA"
  ( cd "$SRC" && cp -f -- *.py *.sh *.yml *.yaml *.md *.txt "$VIA/" 2>/dev/null || true )
  # DGX 端看到的路徑：預設把開頭的 /home 去掉（DGM /home/datalake_* ↔ DGX /datalake_*）
  view="$DGX_VIEW"
  [ -n "$view" ] || view="${VIA#/home}"
  echo
  echo "== 已複製到中繼：$VIA =="
  echo "現在到 DGX 上執行這一行，把檔案放進 staging："
  echo "  mkdir -p ~/dgx_stage/inhouse_af && cp -f $view/*.py $view/*.sh $view/*.yml $view/*.txt ~/dgx_stage/inhouse_af/"
  exit 0
fi

# ---- 模式 B：直接 ssh/scp 到 DGX ---------------------------------------
# ---- 2. 確保 DGX 目的資料夾存在 ----------------------------------------
echo "[2/3] 在 DGX 建立目的資料夾（若不存在）..."
ssh "$DGX" "mkdir -p $DEST"
echo

# ---- 3. 複製程式碼（rsync 優先，否則 scp 指定副檔名）--------------------
echo "[3/3] 複製 scripts/inhouse_af/ 程式碼到 DGX ..."
if command -v rsync >/dev/null 2>&1; then
  # -a 保留權限/時間；排除 __pycache__ 與任何暫存 / 資料檔
  rsync -av \
    --exclude '__pycache__' --exclude '*.pyc' \
    --include '*/' \
    --include '*.py' --include '*.sh' --include '*.yml' --include '*.yaml' \
    --include '*.md' --include '*.txt' \
    --exclude '*' \
    "$SRC/" "$DGX:$DEST/"
else
  echo "      （沒有 rsync，改用 scp 指定副檔名）"
  # shellcheck disable=SC2046
  scp $(cd "$SRC" && ls *.py *.sh *.yml *.yaml *.md *.txt 2>/dev/null | sed "s#^#$SRC/#") \
    "$DGX:$DEST/"
fi

echo
echo "== 完成 =="
echo "DGX 上驗證："
echo "  ssh $DGX 'ls -la $DEST | grep -E \"compare_inhouse_af|validate_gvcf|preflight_glnexus|make_glnexus_list|known_bad\"'"
