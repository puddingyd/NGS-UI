#!/usr/bin/env bash
# Reorganize the on-server layout from
#     /home/n102968/{NGS-UI, vcf, biotools}
# to
#     /home/n102968/NGS_UI/{NGS-UI, vcf, biotools, tertiary_output, data, _index.json}
#
# Run AFTER stopping the systemd unit (sudo systemctl stop ngs-ui).
# Re-run is safe: each step skips when its target already exists.
#
#   bash scripts/migrate_layout.sh
#
# Override defaults with env vars:
#   OLD_HOME=/home/foo  NEW_HOME=/home/foo/NGS_UI  bash scripts/migrate_layout.sh
set -euo pipefail

OLD_HOME="${OLD_HOME:-/home/n102968}"
NEW_HOME="${NEW_HOME:-/home/n102968/NGS_UI}"

OLD_REPO="$OLD_HOME/NGS-UI"
NEW_REPO="$NEW_HOME/NGS-UI"

echo "==> Migration plan"
echo "    OLD_HOME  : $OLD_HOME"
echo "    NEW_HOME  : $NEW_HOME"
echo "    repo      : $OLD_REPO  ->  $NEW_REPO"

if [[ ! -d "$OLD_REPO" ]]; then
  echo "!! $OLD_REPO not found; aborting."
  exit 1
fi

# Refuse to run if the server is still up — files would move under it.
if pgrep -af "uvicorn.*app.main" >/dev/null; then
  echo "!! uvicorn is still running. Stop it first: sudo systemctl stop ngs-ui"
  exit 1
fi

mkdir -p "$NEW_HOME"

move_into_new_home() {
  local src="$1" dst_name="$2"
  local dst="$NEW_HOME/$dst_name"
  if [[ -e "$dst" ]]; then
    echo "skip $dst_name: $dst already exists"
    return
  fi
  if [[ ! -e "$src" ]]; then
    echo "skip $dst_name: source $src missing"
    return
  fi
  echo "mv $src -> $dst"
  mv "$src" "$dst"
}

# 1. Move the repo first; afterwards $NEW_REPO is the canonical checkout.
move_into_new_home "$OLD_REPO" "NGS-UI"

# 2. Move data dirs that previously lived next to the repo.
move_into_new_home "$OLD_HOME/vcf"      "vcf"
move_into_new_home "$OLD_HOME/biotools" "biotools"

# 3. Lift patient data and runtime state out of the repo.
move_into_new_home "$NEW_REPO/tertiary_output" "tertiary_output"
move_into_new_home "$NEW_REPO/data"            "data"

# 4. Lift _index.json out of tertiary_output (one level up).
if [[ -f "$NEW_HOME/tertiary_output/_index.json" && ! -f "$NEW_HOME/_index.json" ]]; then
  echo "mv tertiary_output/_index.json -> _index.json"
  mv "$NEW_HOME/tertiary_output/_index.json" "$NEW_HOME/_index.json"
fi

# 5. Rewrite vcf_path in every sample_metadata.json so existing samples
#    keep pointing at a real file after the vcf/ move.
if [[ -d "$NEW_HOME/tertiary_output" ]]; then
  python3 "$NEW_REPO/scripts/rewrite_vcf_paths.py" \
    --root "$NEW_HOME/tertiary_output" \
    --old  "$OLD_HOME/vcf" \
    --new  "$NEW_HOME/vcf"
fi

# 6. Drop a systemd unit that points at the new layout. We write to
#    /etc/systemd/system/ngs-ui.service via sudo only when the user is
#    root or has passwordless sudo; otherwise we just print the unit
#    contents and let the operator install it manually.
UNIT_PATH="/etc/systemd/system/ngs-ui.service"
read -r -d '' UNIT_CONTENT <<EOF || true
[Unit]
Description=NGS-UI (FastAPI + uvicorn)
After=network.target redis.service

[Service]
Type=simple
User=n102968
WorkingDirectory=$NEW_REPO
Environment=PYTHONPATH=$NEW_REPO/backend
Environment=NGS_UI_HOME=$NEW_HOME
ExecStart=/usr/bin/env python3 -m uvicorn app.main:app --host 0.0.0.0 --port 8765
Restart=on-failure

[Install]
WantedBy=multi-user.target
EOF

echo
echo "==> systemd unit"
if [[ $EUID -eq 0 ]]; then
  printf '%s\n' "$UNIT_CONTENT" >"$UNIT_PATH"
  echo "wrote $UNIT_PATH"
  systemctl daemon-reload
  echo "ran systemctl daemon-reload"
else
  echo "(not root) install this unit manually:"
  echo "----- $UNIT_PATH -----"
  printf '%s\n' "$UNIT_CONTENT"
  echo "----------------------"
  echo "Then: sudo systemctl daemon-reload && sudo systemctl restart ngs-ui"
fi

echo
echo "Migration done. Verify:"
echo "  ls $NEW_HOME"
echo "  ls $NEW_REPO"
echo "  NGS_UI_HOME=$NEW_HOME PYTHONPATH=$NEW_REPO/backend python3 -c 'from app.config import NGS_UI_HOME, TERTIARY_OUTPUT_ROOT, EXOMISER_HOME, LIRICAL_HOME; print(NGS_UI_HOME); print(TERTIARY_OUTPUT_ROOT); print(EXOMISER_HOME); print(LIRICAL_HOME)'"
