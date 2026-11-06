#!/usr/bin/env bash
# Run by the operator on the server after updating/restarting ngs-ui.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
UNIT_DIR="$SCRIPT_DIR/../deploy/systemd"
for unit in ngs-ui-file-index-update.service ngs-ui-file-index-update.timer; do
  sudo install -m 0644 "$UNIT_DIR/$unit" "/etc/systemd/system/$unit"
done
sudo systemctl daemon-reload
# A verified complete generation must exist before disabling either old timer.
sudo systemctl start ngs-ui-file-index-update.service
sudo systemctl enable --now ngs-ui-file-index-update.timer
for old in ngs-ui-secondary-index-update.timer ngs-ui-tertiary-index-update.timer; do
  if systemctl cat "$old" >/dev/null 2>&1; then
    sudo systemctl disable --now "$old"
  fi
done
systemctl status --no-pager ngs-ui-file-index-update.timer
