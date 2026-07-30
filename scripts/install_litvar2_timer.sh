#!/usr/bin/env bash
# Install the monthly LitVar2 updater on the production NGS-UI host.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
UNIT_DIR="$REPO_ROOT/deploy/systemd"

sudo install -m 0644 \
  "$UNIT_DIR/ngs-ui-litvar2-update.service" \
  /etc/systemd/system/ngs-ui-litvar2-update.service
sudo install -m 0644 \
  "$UNIT_DIR/ngs-ui-litvar2-update.timer" \
  /etc/systemd/system/ngs-ui-litvar2-update.timer
sudo systemctl daemon-reload
sudo systemctl enable --now ngs-ui-litvar2-update.timer
systemctl status --no-pager ngs-ui-litvar2-update.timer
