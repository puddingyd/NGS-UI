#!/usr/bin/env bash
# Install the weekly ClinVar updater on the production NGS-UI host.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
UNIT_DIR="$REPO_ROOT/deploy/systemd"

sudo install -m 0644 \
  "$UNIT_DIR/ngs-ui-clinvar-update.service" \
  /etc/systemd/system/ngs-ui-clinvar-update.service
sudo install -m 0644 \
  "$UNIT_DIR/ngs-ui-clinvar-update.timer" \
  /etc/systemd/system/ngs-ui-clinvar-update.timer
sudo systemctl daemon-reload
sudo systemctl enable --now ngs-ui-clinvar-update.timer
sudo systemctl start --no-block ngs-ui-clinvar-update.service
systemctl status --no-pager ngs-ui-clinvar-update.timer
