#!/usr/bin/env bash
# Install the daily tertiary-analysis VCF index updater on the production host.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
UNIT_DIR="$REPO_ROOT/deploy/systemd"

sudo install -m 0644 \
  "$UNIT_DIR/ngs-ui-tertiary-index-update.service" \
  /etc/systemd/system/ngs-ui-tertiary-index-update.service
sudo install -m 0644 \
  "$UNIT_DIR/ngs-ui-tertiary-index-update.timer" \
  /etc/systemd/system/ngs-ui-tertiary-index-update.timer
sudo systemctl daemon-reload
sudo systemctl enable --now ngs-ui-tertiary-index-update.timer
systemctl status --no-pager ngs-ui-tertiary-index-update.timer
