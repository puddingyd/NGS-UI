#!/usr/bin/env bash
# Compatibility installer: migrate to the single shared daily timer.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exec bash "$SCRIPT_DIR/install_file_index_timer.sh"
