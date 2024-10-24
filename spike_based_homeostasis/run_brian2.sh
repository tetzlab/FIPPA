#!/usr/bin/env bash

set -euo pipefail

# If necessary, set up environment for Brian 2, e.g., via
# source ~/brian2/venv/bin/activate

./brian2_homeostasis.py "$@"