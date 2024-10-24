#!/usr/bin/env bash

set -euo pipefail

# If necessary, set up environment for Arbor, e.g., via
source set_arbor_env

./arbor_homeostasis.py "$@"
