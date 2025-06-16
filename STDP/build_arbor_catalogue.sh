#!/usr/bin/env bash

set -euo pipefail

# Set up environment for Arbor, e.g.,
# export PATH=$(readlink -f ~/local/bin):$PATH
# export PYTHONPATH=$(readlink -f ~/local/lib/python3.8/site-packages):${PYTHONPATH-}
# export LD_LIBRARY_PATH=$(readlink -f ~/local/lib):${LD_LIBRARY_PATH-}
# or,
source set_arbor_env

arbor-build-catalogue "$@"
