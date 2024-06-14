#!/usr/bin/env bash

set -euo pipefail

if [[ $# -ne 1 ]]; then
	echo "Missing configuration file. Usage: $(basename "$0") <config_file>"
	#exit 2
fi

config="$1"

# Set up environment for Arbor, e.g.,
# export PATH=$(readlink -f ~/local/bin):$PATH
# export PYTHONPATH=$(readlink -f ~/local/lib/python3.8/site-packages):${PYTHONPATH-}
# export LD_LIBRARY_PATH=$(readlink -f ~/local/lib):${LD_LIBRARY_PATH-}
# or,
# source set_arbor_env

./arbor_homeostasis.py "$config"
