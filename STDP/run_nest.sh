#!/usr/bin/env bash

set -euo pipefail

if [[ $# -ne 1 ]]; then
	echo "Missing configuration variant. Usage: $(basename "$0") <variant>"
	exit 2
fi

variant="$1"

# Using conda/miniforge environment for NEST (needs to be set up beforehand)
#conda activate nest_env

conda run -n nest_env ./nest_stdp_lif.py "$variant"
conda run -n nest_env ./nest_stdp_classical.py "$variant"

#conda deactivate
