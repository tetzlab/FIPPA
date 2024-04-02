#!/usr/bin/env bash

set -euo pipefail

if [[ $# -ne 1 ]]; then
	echo "Missing configuration variant. Usage: $(basename "$0") <variant>"
	exit 2
fi

variant="$1"

# Set up environment for Brian2, e.g.,
# source ~/brian2/venv/bin/activate

./brian2_stdp_lif.py "$variant"
./brian2_stdp_classical.py "$variant"