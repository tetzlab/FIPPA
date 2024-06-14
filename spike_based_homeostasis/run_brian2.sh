#!/usr/bin/env bash

set -euo pipefail

if [[ $# -ne 1 ]]; then
	echo "Missing configuration file. Usage: $(basename "$0") <config_file>"
	exit 2
fi

config="$1"

# Set up environment for Brian2, e.g.,
# source ~/brian2/venv/bin/activate

./brian2_homeostasis.py "$config"