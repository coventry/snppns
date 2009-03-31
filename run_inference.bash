#!/usr/bin/env bash

# Fail on error
set -e
set -o pipefail

source `dirname $0`/config.bash

$PYTHON $reldir/python/inference.py "$@"
