#!/usr/bin/env bash

# Fail on error
set -e
set -o pipefail

source `dirname $0`/config.bash

$reldir/build_dirichlet.bash
$reldir/build_phred.bash
