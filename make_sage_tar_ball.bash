#!/usr/bin/env bash

# This makes a tarball of sage in a known location, for use by
# sage.bash

# Fail on error
set -e
set -o pipefail

source `dirname $0`/config.bash

cd $reldir
ln -sf "`dirname $LOCALSAGE`" sage
tar zcfh sage.tgz sage

