#!/usr/bin/env bash

# Fail on error
set -e
set -o pipefail

source `dirname $0`/config.bash

rm -rf $phredsourcedir
cp -r $originalphredsourcedir $phredsourcedir

# Check that the file to be patched matches
if [ -z "`md5sum $phredsourcedir/setQual.c | grep cc47b60a20f917471b747b10014a0d60`" ] ; then
    echo "setQual.c doesn't match!"
    /bin/false
fi

# Check that the version matches
if [ -z "`fgrep 0.020425.c $phredsourcedir/setQual.c`" ]; then
    echo "phred version doesn't match!"
    /bin/false
fi

# Do the patch
patch $phredsourcedir/setQual.c -i $reldir/diff

# Build the patched phred
cd $phredsourcedir
make clean; make
