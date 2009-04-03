#!/usr/bin/env bash

# This takes a copy of sage to the local scratch space and runs it, to
# minimize load on the NFS.

# Fail on error 
set -e
set -o pipefail

source "`dirname $0`/config.bash"

sagedir=${scratchdir}/sage

if ! [ -e $sagedir ]; then 
    # Copy sage from NFS and unpack it.
    # ...put it in a temp dir, to minimize risk of race conditions
    tmpdir="`make_temp_scratch_dir`" 
    initial_dir=$PWD
    cd "$tmpdir"
    tar zxf ${reldir}/sage.tgz
    # Run sage so that it can adjust to its new location.  Have it
    # wait a random duration, to further reduce the risk of race
    # conditions.
    sage/sage -python -c 'import time, random; time.sleep(random.uniform(0, 5))'
    # Check again for a copy before linking to canonical location.
    # This is a very brief operation, so it should minimize race
    # conditions.
    if ! [ -e $sagedir ] ; then
	ln -sf ${tmpdir}/sage $scratchdir
    fi
    cd $initial_dir
fi

# Run sage 
${sagedir}/sage -python "$@"
