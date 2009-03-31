#!/usr/bin/env bash

# Fail on error
set -e
set -o pipefail

source `dirname $0`/config.bash

# Unpack and build Pyrex
tar zxf ${reldir}/Pyrex*.tar.gz
export firstpath=`tar ztf ${reldir}/Pyrex*.tar.gz | head -1`
cd `dirname $firstpath`
if [ "`basename $PYTHON`" == sage.bash ] ; then
    PYTHON=$LOCALSAGE -python
fi
$PYTHON setup.py install --prefix=$reldir

cd $reldir/python
cp -r ${reldir}/lib/python2.5/site-packages/Pyrex $reldir
PYTHONPATH=$reldir ${reldir}/bin/pyrexc _Dirichlet.pyx
rm -rf build
$PYTHON setup.py build_ext
# Copy the shared object file into the python/ directory
cp `find build/ -name _Dirichlet.so` .
