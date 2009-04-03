set -e
set -o pipefail

source `dirname $0`/config.bash

cd $reldir
dirpath="SnppnS-`cat VERSION`"
mkdir $dirpath
for p in `cat MANIFEST` ; do
    d=`dirname $p`
    mkdir -p $dirpath/$d
    cp $p $dirpath/$d/
done

tar zcf ${dirpath}.tgz $dirpath
rm -rf $dirpath
scp $dirpath.tgz cvp:www/
