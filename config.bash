set -e
set -o pipefail

# Directory containing the release
export reldir=$(readlink -f "`dirname $0`")

# Common local scratch space on the cluster nodes
scratchdir=/state/partition1

# Directory containing the build of the patched phred
originalphredsourcedir=~/src/phred/source

# Path to the blat binary
BLAT=~/bin/i386/blat

# Path to local sage executable (if using sage)
LOCALSAGE=~/src/sage/sage

# Path to the python executable
PYTHON=${reldir}/sage.bash # /usr/local/bin/python

export PYTHONPATH=${PYTHONPATH}:${reldir}/python
phredsourcedir=${reldir}/phred

function make_temp_scratch_dir() {
    # Working directory unique to current process
    tmpdir=${scratchdir}/SnppnStmp.$$
    mkdir -p $tmpdir
    echo $tmpdir
}
