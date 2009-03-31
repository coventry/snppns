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
LOCALSAGEDIR=/code/production/SnppnS/sage-3.4
LOCALSAGE=${LOCALSAGEDIR}/sage

# Location of phred parameter file.
# If not set, will default as follows:
# * the current value of the PHRED_PARAMETER_FILE environment variable
# * the old file that is bundled with phred 0.020425.c
# export PHRED_PARAMETER_FILE=~/phredpar.dat

# Path to directory containing BLAS and LAPACK.
# May not be required for some deployments of scipy.
# If using sage, just uncomment this line.
# export LD_LIBRARY_PATH=${LOCALSAGEDIR}/local/lib

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
