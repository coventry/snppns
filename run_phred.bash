#!/usr/bin/env bash

# This runs phred on a tarball of trace files, finds the longest
# alignment of the phred basecalls to the reference sequence (using
# the call to aggregate.py), and determines the most informative peak
# observation in each patient at each site for both strands.

# Fail on error
set -e
set -o pipefail

source `dirname $0`/config.bash

# Canonicalize the command-line argument paths
tarpath=$(readlink -f $1)
refpath=$(readlink -f $2)
datpath=$(readlink -f $3)

# Use default PHRED_PARAMETER_FILE if it is missing
if [ -z "$PHRED_PARAMETER_FILE" ]; then
    export PHRED_PARAMETER_FILE=./phredpar.dat
fi

# Working directory unique to current process
cd "`make_temp_scratch_dir`"

# Copy the phred and tar files to local scratch
cp ${phredsourcedir}/phredpar.dat .
cp ${phredsourcedir}/phred .
cp $tarpath .

# Unpack the tar file
basetar="`basename $tarpath`"
tar xf $basetar
rm $basetar

# Pull out the subject/trace map generated by make_phred_jobs.py.
# Used by choose_peaks.py.
mv traces/submap.pickle .

# Run phred
./phred -id traces -sa sequences.fa -qa qualities.fa > scores.txt

# Get the alignments of the phred calls to the refseqs
$BLAT $refpath sequences.fa correspondances.psl

# Choose the longest of the alignments, for each trace blat aligned
$PYTHON $reldir/python/aggregate.py

# Choose the most informative tracefile peaks
$PYTHON $reldir/python/choose_peaks.py $refpath

# Save the results
cp final_peaks.pickle  ${datpath}/${basetar}.results

