#!/usr/bin/env python
"""Break the tar files up into cluster jobs, arranged by patient."""

import sys, os, itertools, tempfile, shutil, cPickle
from collections import defaultdict as ddict

trace_suffix = '.scf.gz'
tracepath, mappath, numjobs = sys.argv[1:]
numjobs = int(numjobs)

# Map from trace IDs to patient IDs
patient_map = {}
for line in open(mappath):
    tid, patid = line.split()
    if patient_map.get(tid, patid) != patid:
        raise RuntimeError, (
            'Inconsistent patient map for %s: ' +
            ('%s and %s' % (patid, patient_map[tid])))
    patient_map[tid] = patid
patients_seen = set()
# Map from trace IDs to file paths
tracepath_map = {}
# Map from patients to lists of traces
trace_map = ddict(set)

# Find all the trace paths, and gather them by patient
for path in os.popen('find %s -name \*%s'%(tracepath, trace_suffix)):
    path = path.strip()
    tid = os.path.basename(path)[:-len(trace_suffix)]
    if tid in patient_map:
        if tid in tracepath_map:
            raise RuntimeError, (
                ('Trace path basenames should be unique: '
                 '%s and %s conflict.' )%
                (path, tracepath_map[tid]))
        tracepath_map[tid] = path
        pid = patient_map[tid]
        trace_map[pid].add(tid)
        patients_seen.add(pid)


patients_seen = list(patients_seen)
# Number of patients per job
numpats = (len(patients_seen)/numjobs)
# Make sure it's *at most* numjobs jobs
if (len(patients_seen) % numjobs) != 0:
    numpats += 1

# Group the patients into jobs
jobgroups = [patients_seen[i:i+numpats]
             for i in range(0, len(patients_seen)+1, numpats)]
# Could be an empty list on the end; get rid of it
jobgroups = filter(None, jobgroups)

for jobnum, patients in enumerate(jobgroups):
    # Make a symlink farm to the traces for this job
    tempdir = tempfile.mkdtemp(dir=os.path.abspath('.'))
    tracedir = os.path.join(tempdir, 'traces')
    os.mkdir(tracedir)
    submap = dict((patid, trace_map[patid]) for patid in patients)
    for tid in itertools.chain(*(submap.values())):
        tpath = os.path.abspath(tracepath_map[tid])
        dpath = os.path.join(tracedir, os.path.basename(tpath))
        os.symlink(tpath, dpath)
    # Also take a copy of this portion of the subject/trace map.
    cPickle.dump(submap, open(os.path.join(tracedir, 'submap.pickle'), 'w'))
    # Make a  tarball of it  "c" to create  the tar file, "f"  for the
    # file  system,  "h"  to   dereference  the  symlinks.   No  sense
    # compressing; the trace files are already compressed.
    tarpath = 'job.%i.tar' % jobnum
    assert not os.system('rm -f ' + tarpath)
    tarcmd = 'cd %s ; tar cfh ../%s traces' % (tempdir, tarpath)
    assert not os.system(tarcmd)
    # Tear down the symlink farm
    assert not os.system('rm -rf ' + tempdir)
    
