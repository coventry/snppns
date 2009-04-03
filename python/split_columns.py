import sys, cPickle, scipy

# Get the number of columns per inference job
numcols = int(sys.argv[1])

# Load in the results of the phred runs
results = [cPickle.load(open(rpath)) for rpath in sys.argv[2:]]

# Get the lengths of the reference sequences
lengths = dict((n, len(peaks['+'][0]))
               for n, peaks in results[0][1].items())

# Get a list of the subjects
subjects = set()
for csubs, peaks in results:
    csubs = set(csubs)
    if csubs.intersection(subjects):
        raise RuntimeError, 'Repeated subject!'
    subjects.update(csubs)
    assert peaks.keys() == lengths.keys()
numsubs = len(subjects)
ordered_subjects = sorted(subjects)
subject_order = dict((s, i) for i, s in enumerate(ordered_subjects))

# Make peak arrays for gathering the results together.
peaks = dict(
    (n, dict((strand, scipy.zeros([numsubs, l],'B'))
             for strand in '+-'))
    for n, l in lengths.items())

# Put all of the peaks into these arrays
for csubs, cpeaks in results:
    for subidx, csub in enumerate(csubs):
        fidx = subject_order[csub]
        for refseq in lengths:
            for strand in '+-':
                subpeaks = cpeaks[refseq][strand][subidx]
                peaks[refseq][strand][fidx] = subpeaks

# Break the arrays up into columns for group inference
for refseq, length in lengths.items():
    tables = peaks[refseq]
    for startidx in range(0, length, numcols):
        endidx = startidx+numcols
        portion = dict((strand, table[:,startidx:endidx])
                       for strand, table in tables.items())
        savepath = 'cols.%s.%i-%i.pickle' % (refseq, startidx, endidx)
        cPickle.dump(portion, open(savepath, 'w'))

# Dump the subject list
out = open('subject_list.txt', 'w')
for subname in ordered_subjects:
    print >> out, subname
out.close()

