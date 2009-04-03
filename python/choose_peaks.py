import itertools, scipy.stats, sys, os, cPickle
from itertools import izip
from Fasta_file import Fasta_file
import SNPPosterior, Peak
import BlatMap

# Get the map of these subjects to the tracefiles
subject_map = cPickle.load(open('submap.pickle'))
# Specify an ordering of the subjects
_subject_order = sorted(subject_map)
subject_order = dict((p, i) for i, p in enumerate(_subject_order))
trace_subject_map = dict(
    (tid, pid) for pid, traces in subject_map.items()
    for tid in traces)

refseqs = dict(Fasta_file(sys.argv[1]))
assert all(n.startswith('>') for n in refseqs)
refseqs = dict((n[1:].split()[0], s) for n, s in refseqs.items())

# For storing the peaks in
peaks = dict(
    (n, dict((strand, scipy.zeros((len(subject_order), len(s)), 'B'))
             for strand in '+-'))
    for n, s in refseqs.items())
# Indicate that all cells are currently empty
for tablepair in peaks.values():
    for table in tablepair.values():
        table.fill(255)

class UnindexedError(RuntimeError):
    pass

class TraceInfo:

    def __init__(self, seq, scores, adjscores, bmap):
        self.seq = seq
        self.scores = map(int, scores.split())
        self.bmap = bmap
        assert bmap.tid.endswith('.scf')
        self.tid = bmap.tid[:-4]
        if self.tid not in trace_subject_map:
            raise UnindexedError
        self.strand = bmap.strand
        bmap.make_map()
        adjscores = adjscores.split()
        # Split into groups of five.  Take a list here, as we use it
        # twice when we check the lengths.
        assert (len(adjscores) % 5) == 0
        adjscores = list(izip(*(5*[iter(adjscores)])))
        # Compute the double-peak calls first
        doublesites = set()
        self.peaks = {}
        for peak in adjscores:
            if peak[1] == peak[3]:
                # Not interested in fake double peaks
                continue
            pos = int(peak[0])
            doublesites.add(pos)
            if pos not in bmap.map:
                # Only interested in positions we can confidently
                # assign.
                continue
            peaks = sorted([peak[1:3], peak[3:]])
            peak = (''.join([p[0] for p in peaks]),
                    tuple([Peak.binlocs[int(p[1])]
                           for p in peaks]))
            self.peaks[bmap.map[pos]] = SNPPosterior.rawscores[peak]
        # Compute the single-peaks for the other sites
        self.peaks.update(dict([
            (bmap.map[p], SNPPosterior.rawscores[c,Peak.binlocs[s]])
            for p, (c, s) in enumerate(izip(self.seq, self.scores))
            if ((p in bmap.map) and (c in 'ACGT') and
                # Don't do sites where there was a double peak.
                (p not in doublesites))]))
        self.subject = trace_subject_map[self.tid]

    def fill_entry(self):
        patidx = subject_order[self.subject]
        row = peaks[self.bmap.refseq][self.strand][patidx]
        for pos, peak in self.peaks.items():
            row[pos] = SNPPosterior.choose_best(peak, row[pos])

class TruncatedPhredFilesError(RuntimeError):
    pass

class TraceIter:

    def __init__(self):
        self.streams = [open('final_correspondances.psl'),
                        Fasta_file('final_sequences.fa'),
                        Fasta_file('final_qualities.fa'),
                        Fasta_file('final_scores.fa')]
        self.iter = izip(*self.streams)

    def __iter__(self):
        return self

    def _next(self):
        try:
            return self.iter.next()
        except StopIteration:
            # One of the streams  has stopped.  Check that all streams
            # ended at the same time
            for stream in self.streams:
                try:
                    stream.next()
                except StopIteration:
                    # This stream has stopped.  Good.
                    pass
                else:
                    raise TruncatedPhredFilesError
            else:
                raise
            raise RuntimeError, 'Should never get  here'

    def next(self):
        while 1:
            try:
                (algline, (seqname, seq), (scorename, scores),
                 (adjscorename, adjscores)) = self._next()
                bmap = BlatMap.Map.parse_blat_line(algline, len(seq))
                assert bmap.tid==seqname[1:]==scorename[1:]==adjscorename[1:]
                return TraceInfo(seq, scores, adjscores, bmap)
            except UnindexedError:
                print 'missed', bmap.tid
                pass

if __name__ == '__main__':
    for tcount, trace in enumerate(TraceIter()):
        trace.fill_entry()
    cPickle.dump((_subject_order, peaks),
                 open('final_peaks.pickle', 'w'))
