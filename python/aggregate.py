import sys, os

import itertools, os
from Fasta_file import Fasta_file

seqs = open('final_sequences.fa', 'w')
maps = open('final_correspondances.psl', 'w')
quas = open('final_qualities.fa', 'w')
scos = open('final_scores.fa', 'w')

class ScoresIterator:

    def __init__(self, path):
        self.file = open(path)
        # Get the first name in the file, to start the iteration
        self.cname = self.name(self.file.next())
        self.scores = dict(self)

    def __getitem__(self, name):
        return self.scores.get(name, 'missing')

    def has_key(self, name):
        return True

    def __iter__(self):
        return self

    def name(self, path):
        name = path.strip().split('/')[-1]
        assert name.endswith('.gz')
        return os.path.basename(name)[:-3]

    def nextp(self, line):
        return line.startswith('  traces')

    def next(self):
        if self.cname is None:
            raise StopIteration
        cname, scores = self.cname, []
        for line in self.file:
            if self.nextp(line):
                self.cname = self.name(line)
                break
            scores.append(line.strip())
        else:
            # This was the last entry  in the file; should stop on the
            # next iteration.
            self.cname = None
        return cname, ' '.join(scores)

class QualsIterator(ScoresIterator):

    def __getitem__(self, name):
        return self.scores[name]

    def has_key(self, name):
        return self.scores.has_key(name)

    def name(self, title):
        return title.split()[0][1:]

    def nextp(self, line):
        return line.startswith('>')

cseqs = dict((n.split()[0][1:], s) for n, s in Fasta_file(
    'sequences.fa'))
cquas = QualsIterator('qualities.fa')
cscos = ScoresIterator('scores.txt').scores
cmaps = open('correspondances.psl')
# Skip header
for line in cmaps:
    if set(line.strip()) == set('-'):
        break
alggroups = itertools.groupby(cmaps, lambda m: m.split()[9])
for tid, group in alggroups:
    longest = max(
        group,
        key=lambda m: sum(map(int, m.split()[-3].split(',')[:-1])))
    print >> maps, longest.rstrip()
    for ff, output in ((cseqs, seqs), (cquas, quas), (cscos,scos)):
        if not ff.has_key(tid):
            raise 'Missing sequence!'
        print >> output, '>%s\n%s' % (tid, ff[tid])

for sfile in (seqs, quas, scos, maps):
    sfile.close()
