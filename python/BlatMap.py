class IndexID:

    def __getitem__(self, x):
        if type(x) == slice:
            assert x.step is None
            return xrange(x.start, x.stop)
        return x

class Map:

    def __init__(self, tid, refseq, blockcount, strand, lengths,
                 qstarts, tstarts, transform, seqlen):
        self.tid, self.refseq, self.strand = tid, refseq, strand
        for elt in lengths, qstarts, tstarts:
            assert elt.endswith(',') and (len(elt.split(',')) ==
                                          blockcount+1)
        self.lengths = map(int, lengths.split(',')[:-1])
        self.qstarts = map(int, qstarts.split(',')[:-1])
        self.tstarts = map(int, tstarts.split(',')[:-1])
        self.bdy_coords = zip(self.lengths, self.qstarts, self.tstarts)
        self.transform = transform
        self.seqlen = seqlen

    def length(self):
        return sum(self.lengths)

    def make_map(self):
        if hasattr(self, 'map'):
            return self.map
        self.map = {}
        seqlen = self.seqlen
        for length, qstart, tstart in self.bdy_coords:
            if self.strand == '-':
                adjqstart = seqlen - qstart - 1
                qend = adjqstart - length
                direction = -1
            else:
                adjqstart = qstart
                qend = qstart + length
                direction = 1
            self.map.update(dict(zip(
                range(adjqstart, qend, direction),
                self.transform[tstart:tstart+length])))

    @classmethod
    def parse_blat_line(cls, line, seqlen, gencoordmap=IndexID()):
        line = line.split()
        strand, tid, refseq, blockcount,lengths,qstarts,tstarts = (
            line[8],line[9],line[13],int(line[17]))+tuple(line[18:21])
        return cls(tid, refseq, blockcount, strand, lengths, qstarts,
                   tstarts, gencoordmap, seqlen)

