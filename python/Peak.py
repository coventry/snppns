import scipy

# Construct the score bins
bins = scipy.array([0]+list(reversed(range(60, 10, -10))))
binlocs = scipy.zeros(69)
for binidx, (start, end) in enumerate(zip(bins, list(bins[1:])+[len(binlocs)])):
    for score in range(start, end):
        binlocs[score] = binidx
binlocs[-1] = len(bins)-1

# Possible double peaks
hetpeaks = ['AC', 'AG', 'AT', 'CG', 'CT', 'GT']

class Peak:

    '''Bytewise-storage and reconstitution of peak data'''

    def __init__(self, **kwargs):
        assert len(kwargs) == 1
        if 'raw' in kwargs:
            # Raw peak observation data, build from that
            raw_score = kwargs['raw']
            if len(raw_obs) == 2:
                # It was a single peak.
                nuke, score = raw_obs
                if nuke not in 'ACGT':
                    raise ValueError
                self.double = False
                self.peaks = nuke
                self.bin = binlocs[score]
            else:
                # It was a double peak
                nuke1, score1, nuke2, score2 = raw_obs
                if nuke1 == nuke2:
                    raise ValueError
                self.double = True
                self.peaks = ''.join(sorted([nuke1, nuke2]))
                if not set(self.peaks).issubset('ACGT'):
                    raise ValueError
                self.bin = tuple(
                    binlocs[t[1]] for t in
                    sorted([[nuke1,score1],[nuke2,score2]]))
        else:
            # Bytewise data;  translate into raw  peak observation and
            # build from that.
            byte = kwargs['byte']
            if int(byte) != byte:
                raise ValueError
            byte = int(byte)
            if byte >= 216:
                # It was a single peak
                self.double = False
                byte -= 216
                if byte >= 24:
                    raise ValueError
                self.peaks = 'ACGT'[byte/6]
                self.bin = byte % 6
            else:
                # It was a double peak
                self.double = True
                bin2 = byte % 6
                byte /= 6
                bin1 = byte % 6
                byte /= 6
                self.peaks = hetpeaks[byte]
                self.bin = (bin1, bin2)

    def pack(self):
        if self.double:
            assert 0 <= min(self.bin) <= max(self.bin) < 6
            return (36*hetpeaks.index(self.peaks) + 6*self.bin[0] +
                    self.bin[1])
        else:
            assert len(self.peaks) == 1
            nukeidx = 'ACGT'.find(self.peaks)
            assert nukeidx != -1
            return 216+6*nukeidx+self.bin
