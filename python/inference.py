import sys, cPickle, os, scipy, time
from collections import defaultdict as ddict
import MCMC, Dirichlet, SNPPrior

# Flat prior
hyperprior = Dirichlet.Dirichlet(10*[1])

dict_revcalls = dict(zip(SNPPrior.pcalls, SNPPrior.revcalls))
callidxs = dict((g, i) for i, g in enumerate(SNPPrior.pcalls))

class Sample:

    def __init__(self, gdist, peak_counts, site_dists):
        self.gdist, self.peak_counts = gdist, peak_counts
        self.site_dists = site_dists

    def posterior(self):
        return dict(
            (strand, SNPPrior.PeakModel(
            SNPPrior.call_prior.posterior_mean(counts)))
            for strand, counts in self.peak_counts.items())

    @classmethod
    def zero_counts(cls):
        return dict((strand, SNPPrior.counts()) for strand in '+-')

    def new_sample(self, bytes, independent_sites):
        # Posterior  sample given previous  peak/call counts  (or zero
        # counts, if it's the initial sample.)
        peak_posterior = self.posterior()
        # New peak/call counts
        counts = self.zero_counts()
        # New call counts
        call_counts = scipy.zeros(10)
        site_dists = []
        total_observations = dict((s, 0) for s in '+-')
        for byteidx, cbytes in enumerate(bytes):
            # Get the posterior distributions given each observed peak
            # at this site.
            peaks, dists = zip(*[peak_posterior[s].call_dist(b,s)
                                 for b, s in cbytes])
            # Sample a posterior call
            site_dist = self.gdist*scipy.product(dists, axis=0)
            assert sum(site_dist)
            site_dist = site_dist/sum(site_dist)
            site_dists.append(site_dist)
            call = SNPPrior.sample_call(site_dist)
            strandcall = {'+': call, '-': dict_revcalls[call]}
            # Count the peak/call pairs.
            for peak, (byte, strand) in zip(peaks, cbytes):
                # Undo the reverse complement, if this observation
                # was on the reverse strand.
                ccall = strandcall[strand]
                ccounts = counts[strand]
                # Flat peak/call count
                ccounts[ccall][SNPPrior.peakidxs[peak.peaks]] += 1
                # Call->(peak,score) count
                ccounts['%s-%s' % (ccall, peak.peaks)][peak.bin]+=1
                total_observations[strand] += 1
            # Count of imputed calls
            call_counts[callidxs[call]] += 1
        # Blur the conditional distribution, to speed convergence.
        for strand, ccounts in counts.items():
            for val in ccounts.values():
                val /= total_observations[strand]
                val *= min(total_observations[strand], 100)
        gdist = hyperprior.posterior(call_counts).sample()

        return Sample(gdist, counts, site_dists)

class InitialSample(Sample):

    initial_model = SNPPrior.PeakModel(SNPPrior.call_prior.mean())

    def __init__(self):
        self.gdist = hyperprior.sample()

    def posterior(self):
        # Initial samples  start with  the means of  the Peak/Genotype
        # prior  distributions, to give  them a  hint about  the right
        # direction in which to move.
        return dict((strand, self.initial_model)
                    for strand in '+-')


class PriorSampler:

    def __init__(self, bytes, independent_sites):
        self.bytes, self.independent_sites = bytes, independent_sites
        self.dists = []

    def sample(self, sample=None):
        if sample is None:
            sample = InitialSample()
        new_sample = sample.new_sample(self.bytes,
                                       self.independent_sites)
        self.dists.append(scipy.array(new_sample.site_dists))
        return new_sample

    def compute_estimand(self, sample):
        return sample.gdist

    def average_site_dists(self):
        return sum(self.dists)/len(self.dists)

def main(tables, startidx=None, endidx=None,
         savepath=None, time_allowed=3600):
    sys.path.append('/home/coventry/repository')
    results = []
    length = tables['-'].shape[1]
    for cidx in range(length)[startidx:endidx]:
        print 'column', cidx
        distlen = tables['-'].shape[0]
        lowest_byte = min(scipy.array(
            [t[:,cidx] for t in tables.values()]).flatten())
        if lowest_byte == 255:
            # There were NO calls for  any patient in this column.  No
            # work to do.
            results.append((None, distlen*[None]))
            continue
        # Get rid of any sites where there was no call made.
        byteidxs, bytes = zip(*[
            (i, [(byte,strand) for strand,byte in zip('+-',bytes)
                 if byte < 255])
            for i, bytes in enumerate(zip(
            *[t[:,cidx] for t in tables.values()]))
            if min(bytes) < 255])
        if not bytes:
            # There were NO calls for  any patient in this column.  No
            # work to do.
            results.append((None, None, distlen*[None]))
            continue
        independent_sites = ddict(list)
        for strand, table in tables.items():
            nbds = set()
            for bi, pi in enumerate(byteidxs):
                nbd = tuple(table[pi,cidx-1:cidx+2] )
                if nbd not in nbds:
                    independent_sites[strand].append(bi)
                else:
                    nbds.add(nbd)
        sampler = PriorSampler(bytes, independent_sites)
        series = MCMC.MCMC(sampler)
        # Give up seeking convergence after time_allowed seconds
        stop_time = time.time()+time_allowed
        mean = series.estimate_mean(stop_time=stop_time)
        fulldists = [None]*distlen
        for byteidx, dist in zip(byteidxs,
                                 sampler.average_site_dists()):
            fulldists[byteidx] = dist
        results.append((mean, fulldists))
    cPickle.dump(results, open(savepath, 'w'))

if __name__ == '__main__':
    colpath, datadir = sys.argv[1:]
    tables = cPickle.load(open(colpath))
    savepath = os.path.join(
        datadir, os.path.basename(colpath) + '.results')
    main(tables, savepath=savepath)

