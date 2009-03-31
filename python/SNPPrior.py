import scipy, scipy.stats, cPickle, os
import Dirichlet
from sequences import reverse_complement
from datadir import datadir
from SNPPosterior import (pcalls, ppeaks, peakidxs, revcalls, Peak,
                          bytescores)

ppath = os.path.join(datadir, 'best_mixture_priors.pickle')
priors = cPickle.load(open(ppath))

revorder = [pcalls.index(''.join(sorted(map(reverse_complement, g))))
            for g in pcalls]

class CallDist:

    def __init__(self, scoredists, typedists):
        self.scoredists, self.typedists = scoredists, typedists

class CallPrior(CallDist):

    def mean(self):
        return CallDist(*[
            dict((p, d.mean()) for p, d in dists.items())
            for dists in [self.scoredists, self.typedists]])

    def sample(self):
        return CallDist(*[
            dict((p, d.sample()) for p, d in dists.items())
            for dists in [self.scoredists, self.typedists]])

    def posterior_sample(self, counts):
        return CallDist(*[
            dict((p, d.posterior(counts[p].flatten()).sample())
                 for p, d in dists.items())
            for dists in [self.scoredists, self.typedists]])

    def posterior_mean(self, counts):
        return CallDist(*[
            dict((p, d.posterior(counts[p].flatten()).mean())
                 for p, d in dists.items())
            for dists in [self.scoredists, self.typedists]])

scoredists, typedists = [dict(
    (p.replace(prefix, ''), Dirichlet.DirichletMixture(*reversed(d)))
    for p, d in priors.items() if p.startswith(prefix))
                         for prefix in 'score- type-'.split()]

def counts():
    scorenames = ['%s-%s' % (call, peak)
                  for call in pcalls for peak in ppeaks]
    return dict(
        (p, (scipy.zeros(10) if '-' not in p else
             scipy.zeros((6,)*len(p.split('-')[1]))))
        for p in (pcalls + scorenames))

call_prior = CallPrior(scoredists, typedists)

class PeakModel:

    def __init__(self, dist):
        self.posteriors = {}
        for peakidx, peak in enumerate(ppeaks):
            dim = len(peak)*(6,)
            self.posteriors[peak] = scipy.zeros(
                (10,)+dim, 'f')
            agnostic = scipy.ones(dim, 'f')/scipy.product(dim)
            for gi, g in enumerate(pcalls):
                cslice = (gi,)+len(dim)*(slice(None),)
                self.posteriors[peak][cslice] = (
                    dist.typedists[g][peakidx]*
                    dist.scoredists.get('%s-%s' % (g, peak),
                                        agnostic).reshape(dim))
            self.posteriors[peak] /= sum(self.posteriors[peak])
        self.memo = {}

    def call_dist(self, byte, strand):
        peak = bytescores[byte]
        if (byte, strand) not in self.memo:
            posteriors = self.posteriors[peak.peaks][
                (slice(None),) + (peak.bin if type(peak.bin)==tuple
                                  else (peak.bin,))]
            if strand == '-':
                # Take the reverse-complement of the distribution.
                posteriors = scipy.choose(revorder, posteriors)
            self.memo[byte, strand] = posteriors
        return peak, self.memo[byte, strand]

def sample_call(dist):
    return pcalls[Dirichlet.sample(dist)]
