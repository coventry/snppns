import random, scipy.stats, scipy.integrate
import _Dirichlet
sample = _Dirichlet.sample

class Dirichlet(_Dirichlet.Dirichlet):

    def sample(self, scale=1):
        csample = [random.gammavariate(a, 1) for a in self.alphas]
        for aidx, (sampprob, alpha) in enumerate(zip(csample, self.alphas)):
            while sampprob == 0:
                sampprob = random.gammavariate(alpha, 1)
            csample[aidx] = sampprob
        return scale*scipy.array(csample)/sum(csample)

class DirichletMixture:

    def __init__(self, dirichlets, weights):
        assert len(dirichlets) == len(weights)
        # assert all(isinstance(d, Dirichlet) for d in dirichlets)
        assert len(set(d.numevents for d in dirichlets)) == 1
        # Check that weights is a list of numeric values by upcasting
        # them to floats.
        weights = scipy.array(weights, dtype=float)
        self.dirichlets, self.weights = dirichlets, weights
        self.numevents = self.dirichlets[0].numevents
        self.pairs = zip(self.dirichlets, self.weights)

    def posterior(self, counts):
        newdirichlets = [d.posterior(counts) for d in self.dirichlets]
        logweights = scipy.array(
            [scipy.log(p)+d.lcount_prob(counts)
             for d, p in self.pairs])
        weights = scipy.exp(logweights-max(logweights))
        neweights = weights/sum(weights)
        return DirichletMixture(newdirichlets, neweights)

    def sample(self):
        return self.dirichlets[statistics.sample(self.weights)
                               ].sample()

    def mean(self):
        return sum(w*d.mean() for w, d in zip(
            self.weights, self.dirichlets))/sum(self.weights)


