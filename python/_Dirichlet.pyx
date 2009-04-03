import scipy, random
from scipy.special import gammaln, gamma

cdef class Dirichlet:

    cdef readonly alphas, invalphas, lgalphas
    cdef readonly double norm,lnorm,sumalphas,lgsumalphas,lsumgalphas
    cdef readonly int numevents

    def __init__(self, alphas):
        
        self.alphas = scipy.array(alphas)
        if min(alphas) < 0:
            raise ValueError, 'Counts must be positive'
        self.invalphas = 1./self.alphas
        self.sumalphas = sum(alphas)
        self.lgsumalphas = gammaln(self.sumalphas)
        self.lgalphas = map(gammaln, alphas)
        self.lsumgalphas = sum(self.lgalphas)
        self.lnorm = (self.lgsumalphas - sum(self.lgalphas))
        self.norm = scipy.exp(self.lnorm)
        self.numevents = len(alphas)

    def _sample(self, scale=1):
        csample = scipy.zeros([self.numevents], 'd')
        for idx in range( self.numevents):
            while 1:
                p = random.gammavariate(self.alphas[idx], 1)
                if p != 0:
                    csample[idx] = p
                    break
        return scale*csample/sum(csample)

    def sample(self, scale=1):
        raise UnimplementedError
        csample = foo
        while csample.count(0):
            zeroidx = csample.index(0)
            csample[zeroidx] = random.expovariate(
                self.invalphas[zeroidx], 1)
        return scale*scipy.array(csample)/sum(csample)

    def lpdf(self, P):
        assert len(P) == self.numevents
        rv = self.lnorm
        for p, c in zip(P, self.alphas):
            rv = rv + (c-1)*log(p)
        return rv
        
    def pdf(self, P):
        return scipy.exp(self.lpdf(P))

    def posterior(self, counts):
        return self.__class__(self.alphas+counts)

    def __repr__(self):
        return '%s(%s)' % (self.__class__.__name__,
                           list(self.alphas))

    def lcount_prob(self, counts):
        """Return the log probability of drawing this set of counts in
        a sample of size sum(counts)."""
        assert len(counts) == self.numevents
        n = sum(counts)
        rv = gammaln(n+1)+self.lgsumalphas-gammaln(n+self.sumalphas)
        rv = rv - self.lsumgalphas
        for c, alpha in zip(counts, self.alphas):
            rv = rv + gammaln(c+alpha) - gammaln(c+1)
        return rv

    def count_prob(self, counts):
        """Return the probability  of drawing this set of  counts in a
        sample of size sum(counts)."""
        return scipy.exp(self.lcount_prob(counts))

    def multicount_lprob(self, counts):
        """Return the log probability of independently drawing each of
        these counts, given their totals"""
        rv = 0
        for ccounts in counts:
            rv = rv + self.lcount_prob(ccounts)
        return rv

    def mean(self):
        return self.alphas/self.sumalphas

    def __reduce__(self):
        return self.__class__, self.alphas

def sample(distribution):
    cdef double target, total
    cdef int i, length
    length = len(distribution)
    csum = sum(distribution)
    assert csum > 0
    target = random.uniform(0, csum)
    # It appears  to be important that  total actually be  cdef'd to a
    # float, or  that it  be assigned  a zero float,  rather than  a 0
    # integer.  Very occasional errors were occurring, otherwise.
    total = 0.
    for i from 0 <= i < length:
        total = total + distribution[i]
        if total >= target:
            return i
    print 'sample failed.  Values of distribution and target were',  distribution, target
    raise RuntimeError, 'Should never get here'

