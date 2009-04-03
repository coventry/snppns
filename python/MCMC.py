import itertools, time, scipy

def convergence_diagnostic(estimands):
    """Compute  the within- and  between-variances for  the estimands,
    and use this  to estimate the degree of  convergence so far.  From
    Gelman et al., p. 296."""
    n, m = float(len(estimands[0])), float(len(estimands))
    assert set(map(len, estimands)) == set([n]), 'Must be same lengths.'
    avgs = [sum(sample)/n for sample in estimands]
    totalavg = sum(avgs)/m # This step assumes equal length
    between = sum([(avg-totalavg)**2 for avg in avgs])
    between *= n/(m-1) # This step, too.
    within = sum(sum((est-avg)**2 for est in sample)/(n-1)
                  for avg, sample in zip(avgs, estimands))/m
    var = ((n-1)*within + between)/n
    comparisons = [p for p in (var/within).flatten()
                   if not scipy.isnan(p)]
    if comparisons:
        return max(comparisons)**0.5
    else:
        # All estimands were equal, across all the samples.  That's
        # pretty good convergence.  So just return 1.
        return 1

class MCMC:

    rhat_threshold = 1.1

    def __init__(self, sampler, numseries=5):
        """`sampler'  should  provide two  methods:  `sample' ,  which
        takes a sample  returns another, and `compute_estimand', which
        takes  a sample  and  returns the  associated estimand.   (And
        possibly  `get_samples'  which takes  a  set  of samples,  and
        returns the gibbs samples for all of them."""
        self.sampler = sampler
        self.samples = [sampler.sample()
                        for dummy in itertools.repeat(None,numseries)]
        self.estimands = [
            [] for dummy in itertools.repeat(None, numseries)]
        self.diagnostic_is_current = False
        self.min_series_length = 20
        self._Rhat = 1e6

    def Gibbs_iteration(self):
        """Sample a new set of rate matrices by Gibbs """
        if hasattr(self.sampler, 'get_samples'):
            self.samples = self.sampler.get_samples(self.samples)
        else:
            self.samples=[self.sampler.sample(sample)
                          for sample in self.samples]
        self.diagnostic_is_current = False
        assert len(self.samples) == len(self.estimands)
        for sample, estimand_list in zip(self.samples, self.estimands):
            estimand_list.append(self.sampler.compute_estimand(sample))

    def Rhat(self):
        """Return the estimand series convergence diagnostics."""
        if not self.diagnostic_is_current:
            # Chop off the first 5%  of each series when computing the
            # diagnostic, for "burn-in"
            startidx = int(0.05*len(self.estimands[0]))
            estimands=[series[startidx:] for series in self.estimands]
            self._Rhat = convergence_diagnostic(estimands)
            print self._Rhat
            self.diagnostic_is_current = True
        return self._Rhat

    def finished(self):
        return ((len(self.estimands[0]) >= 2) and
                (self.Rhat() < self.rhat_threshold) and
                (len(self.estimands[0]) >= self.min_series_length))

    def estimated_mean(self):
        """Return the mean of the estimands"""
        estimands=scipy.array(list(itertools.chain(*self.estimands)))
        return sum(estimands)/float(len(estimands))

    def estimate_mean(self, stop_time=1e200, save=lambda *x: None,
                      min_series_length=None):
        """Return the estimated mean posterior stationary distribution
        for this column, once we've sampled enough rate matrices to be
        sure of reasonable convergence. """
        if min_series_length is not None:
            self.min_series_length = min_series_length
        while (not self.finished()) and (time.time() < stop_time):
            assert ((len(self.estimands[0]) <= 3) or
                    (not scipy.isnan(self.Rhat())))
            if 0:
                import hotshot, hotshot.stats
                prof = hotshot.Profile("gibbs prof")
                prof.runcall(self.Gibbs_iteration)
                stats = hotshot.stats.load('gibbs prof')
                stats.strip_dirs()
                stats.sort_stats('time', 'calls')
                stats.print_stats(20)
            self.Gibbs_iteration()
            save(self.samples)
        return self.estimated_mean()

