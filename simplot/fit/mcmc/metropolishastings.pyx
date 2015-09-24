# cython: profile=False

import numpy as np
cimport numpy as np

from libc.math cimport exp

from simplot.mc.montecarlo import generate_timed as _gen_timed
from simplot.mc.montecarlo import generate_events as _gen_events

################################################################################

class McMcSetupError(Exception):
    pass

################################################################################

cdef class MetropolisHastingsAlgorithm(object):
    cdef list _parameter_names
    cdef np.ndarray _start
    cdef np.ndarray _low
    cdef np.ndarray _high
    cdef np.ndarray _theta
    cdef object _function
    cdef int _total
    cdef _success
    cdef int _num_parameters
    cdef object random
    cdef object _proposal
    cdef object _random
    cdef double _likelihood
    def __init__(self, parameter_names, function, proposal, start, parameter_range, seed=None):
        self._parameter_names = parameter_names
        self._start = np.array(start)
        self._low = np.array([x[0] for x in parameter_range])
        self._high = np.array([x[1] for x in parameter_range])
        self._num_parameters = len(self._start)
        self._theta = self._start
        self._function = function
        self._proposal = proposal
        self._total = 0
        self._success = 0
        self._random = np.random.RandomState(seed)
        #compute the intial likelihood values
        self._likelihood = self._function(self._theta) 
        return

    @property
    def parameter_names(self):
        return self._parameter_names

    def __call__(self):
        return self.generate()

    def generate(self):
        theta0 = self._theta
        random = self._random
        cdef double pt, jt, pt0, jt0, r, prob, likelihood;
        while True:
            thetaProposal = self._proposal.generate(theta0)
            if self._check_range(thetaProposal):
                pt = self._function(thetaProposal)
                jt = self._proposal.logDensity(thetaProposal, theta0)
                pt0 = self._function(theta0)
                jt0 = self._proposal.logDensity(theta0, thetaProposal)
                r = pt - jt - pt0 + jt0
                prob = exp(r)
                if prob > 1.0:
                    prob = 1.0
                if random.uniform() < prob:
                    thetaNext = np.copy(thetaProposal)
                    likelihood = pt
                    self._success += 1
                else:
                    thetaNext = theta0
                    likelihood = pt0
                self._theta = thetaNext
                self._likelihood = likelihood
                self._total += 1
                break
        return self._theta, self._likelihood

    cdef _check_range(self, np.ndarray[dtype=double, ndim=1] proposal):
        cdef int N = self._num_parameters
        cdef np.ndarray[dtype=double, ndim=1] low = self._low
        cdef np.ndarray[dtype=double, ndim=1] high = self._high
        for ii in xrange(N):
            if not (low[ii] <= proposal[ii] <= high[ii]):
                return False
        return True

    def efficiency(self):
        eff = 0.0
        if self._total > 0:
            eff = float(self._success) / float(self._total)
        return eff

################################################################################

cdef class AdaptiveMetropolisHastingsAlgorithm(MetropolisHastingsAlgorithm):
    cdef int _period
    cdef int _converged
    cdef int _save_all
    cdef np.ndarray _data
    def __init__(self, *args, **kwargs):
        self._period = 10000
        self._converged = False
        self._save_all = False
        super(AdaptiveMetropolisHastingsAlgorithm, self).__init__(*args, **kwargs)
        self._data = np.zeros(shape=(self._period, len(self._parameter_names)))

    def generate(self):
        if self._save_all:
            return self._generate_all()
        else:
            return self._generate_only_converged()

    def _generate_only_converged(self):
        while not self._converged:
            self._adaptivegenerate()
        return MetropolisHastingsAlgorithm.generate(self)

    def _generate_all(self):
        if not self._converged:
            return self._adaptivegenerate()
        else:
            #converged, behave as normal
            return MetropolisHastingsAlgorithm.generate(self)

    def _adaptivegenerate(self):
        theta, likelihood = MetropolisHastingsAlgorithm.generate(self)
        self._data[self._total-1] = theta # total - 1 as total is incremented in generate
        #not converged, try to adapt
        if self._total % self._period == 0:
            if self._total > 0:
                self._converged = self._proposal.adapt(self.efficiency(), self._data)
                #reset efficiency counters
                self._success = 0
                self._total = 0
        return theta, likelihood

def generate_timed(toymc, maxseconds, maxevents=None, name=None):
    for x in _gen_timed(toymc=toymc, maxseconds=maxseconds, maxevents=maxevents, name=name):
        yield x
    if name is not None:
        _print_eff_str(toymc)
    return

def generate_events(toymc, n, name=None):
    for x in _gen_events(toymc=toymc, n=n, name=name):
        yield x
    if name is not None:
        _print_eff_str(toymc)
    return

def _print_eff_str(toymc):
    try:
        eff = toymc.efficiency()
        effmsg = "%.0f%% (%.2e)" % ((eff * 100.0), eff)
    except AttributeError:
        effmsg = "unknown"
    print "efficiency = ", effmsg
    return
