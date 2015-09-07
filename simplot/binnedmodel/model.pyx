# cython: profile=True

#from sparsehist import SparseArray
from simplot.sparsehist.sparsehist cimport SparseArray
from simplot.sparsehist.sparsehist cimport std_map
import numpy as np
cimport numpy as np

import cython
from cython.operator cimport dereference, preincrement

import itertools
import StringIO

#from rootglobes import crootglobes

from libcpp.vector cimport vector
from libc.stdint cimport uint64_t

ctypedef std_map[uint64_t, double].iterator SparseArrayIterator

DEF NO_DET_DIM = 9999

################################################################################

cdef class BinnedModel:
    cdef SparseArray _N_sel;
    cdef _flux_weights;
    cdef _xsec_weights;
    cdef list _parnames;
    cdef vector[uint64_t] _obs;
    def __init__(self, parnames, N_sel, obs, flux_weights=None, xsec_weights=None):
        self._parnames = parnames
        self._obs = obs
        #self._shape = N_sel.array().shape()
        self._N_sel = N_sel.array()
        if flux_weights is None:
            flux_weights = lambda x: _identity(self._N_sel.shape())
        self._flux_weights = flux_weights
        if xsec_weights is None:
            xsec_weights = lambda x: _identity(self._N_sel.shape())
        self._xsec_weights = xsec_weights
        return

    def __call__(self, pars):
        return self.eval(pars)

    cdef eval(self, pars):
        return self._xsec_weights(pars) * (self._flux_weights(pars) * self._N_sel)

    def observable(self, pars):
        return self.eval(pars).project(self._obs)

    def parameter_names(self):
        return self._parnames

################################################################################

cdef class BinnedModelWithOscillation:
    cdef vector[uint64_t] _shape;
    #cdef SparseArray N_sel;
    cdef SparseArray _eff;
    cdef SparseArray N_nosel;
    cdef vector[uint64_t] _obs;
    cdef uint64_t _flav_dimension;
    cdef uint64_t _enu_dimension;
    cdef uint64_t _det_dimension;
    cdef vector[uint64_t] _otherflav;
    cdef _prob;
    cdef _flux_weights;
    cdef _xsec_weights;
    cdef list _parnames;

    def __init__(self, parnames, N_sel, N_nosel, obs, enudim, flavdim, detdim, detdist, flux_weights=None, xsec_weights=None):
        self._parnames = parnames
        self._shape = N_sel.array().shape()
        self._eff = N_sel.array() / N_nosel.array()
        #self.N_sel = N_sel.array()
        self.N_nosel = N_nosel.array()
        self._obs = obs
        self._flav_dimension = flavdim
        self._enu_dimension = enudim
        if detdim is None:
            detdim = NO_DET_DIM
        self._det_dimension = detdim
        self._otherflav = [1,0,3,2]
        enubinning = N_sel.binning()[enudim]
        self._prob = ProbabilityCache(parnames, enubinning, detdist)
        if flux_weights is None:
            flux_weights = lambda x: _identity(self._shape)
        self._flux_weights = flux_weights
        if xsec_weights is None:
            xsec_weights = lambda x: _identity(self._shape)
        self._xsec_weights = xsec_weights
        return

    def __call__(self, pars):
        return self.eval(pars)

    cdef eval(self, pars):
        return self._xsec_weights(pars) * self._eff * self._osc_flav_rotation(pars, self._flux_weights(pars) * self.N_nosel)
        # avoid unneccessary new copies
        #cdef SparseArray r = self._flux_weights(pars) * self.N_nosel
        #r = self._osc_flav_rotation(pars, self._flux_weights(pars) * self.N_nosel)
        #r *= self._eff
        #r *= self._xsec_weights(pars)
        #return r

    @cython.boundscheck(False)
    cdef SparseArray _osc_flav_rotation(self, pars, SparseArray arr):
        self._prob.update(pars)
        cdef np.ndarray[double, ndim=4] posc = self._prob.array
        cdef uint64_t ienu = self._enu_dimension
        cdef uint64_t idet = self._det_dimension
        cdef uint64_t iflav = self._flav_dimension
        cdef vector[uint64_t] otherflav = self._otherflav
        cdef SparseArray result = SparseArray(self._shape)
        cdef Py_ssize_t enu, det, flav_i, flav_j;
        cdef double pdis, papp;
        cdef SparseArrayIterator it = arr._data.begin()
        cdef SparseArrayIterator end = arr._data.end()
        cdef uint64_t key
        cdef double value, othervalue
        cdef vector[uint64_t] index
        while it != end:
            key = dereference(it).first
            index = arr.decodekey(key)
            value = dereference(it).second
        #for index, value in arr:
            #determine oscillation probability
            enu = index[ienu]
            if idet == NO_DET_DIM:
                det = 0
            else:
                det = index[idet]
            flav_j = index[iflav]
            flav_i = otherflav[flav_j]
            #pdis = posc[enu][det][flav_j][flav_j]
            pdis = posc[enu, det, flav_j, flav_j]
            papp = posc[enu, det, flav_i, flav_j]
            #get N_otherflav
            index[iflav] = flav_i 
            othervalue = arr.get(index)
            index[iflav] = flav_j # reset index
            nosc = (pdis * value) + (papp * othervalue)
            result.set(index, nosc)
            #increment iterator
            preincrement(it)
        return result

    def observable(self, pars):
        return self.eval(pars).project(self._obs)

    def parameter_names(self):
        return self._parnames

    def __str__(self):
        sio = StringIO.StringIO()
        print >>sio, "BinnedModelWithOscillationModel(%s pars, %.2e bins, %.2e max bins)" % (len(self._parnames), len(self.N_nosel), self.N_nosel.max_size())
        indent = "    "
        print >>sio, indent, "shape=", str(self._shape)
        print >>sio, indent, "obs=%s, flavdim=%s, enudim=%s, detdim=%s"
        print >>sio, "sel=%.2e, nosel=%.2e" % ((self._eff*self.N_nosel).sum(), self.N_nosel.sum())
        for ipar, par in enumerate(self.parameter_names()):
            print >>sio, indent, "par%02.0f : %s" % (ipar, par)
        return sio.getvalue()

################################################################################

class ObservableRateVector:
    """Wrapper for model that only returns the observable rate vector."""
    def __init__(self, model):
        self._model = model
    def __call__(self, pars):
        return self._model.observable(pars).flatten()
    def observable(self, pars):
        return self(pars)
    def eval(self, pars):
        return self(pars)
    @property
    def parameter_names(self):
        return self._model.parameter_names()

################################################################################

class SumObservableRateVector:
    """Wrapper for model that only returns the observable rate vector."""
    def __init__(self, model1, model2, frac1):
        self._v1 = ObservableRateVector(model1)
        self._v2 = ObservableRateVector(model2)
        self._frac = frac1
    def __call__(self, pars):
        v1 = np.copy(self._v1(pars))
        v2 = np.copy(self._v2(pars))
        frac = self._frac
        ret = (v1 * frac) + (v2 * (1.0-frac))
        #print "DEBUG SumObservableRateVector", frac, sum(v1), sum(v2)
        #print "DEBUG SumObservableRateVector", sum(ret)
        #raw_input("wait")
        return ret
    def observable(self, pars):
        return self(pars)
    def eval(self, pars):
        return self(pars)
    @property
    def parameter_names(self):
        return self._v1.parameter_names

################################################################################

class ConcatenateObservableRateVector:
    """Wrapper for model that only returns the observable rate vector."""
    def __init__(self, ratevectors):
        self._v = ratevectors
    def __call__(self, pars):
        v = [v(pars) for v in self._v]
        #print "DEBUG ConcatenateObservableRateVector", v
        return np.concatenate(v)
    def observable(self, pars):
        return self(pars)
    def eval(self, pars):
        return self(pars)
    @property
    def parameter_names(self):
        return self._v[0].parameter_names


################################################################################

cdef class ProbabilityCache:
    cdef np.ndarray _enuarray
    cdef _prob
    cdef _detdist
    cdef public np.ndarray array;
    cdef uint64_t _theta12;
    cdef uint64_t _theta23;
    cdef uint64_t _theta13;
    cdef uint64_t _deltacp;
    cdef uint64_t _sdm;
    cdef uint64_t _ldm;
    def __init__(self, parnames, enubinning, detdist, probabilitycalc=None):
        self._parse_parameter_names(parnames)
        self._prob = probabilitycalc #crootglobes.Probability()
        enubinning = np.array(enubinning)
        self._enuarray = (enubinning[1:] + enubinning[:-1])/2.0
        self._detdist = detdist
        enudim = len(enubinning) - 1
        detdim = len(detdist)
        self.array = np.ones(shape=(enudim, detdim, 4, 4), dtype=float)
        for enubin, detbin, flav_i, flav_j in itertools.product(xrange(enudim), xrange(detdim), xrange(4), xrange(4)):
            if flav_i == flav_j:
                self.array[enubin, detbin, flav_i, flav_j] = 1.0
            else:
                self.array[enubin, detbin, flav_i, flav_j] = 0.0

    def _parse_parameter_names(self, parnames):
        theta12 = None
        theta23 = None
        theta13 = None
        deltacp = None
        sdm = None
        ldm = None
        for index, p in enumerate(parnames):
            if p == "theta12":
                theta12 = index
            elif p == "theta23":
                theta23 = index
            elif p == "theta13":
                theta13 = index
            elif p == "deltacp":
                deltacp = index
            elif p == "sdm":
                sdm = index
            elif p == "ldm":
                ldm = index
        if any(p is None for p in [theta12, theta23, theta13, deltacp, sdm, ldm]):
            raise Exception("ProbabilityCache cannot find all oscillation parameters", parnames)
        self._theta12, self._theta23, self._theta13, self._deltacp, self._sdm, self._ldm = theta12, theta23, theta13, deltacp, sdm, ldm
        return
        

    def update(self, pars):
        if self._prob:
            self._prob.setAll(pars[self._theta12], pars[self._theta23], pars[self._theta13], pars[self._deltacp], pars[self._sdm], pars[self._ldm])
            self._prob.update()
            for detbin, detdist in enumerate(self._detdist):
                if detdist == 0:
                    continue
                self._prob.setBaseline(detdist)
                for flav_i, flav_j, flav_init, flav_final, cp in [#appearance
                                                              (0, 0, 2, 2, 1),
                                                              (1, 1, 1, 1, 1),
                                                              (2, 2, 2, 2, -1),
                                                              (3, 3, 1, 1, -1),
                                                              #disappearance
                                                              (0, 1, 2, 1, 1),
                                                              (1, 0, 1, 2, 1),
                                                              (2, 3, 2, 1, -1),
                                                              (3, 2, 1, 2, -1),
                                                              
                ]:
                    for enubin, enu in enumerate(self._enuarray):
                        p = self._prob.getVacuumProbability(flav_init, flav_final, enu, cp)
                        self.array[enubin][detbin][flav_i][flav_j] = p
        return

################################################################################

def _scalar(n, shape):
    s = [0 for s in shape]
    arr = SparseArray(s)
    arr.set(s, n)
    return arr

def _identity(shape):
    return _scalar(1, shape)

def _zero(self, shape):
    return _scalar(0, shape)

################################################################################
