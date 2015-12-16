# cython: profile=False

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
from libc.math cimport asin, sqrt

ctypedef std_map[uint64_t, double].iterator SparseArrayIterator

DEF NO_DET_DIM = 9999

DEF _CODE_THETA = 0
DEF _CODE_SINSQ2THETA = 1
DEF _CODE_SINSQTHETA = 2

################################################################################

class OscParMode:
    THETA = "theta"
    SINSQ2THETA = "sinsq2theta"
    SINSQTHETA = "sinsqtheta"

    @classmethod
    def toint(cls, mode):
        result = None
        if mode == cls.THETA:
            result = _CODE_THETA
        elif mode == cls.SINSQ2THETA:
            result = _CODE_SINSQ2THETA
        elif mode == cls.SINSQTHETA:
            result = _CODE_SINSQTHETA
        return result

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
    cdef SparseArray _N_sel;
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
    cdef _osc_flux_weights;
    cdef list _parnames;

    def __init__(self, parnames, N_sel, N_nosel, obs, enudim, flavdim, detdim, detdist, flux_weights=None, xsec_weights=None, probabilitycalc=None, oscparmode=OscParMode.SINSQTHETA):
        self._parnames = parnames
        self._shape = N_sel.array().shape()
        self._eff = N_sel.array() / N_nosel.array()
        self._N_sel = N_sel.array()
        self.N_nosel = N_nosel.array()
        self._obs = obs
        self._flav_dimension = flavdim
        self._enu_dimension = enudim
        if detdim is None:
            detdim = NO_DET_DIM
        self._det_dimension = detdim
        self._otherflav = [1,0,3,2]
        enubinning = N_sel.binning()[enudim]
        self._prob = ProbabilityCache(parnames, enubinning, detdist, probabilitycalc=probabilitycalc, oscparmode=oscparmode)
        if flux_weights is None:
            flux_weights = lambda x: _identity(self._shape)
        self._flux_weights = flux_weights
        if xsec_weights is None:
            xsec_weights = lambda x: _identity(self._shape)
        self._xsec_weights = xsec_weights
        self._osc_flux_weights = OscFluxWeights(N_nosel, enudim, flavdim, detdim, self._prob)
        return

    def __call__(self, pars):
        return self.eval(pars)

    cdef eval(self, pars):
        return self._xsec_weights(pars) * self._eff * self._osc_flav_rotation(pars, self._flux_weights(pars) * self.N_nosel)
        #return self._xsec_weights(pars) * (self._eff * (self._osc_flux_weights(pars) * (self._flux_weights(pars) * self.N_nosel)))
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
    cdef np.ndarray _detdist
    cdef public np.ndarray array;
    cdef uint64_t _theta12;
    cdef uint64_t _theta23;
    cdef uint64_t _theta13;
    cdef uint64_t _deltacp;
    cdef uint64_t _sdm;
    cdef uint64_t _ldm;

    cdef int _oscparmode;

    cdef double _previous_theta12;
    cdef double _previous_theta23;
    cdef double _previous_theta13;
    cdef double _previous_deltacp;
    cdef double _previous_sdm;
    cdef double _previous_ldm;
    cdef np.ndarray _flav_map;

    def __init__(self, parnames, enubinning, detdist, probabilitycalc=None, oscparmode=OscParMode.SINSQTHETA):
        self._flav_map = np.array([#appearance
                                                              (0, 0, 2, 2, 1),
                                                              (1, 1, 1, 1, 1),
                                                              (2, 2, 2, 2, -1),
                                                              (3, 3, 1, 1, -1),
                                                              #disappearance
                                                              (0, 1, 2, 1, 1),
                                                              (1, 0, 1, 2, 1),
                                                              (2, 3, 2, 1, -1),
                                                              (3, 2, 1, 2, -1),
                                                              
            ], dtype=np.intc)
        self._previous_theta12 = 0.0
        self._previous_theta23 = 0.0
        self._previous_theta13 = 0.0
        self._previous_deltacp = 0.0
        self._previous_sdm = 0.0
        self._previous_ldm = 0.0
        self._oscparmode = OscParMode.toint(oscparmode)
        if probabilitycalc is None and any([d>0 for d in detdist]):
            raise Exception("invalid probability calculator", probabilitycalc, detdist)
        self._parse_parameter_names(parnames)
        self._prob = probabilitycalc #crootglobes.Probability()
        enubinning = np.array(enubinning)
        self._enuarray = (enubinning[1:] + enubinning[:-1])/2.0
        self._detdist = np.array(detdist, dtype=np.intc)
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
            if p == "theta12" or p == "sinsq2theta12" or p == "sinsqtheta12":
                theta12 = index
            elif p == "theta23" or p == "sinsq2theta23" or p == "sinsqtheta23":
                theta23 = index
            elif p == "theta13" or p == "sinsq2theta13" or p == "sinsqtheta13":
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

    cdef _haschanged(self, np.ndarray[double, ndim=1] pars):
        cdef double theta12 = pars[self._theta12]
        cdef double theta23 = pars[self._theta23]
        cdef double theta13 = pars[self._theta13]
        cdef double deltacp = pars[self._deltacp]
        cdef double sdm = pars[self._sdm]
        cdef double ldm = pars[self._ldm]
        cdef int allsame = (theta12 == self._previous_theta12) \
                          and (theta23 == self._previous_theta23) \
                          and (theta13 == self._previous_theta13) \
                          and (deltacp == self._previous_deltacp) \
                          and (ldm == self._previous_ldm) \
                          and (sdm == self._previous_sdm)
        self._previous_theta12 = theta12
        self._previous_theta13 = theta13
        self._previous_theta23 = theta23
        self._previous_deltacp = deltacp
        self._previous_ldm = ldm
        self._previous_sdm = sdm
        return not allsame

    def update(self, np.ndarray[double, ndim=1] pars):
        return self._update(pars)

    cdef _update(self, np.ndarray[double, ndim=1] pars):
        cdef double theta12, theta23, theta13, deltacp, sdm, ldm
        cdef int oscparmode
        if self._prob and self._haschanged(pars):
            #print "DEBUG setting", pars[self._theta12], pars[self._theta23], pars[self._theta13], pars[self._deltacp], pars[self._sdm], pars[self._ldm]
            theta12 = pars[self._theta12]
            theta23 = pars[self._theta23]
            theta13 = pars[self._theta13]
            deltacp = pars[self._deltacp]
            sdm = pars[self._sdm]
            ldm = pars[self._ldm]
            oscparmode = self._oscparmode
            if oscparmode == _CODE_SINSQTHETA:
                theta12 = invsinsqtheta(theta12)
                theta23 = invsinsqtheta(theta23)
                theta13 = invsinsqtheta(theta13)
            elif oscparmode == _CODE_SINSQ2THETA:
                theta12 = invsinsq2theta(theta12)
                theta23 = invsinsq2theta(theta23)
                theta13 = invsinsq2theta(theta13)
            self._prob.setAll(theta12, theta23, theta13, deltacp, sdm, ldm)
            self._prob.update()
            self._fillcache()

    cdef _fillcache(self):
        #get inputs
        prob = self._prob
        cdef np.ndarray[double, ndim=1] enuarray = self._enuarray;
        cdef np.ndarray[double, ndim=4] array = self.array;
        cdef np.ndarray[int, ndim=2] flavmap = self._flav_map;
        cdef np.ndarray[int, ndim=1] detdistarray = self._detdist;
        #temporary variables
        cdef int detbin, detdist, enubin;
        cdef int flav_i, flav_j, flav_init, flav_final, cp;
        cdef double enu, p;
        #iterate over detector bins
        for detbin in xrange(detdistarray.shape[0]):
            detdist = detdistarray[detbin]
            if detdist == 0:
                continue
            #update baseline
            prob.setBaseline(detdist)
            #for flav_i, flav_j, flav_init, flav_final, cp in [#appearance
            #                                                  (0, 0, 2, 2, 1),
            #                                                  (1, 1, 1, 1, 1),
            #                                                  (2, 2, 2, 2, -1),
            #                                                  (3, 3, 1, 1, -1),
            #                                                  #disappearance
            #                                                  (0, 1, 2, 1, 1),
            #                                                  (1, 0, 1, 2, 1),
            #                                                  (2, 3, 2, 1, -1),
            #                                                  (3, 2, 1, 2, -1),
            #                                                  
            #]:
            for flavrow in xrange(flavmap.shape[0]):
                flav_i = flavmap[flavrow, 0]
                flav_j = flavmap[flavrow, 1]
                flav_init = flavmap[flavrow, 2]
                flav_final = flavmap[flavrow, 3]
                cp = flavmap[flavrow, 4]
                #iterate over enubins
                for enubin in xrange(enuarray.shape[0]):
                    enu = enuarray[enubin]
                    p = prob.getVacuumProbability(flav_init, flav_final, enu, cp)
                    #array[enubin][detbin][flav_i][flav_j] = p
                    array[enubin,detbin,flav_i,flav_j] = p
        return

cdef double invsinsqtheta(double x):
    if x < 0.0:
        x = abs(x)
    if x > 1.0:
        x = 1.0
    return asin(sqrt(x))

cdef double invsinsq2theta(double x):
    if x < 0.0:
        x = abs(x)
    if x > 1.0:
        x = 1.0
    return asin(sqrt(x)) / 2.0

################################################################################

cdef class OscFluxWeights:

    cdef np.ndarray _weights;
    cdef np.ndarray _nominal;
    cdef _prob;
    cdef SparseArray _sparse_weights;
    cdef uint64_t _enudim;
    cdef uint64_t _flavdim;
    cdef uint64_t _detdim;
    cdef np.ndarray _otherflav;

    def __init__(self, N_nosel, enudim, flavdim, detdim, prob):
        #determine output shape
        shape = list(N_nosel.array().shape())
        #setup cache arrays
        if detdim == NO_DET_DIM:
            ndet = 1
        else:
            ndet = shape[detdim]
        self._nominal = np.zeros(shape=(shape[enudim], shape[flavdim], ndet), dtype=float)
        self._init_nominal(N_nosel, self._nominal, enudim, flavdim, detdim)
        self._weights = np.copy(self._nominal)
        self._prob = prob
        self._enudim = enudim
        self._flavdim = flavdim
        self._detdim = detdim
        self._otherflav = np.array([1, 0, 3, 2], dtype=np.uint64)
        #setup output array
        for i in xrange(len(shape)):
            if not ((i == enudim) or (i == flavdim) or (i == detdim)):
                shape[i] = 0
        self._sparse_weights = SparseArray(shape)
        self._init_sparse_weights(self._sparse_weights)

    def _init_nominal(self, N_nosel, nominal, enudim, flavdim, detdim):
        if detdim == NO_DET_DIM:
            projdim = (enudim, flavdim)
            proj = N_nosel.array().project(projdim)
            shape = proj.shape()
            for ienu, iflav in itertools.product(xrange(shape[0]), xrange(shape[1])):
                nominal[ienu, iflav, 0] = proj[ienu, iflav, 0]
        else:
            projdim = (enudim, flavdim, detdim)
            proj = N_nosel.array().project(projdim)
            shape = proj.shape()
            for ienu, iflav, idet in itertools.product(xrange(shape[0]), xrange(shape[1]), xrange(shape[2])):
                nominal[ienu, iflav, idet] = proj[ienu, iflav, idet]
        return

    def _init_sparse_weights(self, arr):
        shape = arr.shape()
        ranges = [xrange(max(s, 1)) for s in shape]
        for index in itertools.product(*ranges):
            arr[index] = 1.0
        return

    def __call__(self, pars):
        self._update_weights_array(pars)
        self._update_sparse_array()
        return self._sparse_weights

    @cython.boundscheck(True)
    cdef _update_weights_array(self, pars):
        cdef np.ndarray[double, ndim=3] nominal = self._nominal
        cdef np.ndarray[double, ndim=3] weights = self._weights
        cdef np.ndarray[uint64_t, ndim=1] otherflav = self._otherflav;
        self._prob.update(pars)
        cdef np.ndarray[double, ndim=4] posc = self._prob.array
        cdef int numenubins = nominal.shape[0]
        cdef int numflavbins = nominal.shape[1]
        cdef int numdetbins = nominal.shape[2]
        #declare variables for loop
        cdef double pdis, papp, nom, nom_i, nom_j, w;
        cdef Py_ssize_t ienu, iflav, idet, flav_i, flav_j;
        for ienu in xrange(numenubins):
            for iflav in xrange(numflavbins):
                for idet in xrange(numdetbins):
                    flav_j = iflav
                    flav_i = otherflav[iflav]
                    pdis = posc[ienu, idet, flav_j, flav_j]
                    papp = posc[ienu, idet, flav_i, flav_j]
                    nom_j = nominal[ienu, flav_j, idet]
                    nom_i = nominal[ienu, flav_i, idet]
                    nosc_j = (pdis * nom_j) + (papp * nom_i)
                    if nom_j != 0.0:
                        w = nosc_j / nom_j
                    else:
                        w = 0.0
                    weights[ienu, flav_j, idet] = w
        return

    cdef _update_sparse_array(self):
        cdef np.ndarray[double, ndim=3] weights = self._weights
        cdef SparseArray arr = self._sparse_weights
        cdef SparseArrayIterator it = arr._data.begin()
        cdef SparseArrayIterator end = arr._data.end()
        cdef double value;
        cdef uint64_t enudim = self._enudim
        cdef uint64_t flavdim = self._flavdim
        cdef uint64_t detdim = self._detdim
        cdef uint64_t key, ienu, iflav, idet;
        while it != end:
            key = dereference(it).first
            index = arr.decodekey(key)
            ienu = index[enudim]
            iflav = index[flavdim]
            if detdim == NO_DET_DIM:
                idet = 0
            else:
                idet = index[detdim]
            value = weights[ienu, iflav, idet]
            dereference(it).second = value
            preincrement(it)
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
