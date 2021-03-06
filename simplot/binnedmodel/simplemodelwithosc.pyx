import numpy as np
cimport numpy as np

from simplot.mc.statistics import safedivide
from simplot.binnedmodel.model import ProbabilityCache, OscParMode

DEF _DIM_ENU = 0
DEF _DIM_NUPDG = 1
DEF _DIM_RECO = 2
DEF _NUM_OSC_PARS = 6

cdef class SimpleBinnedModelWithOscillation:
    cdef list _parnames;
    cdef np.ndarray _eff;
    cdef np.ndarray N_nosel;
    cdef np.ndarray _N_nosel_projection;
    cdef np.ndarray _otherflav;
    cdef object _prob;
    cdef np.ndarray _cache3D;
    cdef np.ndarray _cache1D;
    cdef Py_ssize_t _num_enu_bins;
    cdef Py_ssize_t _num_reco_bins;

    def __init__(self, parnames, N_sel, N_nosel, enubinning, detdist, probabilitycalc=None, oscparmode=OscParMode.SINSQTHETA):
        self._num_enu_bins = N_sel.shape[_DIM_ENU]
        self._num_reco_bins = N_sel.shape[_DIM_RECO]
        #check shape
        for arr in [N_sel, N_nosel]:
            if not arr.shape[_DIM_ENU] == len(enubinning)-1:
                raise ValueError("Input array and enubinning array have different sizes.")
            if not arr.shape[_DIM_NUPDG] == 4:
                raise ValueError("Input array has wrong number of nupdg bins.")
            if not arr.shape[_DIM_RECO] > 0:
                raise ValueError("Input array has wrong number of reco bins.")
        self._parnames = parnames
        self._eff = safedivide(N_sel, N_nosel)
        self.N_nosel = N_nosel
        self._N_nosel_projection = np.sum(N_nosel, axis=2)
        self._otherflav = np.array([1,0,3,2], dtype=int)
        self._prob = ProbabilityCache(parnames, enubinning, [detdist], probabilitycalc=probabilitycalc, oscparmode=oscparmode)
        self._cache3D = np.copy(N_sel)
        self._cache1D = np.copy(np.sum(N_sel, axis=(_DIM_NUPDG, _DIM_ENU)))
        return

    def __call__(self, pars):
        return self.eval(pars)

    cdef np.ndarray[np.float64_t, ndim=1] eval(self, np.ndarray[np.float64_t, ndim=1] pars):
        self._updateprediction(pars)
        return np.multiply(self._syst_weights(pars), self._cache1D)

    cdef np.ndarray[np.float64_t, ndim=1] _syst_weights(self, np.ndarray[np.float64_t, ndim=1] pars):
        return pars[_NUM_OSC_PARS:]

    cdef void _updateprediction(self, pars):
        self._osc_flav_rotation(pars)
        np.multiply(self._eff, self._cache3D, out=self._cache3D)
        np.sum(self._cache3D, axis=(_DIM_NUPDG, _DIM_ENU), out=self._cache1D)
        return

#    @cython.boundscheck(False)
    cdef void _osc_flav_rotation(self, pars):
        #get inputs
        cdef np.ndarray[double, ndim=4] posc = self._prob.array
        cdef np.ndarray[np.float64_t, ndim=3] nominal = self.N_nosel
        cdef np.ndarray[np.float64_t, ndim=2] nominal_projection = self._N_nosel_projection
        cdef np.ndarray[np.float64_t, ndim=3] result = self._cache3D
        cdef Py_ssize_t Nenubins = self._num_enu_bins
        cdef Py_ssize_t Nrecobins = self._num_reco_bins
        cdef np.ndarray[dtype=Py_ssize_t, ndim=1] otherflav = self._otherflav
        #update oscillation probabilities
        self._prob.update(pars)
        #calculate
        cdef Py_ssize_t flav_j, flav_i, ienu, ireco
        cdef double pdis, papp, value, othervalue, nosc
        for flav_j in xrange(4):
            for ienu in xrange(Nenubins):
                flav_i = otherflav[flav_j]
                pdis = posc[ienu, 0, flav_j, flav_j]
                papp = posc[ienu, 0, flav_i, flav_j]
                value = nominal_projection[ienu, flav_j]
                othervalue = nominal_projection[ienu, flav_i]
                nosc = (pdis * value) + (papp * othervalue)
                weight = 0.0
                if value != 0.0:
                    weight = nosc / value
                for ireco in xrange(Nrecobins):
                    result[ienu, flav_j, ireco] = weight * nominal[ienu, flav_j, ireco]
        return

    @property
    def parameter_names(self):
        return self._parnames
