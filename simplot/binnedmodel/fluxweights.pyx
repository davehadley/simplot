# cython: profile=True
from simplot.sparsehist.sparsehist cimport SparseArray

from libcpp.vector cimport vector
from libc.stdint cimport uint64_t

import re

_DEFAULT_FLAV_BINMAP = {"numu" : 0, "nue" : 1, "antinumu" : 2, "antinue" : 3}

cdef class FluxWeights:
    cdef SparseArray _arr
    cdef vector[uint64_t] _parindex
    cdef vector[uint64_t] _keys

    def __cinit__(self, parnames, shape, enudim, flavdim, detdim, beammodedim, parametermap=None, flavbinmap=_DEFAULT_FLAV_BINMAP):
        #setup array
        fluxshape = [0 for s in shape]
        fluxshape[enudim] = shape[enudim]
        fluxshape[flavdim] = shape[flavdim]
        fluxshape[detdim] = shape[detdim]
        fluxshape[beammodedim] = shape[beammodedim]
        self._arr = SparseArray(fluxshape)
        # map indices to parameters
        if parametermap is None:
            keys, parindex = self._default_parameter_map(shape, fluxshape, enudim, flavdim, detdim, beammodedim, parnames, flavbinmap)
        else:
            keys, parindex = self._build_parameter_map(shape, fluxshape, enudim, flavdim, detdim, beammodedim, parnames, parametermap, flavbinmap)
        self._keys = keys
        self._parindex = parindex

    def _default_parameter_map(self, shape, fluxshape, enudim, flavdim, detdim, beammodedim, parnames, flavbinmap):
        keys = []
        parindex = []
        found = set()
        for parindex, par in enumerate(parnames):
            match = re.match("f_(.*)_(.*)_(.*)_(.*)", par)
            if match:
                found.add(par)
                detbin, beambin, flav, enubin = match.groups()
                index = list(fluxshape)
                index[enudim] = int(enubin)
                index[detdim] = int(detbin)
                index[flavdim] = int(flavbinmap[flav])
                index[beammodedim] = int(beambin)
                key = self._arr.key(index)
                keys.append(key)
                parindex.append(parindex)
        #check all parameters have been found
        expected = shape[enudim] * shape[flavdim] * shape[detdim] * shape[beammodedim]
        #expected = 1
        #for s in shape:
        #    if not s == 0:
        #        expected *= s
        if not len(found) == expected:
            raise Exception("Wrong number of flux parameters found.", len(found), expected)
        return keys, parindex

    def _build_parameter_map(self, shape, fluxshape, enudim, flavdim, detdim, beammodedim, parnames, parametermap, flavbinmap):
        keys = []
        parindex = []
        found = set()
        for par, range_ in parametermap.iteritems():
            #find parameter index
            pindex = parnames.index(par)
            #find all bins in range
            for detbin, beambin, flav, enubin in range_:
                index = list(fluxshape)
                index[enudim] = int(enubin)
                index[detdim] = int(detbin)
                index[flavdim] = int(flavbinmap[flav])
                index[beammodedim] = int(beambin)
                key = self._arr.key(index)
                keys.append(key)
                parindex.append(pindex)
                found.add(tuple(index))
        #check all parameters have been found
        expected = shape[enudim] * shape[flavdim] * shape[detdim] * shape[beammodedim]
        if not len(found) == expected:
            raise Exception("Wrong number of flux parameters found.", len(found), expected)
        return keys, parindex

    def __call__(self, pars):
        _update(self, pars)
        return self._arr

cdef void _update(FluxWeights self, vector[double]& pars):
        cdef int ii
        cdef uint64_t key
        cdef double value
        for ii in xrange(self._keys.size()):
            key = self._keys[ii]
            parindex = self._parindex[ii]
            value = pars[parindex]
            self._arr._data[key] = value
        return
