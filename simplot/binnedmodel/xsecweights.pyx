# cython: profile=True

import numpy as np
cimport numpy as np

from simplot.sparsehist.sparsehist cimport SparseArray, SparseArrayIterator, array_bisect_right, sparse_array_interpolation

from libc.stdint cimport uint64_t
from libcpp.vector cimport vector

import itertools
from bisect import bisect_right

cimport cython
from cython.operator cimport preincrement, dereference
################################################################################

class XsecWeights:
    def __init__(self, nosel, weightcalc):
        shape = nosel.shape()
        self._shape = shape
        #self._arr = arr
        self._nosel = nosel
        self._xseccalc = weightcalc

    def __call__(self, pars):
        return self._eval(pars)

    def _ones(self):
        cdef SparseArray arr = SparseArray(self._shape)
        cdef SparseArray nosel = self._nosel
        cdef SparseArrayIterator it = nosel._data.begin();
        cdef SparseArrayIterator end = nosel._data.end();
        while it != end:
            arr._data[dereference(it).first] = 1.0
            preincrement(it)
        return arr

    def _eval(self, pars):
        #start off with array of ones
        arr = self._ones()
        #product of weights
        other = self._update(pars)
        return self._product(arr, other)

    def _update(self, pars):
        cdef np.ndarray[object, ndim=1] other = np.zeros(dtype=object, shape=(len(self._xseccalc),));
        for ii, calc in enumerate(self._xseccalc):
            calc.update(pars)
            other[ii] = calc.array()
        return other

#    @cython.boundscheck(False)
#    @cython.boundscheck(None)
#    def _product(self, SparseArray arr, np.ndarray[object, ndim=1] otherarrays):
#        cdef SparseArrayIterator it = arr._data.begin();
#        cdef SparseArrayIterator end = arr._data.end();
#        cdef vector[uint64_t] index;
#        cdef int N = len(otherarrays)
#        cdef SparseArray rhs;
#        cdef uint64_t key;
#        cdef double x = 0.0;
#        while it != end:
#            key = dereference(it).first
#            x = dereference(it).second
#            for ii in xrange(N):
#                rhs = otherarrays[ii]
#                index = rhs.decodekey(key)
#                x *= rhs.get(index)
#            arr.set(index, x)
#            preincrement(it)
#        return arr

    def _product(self, SparseArray arr, np.ndarray[object, ndim=1] otherarrays):
        cdef SparseArray rhs;
        for ii in xrange(len(otherarrays)):
            rhs = otherarrays[ii]
            arr *= rhs
        return arr

################################################################################

class NormWeightCalc:
    def __init__(self, shape, binmap, parname, parameternames):
        weightshape = [0 for s in shape]
        for dimnum in binmap.iterkeys():
            weightshape[dimnum] = shape[dimnum]
        arr = SparseArray(weightshape)
        self._shape = weightshape
        self._arr = arr
        self._binmap = binmap
        self._binkeys = []
        self._binvalues = []
        for k,v in binmap.iteritems():
            self._binkeys.append(k)
            self._binvalues.append(v)
        self._parnum = 0
        if not parname in parameternames:
            raise Exception("parameter not in list of names", parname, parameternames)
        for i, p in enumerate(parameternames):
            if parname == p:
                self._parnum = i
        self._reset()

    def update(self, pars):
        self._set(pars[self._parnum])
        return

    def _reset(self):
        return self._setall(1.0)

    def _setall(self, value):
        arr = self._arr
        ranges = [xrange(s) for s in self._shape if s > 0]
        index = list(self._shape)
        for v in itertools.product(*ranges):
            for ii in xrange(len(v)):
                index[self._binkeys[ii]] = v[ii]
            #print "setting", index, value
            arr[index] = value
        return

    def _set(self, value):
        arr = self._arr
        for index in self._iterindices():
            arr[index] = value
        return

    def _iterindices(self):
        index = list(self._shape)
        for v in itertools.product(*self._binvalues):
            for ii in xrange(len(v)):
                index[self._binkeys[ii]] = v[ii]
            yield index

    def array(self):
        return self._arr

################################################################################

cdef class InterpolatedWeightCalc:
    cdef int _parnum;
    cdef vector[double] _xvec;
    cdef list _yvec;
    cdef SparseArray _arr;
    #SparseArray self._arr;
    cdef str _parname;
    def __init__(self, nominalvalues, parvalues, arrays, parname, parameternames):
        self._check_is_sorted(parvalues)
        self._parnum = self._findparameter(parname, parameternames)
        self._xvec = parvalues
        self._yvec = [arr/nominalvalues for arr in arrays]
        self._arr = None
        self._parname = parname
        #check inputs
        if not len(arrays) == len(self._xvec):
            raise Exception("InterpolatedWeightCalc wrong number of input arrays", parname, len(parvalues), len(arrays))

    def _check_is_sorted(self, values, msg=None):
        l1 = list(values)
        l2 = list(l1)
        l2.sort()
        if not l1 == l2:
            if msg is None:
                msg = "InterpolatedWeightCalc expected sorted values."
            raise Exception(msg, values)

    def __str__(self):
        return "InterpolatedWeightCalc(%02.0f:%s, range=%s)" % (self._parnum, self._parname, ["%.2e"%x for x in self._xvec])

    def _findparameter(self, parname, parameternames):
        parnum = -1
        for i, p in enumerate(parameternames):
            if p == parname:
                parnum = i
        if parnum < 0:
            raise Exception("missing parameter", parname, parameternames)
        return parnum

    def array(self):
        return self._arr

    def update(self, pars):
        #for i in xrange(len(pars)):
        #    print "DEBUG", i, pars[i]
        x = pars[self._parnum]
        self._arr = self.eval(x)

    cdef SparseArray eval(self, double x):
        cdef int last = self._xvec.size() - 1
        if x <= self._xvec[0]:
            return self._yvec[0]
        if x >= self._xvec[last]:
            return self._yvec[last]
        #inside vector
        cdef i = array_bisect_right(self._xvec, x) - 1
        #print "DEBUG", len(yvec), x, xvec, i
        cdef SparseArray y0 = self._yvec[i]
        cdef SparseArray y1 = self._yvec[i+1]
        cdef double x0 = self._xvec[i]
        cdef double x1 = self._xvec[i+1]
        return self._interp(x, x0, x1, y0, y1)

    cdef SparseArray _interp(self, double x, double x0, double x1, SparseArray y0, SparseArray y1):
        cdef double f = (x-x0) / (x1-x0)
        #return f*y1 + (1.0-f)*y0
        return sparse_array_interpolation(f, y0, y1)

################################################################################

cdef class SimpleInterpolatedWeightCalc:
    cdef int _parnum;
    cdef vector[double] _xvec;
    cdef list _yvec;
    cdef vector[double] _arr;
    #SparseArray self._arr;
    cdef str _parname;
    def __init__(self, nominalvalues, parvalues, arrays, parname, parameternames):
        self._check_is_sorted(parvalues)
        self._parnum = self._findparameter(parname, parameternames)
        self._xvec = parvalues
        self._yvec = [arr/nominalvalues for arr in arrays]
        self._arr = None
        self._parname = parname
        #check inputs
        if not len(arrays) == len(self._xvec):
            raise Exception("SimpleInterpolatedWeightCalc wrong number of input arrays", parname, len(parvalues), len(arrays))

    def _check_is_sorted(self, values, msg=None):
        l1 = list(values)
        l2 = list(l1)
        l2.sort()
        if not l1 == l2:
            if msg is None:
                msg = "SimpleInterpolatedWeightCalc expected sorted values."
            raise Exception(msg, values)

    def __str__(self):
        return "SimpleInterpolatedWeightCalc(%02.0f:%s, range=%s)" % (self._parnum, self._parname, ["%.2e"%x for x in self._xvec])

    def _findparameter(self, parname, parameternames):
        parnum = -1
        for i, p in enumerate(parameternames):
            if p == parname:
                parnum = i
        if parnum < 0:
            raise Exception("missing parameter", parname, parameternames)
        return parnum

    def array(self):
        return self._arr

    def update(self, pars):
        #for i in xrange(len(pars)):
        #    print "DEBUG", i, pars[i]
        x = pars[self._parnum]
        self._arr = self.eval(x)

    cdef vector[double] eval(self, double x):
        cdef int last = self._xvec.size() - 1
        if x <= self._xvec[0]:
            return self._yvec[0]
        if x >= self._xvec[last]:
            return self._yvec[last]
        #inside vector
        cdef i = array_bisect_right(self._xvec, x) - 1
        #print "DEBUG", len(yvec), x, xvec, i
        cdef SparseArray y0 = self._yvec[i]
        cdef SparseArray y1 = self._yvec[i+1]
        cdef double x0 = self._xvec[i]
        cdef double x1 = self._xvec[i+1]
        return self._interp(x, x0, x1, y0, y1)

    cdef vector[double] _interp(self, double x, double x0, double x1, vector[double]& y0, vector[double]& y1):
        cdef double f = (x-x0) / (x1-x0)
        cdef vector[double] y = vector[double](y0.size())
        cdef int ii
        for ii in xrange(y.size()):
            y[ii] = f*y1[ii] + (1.0-f)*y0[ii]
        #return f*y1 + (1.0-f)*y0
        return y

################################################################################
