# cython: profile=False
DEF PROFILE_FLAG = False
DEF SHAPE_IS_IDENTICAL = 1
DEF SHAPE_IS_COMPATIBLE = 2

#from libcpp.unordered_map cimport unordered_map
#from libcpp.map cimport map as std_map
from simplot.sparsehist.unordered_map cimport unordered_map as std_map
from libcpp.vector cimport vector

from bisect import bisect_right

import numpy
cimport numpy

cimport cython
from cython.operator cimport preincrement, dereference

from simplot.mplot.histogram import HistogramND, HistogramNDLabel

from unordered_map cimport unordered_map as std_map
from libcpp.vector cimport vector
from libcpp cimport bool
from libc.stdint cimport uint64_t

###############################################################################

cdef class SparseArray:

    #cdef vector[uint64_t] _shape
    #cdef vector[uint64_t] _dimscale
    #cdef std_map[uint64_t, double] _data

    def __cinit__(self, vector[uint64_t] shape, int minsize=10**7):
        self._shape = shape
        self.reserve(minsize)
        cdef uint64_t cumprod = 1
        for s in shape:
            if s > 0:
                self._dimscale.push_back(cumprod)
                cumprod *= s
            else:
                self._dimscale.push_back(0)

    def sum(self):
        cdef SparseArrayIterator it = self._data.begin()
        cdef SparseArrayIterator end = self._data.end()
        cdef double total = 0.0
        while it != end:
            total += dereference(it).second
            preincrement(it)
        return total

    def shape(self):
        return self._shape

    def reserve(self, size):
        #self._data.rehash(size)
        return

    def __str__(self):
        return r"SparseArray(%.2e%%, %.2e/%.2e)" % (100.*self.occupancy(), self.actual_size(), self.max_size())
        #return r"SparseArray(%.2e%%, %.2e/%.2e, buckets=%s, load_factor=%s, max_load_factor=%s)" % (100.*self.occupancy(), self.actual_size(), self.max_size(), self._data.bucket_count(), self._data.load_factor(), self._data.max_load_factor())

    def __len__(self):
        return self.actual_size()

    def __iter__(self):
        for key, value in self._data:
            yield self.decodekey(key), value

    def __setitem__(self, index, double value):
        self.set(index, value)
        return

    def __getitem__(self, index):
        return self.get(index)

    def max_size(self):
        cdef uint64_t size = 1
        for d in self._shape:
            if d > 0:
                size *= d
        return size

    def actual_size(self):
        return len(self._data)

    def occupancy(self):
        return float(self.actual_size())/float(self.max_size())

    cdef double get(self, vector[uint64_t]& index):
        cdef uint64_t key = self.key(index)
        cdef SparseArrayIterator it = self._data.find(key)
        if (it != self._data.end()):
            return dereference(it).second;
        else:
            return 0.0

    cdef void set(self, vector[uint64_t]& index, double value):
        cdef uint64_t key = self.key(index)
        self._data[key] = value
        return

    cdef void add(self, vector[uint64_t]& index, double value):
        cdef uint64_t key = self.key(index)
        self._data[key] += value
        return

    cdef uint64_t key(self, vector[uint64_t]& index):
        self._check_bounds(index)
        cdef uint64_t key = 0
        cdef uint64_t b, s;
        for i in xrange(index.size()):
            b = index[i]
            s = self._dimscale[i]
            key += b*s
        return key

    @cython.cdivision(True)
    cdef vector[uint64_t] decodekey(self, uint64_t key):
        #cdef uint64_t startkey = key
        cdef vector[uint64_t] index;
        index.reserve(self._shape.size())
        cdef uint64_t b, s;
        for i in xrange(self._shape.size()):
            s = self._shape.at(i)
            if s > 0:
                b = key % s
                key = (key - b) / s
            else:
                b = 0
            index.push_back(b)
        #assert self.key(index) == startkey
        return index

    cdef _check_bounds(self, vector[uint64_t]& index):
        cdef uint64_t u, s
        for i in xrange(len(index)):
            u = index[i]
            s = self._shape[i]
            if s > 0 and not (0 <= u < s):
                raise IndexError("%.0fth index out of bounds" % i, self._shape, index,)
        return

    @cython.profile(PROFILE_FLAG)
    def project(self, vector[uint64_t] keep, range_=None):
        newshape = [self._shape[k] for k in keep]
        cdef SparseArray result = SparseArray(newshape)
        #for key, value in self._data:
        #    index = self.decodekey(key)
        #    if range_ is None or self._within_range(index, range_):
        #        index = [index[k] for k in keep]
        #        result.add(index, value)
        cdef SparseArrayIterator it = self._data.begin()
        cdef SparseArrayIterator end = self._data.end()
        cdef uint64_t key;
        cdef vector[uint64_t] index
        cdef vector[uint64_t] newindex = newshape
        cdef uint64_t nkeep = keep.size()
        while it != end:
            key = dereference(it).first
            value = dereference(it).second
            index = self.decodekey(key)
            if range_ is None or self._within_range(index, range_):
                for ii in xrange(nkeep):
                    newindex[ii] = index[keep[ii]]
                result.add(newindex, value)
            preincrement(it)
        return result

    def flatten(self):
        result = numpy.zeros(self.max_size())
        for k, v in self._data:
            result[k] = v
        return result

    def _within_range(self, vector[uint64_t] index, dict range_):
        for k, v in range_.iteritems():
            if not v[0] <= index[k] < v[1]:
                return False
        return  True        

    @cython.profile(PROFILE_FLAG)
    cdef int _check_shape(self, rhs):
        cdef vector[uint64_t] ls = self.shape()
        cdef vector[uint64_t] rs = rhs.shape()
        if vector_content_identical(ls, rs):
            return SHAPE_IS_IDENTICAL
        if not _compatible_shape(ls, rs):
            raise Exception("Arrays have incompatible shape")
        return SHAPE_IS_COMPATIBLE

    @cython.profile(PROFILE_FLAG)
    def __mul__(lhs, rhs):
        # implement result = lhs * rhs
        if isinstance(lhs, float):
            lhs = _makescalar(lhs, rhs.shape())
        #convert to SparseArray
        cdef SparseArray l = lhs
        cdef SparseArray r = rhs
        cdef int mode = l._check_shape(r)
        if mode == SHAPE_IS_IDENTICAL:
            return _multiply_identical_shape_array_with_copy(l, r)
        else:
            return _multiply_array_with_copy(lhs, rhs)

    @cython.profile(PROFILE_FLAG)
    def __add__(SparseArray lhs, SparseArray rhs):
        if not lhs.shape() == rhs.shape():
            raise Exception("cannot __add__, incompatible shape.")
        return _add_array_with_copy(lhs, rhs)

    @cython.profile(PROFILE_FLAG)
    def __sub__(SparseArray lhs, SparseArray rhs):
        if not lhs.shape() == rhs.shape():
            raise Exception("cannot __add__, incompatible shape.")
        return _subtract_array_with_copy(lhs, rhs)

    @cython.profile(PROFILE_FLAG)
    def __imul__(SparseArray self, SparseArray rhs):
        # implement: lhs *= rhs
        cdef int mode = rhs._check_shape(self) # intentionally the opposite order to __mul__
        if mode == SHAPE_IS_IDENTICAL:
            return _multiply_identical_shape_array_inplace(self, rhs)
        else:
            return _multiply_array_inplace(self, rhs)

    @cython.profile(PROFILE_FLAG)
    def __div__(SparseArray self, SparseArray rhs):
        # implement result = lhs / rhs
        self._check_shape(rhs)
        return self._divide_array_with_copy(rhs)

    @cython.profile(PROFILE_FLAG)
    def __idiv__(SparseArray self, SparseArray rhs):
        # implement: lhs /= rhs
        rhs._check_shape(self) # intentionally the opposite order to __mul__
        return self._divide_array_inplace(rhs)
        
    def _divide_array_with_copy(self, SparseArray rhs):
        result = SparseArray(rhs.shape())
        cdef double x = 0.0
        for index, value in rhs:
            if value != 0:
                x = self[index] / value
            else:
                x = 0
            result.set(index, x)
        return result

    def _divide_array_inplace(self, SparseArray rhs):
        cdef double x = 0.0
        for index, value in self:
            if value != 0:
                x = rhs[index] / value
            else:
                x = 0.0
            self.set(index, x)
        return self

    def __reduce__(self):
        dictit = iter(self)
        constructor = SparseArray
        args = (list(self._shape),)
        return (constructor, args, None, None, dictit)

    def clone(self):
        ret = SparseArray(self._shape)
        cdef SparseArrayIterator it = self._data.begin()
        cdef SparseArrayIterator end = self._data.end()
        cdef uint64_t key
        cdef double value
        while it != end:
            key = dereference(it).first
            value = dereference(it).second
            ret._data[key] = value
            preincrement(it)
        return ret

@cython.profile(PROFILE_FLAG)
cdef SparseArray _multiply_array_with_copy(SparseArray lhs, SparseArray rhs):
        cdef SparseArray result = SparseArray(rhs.shape())
        cdef SparseArrayIterator it = rhs._data.begin()
        cdef SparseArrayIterator end = rhs._data.end()
        cdef double x = 0.0
        cdef vector[uint64_t] index
        while it != end:
            key = dereference(it).first
            index = rhs.decodekey(key)
            x = lhs.get(index) * dereference(it).second
            result.set(index, x)
            preincrement(it)
        return result

@cython.profile(PROFILE_FLAG)
cdef SparseArray _multiply_identical_shape_array_with_copy(SparseArray lhs, SparseArray rhs):
        cdef SparseArray result = SparseArray(rhs.shape())
        cdef SparseArrayIterator it = rhs._data.begin()
        cdef SparseArrayIterator end = rhs._data.end()
        cdef double x = 0.0
        cdef vector[uint64_t] index
        while it != end:
            key = dereference(it).first
            x = lhs._data[key] * dereference(it).second
            result._data[key] = x
            preincrement(it)
        return result

@cython.profile(PROFILE_FLAG)
cdef SparseArray _multiply_array_inplace(SparseArray lhs, SparseArray rhs):
        cdef SparseArrayIterator it = lhs._data.begin()
        cdef SparseArrayIterator end = lhs._data.end()
        cdef uint64_t key;
        cdef double x = 0.0
        cdef vector[uint64_t] index
        while it != end:
            key = dereference(it).first
            index = rhs.decodekey(key)
            x = rhs.get(index) * dereference(it).second
            lhs.set(index, x)
            preincrement(it)
        return lhs

@cython.profile(PROFILE_FLAG)
cdef SparseArray _multiply_identical_shape_array_inplace(SparseArray lhs, SparseArray rhs):
        cdef SparseArrayIterator it = lhs._data.begin()
        cdef SparseArrayIterator end = lhs._data.end()
        cdef uint64_t key;
        cdef double x = 0.0
        while it != end:
            key = dereference(it).first
            x = rhs._data[key] * dereference(it).second
            lhs._data[key] = x
            preincrement(it)
        return lhs

@cython.profile(PROFILE_FLAG)
cdef SparseArray _add_array_with_copy(SparseArray lhs, SparseArray rhs):
        cdef SparseArray result = SparseArray(rhs.shape())
        cdef SparseArrayIterator it = rhs._data.begin()
        cdef SparseArrayIterator end = rhs._data.end()
        cdef double x = 0.0
        cdef vector[uint64_t] index
        while it != end:
            key = dereference(it).first
            index = rhs.decodekey(key)
            result.add(index, rhs.get(index))
            preincrement(it)
        it = lhs._data.begin()
        end = lhs._data.end()
        while it != end:
            key = dereference(it).first
            index = lhs.decodekey(key)
            result.add(index, lhs.get(index))
            preincrement(it)
        return result

@cython.profile(PROFILE_FLAG)
cdef SparseArray _subtract_array_with_copy(SparseArray lhs, SparseArray rhs):
        cdef SparseArray result = SparseArray(rhs.shape())
        cdef SparseArrayIterator it = rhs._data.begin()
        cdef SparseArrayIterator end = rhs._data.end()
        cdef double x = 0.0
        cdef vector[uint64_t] index
        while it != end:
            key = dereference(it).first
            index = rhs.decodekey(key)
            result.add(index, -1.0 * rhs.get(index))
            preincrement(it)
        it = lhs._data.begin()
        end = lhs._data.end()
        while it != end:
            key = dereference(it).first
            index = lhs.decodekey(key)
            result.add(index, lhs.get(index))
            preincrement(it)
        return result

def _makescalar(val, shape):
    s = [0 for s in shape]
    arr = SparseArray(s)
    arr.set(s, val)
    return arr

###############################################################################

#cdef class Ones(SparseArray):
#
#    #cdef vector[uint64_t] _shape
#    #cdef vector[uint64_t] _dimscale
#    #cdef std_map[uint64_t, double] _data
#
#    def __cinit__(self, vector[uint64_t] shape, int minsize=10**7):
#        super(Ones, self).__init__(shape, minsize)
#
#    cdef double get(self, vector[uint64_t]& index):
#        cdef uint64_t key = self.key(index)
#        cdef SparseArrayIterator it = self._data.find(key)
#        if it == self._data.end():
#            return 1.0
#        else:
#            return dereference(it).second

###############################################################################
        
cdef class SparseHistogram:
    cdef SparseArray _arr
    cdef vector[vector[double]] _binning
    cdef double _overflow
    cdef _label

    def __init__(self, binning, label=None):
        #convert input label into HistogramNDLabel object
        if label is None:
            label = HistogramNDLabel(binning)
        elif isinstance(label, basestring):
            label = HistogramNDLabel(binning, label=label)
        #set-up member variables.
        cdef vector[uint64_t] shape
        for b in binning:
            shape.push_back(len(b) - 1)
        self._arr = SparseArray(shape)
        self._binning = binning
        self._overflow = 0.0
        self._label = label

    def binning(self):
        return self._binning

    def array(self):
        return self._arr

    def __str__(self):
        return "\n".join([r"SparseHistogram(%.2e%%, %.2e/%.2e)" % (100.*self.occupancy(), self.actual_size(), self.max_size()), 
                          str(self._arr),
                      ])
        

    def __len__(self):
        return len(self._arr)

    def max_size(self):
        return self._arr.max_size()

    def actual_size(self):
        return self._arr.actual_size()

    def occupancy(self):
        return self._arr.occupancy()

    def fill(self, coord, double weight=1.0):
        cdef vector[uint64_t] index = self._findindex(coord)
        if self._isoverflow(index):
            self._overflow += weight
        else:
            self._arr.add(index, weight)
        return

    def eval(self, coord):
        index = self._findindex(coord)
        return self._arr.get(index)

    cdef vector[uint64_t] _findindex(self, vector[double]& coord):
        cdef vector[uint64_t] index
        cdef uint64_t dim
        cdef int i
        cdef double x
        for dim in xrange(coord.size()):
            x = coord[dim]
            #i = bisect_right(self._binning[dim], x) - 1
            i = array_bisect_right(self._binning[dim], x) - 1
            #assert i == (bisect_right(self._binning[dim], x) - 1)
            #if i < 0:
            #    i = self._arr._shape[dim]
            if i < 0:
                #print "DEBUG underflow ", dim, i, x, self._binning[dim]
                i = 0
            if i >= self._arr._shape[dim]:
                #print "DEBUG overflow ", dim, i, x, self._binning[dim]
                i = self._arr._shape[dim] - 1
            index.push_back(i)
        return index

    cdef int _isoverflow(self, vector[uint64_t]& index):
        for ii in xrange(index.size()):
            if index[ii] < 0:
                return True
            if index[ii] >= self._arr._shape[ii]:
                return True
        return False
            

    def project_nd(self, list axes, dict range_=None):
        if range_:
            range_ = self._convert_floatrange_to_binrange(range_)
        result = SparseHistogram([self._binning[k] for k in axes])
        result._arr = self._arr.project(axes, range_=range_)
        return result

    def project_1d(self, unsigned int axis, dict range_=None):
        if range_:
            range_ = self._convert_floatrange_to_binrange(range_)
        binning, values = self._binning.at(axis), self._arr.project((axis,), range_=range_).flatten()
        label = HistogramNDLabel(binning, label = self._label.label,
                                 axislabels=self._label.axislabels[axis],
                                 axisunits=self._label.axisunits[axis],
                                 binlabels=self._label.binlabels[axis],
        )
        return HistogramND(binning, values, label)

    def _convert_floatrange_to_binrange(self, dict range_):
        result = {}
        for axis, v in range_.iteritems():
            low = array_bisect_right(self._binning.at(axis), v[0]) - 1
            high = array_bisect_right(self._binning.at(axis), v[1]) - 1
            result[axis] = (low, high)
        return result

    def split_1d(self, unsigned int axis, unsigned int splitaxis, fmtlabel=None):
        histlist = []
        for ii in xrange(len(self._binning[splitaxis]) - 1):
            #create histogram from sub range
            low = self._binning[splitaxis][ii]
            high = self._binning[splitaxis][ii+1]
            range_ = {splitaxis:(low, high)}
            h = self.project_1d(axis, range_=range_)
            #determine label for histogram
            splitname = self._label.axislabels[splitaxis]
            splitunits = self._label.axisunits[splitaxis]
            binlabel = self._label.getbinlabel(ii, splitaxis)
            if fmtlabel:
                fmtlabel(low, high, splitname, splitunits, binlabel)
                h.label.label = newlabel
            elif binlabel:
                h.label.label = binlabel
            elif splitname:
                newlabel = "%.2g<%s<%.2g" % (low, splitname, high)
                h.label.label = newlabel
            histlist.append(h)
        return histlist

    def __reduce__(self):
        constructor = _unpickle_sparsehistogram
        args = (self._binning, self._arr, self._overflow)
        return (constructor, args, None, None, None)

    def scale(self, float scale):
        cdef SparseArray rhs = self._arr
        cdef SparseArray result = SparseArray(rhs.shape())
        cdef SparseArrayIterator it = rhs._data.begin()
        cdef SparseArrayIterator end = rhs._data.end()
        cdef uint64_t key
        while it != end:
            key = dereference(it).first
            result._data[key] = dereference(it).second * scale
            preincrement(it)
        self._arr = result
        return

    def clone(self):
        ret = SparseHistogram(self._binning, label=self._label)
        ret._overflow = self._overflow
        ret._arr = self._arr.clone()
        return ret

def _unpickle_sparsehistogram(binning, arr, overflow):
    hist = SparseHistogram(binning)
    hist._arr = arr
    hist._overflow = overflow
    return hist


cdef int array_bisect_right(vector[double]& arr, double x):
    cdef int lo = 0
    cdef int hi = arr.size()
    cdef int mid
    while lo < hi:
        mid = (lo + hi) / 2
        if x < arr[mid]:
            hi = mid
        else:
            lo = mid + 1
    return lo

@cython.profile(PROFILE_FLAG)
cdef SparseArray sparse_array_interpolation(double f, SparseArray y0, SparseArray y1):
    if not y1._check_shape(y0) == SHAPE_IS_IDENTICAL:
        raise Exception("ERROR interpolating between arrays with different shapes.")
    cdef SparseArray result = SparseArray(y1.shape())
    cdef SparseArrayIterator it = y1._data.begin()
    cdef SparseArrayIterator end = y1._data.end()
    cdef double v1, v0, v
    cdef vector[uint64_t] index
    cdef uint64_t key
    while it != end:
        key = dereference(it).first
        v1 = dereference(it).second
        v0 = y0._data[key]
        v = f*v1 + (1.0-f)*v0
        result._data[key] = v
        preincrement(it)
    return result

cdef vector_content_identical(vector[uint64_t]& lhs, vector[uint64_t]& rhs):
    if lhs.size() != rhs.size():
        return False
    for ii in xrange(lhs.size()):
        if not lhs[ii] == rhs[ii]:
            return False
    return True

cdef bool _compatible_shape(vector[uint64_t]& ls, vector[uint64_t]& rs):
    if not len(ls) == len(rs):
        return False
    for i in xrange(len(ls)):
        if not (ls[i] == rs[i] or ls[i] == 0):
            return False
    return True
