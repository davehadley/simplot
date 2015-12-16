from unordered_map cimport unordered_map as std_map
from libcpp.vector cimport vector
from libcpp cimport bool
from libc.stdint cimport uint64_t

ctypedef std_map[uint64_t, double] SparseArrayContainer
ctypedef std_map[uint64_t, double].iterator SparseArrayIterator

cdef class SparseArray:
    cdef vector[uint64_t] _shape
    cdef vector[uint64_t] _dimscale
    #cdef unordered_map[uint64_t, double] _data
    cdef SparseArrayContainer _data

    cdef uint64_t key(self, vector[uint64_t]& index);
    cdef vector[uint64_t] decodekey(self, uint64_t key);
    cdef void set(self, vector[uint64_t]& index, double value);
    cdef double get(self, vector[uint64_t]& index);
    cdef void add(self, vector[uint64_t]& index, double value);

    cdef _check_bounds(self, vector[uint64_t]& index);
    cdef int _check_shape(self, rhs);

    #cdef SparseArray _multiply_array_with_copy(self, SparseArray rhs);
    #cdef SparseArray _multiply_array_inplace(self, SparseArray rhs);

#cdef class Ones(SparseArray):
#    pass

cdef int array_bisect_right(vector[double]& arr, double x);
cdef SparseArray sparse_array_interpolation(double f, SparseArray y0, SparseArray y1);

#cdef extern from "<algorithm>" namespace "std" nogil:
#    Iter lower_bound[Iter, T](Iter first, Iter last, const T& value)
