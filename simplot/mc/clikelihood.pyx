import numpy as np
cimport numpy as np
cimport cython

from libc.math cimport log

def gaus_log_density(np.ndarray[double, ndim=1] x, np.ndarray[double, ndim=1] mu, np.ndarray[double, ndim=1] sigma):
    cdef int N = len(mu)
    
    cdef double result = 0
    cdef double chi2
    cdef double s
    for i in xrange(N):
        s = sigma[i]
        if s > 0.0:
            chi2 = (x[i] - mu[i]) / sigma[i]
            expterm = 0.5 * chi2 * chi2
            result -= expterm
    return result

def poisson_log_density(np.ndarray[double, ndim=1] observed, np.ndarray[double, ndim=1] expected):
    return _poisson_log_density(observed, expected)

@cython.boundscheck(False)
cdef _poisson_log_density(np.ndarray[double, ndim=1] observed, np.ndarray[double, ndim=1] expected):
    cdef int N = observed.shape[0]
    cdef double result = 0
    cdef double l, n, e
    cdef int i
    for i in xrange(N):
        n = observed[i]
        e = expected[i]
        if n < 1e-10:
            l = e
        else:
            l = e - n + n*(safelog(n) + safenegativelog(e))
        result -= l
    return result

cdef double safelog(double x):
    cdef double result = 1e10;
    if x > 1e-6:
        result = log(x);
    return result;

cdef double safenegativelog(double x):
    cdef double result = 1e10;
    if x > 1e-6:
        result = -1.0*log(x);
    return result;

