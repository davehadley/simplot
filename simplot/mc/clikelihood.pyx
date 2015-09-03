import numpy as np
cimport numpy as np


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

