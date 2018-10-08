import numpy
import itertools
import operator
import simplot.progress

from simplot.mc.eigendecomp import EigenDecomposition

class Verbosity:
    QUIET = 0
    PRINT_PROGRESS = 1

class Hesse:
    def __init__(self, func, xvec, delta=None, updatedelta=False, verbosity=Verbosity.QUIET, ignore_errors=False):
        self._func = func
        self._xvec = xvec
        if delta is None:
            delta = numpy.ones(shape=(len(xvec),), dtype=float)
        self._delta = delta
        self._matrix_d2dxdy = None
        self._error_matrix = None
        self._error_stddev = None
        self._updatedelta = updatedelta
        self._verbosity = verbosity
        self._ignore_errors = ignore_errors
        self._check_inputs()

    def _check_inputs(self):
        for ii, x in enumerate(self._delta):
            if x == 0.0:
                raise ValueError("Hesse has zero step size for dimension " + str(ii))
        return

    def run(self):
        self._cache_error_matrix()
        return

    def stddev(self):
        if self._error_stddev is None:
            m = self.covariance_matrix()
            diag = self._getdiagonal(m)
            self._error_stddev = numpy.sqrt(diag)
        return self._error_stddev

    def covariance_matrix(self):
        if self._error_matrix is None:
            self._cache_error_matrix()
        return self._error_matrix

    def _cache_error_matrix(self):
        matrix_d2dxdy, inverse = self._calculate_second_partial_derivative_matrix_and_inverse(self._func, self._xvec, self._delta)
        if self._updatedelta:
            #set delta to the diagonal of the matrix
            self._delta = self._getdiagonal(inverse)
            matrix_d2dxdy, inverse = self._calculate_second_partial_derivative_matrix_and_inverse(self._func, self._xvec, self._delta)
        self._matrix_d2dxdy = matrix_d2dxdy
        self._error_matrix = numpy.array(inverse, dtype=float)
        return


    def _calculate_second_partial_derivative_matrix_and_inverse(self, func, xvec, delta):
        matrix = self._calculate_second_partial_derivative_matrix(self._func, self._xvec, self._delta)
        #invert matrix
        try:
            inverse = -numpy.matrix(matrix).I
        except numpy.linalg.LinAlgError:
            if not self._ignore_errors:
                raise Exception("Unable to invert matrix", self._getdiagonal(matrix))
            else:
                matrix = self._fix_matrix(matrix)
                inverse = -numpy.matrix(matrix).I
        return matrix, inverse

    def _getdiagonal(self, m):
        return numpy.array([m[ii, ii] for ii in xrange(len(m))], dtype=float)

    def _iterdim(self, npars):
        for xdim in xrange(npars):
            for ydim in xrange(xdim, npars):
                yield xdim, ydim

    def _calculate_second_partial_derivative_matrix(self, func, pos, delta):
        npars = len(pos)
        result = numpy.zeros(shape=(npars, npars), dtype=float)
        if npars>1:
            numentries = self._n_choose_r(npars, 2) + npars
        else:
            numentries = 1
        iterdim = self._iterdim(npars)
        if self._verbosity >= Verbosity.PRINT_PROGRESS:
            iterdim = simplot.progress.printprogress("Hesse", numentries, iterdim)
        for xdim,ydim in iterdim:
                xdelta = delta[xdim]
                ydelta = delta[ydim]
                testvec = numpy.copy(pos)
                u = numpy.zeros(shape=(3, 3), dtype=float)
                for ii, jj in itertools.product(xrange(3), repeat=2):
                    #skip samples that are not used
                    if xdim == ydim:
                        if not ii == jj:
                            continue
                    else:
                        if ii == 1 or jj == 1:
                            continue
                    #sample likelihood
                    testvec[xdim] = pos[xdim] + (ii-1.0)*xdelta
                    testvec[ydim] = pos[ydim] + (jj-1.0)*ydelta
                    #print ii, jj, testvec
                    u[ii, jj] = func(testvec)
                if xdim == ydim:
                    seconddiff = (u[2][2] - 2.0*u[1][1] +u[0][0]) / xdelta**2
                else:
                    seconddiff = (u[2][2] - u[2][0] - u[0][2] + u[0][0]) / (4.0*xdelta*ydelta)
                #print u
                #print "estimated diff:", seconddiff
                result[xdim][ydim] = seconddiff
                #use symmetry of covariance matrix
                result[ydim][xdim] = result[xdim][ydim]
        return result

    def _n_choose_r(self, n, r):
        #taken from: http://stackoverflow.com/questions/4941753/is-there-a-math-ncr-function-in-python
        r = min(r, n-r)
        if r == 0:
            return 1
        numer = reduce(operator.mul, xrange(n, n - r, -1))
        denom = reduce(operator.mul, xrange(1, r + 1))
        return numer / denom

    def _fix_matrix(self, matrix):
        m = numpy.array(matrix)
        decomp = EigenDecomposition(m)
        A = numpy.copy(decomp.diageigenvalues)
        for ii in xrange(len(A)):
            if A[ii,ii] >= 0.0:
                A[ii,ii] = -1.0e-12 # force small positive eigenvalue
        return numpy.matrix(decomp.transform_matrix_from_eigen_basis(A))
