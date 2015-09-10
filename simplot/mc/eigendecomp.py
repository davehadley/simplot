import itertools

import numpy as np

class EigenDecomposition:
    def __init__(self, A):
        A = np.array(A, copy=True)
        eigenvalues, eigenvectors, diageigenvalues, q, invq = self._decompose(A)
        self.A = A
        self.eigenvectors = eigenvectors
        self.eigenvalues = eigenvalues
        self.diageigenvalues = diageigenvalues
        self._q = q
        self._invq = invq
        self._verify_decomposition()

    def transform_from_eigen_basis(self, x):
        r = x.copy().reshape((len(x), 1))
        r = np.dot(self._q, r)
        r = r.ravel() # convert to flat array
        return r

    def transform_matrix_from_eigen_basis(self, L):
        return np.dot(self._q, np.dot(L, self._invq))

    def transform_to_eigen_basis(self, x):
        r = x.copy().reshape((len(x), 1))
        r = np.dot(self._invq, r)
        r = r.ravel() # convert to flat array
        return r

    def _decompose(self, A):
        eigenvalues,  eigenvectors = np.linalg.eig(A)
        l = np.zeros(shape=A.shape)
        for ii in xrange(len(eigenvalues)):
            l[ii, ii] = eigenvalues[ii]
        q = eigenvectors
        invq = np.linalg.inv(q)
        return eigenvalues, eigenvectors, l, q, invq

    def _verify_decomposition(self, precision=1.e-6):
        A = self.A
        q = self._q
        invq = self._invq
        l = self.diageigenvalues
        #verify we can original matrix back from decomposition
        Aprime = np.dot(q, np.dot(l, invq))
        indices = [xrange(s) for s in A.shape]
        for index in itertools.product(*indices):
            ap = Aprime[index]
            a = A[index]
            if abs(a-ap) > precision:
                raise ValueError("Error in eigen decomposition, cannot get back to input matrix, ")
        return 

