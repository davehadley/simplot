import itertools

import numpy as np

class EigenDecomposition:
    def __init__(self, A):
        A = np.matrix(A)
        eigenvalues, eigenvectors, diageigenvalues, q, invq = self._decompose(A)
        self.A = A
        self.eigenvectors = eigenvectors
        self.eigenvalues = eigenvalues
        self.diageigenvalues = diageigenvalues
        self._q = q
        self._invq = invq
        self._verify_decomposition()

    def transform_from_eigen_basis(self, x):
        x = np.matrix(x)
        r = self._q * x.T
        return np.squeeze(np.asarray(r)) # convert to flat array

    def transform_matrix_from_eigen_basis(self, L):
        L = np.matrix(L)
        return self._q * L * self._invq

    def transform_to_eigen_basis(self, x):
        x = np.matrix(x)
        r = self._invq * x.T
        return np.squeeze(np.asarray(r)) # convert to flat array

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
        Aprime = q * l * invq
        indices = [xrange(s) for s in A.shape]
        for index in itertools.product(*indices):
            ap = Aprime[index]
            a = A[index]
            if abs(a-ap) > precision:
                raise Exception("Error in eigen decomposition, cannot get back to input matrix")
        return 

