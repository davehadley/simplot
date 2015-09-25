import itertools

import numpy as np

class EigenDecomposition:
    def __init__(self, A, dtype=None):
        A = np.array(A, copy=True, dtype=dtype)
        eigenvalues, eigenvectors, inveigenvectors = self._decompose(A)
        self.A = A
        self.eigenvectors = np.array(eigenvectors, dtype=dtype)
        self.eigenvalues = np.array(eigenvalues, dtype=dtype)
        self.diageigenvalues = np.diag(self.eigenvalues)
        self.inveigenvectors = np.array(inveigenvectors, dtype=dtype)
        self._verify_decomposition()

    def transform_from_eigen_basis(self, x):
        r = x.copy().reshape((len(x), 1))
        r = np.dot(self.eigenvectors, r)
        r = r.ravel() # convert to flat array
        return r

    def transform_matrix_from_eigen_basis(self, L):
        return np.dot(self.eigenvectors, np.dot(L, self.inveigenvectors))

    def transform_to_eigen_basis(self, x):
        r = x.copy().reshape((len(x), 1))
        r = np.dot(self.inveigenvectors, r)
        r = r.ravel() # convert to flat array
        return r

    def _decompose(self, A):
        eigenvalues,  eigenvectors = np.linalg.eig(A)
        inveigenvectors = np.linalg.inv(eigenvectors)
        return eigenvalues, eigenvectors, inveigenvectors

    def _verify_decomposition(self, precision=1.e-6):
        A = self.A
        q = self.eigenvectors
        invq = self.inveigenvectors
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

