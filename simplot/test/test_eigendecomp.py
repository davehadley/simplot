import unittest
from simplot.mc.eigendecomp import EigenDecomposition
from collections import namedtuple

import numpy as np

TestMatrix = namedtuple("TestMatrix", ["matrix", "eigenvalues", "eigenvectors"])

class TestEigenDecomposition(unittest.TestCase):

    def test_diagonal(self):
        m = TestMatrix(np.diag(np.arange(5) + 1), np.arange(5) + 1, np.diag(np.ones(5)))
        decomp = EigenDecomposition(m.matrix)
        self._compare_expected(m, decomp)

    def test_good(self):
        m = [[2.0, 1.0,], 
             [1.0, 2.0]
        ]
        egval = [3.0, 1.0]
        sqh = np.sqrt(0.5)
        egvec = [[sqh, -sqh], [sqh, sqh]]
        m = TestMatrix(m, egval, egvec)
        decomp = EigenDecomposition(m.matrix)
        self._compare_expected(m, decomp)

    def test_bad(self):
        #try some singular matrices, expect errors
        egval, egvec = [0.0, 0.0], [[1.0, 0.0], [0.0, 1.0]]
        m0 = [[0.0, 0.0], [0.0, 0.0]] # zero eigenvalues, should work
        m1 = [[0.0, 0.0], [0.0, 1.0]] # zero eigen-values, should work
        m2 = [[0.0, 0.0], [1.0, 0.0]] # not diagonalizable, should raise ValueError
        m3 = [[1.0, 1.0], [0.0, 1.0]] # not diagonalizable, should raise ValueError
        for m, err in zip([m0, m1, m2, m3], [None, None, ValueError, ValueError]):
            m = TestMatrix(m, np.diag(m), egvec)
            if err:
                with self.assertRaises(err):
                    decomp = EigenDecomposition(m.matrix)
            else:
                decomp = EigenDecomposition(m.matrix)
                self._compare_expected(m, decomp)
    
    def _compare_expected(self, m, decomp):
        self._compare_expected_values(m, decomp)
        self._compare_expected_transformation(m, decomp)
        return

    def _compare_expected_values(self, m, decomp):
        #check number of values
        self.assertEquals(len(m.eigenvalues), len(decomp.eigenvalues))
        #check all public members are as expected
        self.assertEquals(len(m.eigenvectors), len(decomp.eigenvectors))
        N = len(m.eigenvalues)
        for l, r in zip(m.eigenvalues, decomp.eigenvalues):
            self.assertEquals(l, r)
        for l, r in zip(m.eigenvalues, decomp.eigenvalues):
            self.assertEquals(l, r)
        self._comparray(decomp.diageigenvalues, np.diag(m.eigenvalues))

    def _compare_expected_transformation(self, m, decomp):
        N = len(m.eigenvalues)
        #check transformation methods
        for ii, (egval, egvec) in enumerate(zip(m.eigenvalues, m.eigenvectors)):
            for size in np.arange(1, 10, dtype=float):
                # create vector of size along ii-axis
                x = np.zeros(N)
                x[ii] = size
                #expect in the eigen basis the eigen vector * size
                y = decomp.transform_to_eigen_basis(x)
                yexp = np.multiply(size, egvec)
                #expect inverse operation to return original vector
                xexp = decomp.transform_from_eigen_basis(yexp)
                self._comparray(y, yexp)
                self._comparray(x, xexp)
        return

    def _comparray(self, a1, a2):
        success = np.allclose(a1, a2)
        msg  = ""
        if not success:
            msg = str(a1) + " != " + str(a2)
        self.assertTrue(success, msg=msg)
        return
                
def main():
    unittest.main()
    return

if __name__ == "__main__":
    main()
