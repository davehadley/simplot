import itertools
import unittest
import numpy as np
from simplot.mc.reducecov import ReducedGaussian

################################################################################

class TestMarginalisation(unittest.TestCase):

    def test_diagonal(self):
        for Nrows in xrange(2, 10):
            mu = np.arange(Nrows)
            cov = np.diag(np.arange(Nrows) + 1)
            parameter_names = list(range(Nrows))
            self._check_result(parameter_names, mu, cov)
        return

    def test_non_diagonal(self):
        for Nrows in xrange(2, 10):
            rowindices = list(range(Nrows))
            mu = np.arange(Nrows)
            cov = np.diag(np.zeros(Nrows, dtype=float))
            for ii, jj in itertools.product(xrange(Nrows), repeat=2):
                cov[ii, jj] = (0.9)**(abs(ii-jj) + 1)
            parameter_names = range(Nrows)
            self._check_result(parameter_names, mu, cov)
        return

    def _check_result(self, parameter_names, mu, cov):
        """Maginalised matrix should just drop the rows and columns."""
        rowindices = list(range(len(mu)))
        for removed in itertools.permutations(rowindices, r=2):
            rg = ReducedGaussian(parameter_names, mu, cov, removed)
            for inew, jnew in itertools.product(xrange(len(rg.parameter_names)), repeat=2):
                iorig = parameter_names.index(rg.parameter_names[inew])
                jorig = parameter_names.index(rg.parameter_names[jnew])
                self.assertAlmostEquals(cov[iorig, jorig], rg.cov[inew, jnew])
                self.assertAlmostEquals(mu[iorig], rg.mu[inew])
                self.assertAlmostEquals(mu[jorig], rg.mu[jnew])
        return

################################################################################

class TestConditional(unittest.TestCase):

    def test_diagonal(self):
        #diagonal has no correlation and therefore conditional pdf is the same 
        #as the marginalised pdf
        for Nrows in xrange(2, 10):
            mu = np.arange(Nrows)
            cov = np.diag(np.arange(Nrows) + 1)
            parameter_names = list(range(Nrows))
            for removed in itertools.permutations(xrange(Nrows), r=2):
                rg = ReducedGaussian(parameter_names, mu, cov, marginalise=removed)
                expectedmu = [m for ii, m in enumerate(mu) if ii not in removed]
                expectedcov = np.diag([cov[ii][ii] for ii in xrange(Nrows) if ii not in removed])
                self._check_result(rg, expectedmu, expectedcov)
        return

    def test_bivariate(self):
        m1 = 1.0
        m2 = 2.0
        s1 = 2.0
        s2 = 4.0
        #some fake observation
        x1 = -1.0
        x2 = 6.0
        for cor in np.linspace(0.0, 0.99, num=10):
            mu = [m1, m2]
            cov = np.array([[s1**2, cor*s1*s2], [cor*s1*s2, s2**2]])
            #if we observed "b", parameter 2
            rg = ReducedGaussian(["a", "b"], mu, cov, conditional={"b":x2})
            expectedmu = [m1 + cor*(s1/s2)*(x2-m2)]
            expectedcov = np.array([[(1.0 - cor**2)*s1**2]])
            self._check_result(rg, expectedmu, expectedcov)
            #if we observe "a", parameter1
            rg = ReducedGaussian(["a", "b"], mu, cov, conditional={"a":x1})
            expectedmu = [m2 + cor*(s2/s1)*(x1-m1)]
            expectedcov = np.array([[(1.0 - cor**2)*s2**2]])
            self._check_result(rg, expectedmu, expectedcov)

        return

    def _check_result(self, reducedgaussian, expectedmu, expectedcov):
        rg = reducedgaussian
        for inew, jnew in itertools.product(xrange(len(rg.parameter_names)), repeat=2):
            self.assertAlmostEquals(expectedcov[inew, jnew], rg.cov[inew, jnew])
            self.assertAlmostEquals(expectedmu[inew], rg.mu[inew])
            self.assertAlmostEquals(expectedmu[jnew], rg.mu[jnew])
        return


################################################################################

def main():
    return unittest.main()

if __name__ == "__main__":
    main()

