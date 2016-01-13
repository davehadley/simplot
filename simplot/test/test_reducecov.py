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

def main():
    return unittest.main()

if __name__ == "__main__":
    main()

