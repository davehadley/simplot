import itertools
import unittest

import numpy as np

from simplot.mc.likelihood import GaussianLikelihood, MultiVariateGaussianLikelihood
from simplot.fit.hesse import Hesse, Verbosity

class TestHesse(unittest.TestCase):
    def test_gaussian(self):
        num_dim = 100
        #build gaussian likelihood
        mu = np.arange(num_dim, dtype=float)
        sigma = np.arange(1, num_dim + 1, dtype=float)
        parameter_names = ["par" + str(m) for m in mu]
        lhd = GaussianLikelihood(parameter_names, mu, sigma)
        self.check_hesse(lhd, mu, sigma)

    def check_hesse(self, lhd, mu, sigma, expectedcov=None):
        if expectedcov is None:
            expectedcov = np.diag(np.power(sigma, 2))
        h = Hesse(lhd, mu, delta=sigma, updatedelta=True)
        h.run()
        #check errors
        for o, e in itertools.izip_longest(h.stddev(), sigma):
            self.assertAlmostEquals(o, e)
        #check covariance
        cov = h.covariance_matrix()
        for ii, jj in itertools.product(xrange(len(mu)), repeat=2):
            expected = expectedcov[ii][jj]
            observed = cov[ii][jj]
            self.assertAlmostEquals(expected, observed)
        return

    def test_multivargaussian(self):
        num_dim = 100
        #build gaussian likelihood
        mu = np.arange(num_dim, dtype=float)
        cov = np.zeros(shape=(num_dim, num_dim), dtype=float)
        for ii, jj in itertools.product(xrange(num_dim), repeat=2):
            cov[ii,jj] = (0.5)**(abs(ii-jj)) * (ii + 1)*(jj + 1)
        sigma = np.sqrt(np.diag(cov))
        parameter_names = ["par" + str(m) for m in mu]
        lhd = MultiVariateGaussianLikelihood(parameter_names, mu, cov)
        self.check_hesse(lhd, mu, sigma, expectedcov=cov)
        return

def main():
    unittest.main()
    return

if __name__ == "__main__":
    main()
