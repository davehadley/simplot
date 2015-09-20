import unittest

import numpy as np

import itertools

from simplot.mc.statistics import calculate_statistics, Mean, StandardDeviation, Covariance
from simplot.mc.likelihood import GaussianLikelihood, MultiVariateGaussianLikelihood
from simplot.fit.mcmc.metropolishastings import MetropolisHastingsAlgorithm
from simplot.fit.mcmc.proposalfunc import GaussianProposalFunction

class TestMcMc(unittest.TestCase):

    def setUp(self):
        self._num_dim = 2

    def test_gaussian(self):
        num_dim = self._num_dim
        mu = np.arange(0, num_dim, dtype=float)
        sigma = np.arange(1, num_dim + 1, dtype=float)
        parameter_names = ["par_" + str(m) for m in mu]
        lhd = GaussianLikelihood(parameter_names, mu, sigma)
        proposal = GaussianProposalFunction(sigma=sigma*0.4)
        start = np.copy(mu)
        mc = MetropolisHastingsAlgorithm(lhd, proposal, start, parameter_range=[(m - 10.0*s, m + 10.0*s) for m, s in zip(mu, sigma)])
        self.check_result(mc, mu, sigma)
        return

    def check_result(self, toymc, mu, sigma, expectedcov=None, npe=10**5, burnin=1000):
        if expectedcov is None:
            expectedcov = np.diag(np.power(sigma, 2.0))
        mean = Mean()
        stddev = StandardDeviation()
        cov = Covariance()
        statistics = [mean, stddev, cov]
        for _ in xrange(burnin):
            toymc() # burn in
        def func(toymc=toymc):
            ret = toymc()
            return ret[0]
        calculate_statistics(func, statistics, npe)
        # compare mean results
        for x, err, expected in zip(mean.eval(), mean.err(), mu):
            self.assertAlmostEquals(x, expected, delta=10.0*err) # use large window as error esitmate too small for MCMC
        # compare stddev results
        for x, err, expected in zip(stddev.eval(), stddev.err(), sigma):
            self.assertAlmostEquals(x, expected, delta=10.0*err)
        # compare covariance results
        covval = cov.eval()
        coverr = cov.err()
        for ii, jj in itertools.product(xrange(len(expectedcov)), repeat=2):
            self.assertAlmostEquals(covval[ii,jj], expectedcov[ii,jj], delta=10.0*coverr[ii,jj])
        return 

def main():
    unittest.main()
    return

if __name__ == "__main__":
    main()
