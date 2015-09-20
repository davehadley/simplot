import unittest

import numpy as np

import itertools

from simplot.mc.statistics import calculate_statistics, Mean, StandardDeviation, Covariance
from simplot.mc.likelihood import GaussianLikelihood, MultiVariateGaussianLikelihood
from simplot.fit.mcmc.metropolishastings import MetropolisHastingsAlgorithm
from simplot.fit.mcmc.proposalfunc import GaussianProposalFunction
from simplot.fit.mcmc.io import write_mc_to_root_file

class TestMcMc(unittest.TestCase):

    def setUp(self):
        self._num_dim = 10

    def build_gaussian(self):
        num_dim = self._num_dim
        mu = np.arange(1, num_dim + 1, dtype=float)
        sigma = np.arange(1, num_dim + 1, dtype=float)
        parameter_names = ["par_" + str(int(m)) for m in mu]
        lhd = GaussianLikelihood(parameter_names, mu, sigma)
        proposal = GaussianProposalFunction(sigma=sigma*0.4)
        start = np.copy(mu)
        mc = MetropolisHastingsAlgorithm(lhd.parameter_names, lhd, proposal, start, parameter_range=[(m - 10.0*s, m + 10.0*s) for m, s in zip(mu, sigma)])
        return mc, mu, sigma

    def test_gaussian(self):
        mc, mu, sigma = self.build_gaussian()
        self.check_result(mc, mu, sigma)
        return

    def test_write(self):
        #build model
        mc, mu, sigma = self.build_gaussian()
        #try various options for writing to a file
        write_mc_to_root_file(mc, "test_fit_mcmc_gaussian_1.root", "mc", 10**5)
        write_mc_to_root_file(mc, "test_fit_mcmc_gaussian_2.root", "mc", None, seconds=10)
        write_mc_to_root_file(mc, "test_fit_mcmc_gaussian_3.root", "mc", 10**5, burnin=1000)
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
        # use large window as error esitmate too small for MCMC
        delta = 0.1
        # compare mean results
        for x, err, expected in zip(mean.eval(), mean.err(), mu):
            self.assertAlmostEquals(x, expected, delta=delta*expected) 
        # compare stddev results
        for x, err, expected in zip(stddev.eval(), stddev.err(), sigma):
            self.assertAlmostEquals(x, expected, delta=delta*expected)
        # compare covariance results
        covval = cov.eval()
        coverr = cov.err()
        for ii, jj in itertools.product(xrange(len(expectedcov)), repeat=2):
            self.assertAlmostEquals(covval[ii,jj], expectedcov[ii,jj], delta=delta*max(expectedcov[ii,ii], expectedcov[jj, jj]))
        return 

def main():
    unittest.main()
    #test = TestMcMc("test_gaussian")
    #test.setUp()
    #test.test_gaussian()
    return

if __name__ == "__main__":
    from simplot.profile import profile_func
    profile_func(main)
