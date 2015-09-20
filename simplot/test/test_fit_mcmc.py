import unittest

import numpy as np

from simplot.mc.statistics import calculate_statistics, Mean, StandardDeviation, Covariance
from simplot.mc.likelihood import GaussianLikelihood, MultiVariateGaussianLikelihood
from simplot.fit.mcmc.metropolishastings import MetropolisHastingsMC, MetropolisHastingsAlgorithm, ListDataset
from simplot.fit.mcmc.proposalfunc import GaussianProposalFunction

class TestMcMc(unittest.TestCase):
    def test_gaussian(self):
        mu = np.arange(0, num_dim)
        sigma = np.arange(1, num_dim + 1)
        parameter_names = ["par_" + m for m in mu]
        lhd = GaussianLikelihood(parameter_names, mu, sigma)
        proposal = GaussianProposalFunction(sigma=sigma)
        start = np.copy(mu)
        alg = MetropolisHastingsAlgorithm(lhd, proposal, start, parameter_range=[(m - 10.0*s, m + 10.0*s) for m, s in zip(mu, sigma)])
        dataset = ListDataset()
        mc = MetropolisHastingsMC(alg, dataset)
        self.check_result(mc, mu, sigma)
        return

    def check_result(toymc, mu, sigma, expectedcov, npe=10**5):
        mean = Mean()
        stddev = StandardDeviation()
        cov = Covariance()
        statistics = [mean, stddev, cov]
        for _ in toymc:
            toymc() # burn in
        calculate_statistics(mc, statistics, npe)
        # compare mean results
        for x, err, expected in zip(mean.eval(), mean.err(), mu):
            self.assertAlmostEquals(x, expected, delta=5.0*err)
        # compare stddev results
        for x, err, expected in zip(stddev.eval(), stddev.err(), sigma):
            self.assertAlmostEquals(x, expected, delta=5.0*err)
        # compare covariance results
        covval = cov.eval()
        coverr = cov.err()
        for ii, jj in itertools.product(xrange(len(expectedcov)), repeat=2):
            self.assertAlmostEquals(covval[ii,jj], expectedcov[ii,jj], delta=5.0*coverr[ii,jj])
        return 

def main():
    unittest.main()
    return

if __name__ == "__main__":
    main()
