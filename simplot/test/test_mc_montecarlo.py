
import itertools
import math
import unittest

import numpy as np

from simplot.mc.generators import GaussianGenerator
#from simplot.mc.Likelihood import Minus2LnLikelihood, GaussianLikelihood
from simplot.mc.montecarlo import ToyMC
from simplot.mc.statistics import calculate_statistics_from_toymc, Mean, StandardDeviation, FractionalStandardDeviation, Covariance

class TestMonteCarlo(unittest.TestCase):
    def setUp(self):
        names = ["a", "b"]
        def simplemodel(pars):
            return np.copy(pars)
        simplemodel.parameter_names = names
        self.model = simplemodel
        self.mu = [1.0, 2.0]
        self.sigma = [1.0, 2.0]
        self.gen = GaussianGenerator(names, self.mu, self.sigma, seed=12912)

    def test_toymc_stat(self):
        npe = 1000
        toymc = ToyMC(self.model, self.gen)
        for _ in xrange(npe):
            mc()
        return

    def test_toymc_stat(self):
        npe = 10**4
        toymc = ToyMC(self.model, self.gen)
        statistics = [Mean(), StandardDeviation(), FractionalStandardDeviation()]
        calculate_statistics_from_toymc(toymc, statistics, npe)
        expected = [self.mu, self.sigma, np.divide(self.sigma,self.mu), ]
        for s, exp in itertools.izip_longest(statistics, expected):
            val = s.eval()
            error = s.err()
            for v, er, ex in itertools.izip_longest(val, error, exp):
                self.assertAlmostEquals(v, ex, delta=3.0*er)
        return

    def test_toymc_covariance(self):
        npe = 10**4
        toymc = ToyMC(self.model, self.gen)
        expected = np.diag(np.power(self.sigma, 2))
        cov = calculate_statistics_from_toymc(toymc, Covariance(), npe)[0]
        cval = cov.eval()
        cerr = cov.err()
        for ii, jj in itertools.product(xrange(len(self.mu)), repeat=2):
            self.assertAlmostEquals(cval[ii, jj], expected[ii][jj], delta=3.0*cerr[ii, jj])
        return
        

def main():
    unittest.main()

if __name__ == "__main__":
    main()


