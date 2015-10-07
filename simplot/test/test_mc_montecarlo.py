
import itertools
import math
import unittest

import numpy as np

from simplot.mc.generators import GaussianGenerator
#from simplot.mc.Likelihood import Minus2LnLikelihood, GaussianLikelihood
from simplot.mc.montecarlo import ToyMC, MonteCarloException, MonteCarloParameterMismatch
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

    def test_toymc_exceptions(self):
        with self.assertRaises(MonteCarloParameterMismatch):
            ToyMC(self.model, GaussianGenerator(list(reversed(self.model.parameter_names)), self.mu, self.sigma))
        with self.assertRaises(MonteCarloParameterMismatch):
            ToyMC(self.model, GaussianGenerator(["x", "y"], self.mu, self.sigma))
        with self.assertRaises(MonteCarloParameterMismatch):
            ToyMC(self.model, GaussianGenerator(self.model.parameter_names + ["z"], self.mu + [3.0], self.sigma + [3.0]))

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

#    def test_toymc_clone(self):
#        npe = 10**4
#        toymc = ToyMC(self.model, self.gen)
#        statistics1 = [Mean(), StandardDeviation(), FractionalStandardDeviation()]
#        calculate_statistics_from_toymc(toymc, statistics1, npe)
#        toymc = toymc.clone()
#        statistics2 = [Mean(), StandardDeviation(), FractionalStandardDeviation()]
#        calculate_statistics_from_toymc(toymc, statistics2, npe)
#        for s1, s2 in itertools.izip_longest(statistics1, statistics2):
#            val1, val2 = s1.eval(), s2.eval()
#            err1, err2 = s1.err(), s2.err()
#            for v1, v2, e1, e2 in itertools.izip_longest(val1, val2, err1, err2):
#                delta = 5.0 * np.sqrt(e1**2 + e2**2)
#                self.assertAlmostEquals(v1, v2, delta=delta)
#        return

    def test_toymcexperiment_clone(self):
        toymc = ToyMC(self.model, self.gen)
        exp1 = toymc()
        exp2 = exp1.clone()
        self.assertTrue(np.array_equal(exp1.vec, exp2.vec))
        self.assertTrue(np.array_equal(exp1.pars, exp2.pars))
        self.assertFalse(exp1.vec is exp2.vec)
        self.assertFalse(exp1.pars is exp2.pars)
        return

    def test_toymcexperiment_infostring(self):
        toymc = ToyMC(self.model, self.gen)
        str(toymc) # just check no exceptions are raised
        return

def main():
    unittest.main()

if __name__ == "__main__":
    main()


