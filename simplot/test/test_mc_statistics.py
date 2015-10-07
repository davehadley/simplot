import itertools
import unittest

import numpy as np

from simplot.mc.generators import MultiVariateGaussianGenerator, GaussianGenerator
from simplot.mc.statistics import *
from simplot.nostdout import nostdout

class TestStatistics(unittest.TestCase):

    def test_safedivide(self):
        lhs = np.array([0.0, 2.0, 2.0, 1.0, 0.0, 1.0])
        rhs = np.array([1.0, 1.0, 2.0, 0.0, 0.0, 1.0])
        expected = np.array([0.0, 2.0, 1.0, 0.0, 0.0, 1.0])
        result = safedivide(lhs, rhs)
        self.assertTrue(expected.shape==result.shape)
        for x, exp in zip(result, expected):
            self.assertEquals(x, exp)
        return

    def test_mean_and_stddev(self):
        expectedmu = np.arange(0.0, 10.0)
        names = ["par"+str(x) for x in expectedmu]
        expectedsigma = np.arange(0.0, 10.0)
        mu = Mean()
        gen = GaussianGenerator(names, expectedmu, expectedsigma, seed=1290)
        statistics = [Mean(), StandardDeviation(), FractionalStandardDeviation()]
        with nostdout():
            calculate_statistics(gen, statistics, 10**4, name="test")
        for s, expected in itertools.izip_longest(statistics, [expectedmu, expectedsigma, safedivide(expectedsigma, expectedmu)]):
            for x, err, exp in itertools.izip_longest(s.eval(), s.err(), expected):
                self.assertAlmostEquals(x, exp, delta=5.0*err)
        return

    def test_covariance(self):
        expectedmu = np.array([2.0, 10.0])
        expectedcov = np.array([[1.0, -0.5], [-0.5, 1.0]])
        expectedfraccov = np.copy(expectedcov)
        for ii, jj  in itertools.product(xrange(len(expectedfraccov)), repeat=2):
            expectedfraccov[ii,jj] /= expectedmu[ii]*expectedmu[jj]
        names = ["a","b"]
        gen = MultiVariateGaussianGenerator(names, expectedmu, expectedcov, seed=1290)
        statistics = [Covariance(), Covariance(fractional=True)]
        calculate_statistics(gen, statistics, 10**4)
        for s, expected in itertools.izip_longest(statistics, [expectedcov, expectedfraccov]):
            x, err, exp = s.eval(), s.err(), expected
            for ii, jj in itertools.product(xrange(len(expectedfraccov)), repeat=2):
                self.assertAlmostEquals(x[ii,jj], exp[ii,jj], delta=5.0*err[ii,jj])
        return

    def test_roothistogram(self):
        names = ["a", "b"]
        expectedmu = np.array([2.0, 4.0])
        expectedsigma = np.array([1.0, 2.0])
        gen = GaussianGenerator(names, expectedmu, expectedsigma)
        #try 3 ways of constructing object to ensure all code paths are tested
        hist1 = RootHistogram("testhist", "testhist")
        hist2 = RootHistogram("testhist", "testhist", range=[-5.0, 5.0])
        hist3 = RootHistogram("testhist", "testhist", range=[(-5.0, 5.0), (-10., 10.)])
        mean = Mean()
        sigma = StandardDeviation()
        calculate_statistics(gen, [hist1, hist2, hist3, mean, sigma], 10**4)
        for h1, h2, h3, mu, sigma in itertools.izip_longest(hist1.eval(), hist2.eval(), hist2.eval(), mean.eval(), sigma.eval()):
            self.assertAlmostEqual(h1.GetMean(), mu)
            self.assertAlmostEqual(h1.GetRMS(), sigma)
            self.assertAlmostEqual(h2.GetMean(), mu)
            self.assertAlmostEqual(h2.GetRMS(), sigma)
            self.assertAlmostEqual(h3.GetMean(), mu)
            self.assertAlmostEqual(h3.GetRMS(), sigma)
        return

def main():
    unittest.main()
    return

if __name__ == "__main__":
    main()
