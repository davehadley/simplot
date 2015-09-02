import math
import unittest

import numpy as np

from simplot.mc.generators import GaussianGenerator, MultiVariateGaussianGenerator, ConstantGenerator, GeneratorList

class TestGenerators(unittest.TestCase):

    def setUp(self):
        self.mu = range(-5, 5)
        self.sigma = range(0, 10)
        self.names = ["par_"+str(m) for m in self.mu]


    def test_gaussian(self):
        mu = self.mu
        sigma = self.sigma
        names = self.names
        gen = GaussianGenerator(names, mu, sigma, seed=19021)
        npe = 10**4
        precision = 3.0*1.0/math.sqrt(npe)
        data = np.array([gen() for _ in xrange(npe)])
        self._checkmean(data)
        self._checkstddev(data)

    def test_gaussian_withfixed(self):
        mu = self.mu
        sigma = self.sigma
        names = self.names
        gen = GaussianGenerator(names, mu, sigma, seed=19021)
        fixed = {n:ii for ii, n in enumerate(self.names[-5:])}
        fixedmu = list(self.mu)
        fixedsigma = list(self.sigma)
        for n, v in fixed.iteritems():
            fixedmu[names.index(n)] = v
            fixedsigma[names.index(n)] = 0.0
        gen.setfixed(fixed)
        npe = 10**4
        precision = 3.0*1.0/math.sqrt(npe)
        data = np.array([gen() for _ in xrange(npe)])
        self._checkmean(data, mu=fixedmu, sigma=fixedsigma)
        self._checkstddev(data, sigma=fixedsigma)

    def _checkmean(self, data, mu=None, sigma=None, precision=None):
        if mu is None:
            mu = self.mu
        if sigma is None:
            sigma = self.sigma
        if precision is None:
            precision = 3.0*1.0/math.sqrt(len(data))
        #check mean
        for ii, m in enumerate(mu):
            self.assertAlmostEqual(np.mean(data[:, ii]), m, delta=sigma[ii]*precision)

    def _checkstddev(self, data, sigma=None, precision=None):
        if sigma is None:
            sigma = self.sigma
        if precision is None:
            precision = 3.0*1.0/math.sqrt(len(data))
        #check standard deviation
        for ii, s in enumerate(sigma):
            self.assertAlmostEqual(np.std(data[:, ii]), s, delta=s*precision*math.sqrt(0.5))

    def test_constant(self):
        mu = range(-5, 5)
        names = ["par_"+str(mu) for m in mu]
        gen = ConstantGenerator(names, mu)
        npe = 10**4
        data = np.array([gen() for _ in xrange(npe)])
        for ii, m in enumerate(mu):
            for jj in xrange(len(data)):
                self.assertEquals(data[jj,ii], mu[ii])
        return


def main():
    unittest.main()

if __name__ == "__main__":
    main()

            
