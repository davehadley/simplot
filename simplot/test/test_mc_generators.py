import itertools
import math
import unittest

import numpy as np

from simplot.mc.generators import UniformGenerator, GaussianGenerator, MultiVariateGaussianGenerator, ConstantGenerator, GeneratorList, Generator

class TestGenerators(unittest.TestCase):

    def setUp(self):
        self.mu = np.arange(-5, 5, dtype=float)
        self.sigma = np.arange(0, 10, dtype=float)
        self.names = ["par_"+str(m) for m in self.mu]
        self.cov = np.zeros(shape=(len(self.mu), len(self.mu)))
        for ii, jj in itertools.product(xrange(len(self.mu)), repeat=2):
            cor = 0.9
            self.cov[ii,jj] = cor**(abs(ii-jj)) * self.sigma[ii]*self.sigma[jj]
        return

    def test_generatorexception(self):
        g = Generator(["a","b"], [0.0, 0.0])
        with self.assertRaises(NotImplementedError):
            g() # not implemented method
        with self.assertRaises(NotImplementedError):
            g.getsigma("a") # not implemented method
        with self.assertRaises(ValueError):
            g.setfixed({"x":0.0}) # fix non-existant parameter
        with self.assertRaises(ValueError):
            Generator(["a","b"], [0.0, 0.0, 0.0]) # wrong num parameters
        with self.assertRaises(ValueError):
            GaussianGenerator(["a","b"], [0.0, 0.0], [1.0, 1.0, 1.0]) # wrong num parameters (gaus)
        with self.assertRaises(ValueError):
            MultiVariateGaussianGenerator(["a","b"], [0.0, 0.0], [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]) # wrong num parameters (gaus)
        with self.assertRaises(ValueError):
            Generator(["a","a"], [0.0, 0.0]) # duplicate par name

    def test_generatorlistexception(self):
        gen1 = GaussianGenerator(["a", "b"], [1.0, 1.0], [1.0, 1.0])
        gen2 = GaussianGenerator(["c", "d"], [1.0, 1.0], [1.0, 1.0])
        lgen = GeneratorList(gen1, gen2)
        with self.assertRaises(ValueError):
            lgen.getmu("x")
        with self.assertRaises(ValueError):
            lgen.getsigma("x")
        with self.assertRaises(ValueError):
           lgen.getcovariance("a", "x")
        

    def test_fixedexcept_parameters(self):
        npe = 10**4
        mu = self.mu
        sigma = self.sigma
        names = self.names
        gen = GaussianGenerator(names, mu, sigma, seed=19021)
        for low, high in [(3, 7), (0, 10), (1, 3)]:
            fixedval = 0.5
            fixed = {n:0.5 for n in names[low:high]}
            gen.fixallexcept(fixed)
            data = np.array([gen() for _ in xrange(npe)])
            expectedmu = np.copy(mu)
            expectedsigma = np.copy(sigma)
            expectedsigma[:low] = 0.0
            expectedsigma[high:] = 0.0
            self._checkmean(data, mu=expectedmu)
            self._checkstddev(data, sigma=expectedsigma)
        return

    def test_fixed_parameters(self):
        npe = 10**4
        mu = self.mu
        sigma = self.sigma
        names = self.names
        gen = GaussianGenerator(names, mu, sigma, seed=19021)
        for low, high in [(3, 7), (0, 10), (1, 3)]:
            fixedval = 0.5
            fixed = {n:0.5 for n in names[low:high]}
            gen.setfixed(fixed)
            data = np.array([gen() for _ in xrange(npe)])
            expectedmu = np.copy(mu)
            expectedmu[low:high] = fixedval
            expectedsigma = np.copy(sigma)
            expectedsigma[low:high] = 0.0
            self._checkmean(data, mu=expectedmu)
            self._checkstddev(data, sigma=expectedsigma)
        #check switching if off we get original behaviour
        gen.setfixed(None)
        data = np.array([gen() for _ in xrange(npe)])
        self._checkmean(data, mu=mu)
        self._checkstddev(data, sigma=sigma)
        return


    def test_gaussian(self):
        mu = self.mu
        sigma = self.sigma
        names = self.names
        gen = GaussianGenerator(names, mu, sigma, seed=19021)
        npe = 10**4
        precision = 5.0*1.0/math.sqrt(npe)
        data = np.array([gen() for _ in xrange(npe)])
        self._checkmean(data, gen=gen)
        self._checkstddev(data, gen=gen)
        self._checkcovariance(data, cov=np.diag(np.power(sigma,2)), gen=gen)

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
        precision = 5.0*1.0/math.sqrt(npe)
        data = np.array([gen() for _ in xrange(npe)])
        self._checkmean(data, mu=fixedmu, sigma=fixedsigma, gen=gen)
        self._checkstddev(data, sigma=fixedsigma)
        self._checkcovariance(data, cov=np.diag(np.power(fixedsigma,2)))

    def test_constant(self):
        mu = range(-5, 5)
        names = ["par_"+str(m) for m in mu]
        gen = ConstantGenerator(names, mu)
        npe = 10**4
        data = np.array([gen() for _ in xrange(npe)])
        for ii, m in enumerate(mu):
            for jj in xrange(len(data)):
                self.assertEquals(data[jj,ii], mu[ii])
        sigma = np.zeros(len(mu))
        self._checkmean(data, mu=mu, sigma=sigma, gen=gen)
        self._checkstddev(data, sigma=sigma, gen=gen)
        self._checkcovariance(data, cov=np.diag(np.power(sigma,2)), gen=gen)
        return

    def test_uniform(self):
        mu = np.arange(0, 10, dtype=float)
        low = np.array([m - 5.0 for m in mu], dtype=float)
        high = np.array([m + 5.0 for m in mu], dtype=float)
        sigma = (high - low) / np.sqrt(12)
        names = ["par"+str(m) for m in mu]
        gen = UniformGenerator(names, mu, range_=zip(low, high))
        npe = 10**4
        data = np.array([gen() for _ in xrange(npe)])
        self._checkmean(data, mu=mu, sigma=sigma, gen=gen)
        self._checkstddev(data, sigma=sigma, gen=gen)
        self._checkcovariance(data, cov=np.diag(np.power(sigma,2)), gen=gen)
        return

    def test_multivargaus(self):
        mu = self.mu
        cov = self.cov
        names = self.names
        gen = MultiVariateGaussianGenerator(names, mu, cov, seed=19021)
        npe = 10**4
        #precision = 10.0*1.0/math.sqrt(npe)
        data = np.array([gen() for _ in xrange(npe)])
        self._checkmean(data, gen=gen)
        self._checkstddev(data, gen=gen)
        self._checkcovariance(data, gen=gen)
        return

    def test_generatorlist(self):
        #make a generator list from every kind of generator we have
        mu = self.mu
        sigma = self.sigma
        cov = self.cov
        names = self.names
        gnames = ["gaus"+n for n in names]
        cnames = ["const"+n for n in names]
        mnames = ["multigaus"+n for n in names]
        gengaus = GaussianGenerator(gnames, mu, sigma, seed=12021)
        genconst = ConstantGenerator(cnames, mu)
        genmulti = MultiVariateGaussianGenerator(mnames, mu, cov, seed=2314)
        gen = GeneratorList(gengaus, genconst, genmulti)
        #generate data
        npe = 10**4
        data = np.array([gen() for _ in xrange(npe)])
        #compare to expected
        exmu = np.concatenate([mu, mu, mu])
        exsigma = np.concatenate([sigma, np.zeros(len(mu)), sigma])
        excov = np.diag(np.power(exsigma, 2))
        excov[len(mu)*2:,len(mu)*2:] = cov
        #do checks
        self._checkmean(data, mu=exmu, sigma=exsigma, gen=gen)
        self._checkstddev(data, sigma=exsigma, gen=gen)
        self._checkcovariance(data, gen=gen, cov=excov)
        return

    def _checkmean(self, data, mu=None, sigma=None, precision=None, gen=None):
        if mu is None:
            mu = self.mu
        if sigma is None:
            sigma = self.sigma
        if precision is None:
            precision = 5.0*1.0/math.sqrt(len(data))
        #check mean
        for ii, m in enumerate(mu):
            dm = np.mean(data[:, ii])
            self.assertAlmostEqual(dm, m, delta=sigma[ii]*precision)
            if gen:
                gm = gen.getmu(gen.parameter_names[ii])
                self.assertAlmostEqual(dm, gm, delta=sigma[ii]*precision)
                self.assertEqual(gm, m)

    def _checkstddev(self, data, sigma=None, precision=None, gen=None):
        if sigma is None:
            sigma = self.sigma
        if precision is None:
            precision = 5.0*1.0/math.sqrt(len(data))
        #check standard deviation
        for ii, s in enumerate(sigma):
            ds = np.std(data[:, ii])
            self.assertAlmostEqual(ds, s, delta=s*precision*math.sqrt(0.5))
            if gen:
                gs = gen.getsigma(gen.parameter_names[ii])
                self.assertAlmostEqual(ds, gs, delta=s*precision*math.sqrt(0.5))
                self.assertEqual(gs, s)

    def _checkcovariance(self, data, cov=None, precision=None, gen=None):
        if cov is None:
            cov = self.cov
        if precision is None:
            precision = 5.0*1.0/math.sqrt(len(data))
        datacov = np.cov(data, rowvar=0)
        self.assertEquals(cov.shape, datacov.shape)
        for ii, jj in itertools.product(xrange(datacov.shape[0]), repeat=2):
            self.assertAlmostEqual(cov[ii,jj], datacov[ii][jj], delta=cov[ii,ii]*cov[jj,jj]*precision)
            if gen:
                gs = gen.getcovariance(gen.parameter_names[ii], gen.parameter_names[jj])
                self.assertAlmostEqual(gs, datacov[ii][jj], delta=cov[ii,ii]*cov[jj,jj]*precision)
                self.assertEqual(gs, cov[ii,jj])
        return

def main():
    unittest.main()

if __name__ == "__main__":
    main()

            
