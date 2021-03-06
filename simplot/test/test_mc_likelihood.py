import itertools
import math
import unittest

import numpy as np

from simplot.mc.likelihood import Minus2LnLikehood, ConstantLikelihood, \
    MultiVariateGaussianLikelihood, GaussianLikelihood, SumLikelihood, \
    CombinedLikelihood, EventRateLikelihood, LikelihoodParametersMismatch, \
    LikelihoodException, Likelihood

class TestLikelihood(unittest.TestCase):

    def test_likelihood_exception(self):
        abstract = Likelihood(["a"])
        with self.assertRaises(NotImplementedError):
            abstract([0.0])
        return

    def test_gaussianlikelihood_exception(self):
        gaus = GaussianLikelihood(["a", "b"], [0.0, 1.0], [1.0, 1.0])
        with self.assertRaises(LikelihoodParametersMismatch):
            gaus([0.0])
        with self.assertRaises(ValueError):
            GaussianLikelihood(["a", "b"], [0.0, 1.0], [1.0, 1.0, 1.])
        with self.assertRaises(ValueError):
            GaussianLikelihood(["a"], [0.0, 1.0], [1.0, 1.0])
        with self.assertRaises(ValueError):
            GaussianLikelihood(["a"], [[0.0, 1.0]], [1.0, 1.0])
        return

    def test_multivargauslikelihood_exception(self):
        with self.assertRaises(ValueError):
            MultiVariateGaussianLikelihood(["a", "b"], [1.0, 1.0, 1.0], np.diag([1.0, 2.0]))
        with self.assertRaises(ValueError):
            MultiVariateGaussianLikelihood(["a", "b"], [1.0, 1.0], [[1.0, 1.0, 1.0],[1.0, 1.0, 1.0]])
        with self.assertRaises(ValueError):
            MultiVariateGaussianLikelihood(["a", "b"], [[0.0, 1.0], [0.0, 1.0]], [[1.0, 0.0],[0.0,1.0]])
        return

    def test_duplicate_parameters(self):
        parnames = ["par", "par"]
        mu = [1.0, 1.0]
        sigma = [1.0, 1.0]
        settings = [(ConstantLikelihood, (parnames, 1.0)),
                    (GaussianLikelihood, (parnames, mu, sigma)),
                    (MultiVariateGaussianLikelihood, (parnames, mu, np.diag(np.power(sigma,2)))),
        ]
        for f, args in settings:
            with self.assertRaises(LikelihoodParametersMismatch):
                f(*args)

    def test_minus2lnlikelihood(self):
        for y in range(-10, 10):
            f = ConstantLikelihood(["parameter"], y)
            lhd = Minus2LnLikehood(f)
            self.assertEquals(lhd([0.0]), -2 * y)
        return

    def test_gaussian(self):
        mu = range(-5, 5)
        parnames = [str(m) for m in mu]
        sigma = range(0, 10)
        lhd = Minus2LnLikehood(GaussianLikelihood(parnames, mu, sigma))
        pars = np.array(mu, dtype=float)
        for i in xrange(len(mu)):
            for nsigma in xrange(-5, 6):
                pars[i] = mu[i] + nsigma*sigma[i]
                l = lhd(pars)
                if sigma[i]:
                    expected = nsigma**2
                else:
                    expected = 0.0
                self.assertEquals(expected, l)
            pars[i] = mu[i]
        return

    def test_multivargaus_diag(self):
        mu = range(-5, 5)
        parnames = [str(m) for m in mu]
        sigma = np.arange(1, 11, dtype=float)
        cov = np.diag(np.power(sigma, 2))
        lhd = Minus2LnLikehood(MultiVariateGaussianLikelihood(parnames, mu, cov))
        pars = np.array(mu, dtype=float)
        for i in xrange(len(mu)):
            for nsigma in xrange(-5, 6):
                pars[i] = mu[i] + nsigma*sigma[i]
                l = lhd(pars)
                if sigma[i]:
                    expected = nsigma**2
                else:
                    expected = 0.0
                self.assertAlmostEquals(expected, l)
            pars[i] = mu[i]
        return

    def test_sum(self):
        mu = range(-5, 5)
        parnames = [str(m) for m in mu]
        sigma = range(0, 10)
        funcs = []
        for ii in xrange(5):
            funcs.append(Minus2LnLikehood(GaussianLikelihood(parnames, mu, sigma)))
        lhd = SumLikelihood(funcs)
        pars = np.array(mu, dtype=float)
        for i in xrange(len(mu)):
            for nsigma in xrange(-5, 6):
                pars[i] = mu[i] + nsigma*sigma[i]
                l = lhd(pars)
                if sigma[i]:
                    expected = len(funcs) * nsigma**2
                else:
                    expected = 0.0
                self.assertEquals(expected, l)
            pars[i] = mu[i]
        return

    def test_sum_throws(self):
        funcs = [ConstantLikelihood(["par1"], 0.0),
                 ConstantLikelihood(["par2"], 0.0),
                 ]
        with self.assertRaises(LikelihoodParametersMismatch):
            SumLikelihood(funcs)

    def test_combined(self):
        mu = range(-5, 5)
        parnames = [str(m) for m in mu]
        sigma = range(0, 10)
        funcs = []
        N = 2
        for ii in xrange(0, len(mu), N):
            funcs.append(Minus2LnLikehood(GaussianLikelihood(parnames[ii:ii+N], mu[ii:ii+N], sigma[ii:ii+N])))
        lhd = CombinedLikelihood(funcs)
        pars = np.array(mu, dtype=float)
        for i in xrange(len(mu)):
            for nsigma in xrange(-5, 6):
                pars[i] = mu[i] + nsigma*sigma[i]
                l = lhd(pars)
                if sigma[i]:
                    expected = nsigma**2
                else:
                    expected = 0.0
                self.assertEquals(expected, l)
            pars[i] = mu[i]
        return

    def test_eventratelikelihood(self):
        mu = np.linspace(1e9, 1e10) # use large numbers to be very close to gaussian
        names = ["par" + str(int(m)) for m in mu]
        def model(x):
            return np.array(x)
        model.parameter_names = names
        lhd = Minus2LnLikehood(EventRateLikelihood(model, mu))
        pars = np.array(mu, dtype=float)
        for i in xrange(len(mu)):
            for nsigma in xrange(-3, 4):
                pars[i] = mu[i] + nsigma*np.sqrt(mu[i])
                l = lhd(pars)
                expected = nsigma**2 # gaussian approximation (not exact)
                self.assertAlmostEquals(expected, l, delta=1e-3)
            pars[i] = mu[i]
        return



def main():
    unittest.main()

if __name__ == "__main__":
    main()

