import numpy as np

from clikelihood import gaus_log_density, poisson_log_density

################################################################################

class LikelihoodException(Exception):
    pass

class LikelihoodParametersMismatch(LikelihoodException):
    pass

################################################################################

class Likelihood(object):
    def __init__(self, parameter_names):
        self.parameter_names = parameter_names
        self._npars = len(self.parameter_names)
        if not len(self.parameter_names) == len(set(self.parameter_names)):
            duplicates = [(i,p) for p in enumerate(self.parameter_names) if self.parameter_names.count(p) > 1]
            raise LikelihoodParametersMismatch("duplicate parameter names", duplicates)

    def __call__(self, x):
        raise NotImplementedError("Sub-class must override this method.")

    def _checksize(self, x):
        if not len(x) == self._npars:
            raise LikelihoodParametersMismatch("called with wrong number of parameters num args(%s) != num pars(%s)" % (len(x), self._npars))
        

################################################################################

class Minus2LnLikehood(Likelihood):
    def __init__(self, f):
        super(Minus2LnLikehood, self).__init__(f.parameter_names)
        self._func = f
    def __call__(self, x):
        return -2.0*self._func(x)

class ConstantLikelihood(Likelihood):
    def __init__(self, parameter_names, value):
        super(ConstantLikelihood, self).__init__(parameter_names)
        self._value = float(value)
    def __call__(self, x):
        self._checksize(x)
        return self._value

################################################################################

class GaussianLikelihood(Likelihood):
    def __init__(self, parameter_names, mu, sigma, startindex=0):
        super(GaussianLikelihood, self).__init__(parameter_names)
        self._mu = np.array(mu, dtype=float, copy=True)
        self._sigma = np.array(sigma, dtype=float, copy=True)
        self._verify()

    def _verify(self):
        if not (len(self._mu.shape) == len(self._sigma.shape) == 1):
            raise ValueError("GaussianLikelihood mu/sigma have the wrong shape", self._mu.shape == self._sigma.shape)
        if not len(self._mu) == len(self._sigma) == self._npars:
            raise ValueError("GaussianLikelihood given mu and sigma of different length", self._mu, self._sigma)

    def __call__(self, x):
        self._checksize(x)
        mu = self._mu
        sigma = self._sigma
        return gaus_log_density(x, mu, sigma)

################################################################################

class MultiVariateGaussianLikelihood(Likelihood):
    def __init__(self, parameter_names, mu, cov):
        super(MultiVariateGaussianLikelihood, self).__init__(parameter_names)
        self._mu = np.array(mu, dtype=float, copy=True)
        cov = np.array(cov, copy=True)
        self._verify(cov) # check inputs before trying matrix inversion.
        self._invcov = np.linalg.inv(cov)


    def _verify(self, cov):
        if len(self._mu.shape) != 1:
            raise ValueError("MultiVariateGaussianLikelihood given mu with the wrong shape", self._mu.shape == cov.shape)
        if not len(self._mu) == len(cov) == self._npars:
            raise ValueError("MultiVariateGaussianLikelihood given mu and cov of different length", len(self._mu), len(cov), self._npars)
        #check covariance is square (not needed as np.linalg.inv should raise a LinAlgError)
        shape = cov.shape
        if len(shape)!=2 or shape[0] != shape[1]:
            raise ValueError("MultiVariateGaussianLikelihood given cov with the wrong shape", shape)
        
    def __call__(self, x):
        self._checksize(x)
        xminmu = x - self._mu
        xminmu = xminmu.reshape((1, self._npars))
        return -0.5 * np.dot(xminmu, np.dot(self._invcov, xminmu.T))[0,0]

################################################################################

class SumLikelihood(Likelihood):
    def __init__(self, funcs):
        parameter_names = self._checkparnames(funcs)
        super(SumLikelihood, self).__init__(parameter_names)
        self._funcs = funcs

    def __call__(self, x):
        self._checksize(x)
        return sum(f(x) for f in self._funcs)

    def _checkparnames(self, funcs):
        parameter_names = None
        for f in funcs:
            if parameter_names is None:
                parameter_names = f.parameter_names
            if not parameter_names == f.parameter_names:
                raise LikelihoodParametersMismatch("SumLikelihood given functions with mis-matched parameter names", parameter_names, f.parameter_names, str(f))
        return parameter_names

################################################################################

class CombinedLikelihood(Likelihood):
    def __init__(self, funcs):
        parameter_names = []
        self._funcs = []
        start = 0
        stop = 0
        for f in funcs:
            n = f.parameter_names
            parameter_names.extend(n)
            stop = start + len(n)
            self._funcs.append((f, start, stop))
            start = stop
        super(CombinedLikelihood, self).__init__(parameter_names)
    
    def __call__(self, x):
        self._checksize(x)
        return sum(f(x[start:stop]) for f, start, stop in self._funcs)

################################################################################

class EventRateLikelihood(Likelihood):
    def __init__(self, model, data):
        parameter_names = model.parameter_names
        self._model = model
        self._observed = np.copy(data)
        super(EventRateLikelihood, self).__init__(parameter_names)

    def __call__(self, x):
        self._checksize(x)
        observed = self._observed
        expected = self._model(x)
        return poisson_log_density(observed, expected)

################################################################################

class EventRateLikelihoodWithScale(Likelihood):
    def __init__(self, model, data, scale):
        parameter_names = model.parameter_names
        self._scale = float(scale)
        self._model = model
        self._observed = np.copy(data) * scale
        super(EventRateLikelihoodWithScale, self).__init__(parameter_names)

    def __call__(self, x):
        self._checksize(x)
        observed = self._observed #observed is already scaled
        expected = self._model(x) * self._scale
        return poisson_log_density(observed, expected)
