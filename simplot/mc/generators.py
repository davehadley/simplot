
import math

import numpy as np
import scipy.stats

###############################################################################

class GeneratorException(Exception):
    pass

###############################################################################

class Generator(object):
    """Interface for generators. The generators must implement _generate()."""
    def __init__(self, parameter_names, start_values):
        self.parameter_names = list(parameter_names)
        self.start_values = np.array(start_values)
        if not len(self.parameter_names) == len(self.start_values):
            raise GeneratorException("ERROR: different number of parameter names and starting values given", len(self.parameter_names), len(self.start_values))
        self._fixed = {}

    def __call__(self):
        v = self._generate()
        #replace fixed parameters with their set values
        for index, val in self._fixed.iteritems():
            v[index] = val
        return v

    def _generate(self):
        raise GeneratorException("ERROR: child class must implement _generate method.")
        
    def fixallexcept(self, varied):
        fixed = set(self.parameter_names)
        if varied is not None:
            for nf in varied:
                fixed.remove(nf)
        return self.setfixed(fixed)

    def setfixed(self, fixed):
        '''Fixed is either a list of parameter names, or a dictionary mapping parameter names to their fixed values.'''
        if fixed is None:
            fixed = {}
        try:
            iteritems = fixed.iteritems()
        except AttributeError:
            fixed = dict.fromkeys(fixed, None)
            iteritems = fixed.iteritems()
        for k, v in iteritems:
            try:
                index = self.parameter_names.index(k)
            except ValueError:
                raise GeneratorException("ERROR: generator has no parameter with name", k)
            if v is None:
                v = self.start_values[index]
            self._fixed[index] = v
        return

    def getmu(self, parname):
        index = self.parameter_names.index(parname)
        return self.start_values[index]


###############################################################################

class GaussianGenerator(Generator):
    def __init__(self, parameter_names, mu, sigma, seed=1921):
        super(GaussianGenerator, self).__init__(parameter_names, start_values=mu)
        self._mu = np.array(mu)
        self._sigma = np.array(sigma)
        self._rng = scipy.stats.norm(loc=self._mu, scale=self._sigma)
        self._verify()

    def _generate(self):
        return self._rng.rvs()

    def getsigma(self, parname):
        index = self.parameter_names.index(parname)
        return self.sigma[index]
    
    def getcovariance(self, par1, par2):
        return 0.0
    
    def _verify(self):
        all_lists = [self.parameter_names, self._mu, self._sigma, self.start_values]
        if not all(len(v) == len(self._mu) for v in all_lists):
            raise Exception("GaussianGenerator initialisation list arguments are not the same length",
                            all_lists,
                            )

###############################################################################

class ConstantGenerator(Generator):
    def __init__(self, parameter_names, mu):
        super(ConstantGenerator, self).__init__(parameter_names, start_values=mu)

    def _generate(self):
        return np.copy(self.start_values)

    def getsigma(self, parname):
        index = self.parameter_names.index(parname)
        return 0.0
    
    def getcovariance(self, par1, par2):
        return 0.0

###############################################################################

class GeneratorList(Generator):
    def __init__(self, *args):
        self._generators = list(args)
        parameter_names = []
        start_values = []
        for gen in args:
            parameter_names.extend(gen.parameter_names)
            start_values.extend(gen.start_values)
        super(GeneratorList, self).__init__(parameter_names, start_values)

    def _generate(self):
        return np.concatenate([gen._generate() for gen in self._generators])

    def getmu(self, parname):
        result = None
        for gen in self._generators:
            if parname in g.parameter_names:
                result = g.getmu(parname)
                return result
        #could not find parameter
        raise ValueError("No known parameter", parname)

    def getsigma(self, parname):
        result = None
        for gen in self._generators:
            if parname in g.parameter_names:
                result = g.getsigma(parname)
                return result
        #could not find parameter
        raise ValueError("No known parameter", parname)

    def getcovariance(self, par1, par2):
        result = None
        for gen in self._generators:
            if par1 in g.parameter_names:
                result = g.getcovariance(par1, par2)
                return result
        #could not find parameter
        raise ValueError("No known parameter", parname)                

###############################################################################

class MultiVariateGaussianGenerator(Generator):
    def __init__(self, parameter_names, mu, cov, seed=1921):
        super(GaussianGenerator, self).__init__(parameter_names, start_values=mu)
        self._mu = np.array(mu)
        self._sigma = np.array(sigma)
        self._rng = scipy.stats.multivariate_normal(mean=self._mu, cov=cov)
        self._verify()

    def _generate(self):
        return self._rng.rvs()

    def getsigma(self, parname):
        index = self.parameter_names.index(parname)
        return self.sigma[index]
    
    def getcovariance(self, par1, par2):
        return 0.0
    
    def _verify(self):
        all_lists = [self._parameter_names, self._mu, self._cov]
        if not all(len(v) == len(mu) for v in all_lists):
            raise Exception("MultiVariateGaussianGenerator initialisation list arguments are not the same length",
                            all_lists,
                            )


