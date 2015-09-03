
import math

import numpy as np
import scipy.stats

from simplot.mc.eigendecomp import EigenDecomposition

###############################################################################

class GeneratorException(Exception):
    pass

###############################################################################

class Generator(object):
    """Interface for generators. The generators must implement _generate()."""
    def __init__(self, parameter_names, start_values):
        self.parameter_names = list(parameter_names)
        self.start_values = np.array(start_values, copy=True)
        self.start_values.setflags(write=False)
        self._fixed = {}
        self._verify_generator()

    def __call__(self):
        v = self._generate()
        #replace fixed parameters with their set values
        for index, val in self._fixed.iteritems():
            v[index] = val
        return v

    def _generate(self):
        raise NotImplemented("ERROR: child class must implement _generate method.")
        
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

    def getcovariance(self, par1, par2):
        if par1 == par2:
            return self.getsigma(par1)**2
        else:
            return 0.0

    def getsigma(self, parname):
        raise NotImplemented("ERROR: child class must implement getsigma method.")

    def _verify_generator(self):
        if not len(self.parameter_names) == len(self.start_values):
            raise GeneratorException("ERROR: different number of parameter names and starting values given", len(self.parameter_names), len(self.start_values))
        if not len(self.parameter_names) == len(set(self.parameter_names)):
            duplicates = [(ii, p) for ii, p in enumerate(self.parameter_names) if self.parameter_names.count(p) > 1]
            raise GeneratorException("ERROR: generator initialised with duplicate parameter names", duplicates)
        return


###############################################################################

class GaussianGenerator(Generator):
    def __init__(self, parameter_names, mu, sigma, seed=None):
        super(GaussianGenerator, self).__init__(parameter_names, start_values=mu)
        self._mu = np.array(mu, copy=True)
        self._sigma = np.array(sigma, copy=True)
        self._mu.setflags(write=False)
        self._sigma.setflags(write=False)
        self._rng = np.random.RandomState(seed=seed)
        self._verify()

    def _generate(self):
        mu = self._mu
        sigma = self._sigma
        return mu + sigma*self._rng.normal(size=len(mu))

    def getsigma(self, parname):
        index = self.parameter_names.index(parname)
        return self._sigma[index]
    
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
            if parname in gen.parameter_names:
                result = gen.getmu(parname)
                return result
        #could not find parameter
        raise ValueError("No known parameter", parname)

    def getsigma(self, parname):
        result = None
        for gen in self._generators:
            if parname in gen.parameter_names:
                result = gen.getsigma(parname)
                return result
        #could not find parameter
        raise ValueError("No known parameter", parname)

    def getcovariance(self, par1, par2):
        foundpar1 = False
        foundpar2 = False
        for gen in self._generators:
            if par1 in gen.parameter_names:
                foundpar1 = True
            if par2 in gen.parameter_names:
                foundpar2 = True
            if par1 in gen.parameter_names and par2 in gen.parameter_names:
                return gen.getcovariance(par1, par2)
        if foundpar1 and foundpar2:
            #parameters are not in the same generator, assume no correlation
            return 0.0
        #could not find parameter
        raise ValueError("No known parameter", par1, par2)


###############################################################################

class MultiVariateGaussianGenerator(Generator):
    def __init__(self, parameter_names, mu, cov, seed=None):
        super(MultiVariateGaussianGenerator, self).__init__(parameter_names, start_values=mu)
        self._mu = np.array(mu, copy=True)
        self._cov = np.array(cov, copy=True)
        self._sigma = np.sqrt(np.diag(self._cov))
        self._decomp = EigenDecomposition(self._cov)
        self._eigensigma = np.sqrt(self._decomp.eigenvalues)
        for arr in [self._mu, self._cov, self._sigma, self._eigensigma]:
            arr.setflags(write=False)
        self._rng = np.random.RandomState(seed=seed)
        self._verify()

    def _generate(self):
        mu = self._mu
        sigma = self._eigensigma
        x = self._rng.normal(size=len(mu))
        x = x * sigma
        x = self._decomp.transform_from_eigen_basis(x)
        x = x + mu
        return x

    def getsigma(self, parname):
        index = self.parameter_names.index(parname)
        return self._sigma[index]
    
    def getcovariance(self, par1, par2):
        i1 = self.parameter_names.index(par1)
        i2 = self.parameter_names.index(par2)
        cov = self._cov[i1,i2]
        return cov
    
    def _verify(self):
        all_lists = [self.parameter_names, self._mu, self._cov]
        if not all(len(v) == len(self._mu) for v in all_lists):
            raise Exception("MultiVariateGaussianGenerator initialisation list arguments are not the same length",
                            all_lists,
                            )


