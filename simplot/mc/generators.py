# cython: profile=True
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
        raise NotImplementedError("ERROR: child class must implement _generate method.")
        
    def fixallexcept(self, varied):
        fixed = set(self.parameter_names)
        if varied is not None:
            for nf in varied:
                fixed.remove(nf)
        return self.setfixed(fixed)

    def setfixed(self, fixed):
        '''Fixed is either a list of parameter names, or a dictionary mapping parameter names to their fixed values.'''
        self._fixed = {} # reset fixed parameters
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
                raise ValueError("ERROR: generator has no parameter with name", k)
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
        raise NotImplementedError("ERROR: child class must implement getsigma method.")

    def _verify_generator(self):
        if not len(self.parameter_names) == len(self.start_values):
            raise ValueError("ERROR: different number of parameter names and starting values given", len(self.parameter_names), len(self.start_values))
        if not len(self.parameter_names) == len(set(self.parameter_names)):
            duplicates = [(ii, p) for ii, p in enumerate(self.parameter_names) if self.parameter_names.count(p) > 1]
            raise ValueError("ERROR: generator initialised with duplicate parameter names", duplicates)
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
        #return mu + sigma*self._rng.normal(size=len(mu))
        x = self._rng.normal(size=len(mu))
        np.add(mu, np.multiply(x, sigma, x), x)
        return x

    def getsigma(self, parname):
        index = self.parameter_names.index(parname)
        return self._sigma[index]
    
    def _verify(self):
        all_lists = [self.parameter_names, self._mu, self._sigma, self.start_values]
        if not all(len(v) == len(self._mu) for v in all_lists):
            raise ValueError("GaussianGenerator initialisation list arguments are not the same length",
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

class UniformGenerator(Generator):
    def __init__(self, parameter_names, value, range_, seed=None):
        super(UniformGenerator, self).__init__(parameter_names, start_values=value)
        self._scale = np.array([x[1]-x[0] for x in range_], dtype=float)
        self._shift = np.array([x[0] for x in range_], dtype=float)
        self._rng = np.random.RandomState(seed=seed)

    def _generate(self):
        scale = self._scale
        shift = self._shift
        x = self._rng.uniform(size=len(scale))
        return (scale*x) + shift

    def getsigma(self, parname):
        index = self.parameter_names.index(parname)
        return self._scale[index] / np.sqrt(12)
    
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
        self._mu = np.array(mu, copy=True, dtype=float)
        self._cov = np.array(cov, copy=True, dtype=float)
        self._sigma = np.sqrt(np.diag(self._cov), dtype=float)
        self._decomp = EigenDecomposition(self._cov)
        eigenvalues = self._removenegative(self._decomp.eigenvalues)
        self._eigensigma = np.sqrt(eigenvalues)
        self._transform = np.array(self._decomp._q)
        for arr in [self._mu, self._cov, self._sigma, self._eigensigma]:
            arr.setflags(write=False)
        self._rng = np.random.RandomState(seed=seed)
        self._verify()

    def _removenegative(self, eigenvalues):
        printwarning = True
        eigenvalues = np.array(self._decomp.eigenvalues)
        #remove negative eigenvalues
        for ii in xrange(len(eigenvalues)):
            if eigenvalues[ii] < 0.0:
                if printwarning:
                    print "WARNING: MultiVariateGaussianGenerator matrix has negative eigenvalues"
                    printwarning = False
                eigenvalues[ii] = 0.0
        return eigenvalues

    def _generate(self):
        mu = self._mu
        sigma = self._eigensigma
        x = self._rng.normal(size=len(mu))
        np.multiply(x, sigma, x)
        x = self._decomp.transform_from_eigen_basis(x)
        np.add(x, mu, x)
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
            raise ValueError("MultiVariateGaussianGenerator initialisation list arguments are not the same length",
                            all_lists,
                            )


