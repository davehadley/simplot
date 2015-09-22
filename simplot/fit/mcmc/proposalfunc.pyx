# cython: profile=False

cimport cython
import numpy as np
cimport numpy as np
import StringIO
import itertools
import math
from simplot.fit.hesse import Hesse
from simplot.fit.hesse import Verbosity as HesseVerbosity
from simplot.fit.mcmc.metropolishastings import McMcSetupError
from simplot.mc.likelihood import MultiVariateGaussianLikelihood
from simplot.mc.generators import MultiVariateGaussianGenerator

from simplot.mc.clikelihood import gaus_log_density

from libc.math cimport sin, sqrt, log, M_PI, asin

###############################################################################

_RANDOM = np.random.RandomState(231321)

cdef class NormalGenerator:
    cdef object _random
    cdef np.ndarray _cache;
    cdef int _N
    cdef int _count
    def __cinit__(self):
        self._random = np.random.RandomState(231321)
        self._N = 1000
        self._gen()

    cdef _gen(self):
        self._cache = self._random.normal(size=self._N)
        self._count = 0
        return

    def normal(self):
        if self._count >= self._N:
            self._gen()
        cdef double r = self._cache[self._count]
        self._count += 1
        return r
        

_NORMAL = NormalGenerator()

cdef double _normal():
    return _NORMAL.normal()

###############################################################################

@cython.boundscheck(False)
cdef void _gaus_generate(np.ndarray[double, ndim=1] result, np.ndarray[double, ndim=1] mu, np.ndarray[double, ndim=1] sigma):
    cdef int N = result.shape[0]
    #cdef np.ndarray[double, ndim=1] normal = _RANDOM.normal(size=N);
    cdef int ii
    for ii in xrange(N):
        result[ii] = mu[ii] + sigma[ii] * _normal()
    return

cdef double _sinsqtheta2gaussian_log_density_i(x, mu, sigma):
    pass

cdef _sinsqtheta2gaussian_generate_i(x, mu, sigma):
    pass

cdef _sinsq2theta(theta):
    return np.power(np.sin(2.0*theta), 2);

cdef double _sinthetagaussian_log_density_i(x, mu, sigma):
        cdef int N = len(x)
        cdef double y, chi2
        cdef double jacobian = 0.0
        cdef double gaus = 0.0
        cdef double result = 0.0
        #compute jacobian
        for ii in xrange(N):
            y = sin(2.0*x[ii])
            y = y*y # take square
            #set limits on y
            if(y < 0.001):
                y = 0.001;
            elif(y > (1.0-0.001)):
                y = (1.0-0.001)
            chi2 = (y - mu[ii]) / sigma[ii]
            gaus -= 0.5 * chi2 * chi2
            jacobian += -log(4.0 * sqrt((1.0 - y)*y));
        result += jacobian
        return result;

cdef _sintheta(x):
   return np.sin(x)

cdef _sinthetagaussian_generate_i(result, mu, sigma):
    cdef int N = len(result)
    cdef np.ndarray[double, ndim=1] normal = _RANDOM.normal(size=N);
    cdef int ii
    for ii in xrange(N):
        result[ii] = mu[ii] + sigma[ii] * normal[ii]
    cdef double y
    for ii in xrange(N):
        #generate a y that is between -1, and 1
        y = -2.0;
        while (y < -1.0 or y > 1.0):
            n = _RANDOM.normal()
            y = mu[ii] + n * sigma[ii]
        y = asin(y);
        #get mirror solutions around pi
        if (_RANDOM.uniform() < 0.5):
            if (y > 0.0):
                y = M_PI - y;
        else:
            y = -M_PI - y;
        result[ii] = y
    return

###############################################################################

class GaussianProposalFunction(object):
    def __init__(self, sigma):
        self.sigma = np.array(sigma, dtype=float)
        self._result = np.zeros(len(self.sigma))
        self._sanityCheckState()

    def _sanityCheckState(self):
        if not np.all(np.isfinite(self.sigma)):
            raise McMcSetupError("sigma list contains NaN", self.infoString())
    
    def generate(self, parameters):
        x = self._result
        _gaus_generate(x, parameters, self.sigma)
        return x
    
    def logDensity(self, x, p):
        return gaus_log_density(x, p, self.sigma)

###############################################################################

class FixedGaussianProposalFunction(GaussianProposalFunction):
    def __init__(self, mu, sigma, parameterRanges=None):
        if parameterRanges is None:
            parameterRanges = [(m - 5.0*s, m + 5.0*s) for m,s in zip(mu, sigma)]
        super(FixedGaussianProposalFunction, self).__init__(parameterRanges=parameterRanges, stepSize=None, sigma=sigma)
        self.mu = np.array(mu, dtype=float)
    def generate(self, parameters):
        x = self._result
        _gaus_generate(x, self.mu, self.sigma)
        return x
    
    def logDensity(self, x, p):
        return gaus_log_density(x, self.mu, self.sigma)

###############################################################################

class GaussianTransformProposalFunction(object):
    """Same as GaussianProposalFunction but the variables can be Gaussian in a 
    different space (eg theta13 is Gaussian in sin^2 theta13)."""
    def __init__(self, parameterRanges, stepSize=0.1, sigma=None, func=None, invfunc=None, indomain=None):
        '''By default the step-size will stepSize*parameterRange.
           Alternatively, a specific sigma can be specified witht he sigma argument.
        
        parameterRanges a list of tuples with min and max values for each parameter.
        priorSigma a list containing the force sigma values for each parameter.
        '''
        self.parameterRanges = parameterRanges
        self.nParameters = len(parameterRanges)
        self.parNList = list(range(self.nParameters))
        if sigma is None:
            sigma = [None] * self.nParameters
        if func is None:
            func = [None] * self.nParameters
        if invfunc is None:
            invfunc = [None] * self.nParameters
        if indomain is None:
            indomain = [None] * self.nParameters
        self.sigma = np.array(sigma, dtype=float)
        self.func = func
        self.invfunc = invfunc
        self.indomain = indomain
        self.stepSize = stepSize
        self._result = np.zeros(self.nParameters)
        self._sanityCheckInput()
        self._autoCalcWidthParameters()

    def _sanityCheckInput(self):
        if not len(self.sigma) == self.nParameters:
            raise McMcSetupError("sigma list must have an entry for each parameter (")
    
    def _autoCalcWidthParameters(self):
        #automatically calculate width parameters
        for i,(vMin,vMax) in enumerate(self.parameterRanges):
            if self.sigma[i] is None:
                s = self.stepSize*(vMax-vMin)
                self.sigma[i] = s
        return
    
    def generate(self, parameters):
        x = self._result
        #convert parameters values to the space in which variables are Gaussian
        p = self._applytransformcopy(parameters)
        while True:
            _gaus_generate(x, p, self.sigma)
            if self._checkdomain(x):
                #convert x values back to the original space
                self._applyinversetransform(x)
                break
        return x
    
    def _checkdomain(self, x):
        #apply transform
        for i in xrange(len(x)):
            f = self.indomain[i]
            if f and not f(x[i]):
                return False
        return True

    def logDensity(self, x, p):
        #convert x values to the space in which variables are Gaussian
        x = self._applytransformcopy(x)
        p = self._applytransformcopy(p)
        return gaus_log_density(x, p, self.sigma)

    def _applytransformcopy(self, parameters):
        parameters = np.copy(parameters)
        #apply transform
        for i in xrange(len(parameters)):
            f = self.func[i]
            if f:
                parameters[i] = f(parameters[i])
        return parameters
    
    def _applyinversetransformcopy(self, parameters):
        #apply transform
        parameters = np.copy(parameters)
        for i in xrange(len(parameters)):
            f = self.invfunc[i]
            if f:
                parameters[i] = f(parameters[i])
        return parameters
    
    def _applytransform(self, parameters):
        #apply transform
        for i in xrange(len(parameters)):
            f = self.func[i]
            if f:
                parameters[i] = f(parameters[i])
        return parameters
    
    def _applyinversetransform(self, parameters):
        #apply transform
        for i in xrange(len(parameters)):
            f = self.invfunc[i]
            if f:
                parameters[i] = f(parameters[i])
        return parameters

    def infoString(self):
        sio = StringIO.StringIO()
        print >>sio, "GaussianProposalFunction(npars={}, stepsize={}".format(len(self.parameterRanges),
                                                                                 self.stepSize
                                                                            )
        for i,(sigma, parrange) in enumerate(zip(self.sigma, self.parameterRanges)):
            print >>sio, "    {} : sigma={}, low={}, high={}".format(i, sigma, parrange[0], parrange[1])
        return sio.getvalue()


###############################################################################

class _MultiVariateGaussianLikelihoodWrapper(MultiVariateGaussianLikelihood):
    def __init__(self, mu, cov):
        parameter_names = ["dummy_par_" + str(i) for i in xrange(len(mu)) ]
        super(_MultiVariateGaussianLikelihoodWrapper, self).__init__(parameter_names, mu, cov)

    def __call__(self, x, mu):
        self._mu = mu # update mean
        return MultiVariateGaussianLikelihood.__call__(self, x)

###############################################################################

class MultiVariateGaussianProposal(object):
    def __init__(self, mu, cov, startindex=0):
        cov = np.copy(cov)
        self._startindex = startindex
        self._stopindex = startindex + len(mu)
        self._lhd = _MultiVariateGaussianLikelihoodWrapper(mu, cov)
        parameter_names = ["par_"+str(i) for i in xrange(len(mu))]
        self._gen = MultiVariateGaussianGenerator(parameter_names, mu=[0.0]*len(mu), cov=cov)

    def logDensity(self, xvec, mu):
        return self._lhd(xvec, mu)

    def generate(self, parameters):
        return self._gen() + parameters[self._startindex:self._stopindex]

###############################################################################

class FixedMultiVariateGaussianProposalFunction(MultiVariateGaussianProposal):
    def __init__(self, mu, cov, startindex=0):
        super(FixedMultiVariateGaussianProposalFunction, self).__init__(cov=cov, mu=mu, startindex=startindex)
        self._lhd = MultiVariateGaussianLikelihood(["dummy_par_" + str(i) for i in xrange(len(mu))], mu, cov)

    def generate(self, parameters):
        x = self._result
        _gaus_generate(x, self._mu, self.sigma)
        return x

    def logDensity(self, xvec, mu):
        return self._lhd(xvec)

    def generate(self, parameters):
        gen = self._gen.next()
        return self._mu + gen

###############################################################################

class CombinedProposalFunction(object):
    def __init__(self, funcs, num_pars):
        self._prop = funcs
        self._fi = []
        start = 0
        for i in xrange(len(funcs)):
            f = funcs[i]
            stop = start + num_pars[i]
            self._fi.append( (f, start, stop) )
            start = stop

    def logDensity(self, xvec, mu):
        return sum([p.logDensity(xvec[start:stop], mu[start:stop]) for p, start, stop in self._fi])

    def generate(self, parameters):
        return np.concatenate([p.generate(parameters[start:stop]) for p, start, stop in self._fi])

    def adapt(self, eff, data=None):
        ret = False
        for f, start, stop in self._fi:
            try:
                ret = ret and f.adapt(eff, data)
            except AttributeError:
                pass # adapt has no been inmplemented
        return ret


###############################################################################

class SinSqTheta2FixedGaussianProposalFunction(object):
    def __init__(self, mu, sigma):
        self._mu = mu
        self._sigma = sigma
        self._result = np.zeros(shape=(1), dtype=float)

    def logDensity(self, xvec, mu):
        return _sinsqtheta2gaussian_log_density_i(xvec[0], self._mu, self._sigma)

    def generate(self, parameters):
        _sinsqtheta2gaussian_generate_i(self._result, self._mu, self._sigma)
        return self._result

###############################################################################

class SinSqTheta2GaussianProposalFunction(object):
    def __init__(self, mu, sigma):
        self._mu = mu
        self._sigma = sigma
        self._result = np.zeros(shape=(1), dtype=float)

    def logDensity(self, xvec, mu):
        return _sinsqtheta2gaussian_log_density_i(xvec[0], _sinsq2theta(mu[0]), self._sigma)

    def generate(self, parameters):
        _sinsqtheta2gaussian_generate_i(self._result, _sinsq2theta(parameters[0]), self._sigma)
        return self._result

###############################################################################

class SinGaussianProposalFunction(object):
    def __init__(self, mu, sigma):
        self._mu = mu
        self._sigma = sigma
        self._result = np.zeros(shape=(1), dtype=float)

    def logDensity(self, xvec, mu):
        return _sinthetagaussian_log_density_i(xvec[0], _sintheta(mu[0]), self._sigma)

    def generate(self, parameters):
        _sinthetagaussian_generate_i(self._result, _sintheta(parameters[0]), self._sigma)
        return self._result

###############################################################################

class SinGaussianFixedProposalFunction(object):
    def __init__(self, mu, sigma):
        self._mu = mu
        self._sigma = sigma
        self._result = np.zeros(shape=(1), dtype=float)

    def logDensity(self, xvec, mu):
        return _sinthetagaussian_log_density_i(xvec[0], _sintheta(self._mu[0]), self._sigma)

    def generate(self, parameters):
        _sinthetagaussian_generate_i(self._result, self._mu, self._sigma)
        return self._result

###############################################################################

class HesseProposalFunction(FixedMultiVariateGaussianProposalFunction):
    def __init__(self, mu, func, initerr, startindex=0):
        cov = self._gethessecovariance(func, mu, initerr)
        super(HesseProposalFunction, self).__init__(cov=cov, mu=mu, startindex=startindex)

    def _gethessecovariance(self, func, mu, initerr):
        hesse = Hesse(func, mu, delta=initerr, verbosity=HesseVerbosity.PRINT_PROGRESS)
        cov = hesse.covariance_matrix()
        return cov
###############################################################################

class ProposalWithSomeParametersFixed(object):
    def __init__(self, proposalfunc, fixedvalues):
        self._fixedvalues = fixedvalues
        self._proposalfunc = proposalfunc

    def logDensity(self, x, p):
        x = np.copy(x)
        p = np.copy(p)
        if not len(x) == len(self._fixedvalues):
            raise Exception("ProposalWithSomeParametersFixed mismatched number of parameters.", len(self._fixedvalues), len(self._fixedvalues))
        for index, val in enumerate(self._fixedvalues):
            if val is not None:
                x[index] = val
                p[index] = val
        return self._proposalfunc.logDensity(x, p)

    def generate(self, parameters):
        result = self._proposalfunc.generate(parameters)
        if not len(result) == len(self._fixedvalues):
            raise Exception("ProposalWithSomeParametersFixed mismatched number of parameters.", len(result), len(self._fixedvalues))
        for index, val in enumerate(self._fixedvalues):
            if val is not None:
                result[index] = val
        return result

    def adapt(self, eff, data=None):
        return self._proposalfunc.adapt(eff, data)

###############################################################################

class SimpleAdaptiveMultiVariateGaussianProposal(object):
    def __init__(self, mu, cov, startindex=0):
        self._mu = mu
        self._cov = cov
        self._startindex = startindex
        self._update()
        self._converged = False
        self._total_scale = 1.0

    def _update(self):
        self._gen = MultiVariateGaussianProposal(self._mu, self._cov, self._startindex)
        return

    def generate(self, parameters):
        return self._gen.generate(parameters)

    def logDensity(self, x, p):
        return self._gen.logDensity(x, p)

    def adapt(self, eff, data=None):
        print "SimpleAdaptiveMultiVariateGaussianProposal eff=%.2f, scale=%.2f" % (eff, self._total_scale)
        if not self._converged:
            if abs(eff - 0.4) < 0.1:
                self._converged = True
                print "SimpleAdaptiveMultiVariateGaussianProposal CONVERGED eff=%.2f, scale=%.2f" % (eff, self._total_scale)
            else:
                cov = self._cov
                scale = 1.0 + 2.0*(eff - 0.4)
                self._total_scale *= scale
                cov *= scale
                self._cov = cov
                self._update()
        return self._converged

###############################################################################

class AdaptiveMultiVariateGaussianProposal(MultiVariateGaussianProposal):
    def __init__(self, *args, **kwargs):
        super(AdaptiveMultiVariateGaussianProposal, self).__init__(*args, **kwargs)
        self._converged = False
        self._total_scale = 1.0
        self._updatedcov = None

    def adapt(self, eff, data):
        print "AdaptiveMultiVariateGaussianProposal eff=%.2f, scale=%.2f" % (eff, self._total_scale)
        if not self._converged:
            if abs(eff - 0.4) < 0.1:
                self._converged = True
                print "CONVERGED"
            else:
                scale = 1.0 + 2.0*(eff - 0.4)
                self._total_scale *= scale
                if eff > 0.1 and self._updatedcov is None:
                    #covariance is likely to be unreliable at low efficiency
                    print "AdaptiveMultiVariateGaussianProposal updating covariance matrix"
                    newcov = self._calccov(data) * self._total_scale
                    self._updatedcov = newcov
                    assert newcov.shape == self._numcov.shape
                else:
                    newcov = self._numcov * scale
                self._numcov = newcov
                self._gen = self._iterrandom(self._numcov)
        return self._converged

    def _calccov(self, data):
        cov = np.cov(data.T)
        return cov

class ProposalWithSomeParametersFixed(object):
    def __init__(self, proposalfunc, fixedvalues):
        self._fixedvalues = fixedvalues
        self._proposalfunc = proposalfunc

    def logDensity(self, x, p):
        x = np.copy(x)
        p = np.copy(p)
        if not len(x) == len(self._fixedvalues):
            raise Exception("ProposalWithSomeParametersFixed mismatched number of parameters.", len(x), len(self._fixedvalues))
        for index, val in enumerate(self._fixedvalues):
            if val is not None:
                x[index] = val
                p[index] = val
        return self._proposalfunc.logDensity(x, p)

    def generate(self, parameters):
        result = self._proposalfunc.generate(parameters)
        if not len(result) == len(self._fixedvalues):
            raise Exception("ProposalWithSomeParametersFixed mismatched number of parameters.", len(result), len(self._fixedvalues))
        for index, val in enumerate(self._fixedvalues):
            if val is not None:
                result[index] = val
        return result

    def adapt(self, eff, data=None):
        return self._proposalfunc.adapt(eff, data)
