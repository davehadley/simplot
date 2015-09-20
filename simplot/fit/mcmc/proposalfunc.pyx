import numpy as np
cimport numpy as np
import StringIO
import itertools
import math
from simplot.fit.hesse import Hesse
from simplot.fit.hesse import Verbosity as HesseVerbosity
from simplot.fit.mcmc.metropolishastings import McMcSetupError
from simplot.mc.likelihood import MultiVariateGaussianLikelihood

from simplot.mc.clikelihood import gaus_log_density

from libc.math cimport sin, sqrt, log, M_PI, asin

###############################################################################

_RANDOM = np.random.RandomState(231321)

cdef _gaus_generate(np.ndarray[float, ndim=1] result, np.ndarray[float, ndim=1] mu, np.ndarray[float, ndim=1] sigma):
    cdef int N = len(result)
    cdef np.ndarray[float, ndim=1] normal = _RANDOM.normal(size=N);
    cdef int ii
    for ii in xrange(N):
        result[ii] = mu[ii] + sigma[ii] * normal[ii]
    return

cdef double _sinsqtheta2gaussian_log_density_i(x, mu, sigma):
    pass

cdef _sinsqtheta2gaussian_generate_i(x, mu, sigma):
    pass

cdef _sinsq2theta(theta):
    return np.power(np.sin(2.0*theta), 2);

cdef double _sinthetagaussian_log_density_i(x, mu, sigma):
        cdef int N = len(x)
        cdef double y
        cdef double jacobian = 0.0
        cdef double result = gaus_log_density(y, mu, sigma)
        #compute jacobian
        for ii in xrange(N):
            y = sin(2.0*x[ii])
            y = y*y # take square
            #set limits on y
            if(y < 0.001):
                y = 0.001;
            elif(y > (1.0-0.001)):
                y = (1.0-0.001)
            jacobian += -log(4.0 * sqrt((1.0 - y)*y));
        result += jacobian
        return result;

cdef _sintheta(x):
   return np.sin(x)

cdef _sinthetagaussian_generate_i(result, mu, sigma):
    cdef int N = len(result)
    cdef np.ndarray[float, ndim=1] normal = _RANDOM.normal(size=N);
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
    def __init__(self, parameterRanges, stepSize=0.4, sigma=None):
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
        self.stepSize = stepSize
        self._result = np.zeros(self.nParameters)
        self.sigma = self._calcWidthParameters(sigma, parameterRanges)
        self._sanityCheckState()
        

    def _sanityCheckState(self):
        if not len(self.sigma) == self.nParameters:
            raise McMcSetupError("sigma list must have an entry for each parameter", self.infoString())
        if not np.all(np.isfinite(self.sigma)):
            raise McMcSetupError("sigma list contains NaN", self.infoString())
        
    
    
    def _calcWidthParameters(self, sigma, parameterRanges):
        #automatically calculate width parameters
        for i,(vMin,vMax) in enumerate(parameterRanges):
            if sigma[i] is None:
                s = self.stepSize*(vMax-vMin)
                sigma[i] = s
        sigma = np.array(sigma, dtype=float)
        return sigma
    
    def generate(self, parameters):
        x = self._result
        _gaus_generate(x, parameters, self.sigma)
        return x
    
    def logDensity(self, x, p):
        return gaus_log_density(x, p, self.sigma)
    
    def infoString(self):
        sio = StringIO.StringIO()
        print >>sio, "GaussianProposalFunction(npars={}, stepsize={}".format(len(self.parameterRanges),
                                                                                 self.stepSize
                                                                            )
        for i,(sigma, parrange) in enumerate(zip(self.sigma, self.parameterRanges)):
            print >>sio, "    {} : sigma={}, low={}, high={}".format(i, sigma, parrange[0], parrange[1])
        return sio.getvalue()

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
        self._npars = cov.GetNrows()
        self._numcov = self._convertmatrix_roottoarray(cov)
        self._startindex = startindex
        self._endindex = cov.GetNrows() + startindex
        self._gen = self._iterrandom(cov)
        self._invcov = cov.Clone().Invert()
        self._lhd = _MultiVariateGaussianLikelihoodWrapper(mu, cov)
        try:
            iter(mu)
        except TypeError:
            #not iterable, use same value for all rows
            self._mu = np.array([mu]*cov.GetNrows(), dtype=float)
        else:
            self._mu = np.array(mu)
        self._zeroes = np.zeros(cov.GetNrows(), dtype=float)

    def _iterrandom(self, cov):
            while True:
                #for speed, generate many events at a time
                size = 10000
                #convert cov to np format
                val = np.random.multivariate_normal(self._zeroes, self._numcov, size)
                for i in xrange(size):
                    yield val[i]

    def logDensity(self, xvec, mu):
        return self._lhd(xvec, mu)

    def generate(self, parameters):
        vec = self._gen.next()
        return vec + parameters


    def infoString(self):
        sio = StringIO.StringIO()
        print >>sio, "MultiVariateGaussianProposalFunction(npars={})".format(self._npars),
        for i in xrange(self._npars):
            sigma = math.sqrt(abs(self._numcov[i][i]))
            mu = math.sqrt(abs(self._mu[i]))
            print >>sio, "    {} : mu={}, sigma={}".format(i, mu, sigma)
        print >>sio, self._numcov
        return sio.getvalue()

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

class SimpleAdaptiveMultiVariateGaussianProposal(MultiVariateGaussianProposal):
    def __init__(self, *args, **kwargs):
        super(SimpleAdaptiveMultiVariateGaussianProposal, self).__init__(*args, **kwargs)
        self._converged = False
        self._total_scale = 1.0

    def adapt(self, eff, data=None):
        print "SimpleAdaptiveMultiVariateGaussianProposal eff=%.2f, scale=%.2f" % (eff, self._total_scale)
        if not self._converged:
            if abs(eff - 0.4) < 0.1:
                self._converged = True
                print "SimpleAdaptiveMultiVariateGaussianProposal CONVERGED eff=%.2f, scale=%.2f" % (eff, self._total_scale)
            else:
                cov = self._numcov
                scale = 1.0 + 2.0*(eff - 0.4)
                self._total_scale *= scale
                cov *= scale
                # if eff < 0.2:
                #     #efficiency is low, decrease covariance size
                #     cov *= scale
                # elif eff > 0.5:
                #     #efficiency is high, increase covariance size
                #     cov *= 1.1
                self._numcov = cov
                #update generator
                self._gen = self._iterrandom(cov)
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
