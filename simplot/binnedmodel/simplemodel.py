from collections import OrderedDict

import numpy as np

from simplot.mc.generators import GeneratorList, MultiVariateGaussianGenerator, GeneratorSubset
from simplot.mc.montecarlo import ToyMC
from simplot.mc.statistics import Covariance, Mean, calculate_statistics_from_toymc
from simplot.cache import cache

from simplot.binnedmodel.xsecweights import SimpleInterpolatedWeightCalc
from simplot.binnedmodel.sample import Sample

_PAR_BIN_FORMAT = "bin%02.0f"

################################################################################

class SimpleMcBuilder(object):

    def build(self, name, toymc, keep=None, cache_name=None, npe=1000):
        self.name = name
        try:
            #assume keep is dict(parnames, splinepoints)
            items = keep.items()
            keep, spline_points = zip(*items)
        except:
            #assume keep is list(parnames)
            spline_points = None
        cov, mean = self._generate_covariance_with_cache(toymc=toymc, keep=keep, npe=npe, cache_name=cache_name)
        splines = self._generate_splines_with_cache(toymc=toymc, keep=keep, spline_points=spline_points, cache_name=cache_name)
        ratevector = self._buildratevector(mean, splines)
        generator = self._buildgenerator(toymc, keep, cov)
        toymc = ToyMC(ratevector, generator)
        return toymc

    def _generate_covariance_with_cache(self, toymc, keep, npe=1000, cache_name=None):    
        def func(self=self, toymc=toymc, keep=keep):
            return self._generate_covariance(toymc, keep)
        if cache_name is not None:
            cov = cache("SimpleMcBuilderCovariance_" + cache_name, func)
        else:
            cov = func()
        return cov

    def _generate_splines_with_cache(self, toymc, keep=None, spline_points=None, cache_name=None):
        def func(self=self, toymc=toymc, keep=keep, spline_points=spline_points):
            return self._generate_splines(toymc, keep, spline_points)
        if cache_name is not None:
            cov = cache("SimpleMcBuilderSplines_" + cache_name, func)
        else:
            cov = func()
        return cov

    def _generate_splines(self, toymc, keep=None, spline_points=None):
        result = OrderedDict()
        if spline_points is None:
            spline_points = self._autosplinepoints(toymc, keep)
        for par, xpoints in spline_points.iteritems():
            y = []
            for x in points:
                index = toymc.generator.parameter_names.index(par)
                pars = np.array(toymc.generator.start_values)
                pars[index] = x
                y.append(np.array(toymc.ratevector(pars)))
            wc = SimpleInterpolatedWeightCalc(x=points, y=y)
            result[par] = wc
        return result

    def _autosplinepoints(self, toymc, keep):
        spline_points = OrderedDict()
        if keep is None:
            keep = []
        for par in keep:
            sigma = toymc.generator.getsigma(par)
            spline_points[par] = [float(ii)*sigma for ii in xrange(-5, 6)]
        return spline_points

    def _generate_covariance(self, toymc, keep, npe=1000):
        cov = Covariance(fractional=True)
        mean = Mean()
        generator = toymc.generator
        if keep is not None:
            generator.setfixed(set(keep))
        name = "generate covariance matrix for " + str(self.name)
        calculate_statistics_from_toymc(toymc, [cov, mean], npe, name=name)
        if keep is not None:
            generator.setfixed(None)
        return cov.eval(), mean.eval()

    def _buildratevector(self, mean, splines):
        model = SimpleModel(mean, splines)
        return model

    def _buildgenerator(self, toymc, keep, cov):
        N = len(cov)
        names = [(_PAR_BIN_FORMAT % ii) for ii in xrange(N)]
        mu = [1.0 for ii in xrange(N)]
        gen = MultiVariateGaussianGenerator(names, mu=mu, cov=cov)
        if keep is not None:
            subset = GeneratorSubset(keep, toymc.generator)
            gen = GeneratorList(subset, mvgaus)
        return gen

################################################################################

class SimpleModel(Sample):
    def __init__(self, nominal, splines):
        self._nominal = np.array(nominal, dtype=float)
        parameter_names = []
        interp = []
        for parname, wc in splines:
            parameter_names.append(parname)
            interp.append(wc)
        self._interp = interp
        self._binweights_start = len(parameter_names)
        self._binweights_end = self._binweights_start + len(self._nominal)
        for ii in xrange(len(self._nominal)):
            parameter_names.append(_PAR_BIN_FORMAT % ii)
        super(SimpleModel, self).__init__(parameter_names)
    
    def __call__(self, x):
        binweights = x[self._binweights_start:self._binweights_end]
        interpolatedweights = self._interpolatedweights(x)
        return interpolatedweights * binweights * self._nominal

    def _interpolatedweights(self, x):
        result = np.ones(len(self._nominal))
        for wc in self._interp:
            result *= wc(x)
        return result

################################################################################

class SimpleModelWithOscillation(SimpleModel):
    pass
