from collections import OrderedDict

import numpy as np
import ROOT

from simplot.mc.generators import GeneratorList, MultiVariateGaussianGenerator, GeneratorSubset
from simplot.mc.montecarlo import ToyMC
from simplot.mc.statistics import Covariance, Mean, calculate_statistics_from_toymc
from simplot.cache import cache

from simplot.binnedmodel.xsecweights import SimpleInterpolatedWeightCalc
from simplot.binnedmodel.sample import Sample, BinnedSample, BinnedSampleWithOscillation, CombinedBinnedSample

from simplot.binnedmodel.simplemodelwithosc import SimpleBinnedModelWithOscillation

_PAR_BIN_FORMAT = "bin%02.0f"

################################################################################

class SimpleMcBuilder(object):

    def build(self, name, toymc, keep=None, cache_name=None, npe=1000, fixed=None):
        self.name = name
        try:
            #assume keep is dict(parnames, splinepoints)
            keys = keep.keys()
            spline_points = OrderedDict()
            for k in keep:
                spline_points[k] = keep[k]
        except AttributeError:
            #assume keep is list(parnames)
            spline_points = None
        cov, mean = self._generate_covariance_with_cache(toymc=toymc, keep=keep, npe=npe, cache_name=cache_name, fixed=fixed)
        splines = self._generate_splines_with_cache(toymc=toymc, nominal=toymc.asimov().vec, keep=keep, spline_points=spline_points, cache_name=cache_name)
        ratevector = self._buildratevector(mean, splines)
        generator = self._buildgenerator(toymc, keep, cov)
        toymc = ToyMC(ratevector, generator)
        return toymc, cov

    def _generate_covariance_with_cache(self, toymc, keep, npe=1000, cache_name=None, fixed=None):    
        def func(self=self, toymc=toymc, keep=keep):
            return self._generate_covariance(toymc, keep, npe=npe, fixed=fixed)
        if cache_name is not None:
            cov = cache("SimpleMcBuilderCovariance_" + cache_name, func)
        else:
            cov = func()
        return cov

    def _generate_splines_with_cache(self, toymc, nominal, keep=None, spline_points=None, cache_name=None):
        def func(self=self, toymc=toymc, nominal=nominal, keep=keep, spline_points=spline_points):
            return self._generate_splines(toymc, nominal, keep, spline_points)
        if cache_name is not None:
            cov = cache("SimpleMcBuilderSplines_" + cache_name, func)
        else:
            cov = func()
        return cov

    def _generate_splines(self, toymc, nominal, keep=None, spline_points=None):
        result = OrderedDict()
        if spline_points is None:
            spline_points = self._autosplinepoints(toymc, keep)
        for par, xpoints in spline_points.iteritems():
            ypoints = []
            for x in xpoints:
                index = toymc.generator.parameter_names.index(par)
                pars = np.array(toymc.generator.start_values)
                pars[index] = x
                ypoints.append(np.array(toymc.ratevector(pars)))
            wc = SimpleInterpolatedWeightCalc(nominalvalues=nominal, parvalues=xpoints, arrays=ypoints, parname=par, parameternames=keep)
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

    def _generate_covariance(self, toymc, keep, npe=1000, fixed=None):
        cov = Covariance(fractional=True)
        mean = Mean()
        generator = toymc.generator
        if keep is not None and fixed is not None:
            generator.setfixed(set(keep) + set(fixed))
        elif keep is not None:
            generator.setfixed(set(keep))
        elif fixed is not None:
            generator.setfixed(set(fixed))
        name = "generate covariance matrix for " + str(self.name)
        calculate_statistics_from_toymc(toymc, [cov, mean], npe=npe, name=name)
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
            gen = GeneratorList(subset, gen)
        return gen

################################################################################

class SimpleModel(Sample):
    def __init__(self, nominal, splines, binoffset=0):
        self._nominal = np.array(nominal, dtype=float)
        parameter_names = []
        interp = []
        for parname, wc in splines.iteritems():
            parameter_names.append(parname)
            interp.append(wc)
        self._interp = interp
        self._binweights_start = len(parameter_names)
        self._binweights_end = self._binweights_start + len(self._nominal)
        for ii in xrange(len(self._nominal)):
            parameter_names.append(_PAR_BIN_FORMAT % (ii+binoffset))
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

class SimpleMcWithOscillationBuilder(SimpleMcBuilder):
    def build(self, name, toymc, sample, cache_name=None, npe=1000, probabilitycalc=None, fixed=None):
        self.name = name
        oscpars = toymc.generator.parameter_names[:6]
        cov, mean = self._generate_covariance_with_cache(toymc=toymc, keep=oscpars, npe=npe, cache_name=cache_name, fixed=fixed)
        generator = self._buildgenerator(toymc, oscpars, cov)
        ratevector = self._buildratevector(oscpars, toymc, sample, probabilitycalc=probabilitycalc)
        toymc = ToyMC(ratevector, generator)
        return toymc, cov

        
    def _buildratevector(self, oscpars, toymc, sample, binoffset=0, probabilitycalc=None):
        N_sel, N_nosel, observables, enubinning, dim_enu, dim_nupdg, detdist = self._determine_properties(sample)
        N_sel = self._transform_array(N_sel, observables, enubinning, dim_enu, dim_nupdg)
        N_nosel = self._transform_array(N_nosel, observables, enubinning, dim_enu, dim_nupdg)
        parnames = oscpars + [_PAR_BIN_FORMAT % (ii+binoffset) for ii in xrange(N_sel.shape[2])]
        if probabilitycalc is None:
            probabilitycalc = ROOT.crootprob3pp.Probability()
        ratevector = SimpleBinnedModelWithOscillation(parnames, N_sel, N_nosel, enubinning, detdist, probabilitycalc=probabilitycalc)
        return ratevector

    def _transform_array(self, arr, observables, enubinning, dim_enu, dim_flav):
        keep = observables
        enuvec = []
        for ienu in xrange(len(enubinning) - 1):
            flavvec = []
            for iflav in xrange(4):
                range_ = {dim_enu:(ienu, ienu+1), dim_flav:(iflav, iflav+1)}
                r = list(arr.project(keep, range_=range_).flatten())
                flavvec.append(r)
            enuvec.append(flavvec)
        return np.array(enuvec)

    def _determine_properties(self, sample):
        dim_enu = sample.axisnames.index(sample._enu_axis_name)
        dim_nupdg = sample.axisnames.index(sample._flav_axis_name)
        enubinning = sample.binedges[dim_enu]
        observables = [sample.axisnames.index(o) for o in sample.observables]
        N_sel = sample.N_sel.array()
        N_nosel = sample.N_nosel.array()
        detdist = sample._distance
        return N_sel, N_nosel, observables, enubinning, dim_enu, dim_nupdg, detdist

################################################################################

class SimpleCombinedMcWithOscillationBuilder(SimpleMcWithOscillationBuilder):

    def build(self, name, toymc, sample, cache_name=None, npe=1000, probabilitycalc=None, fixed=None):
        self.name = name
        oscpars = toymc.generator.parameter_names[:6]
        cov, mean = self._generate_covariance_with_cache(toymc=toymc, keep=oscpars, npe=npe, cache_name=cache_name, fixed=fixed)
        generator = self._buildgenerator(toymc, oscpars, cov)
        ratevector = self._buildratevector(oscpars, toymc, sample, generator, probabilitycalc=probabilitycalc)
        toymc = ToyMC(ratevector, generator)
        return toymc, cov
        
    def _buildratevector(self, oscpars, toymc, samples, generator, probabilitycalc=None):
        converted = []
        binoffset = 0
        for s in samples:
            s = self._convert_sample(oscpars, toymc, s, binoffset=binoffset, probabilitycalc=probabilitycalc)
            binoffset += len(s.parameter_names)
            converted.append(s)
        return CombinedBinnedSample(converted, generator.parameter_names)

    def _convert_sample(self, oscpars, toymc, sample, binoffset, probabilitycalc=None):
        result = None
        if isinstance(sample, BinnedSampleWithOscillation):
            result = SimpleMcWithOscillationBuilder._buildratevector(self, oscpars, toymc, sample, binoffset=binoffset, probabilitycalc=probabilitycalc)
        elif isinstance(sample, BinnedSample):
            #determine nominal value for this object
            asimovpars = toymc.asimov().pars
            samplepars = [asimovpars[toymc.ratevector.parameter_names.index(n)] for n in  sample.parameter_names]
            nominal = sample(samplepars)
            result = SimpleModel(nominal, {}, binoffset=binoffset)
        else:
            raise Exception("SimpleCombinedMcWithOscillationBuilder does not know how to convert sample of type " + str(type(sample)), sample)
        return result
