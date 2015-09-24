import itertools

from simplot.pdg import PdgNeutrinoOscillationParameters
from simplot.cache import cache
from simplot.mc.montecarlo import MonteCarloParameterMismatch
import simplot.sparsehist.sparsehist
from simplot.sparsehist import SparseHistogram
from simplot.binnedmodel.model import BinnedModel as _BinnedModel
from simplot.binnedmodel.model import BinnedModelWithOscillation as _BinnedModelWithOscillation

import numpy as np

################################################################################

class Sample(object):
    def __init__(self, parameter_names):
        self.parameter_names = parameter_names
    def __call__(self, x):
        raise NotImplementedError("ERROR: child class should override __call__.")

################################################################################

class BinnedSample(Sample):
    def __init__(self, name, binning, observables, data, cache_name=None, systematics=None):
        parameter_names = self._build_parameter_names(systematics)
        super(BinnedSample, self).__init__(parameter_names)
        self.name = name
        self.axisnames = [n for n, _ in binning]
        self.binedges = [np.array(edges, copy=True) for _, edges in binning]
        self.observables = observables
        def func():
            return self._loaddata(data, systematics)
        if cache_name:
            data = cache(cache_name, func)
        else:
            data = func()
        self._model, self.N_sel, self.N_nosel = self._buildmodel(systematics, data, observables)

    def _build_parameter_names(self, systematics):
        parameter_names = []
        if systematics:
            parameter_names += systematics.parameter_names
        return parameter_names

    def _buildmodel(self, systematics, data, observables):
        hist, systhist = data
        observabledim = [self.axisnames.index(p) for p in observables]
        #xsec_weights = self._buildxsecweights(systematics, systhist, hist)
        #flux_weights = self._buildfluxweights(fluxsystematics)
        xsec_weights, flux_weights = None, None
        if systematics:
            xsec_weights, flux_weights = systematics(self.parameter_names, systhist, hist)
        return _BinnedModel(self.parameter_names, hist, observabledim, xsec_weights=xsec_weights, flux_weights=flux_weights), hist, None

    def __call__(self, x):
        if len(x) != len(self.parameter_names):
            raise ValueError("Sample called with wrong number of parameters")
        return self._model.observable(x).flatten()

    def array(self, x):
        return self._model(x)

    def _loaddata(self, data, systematics):
        hist = SparseHistogram(self.binedges)
        if systematics:
            systhist = [[SparseHistogram(self.binedges) for val in values] for syst, values in systematics.spline_parameter_values]
        else:
            systhist = []
        for coord, selweight, systweight in data:
            hist.fill(coord, selweight)
            for isyst in xrange(len(systhist)):
                for ival in xrange(len(systhist[isyst])):
                    systhist[isyst][ival].fill(coord, selweight * systweight[isyst][ival])
        return hist, systhist

################################################################################

class BinnedSampleWithOscillation(BinnedSample):
    def __init__(self, name, binning, observables, data, enuaxis, flavaxis, distance, cache_name=None, systematics=None, probabilitycalc=None):
        self._enu_axis_name = enuaxis
        self._flav_axis_name = flavaxis
        self._distance = distance
        self._probabilitycalc = probabilitycalc
        super(BinnedSampleWithOscillation, self).__init__(name=name, 
                                                          binning=binning, 
                                                          observables=observables,
                                                          data=data, 
                                                          cache_name=cache_name,
                                                          systematics=systematics,
        )

    def _build_parameter_names(self, systematics):
        parameter_names = list(PdgNeutrinoOscillationParameters.ALL_PARS_SINSQ2)
        if systematics:
            parameter_names += systematics.parameter_names
        return parameter_names

    def _buildmodel(self, systematics, data, observables):
        selhist, noselhist, selsysthist, noselsysthist = data
        observabledim = [self.axisnames.index(p) for p in observables]
        enudim = self.axisnames.index(self._enu_axis_name)
        flavdim = self.axisnames.index(self._flav_axis_name)
        xsec_weights, flux_weights = None, None
        if systematics:
            xsec_weights, flux_weights = systematics(self.parameter_names, selsysthist, selhist)
        probabilitycalc = self._probabilitycalc
        if probabilitycalc is None:
            #no user supplied probability calculator, use prob3++
            import simplot.rootprob3pp.lib
            import ROOT
            probabilitycalc = ROOT.crootprob3pp.Probability()
        return _BinnedModelWithOscillation(self.parameter_names, selhist, noselhist, observabledim, enudim, flavdim, None, [self._distance], xsec_weights=xsec_weights, flux_weights=flux_weights, probabilitycalc=probabilitycalc), selhist, noselhist

    def _loaddata(self, data, systematics):
        selhist = SparseHistogram(self.binedges)
        noselhist = SparseHistogram(self.binedges)
        if systematics:
            selsysthist = [[SparseHistogram(self.binedges) for val in values] for syst, values in systematics.spline_parameter_values]
            noselsysthist = [[SparseHistogram(self.binedges) for val in values] for syst, values in systematics.spline_parameter_values]
        else:
            selsysthist = []
            noselsysthist = []
        for coord, selweight, noselweight, systweight in data:
            if selweight != 0:
                selhist.fill(coord, selweight)
            noselhist.fill(coord, noselweight)
            for isyst in xrange(len(selsysthist)):
                for ival in xrange(len(selsysthist[isyst])):
                    if selweight != 0:
                        selsysthist[isyst][ival].fill(coord, selweight * systweight[isyst][ival])
                    noselsysthist[isyst][ival].fill(coord, noselweight * systweight[isyst][ival])
        return selhist, noselhist, selsysthist, noselsysthist

################################################################################

class CombinedBinnedSample(Sample):
    def __init__(self, samples, parameter_order=None, ignoreerrors=False):
        self._samples = samples
        parameter_names, mapping = self._determine_parameter_mapping(samples, parameter_order=parameter_order, ignoreerrors=ignoreerrors)
        self._par_map = mapping
        super(CombinedBinnedSample, self).__init__(parameter_names)

    def sample_parameters(self, pars, samplenum):
        return self._get_args(pars, samplenum)

    def eval_sample(self, pars, samplenum):
        return self._samples[samplenum](self._get_args(pars, samplenum))

    def array_sample(self, pars, samplenum):
        return self._samples[samplenum].array(self._get_args(pars, samplenum))

    def __call__(self, x):
        if len(x) != len(self.parameter_names):
            raise ValueError("Sample called with wrong number of parameters")
        return np.concatenate([s(self._get_args(x, i)) for i, s in enumerate(self._samples)])

    def _get_args(self, x, samplenum):
        x2 = np.fromiter(itertools.imap(x.__getitem__, self._par_map[samplenum]), x.dtype)
        return x2

    def _determine_parameter_mapping(self, samples, parameter_order=None, ignoreerrors=False):
        parameter_names = []
        for s in samples:
            for p in s.parameter_names:
                if p not in parameter_names:
                    parameter_names.append(p)
        if parameter_order:
            if (not set(parameter_names) == set(parameter_order)) and (not ignoreerrors):
                raise Exception("user parameter order does not contain the correct parameters", MonteCarloParameterMismatch.compare_parameters_message(parameter_names, parameter_order))
            parameter_names = list(parameter_order)
        mapping = []
        for s in samples:
            m = np.array([parameter_names.index(p) for p in s.parameter_names], dtype=int)
            mapping.append(m)
        return parameter_names, mapping

    @property
    def samples(self):
        return self._samples

################################################################################

