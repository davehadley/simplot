import itertools

from simplot.pdg import PdgNeutrinoOscillationParameters
from simplot.cache import cache
import simplot.sparsehist.sparsehist
from simplot.sparsehist import SparseHistogram
from simplot.binnedmodel.model import BinnedModel as _BinnedModel
from simplot.binnedmodel.model import BinnedModelWithOscillation as _BinnedModelWithOscillation
from simplot.binnedmodel.xsecweights import XsecWeights, InterpolatedWeightCalc

import numpy as np

################################################################################

class Sample(object):
    def __init__(self, parameter_names):
        self.parameter_names = parameter_names
    def __call__(self, x):
        raise NotImplemented("ERROR: child class should override __call__.")

################################################################################

class BinnedSample(Sample):
    def __init__(self, name, binning, observables, data, cache_name=None, systematics=None, fluxsystematics=None):
        parameter_names = self._build_parameter_names(systematics, fluxsystematics)
        super(BinnedSample, self).__init__(parameter_names)
        self.name = name
        self.axisnames = [n for n, _ in binning]
        self.binedges = [np.array(edges, copy=True) for _, edges in binning]
        def func():
            return self._loaddata(data, systematics)
        if cache_name:
            data = cache(cache_name, func)
        else:
            data = func()
        self._model = self._buildmodel(systematics, data, observables, fluxsystematics)

    def _build_parameter_names(self, systematics, fluxsystematics):
        parameter_names = []
        if systematics:
            for s, _ in systematics:
                parameter_names.append(s)
        if fluxsystematics:
            for s in fluxsystematics.parameter_names:
                parameter_names.append(s)
        return parameter_names

    def _buildmodel(self, systematics, data, observables, fluxsystematics):
        hist, systhist = data
        observabledim = [self.axisnames.index(p) for p in observables]
        xsec_weights = self._buildxsecweights(systematics, systhist, hist)
        flux_weights = self._buildfluxweights(fluxsystematics)
        return _BinnedModel(self.parameter_names, hist, observabledim, xsec_weights=xsec_weights, flux_weights=flux_weights)

    def __call__(self, x):
        return self._model.observable(x).flatten()

    def _loaddata(self, data, systematics):
        hist = SparseHistogram(self.binedges)
        if systematics:
            systhist = [[SparseHistogram(self.binedges) for val in values] for syst, values in systematics]
        else:
            systhist = []
        for coord, selweight, systweight in data:
            hist.fill(coord, selweight)
            for isyst in xrange(len(systhist)):
                for ival in xrange(len(systhist[isyst])):
                    systhist[isyst][ival].fill(coord, selweight * systweight[isyst][ival])
        return hist, systhist

    def _buildxsecweights(self, systematicsvalues, systhist, hist):
        xsecweights = None
        if systematicsvalues is not None:
            wclist = []
            for isyst, (syst, parval) in enumerate(systematicsvalues):
                assert len(parval) == len(systhist[isyst])
                l = zip(parval, systhist[isyst])
                l.sort() # sort by parameter value
                weights = [x[1].array() for x in l] 
                parval = [x[0] for x in l]
                wc = InterpolatedWeightCalc(hist.array(), parval, weights, syst, self.parameter_names)
                wclist.append(wc)
            xsecweights = XsecWeights(hist.array(), wclist)
        return xsecweights

    def _buildfluxweights(self, fluxsystematics):
        fw = None
        if fluxsystematics:
            #we don't know how to build the flux systematics,
            # user must provide a callable that constructs the object
            fw = fluxsystematics(self.parameter_names)
        return fw

################################################################################

class BinnedSampleWithOscillation(BinnedSample):
    def __init__(self, name, binning, observables, data, enuaxis, flavaxis, distance, cache_name=None, systematics=None, fluxsystematics=None, probabilitycalc=None):
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
                                                          fluxsystematics=fluxsystematics,
        )

    def _build_parameter_names(self, systematics, fluxsystematics):
        parameter_names = list(PdgNeutrinoOscillationParameters.ALL_PARS_SINSQ2)
        if systematics:
            for s, _ in systematics:
                parameter_names.append(s)
        if fluxsystematics:
            for s in fluxsystematics.parameter_names:
                parameter_names.append(s)
        return parameter_names

    def _buildmodel(self, systematics, data, observables, fluxsystematics):
        selhist, noselhist, selsysthist, noselsysthist = data
        observabledim = [self.axisnames.index(p) for p in observables]
        enudim = self.axisnames.index(self._enu_axis_name)
        flavdim = self.axisnames.index(self._flav_axis_name)
        xsec_weights = self._buildxsecweights(systematics, selsysthist, selhist)
        flux_weights = self._buildfluxweights(fluxsystematics)
        probabilitycalc = self._probabilitycalc
        return _BinnedModelWithOscillation(self.parameter_names, selhist, noselhist, observabledim, enudim, flavdim, None, [self._distance], xsec_weights=xsec_weights, flux_weights=flux_weights, probabilitycalc=probabilitycalc)

    def _loaddata(self, data, systematics):
        selhist = SparseHistogram(self.binedges)
        noselhist = SparseHistogram(self.binedges)
        if systematics:
            selsysthist = [[SparseHistogram(self.binedges) for val in values] for syst, values in systematics]
            noselsysthist = [[SparseHistogram(self.binedges) for val in values] for syst, values in systematics]
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
    def __init__(self, samples, parameter_order=None):
        self._samples = samples
        parameter_names, mapping = self._determine_parameter_mapping(samples, parameter_order=parameter_order)
        self._par_map = mapping
        super(CombinedBinnedSample, self).__init__(parameter_names)

    def __call__(self, x):
        return np.concatenate([s(self._get_args(x, i)) for i, s in enumerate(self._samples)])

    def _get_args(self, x, samplenum):
        x2 = np.fromiter(itertools.imap(x.__getitem__, self._par_map[samplenum]), x.dtype)
        return x2

    def _determine_parameter_mapping(self, samples, parameter_order=None):
        parameter_names = []
        for s in samples:
            for p in s.parameter_names:
                if p not in parameter_names:
                    parameter_names.append(p)
        if parameter_order:
            if not set(parameter_names) == set(parameter_order):
                raise Exception("user parameter order does not contain the correct parameters")
            parameter_names = list(parameter_order)
        mapping = []
        for s in samples:
            m = np.array([parameter_names.index(p) for p in s.parameter_names], dtype=int)
            mapping.append(m)
        return parameter_names, mapping

################################################################################

