import itertools

from simplot.cache import cache
import simplot.sparsehist.sparsehist
from simplot.sparsehist import SparseHistogram
from simplot.binnedmodel.model import BinnedModel as _BinnedModel
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
    def __init__(self, name, binning, observables, data, cache_name=None, systematics=None):
        parameter_names = []
        if systematics:
            for s, _ in systematics:
                parameter_names.append(s)
        super(BinnedSample, self).__init__(parameter_names)
        self.name = name
        self.axisnames = [n for n, _ in binning]
        self.binedges = [np.array(edges, copy=True) for _, edges in binning]
        def func():
            return self._loaddata(data, systematics)
        if cache_name:
            hist, systhist = cache(cache_name, func)
        else:
            hist, systhist = func()
        observabledim = [self.axisnames.index(p) for p in observables]
        xsec_weights = self._buildxsecweights(systematics, systhist, hist)
        self._model = _BinnedModel(parameter_names, hist, observabledim, xsec_weights=xsec_weights)

    def __call__(self, x):
        return self._model.observable(x).flatten()

    def _loaddata(self, data, systematics):
        hist = SparseHistogram(self.binedges)
        systhist = [[SparseHistogram(self.binedges) for val in values] for syst, values in systematics]
        for coord, selweight, systweight in data:
            hist.fill(coord, selweight)
            for isyst in xrange(len(systhist)):
                for ival in xrange(len(systhist[isyst])):
                    systhist[isyst][ival].fill(coord, systweight[isyst, ival])
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

################################################################################

class BinnedSampleWithOscillation(BinnedSample):
    def __init__(self, name, binning, data, distance, cache_name=None):
        self._distance = distance
        super(BinnedSampleWithOscillation, self, name, binning, data, cache_name)

################################################################################

class CombinedBinnedSample(Sample):
    def __init__(self, samples):
        self._samples = samples
        parameter_names = samples[0].parameter_names
        super(CombinedBinnedSample, self).__init__(parameter_names)

    def __call__(self, x):
        return np.concatenate([s(x) for s in self._samples])

################################################################################

