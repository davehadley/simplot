import numpy as np
import math
import copy

from bisect import bisect_right

###############################################################################

def _countnd(obj):
    count = 0
    try:
        a = obj
        while True:
            a = iter(a).next()
            count += 1
    except Exception:
        pass
    return count

###############################################################################

class HistogramND(object):
    def __init__(self, binning, values=None, label=None):
        self.overflow = 0.0
        self.overflowsumw2 = 0.0
        #check dimensionality of input objects matches the given nd
        if _countnd(binning) == 1:
            binning = (binning,)
        nd = len(binning)
        binning = [get_binning_array(b) for b in binning] #convert to array
        if values is None:
            shape = tuple(len(b)-1 for b in binning)
            values = np.zeros(shape=shape)
        ndv = _countnd(values)
        if not (nd == ndv):
            raise Exception("HistogramND inputs have inconsistent dimensionality", nd, ndv)
        #convert input label into HigtogramNDLabel object
        if label is None:
            label = HistogramNDLabel(binning)
        elif isinstance(label, basestring):
            label = HistogramNDLabel(binning, label=label)
        self.nd = nd
        self.binning = binning
        self.values = values
        self.sumw2 = np.power(values, 2)
        self.label = label

    @property
    def xbinning(self):
        """convenient for plotting code that only plots 1D histograms."""
        return self.binning[0]

    @property
    def ybinning(self):
        """convenient for plotting code that only plots 1D histograms."""
        return self.binning[1]

    @property
    def errors(self):
        return np.sqrt(self.sumw2)

    def fill(self, coord, weight=1.0):
        index = self.binindex(coord)
        try:
            self.values[index] += weight
            self.sumw2[index] += weight**2
        except IndexError:
            self.overflow += weight
            self.overflowsumw2 += weight**2

    def setbin(self, index, value):
        self.values[index] = value

    def getbin(self, index):
        return self.values[index]

    def eval(self, coord):
        index = self.binindex(coord)
        return self.getbin(index)

    def binindex(self, coord):
        if isinstance(coord, float):
            coord = (coord,)
        return tuple(self._findbin(self.binning[i], x) for i,x in enumerate(coord))

    def _findbin(self, array, x):
        return bisect_right(array, x) - 1

    def scale(self, scale):
        self.values *= scale
        self.sumw2 *= scale**2
        self.overflow *= scale
        self.overflowsumw2 *= scale**2
        return

    def sum(self):
        return np.sum(self.values)

    def add(self, hist):
        self.values += hist.values
        self.sumw2 += hist.sumw2
        self.overflow += hist.overflow
        self.overflowsumw2 += hist.overflowsumw2

    def clone(self):
        return copy.deepcopy(self)

###############################################################################

class HistogramNDLabel(object):
    def __init__(self, binning, label=None, axislabels=None, axisunits=None, binlabels=None):
        if _countnd(binning) == 1:
            binning = (binning,)
        nd = len(binning)
        #guarantee labels have required dimensionality
        if axislabels is None:
            axislabels = tuple([None] * (nd + 1))
        if axisunits is None:
            axisunits = tuple([None] * (nd + 1))
        if isinstance(axislabels, basestring):
            axislabels = (axislabels, None)
        if isinstance(axisunits, basestring):
            axisunits = (axisunits, None)
        #special case, axes are labelled except for the histogram value axis
        if len(axislabels) == nd:
            axislabels = tuple(list(axislabels) + [None])
        if len(axisunits) == nd:
            axislabels = tuple(list(axisunits) + [None])
        if not len(axislabels) == nd + 1:
            raise Exception("wrong number of axis labels given", nd, axislabels)
        if not len(axisunits) == nd + 1:
            raise Exception("wrong number of axis units given", nd, axisunits)
        #setup bin labels
        if binlabels is None:
            binlabels = [None] * nd
        else:
            for label in binlabels:
                if isinstance(label, basestring):
                    binlabels = [binlabels]
                    break
        if not len(binlabels) == nd:
            raise Exception("wrong dimensionality for bin labels", nd, len(binlabels))
        #store values
        self.nd = nd
        self.label = label
        self.axislabels = axislabels
        self.axisunits = axisunits
        self.binlabels = binlabels

    def getbinlabel(self, binnum, axis=0):
        ret = None
        labels = self.binlabels[axis]
        if labels is not None:
            ret = labels[binnum]
        return ret            


###############################################################################

def get_binning_array(xbinning):
    try:
        #may be a complex object that provide binedges
        b = xbinning.binedges
    except AttributeError:
        #object does not have binedges attribute, 
        #assume that we were given an array-like of bin edges. 
        b = np.array(xbinning, dtype=float)
    return b


