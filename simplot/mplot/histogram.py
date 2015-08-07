import numpy as np
import math

###############################################################################

class Stats:
    poisson = "poisson"
    gaussian = "gaussian" 

###############################################################################

def _countnd(obj):
        count = 0
        try:
            a = obj
            while True:
                a = a[0]
                count += 1
        except Exception:
            pass
        return count

###############################################################################

class HistogramND(object):
    def __init__(self, binning, values=None, label=None):
        #check dimensionality of input objects matches the given nd
        if _countnd(binning) == 1:
            binning = (binning,)
        nd = len(binning)
        if values is None:
            shape = (len(b)-1 for b in binning)
            values = np.zeroes(shape=shape)
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
        self.label = label

    @property
    def xbinning(self):
        """convenient for plotting code that only plots 1D histograms."""
        return self.binning[0]

    @property
    def ybinning(self):
        """convenient for plotting code that only plots 1D histograms."""
        return self.binning[1]

    def fill(self, coord, weight=1.0):
        self.values[coord] += weight

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

def _get_binning_array(xbinning):
    try:
        b = xbinning.binedges
    except AttributeError:
        #object does not have binedges attribute, 
        #assume that we were given an array-like of bin edges. 
        b = np.array(xbinning, dtype=float)
    return b

###############################################################################

def _get_bin_centres(b):
    return (b[1:] + b[:-1]) / 2.0

###############################################################################

def _get_bin_widths(b):
    return (b[1:] - b[:-1])

###############################################################################

def plot_hist_points(ax, xbinning, y, yerr=None, 
              yerrmode=Stats.gaussian,
              color="black",
              marker="o",
              markersize=10.0,
              linestyle="None",
              show_empty_bins=False,
              *args, **kwargs):
    b = _get_binning_array(xbinning)
    #bin centers
    x = _get_bin_centres(b)
    #"xerr" is 1/2 bin width
    xerr = _get_bin_widths(b) / 2.0
    if yerr is None:
        #if no yerr is given assume sqrt(N)
        if yerrmode == Stats.gaussian:
            #TODO deal with weights correctly
            yerr = np.sqrt(y) 
        #if no yerr is given assume poisson statistics
        elif yerrmode == Stats.poisson:
            #TODO implement Poisson errors
            raise NotImplementedError()
    if not show_empty_bins:
        #filter out bins with zero content
        filtered = filter(lambda v: v[2]!=0.0, zip(x, xerr, y, yerr))
        x, xerr, y, yerr = map(np.array, zip(*filtered))
    #draw the plot
    (plotline, caplines, barlinecols) = ax.errorbar(x, y, yerr=yerr, xerr=xerr,
                                                    color=color,
                                                    marker=marker,
                                                    markersize=markersize,
                                                    linestyle=linestyle,
                                                    *args, **kwargs)
    return (plotline, caplines, barlinecols)

###############################################################################

def plot_hist_bars(ax, xbinning, y,
                   *args, **kwargs):
    #convert the input data for passing to Axes.hist
    binedges = _get_binning_array(xbinning)
    x = _get_bin_centres(binedges)
    weights = y
    return ax.hist(x, bins=binedges, weights=weights,
                   *args, **kwargs)

###############################################################################

def plot_hist_lines(ax, xbinning, y, *args, **kwargs):
    binedges = _get_binning_array(xbinning)
    xline = []
    yline = []
    for ii in xrange(len(binedges)):
            eb = binedges[ii]
            xline.append(eb)
            xline.append(eb)
            try:
                low = y[ii-1]
            except IndexError:
                low = 0.0
            try:
                high = y[ii]
            except IndexError:
                high = 0.0
            yline.append(low)
            yline.append(high)
    return ax.plot(xline, yline, *args, **kwargs)

###############################################################################

def _test_plot_hist():
    from nuOscillation.model import binning
    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    dataset = np.random.normal(size=100)
    xbinning = binning.Binning(10, -5.0, 5.0)
    y, binedges = np.histogram(dataset, xbinning.binedges)
    plot_hist_points(ax, xbinning, y)
    plot_hist_bars(ax, xbinning, y)
    plt.show()
    raw_input("wait")
    return

###############################################################################

def main():
    _test_plot_hist()
    return

###############################################################################

if __name__ == "__main__":
    main()
