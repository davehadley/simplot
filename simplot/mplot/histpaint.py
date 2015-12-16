from collections import OrderedDict
import matplotlib.pyplot as plt
import numpy as np

from simplot.rootplot import drawoptions
from simplot.rootplot import drawtools
from simplot.mplot.histogram import get_binning_array

###############################################################################

class Stats:
    poisson = "poisson"
    gaussian = "gaussian" 

###############################################################################

class HistogramCollectionPainter:
    def __init__(self):
        pass

    def paint(self, hists, *options, **kwargs):
        ax = None
        if "ax" in kwargs:
            ax = kwargs["ax"]
        if not drawtools._checkAllOptionsAreValid(options):
            raise Exception("Invalid option provided", options)
        #clone input histograms
        hclone = OrderedDict()
        for k, v in hists.iteritems():
            hclone[k] = v.clone()
        hists = hclone
        #draw histograms
        hists = self._showoverflow(hists, options)
        hists = self._dividebinwidth(hists, options)
        hists = self._normalisehistograms(hists, options)
        hists = self._buildstack(hists, options)
        hists, fig = self._drawhistograms(hists, options, ax=ax)
        return fig

    def _findoption(self, opttype, options, default=None):
        return drawtools._findOption(opttype, options, default)

    def _normalisehistograms(self, hists, options):
        norm = self._findoption(drawoptions.Normalisation, options, default=drawoptions.Normalisation(drawoptions.Normalisation.noNormalisation))
        treatAsData = self._findoption(drawoptions.TreatAsData, options, default=drawoptions.TreatAsData())
        nt = _NormalisationTool()
        mode = norm.getMode()
        hists = nt.normalise(hists, mode, treatAsData)
        return hists

    def _buildstack(self, hists, options):
        stack = self._findoption(drawoptions.Stack, options, default=drawoptions.Stack(drawoptions.Stack.overlay))
        treatasdata = self._findoption(drawoptions.TreatAsData, options, default=drawoptions.TreatAsData())
        st = _StackTool()
        hists = st.stack(hists, stack, treatasdata)
        return hists

    def _drawhistograms(self, hists, options, ax=None):
        formatter = self._findoption(drawoptions.Stack, options, default=None)
        treatasdata = self._findoption(drawoptions.TreatAsData, options, default=drawoptions.TreatAsData())
        canvassize = self._findoption(drawoptions.CanvasSize, options, default=None)
        axislabels = self._findoption(drawoptions.AxisLabels, options, default=None)
        binlabels = self._findoption(drawoptions.BinLabels, options, default=None)
        axisscale = self._findoption(drawoptions.AxisScale, options, default=None)
        framerange = self._findoption(drawoptions.FrameRange, options, default=None)
        legendposition = self._findoption(drawoptions.LegendPosition, options, default=None)
        dt = _DrawTool()
        return dt.draw(hists, treatasdata, formatter, canvassize, axislabels, binlabels, axisscale, framerange, legendposition, ax=ax)

    def _showoverflow(self, hists, options):
        showoverflow = self._findoption(drawoptions.ShowOverflow, options, default=None)
        if showoverflow is not None and showoverflow.flag:
            for h in hists.itervalues():
                #merge overflow with the last bin
                h.values[-1] += h.overflow
                h.sumw2[-1] += h.overflowsumw2
                h.overflow = 0.0
                h.overflowsumw2 = 0.0
        return hists

    def _dividebinwidth(self, hists, options):
        dbw = self._findoption(drawoptions.DivideByBinWidth, options, default=None)
        if dbw is not None and dbw.flag:
            for h in hists.itervalues():
                bw = _get_bin_widths(h.xbinning)
                h.values = h.values / bw
                h.scale(dbw.scale)
        return hists

###############################################################################

class _NormalisationTool:
    def __init__(self):
        pass

    def normalise(self, hists, mode, treatasdata):
        if mode == drawoptions.Normalisation.noNormalisation:
            #do nothing
            pass
        elif mode == drawoptions.Normalisation.mcToData:
            hists = self._normmctodata(hists, treatasdata)
        elif mode == drawoptions.Normalisation.totalUnitArea:
            hists = self._normtotalarea(hists, treatasdata)
        elif mode == drawoptions.Normalisation.unitArea:
            hists = self._normunitarea(hists)
        return hists

    def _normmctodata(self, hists, treatasdata):
        mc = self._summc(hists, treatasdata)
        data = self._sumdata(hists, treatasdata)
        try:
            scale = data / mc
        except ZeroDivisionError:
            scale = 1.0
        result = {}
        for k, h in hists.iteritems():
            if self._ismc(k, treatasdata):
                h.scale(scale)
        return hists

    def _normtotalarea(self, hists, treatasdata):
        mc = self._summc(hists, treatasdata)
        data = self._sumdata(hists, treatasdata)
        result = {}
        for k, h in hists.iteritems():
            if self._ismc(k, treatasdata):
                norm = mc
            else:
                norm = data
            try:
                scale = 1.0 / norm
            except ZeroDivisionError:
                scale = 1.0
            h.scale(scale)
        return hists        

    def _normunitarea(self, hists):
        for k, h in hists.iteritems():
            try:
                scale = 1.0 / h.sum()
            except ZeroDivisionError:
                scale = 1.0
            h.scale(scale)
        return hists

    def _summc(self, hists, treatasdata):
        return sum(h.sum() for k, h in hists.iteritems() if self._ismc(k, treatasdata))

    def _sumdata(self, hists, treatasdata):
        return sum(h.sum() for k, h in hists.iteritems() if self._isdata(k, treatasdata))

    def _ismc(self, key, treatasdata):
        return not key in treatasdata

    def _isdata(self, key, treatasdata):
        return not self._ismc(key, treatasdata)

###############################################################################

class _StackTool:
    def __init__(self):
        pass

    def stack(self, hists, stackopt, treatasdata):
        if stackopt.getStackMode() == drawoptions.Stack.stackMC:
            hists = self._stackmc(hists, stackopt, treatasdata)
        return hists

    def _stackmc(self, hists, stackopt, treatasdata):
        result = OrderedDict()
        order = self._getorder(hists, stackopt)
        order.reverse() # start from bottom,
        #add MC histograms to object
        for ii in xrange(len(order)): 
            k = order[ii]
            if not k in treatasdata:
                h = hists[k].clone()
                #add histograms underneath
                for jj in xrange(0, ii):
                    if not order[jj] in treatasdata:
                        h.add(hists[order[jj]])
                result[k] = h
        #add data histograms to object
        for k in order: 
            if k in treatasdata:
                result[k] = hists[k]
        return result

    def _getorder(self, hists, stackopt):
        order = []
        #add user defined order
        order.extend(stackopt.getUserOrder())
        keys = hists.keys()
        if stackopt.getOrderMode() == drawoptions.Stack.orderByIntegral:
            keys.sort(key=lambda x: hists[x].sum())
        elif stackopt.getOrderMode() == drawoptions.Stack.orderByInsertion:
            #do nothing, key should already be sorted
            pass
        else:
            #sort alphabetically by default
            keys.sort()
        #add remaining histograms
        for k in keys:
            if k not in order:
                order.append(k)
        return order
        

###############################################################################

class _DrawTool:
    def __init__(self, ):
        pass

    def draw(self, hists, treatasdata, formatter, canvassize, axislabels, binlabels, axisscale, framerange, legendpos, ax=None):
        if ax is None:
            fig, ax = self._build_figure(canvassize)
        else:
            fig = None
        self._draw_hists(ax, hists, treatasdata)
        self._label_axes(ax, hists, axislabels, binlabels)
        self._draw_legend(ax, legendpos)
        self._setscale(ax, axisscale)
        #re-optimise the figure layout once everything has been drawn
        self._setframerange(ax, framerange)
        if fig:
            fig.tight_layout()
        return hists, fig

#    def _draw_legend(self, ax, legendpos):
#        if legendpos is not None:
#            xlow, ylow, xhigh, yhigh = legendpos.calculateLegendLimits()
#            #transform to axes coordinates
#            xl, xh = ax.get_xlim()
#            yl, yh = ax.get_ylim()
#            xlow = xl + float(xh-xl) * xlow
#            xhigh = xl + float(xh-xl) * xhigh
#            ylow = yl + float(yh-yl) * ylow
#            yhigh = yl + float(yh-yl) * yhigh
#            #set legend co-cordinates
#            loc = (xlow, ylow)
#            #no way to manually set size?
#            #we can only set the lower left coordinate
#            ax.legend(loc=loc)
#        else:
#            ax.legend()
#        return

    def _draw_legend(self, ax, legendpos):
        if legendpos is not None:
            xlow, ylow, xhigh, yhigh = legendpos.calculateLegendLimits()
            mid = 0.5
            if xlow > mid and ylow > mid:
                loc = 1 # upper right
            elif xlow < mid and ylow > mid:
                loc = 2 # upper left
            elif xlow < mid and ylow < mid:
                loc = 3 # lower left
            elif xlow > mid and ylow < mid:
                loc = 4 # lower right
            else:
                loc = 0
            ax.legend(loc=loc)
        else:
            ax.legend()
        return


    def _setframerange(self, ax, framerange):
        if framerange is not None:
            xlow, ylow, xhigh, yhigh = framerange.xMin, framerange.yMin, framerange.xMax, framerange.yMax
            start_xlow, start_xhigh = ax.get_xlim()
            start_ylow, start_yhigh = ax.get_ylim()
            auto = drawoptions.FrameRange.auto
            if xlow == auto:
                xlow = start_xlow
            if ylow == auto:
                ylow = start_ylow
            if xhigh == auto:
                xhigh = start_xhigh
            if yhigh == auto:
                yhigh = start_yhigh
            ax.set_xlim(xlow, xhigh)
            ax.set_ylim(ylow, yhigh)
        return

    def _setscale(self, ax, axisscale):
        if axisscale is not None:
            if axisscale.logX:
                ax.set_xscale("log")
            if axisscale.logY:
                ax.set_yscale("log")
        return

    def _build_figure(self, canvassize):
        fig = plt.figure()
        figsize = None
        dpi = fig.get_dpi()
        if canvassize is not None:
            figsize = np.array([canvassize.x, canvassize.y], dtype=float)
            figsizeinches = figsize / dpi
            fig.set_size_inches(figsizeinches)
        ax = fig.add_subplot(111)
        return fig, ax

    def _draw_hists(self, ax, hists, treatasdata):
        for k, h  in hists.iteritems():
            label = k
            if k in treatasdata:
                plot_hist_points(ax, xbinning=h.xbinning, y=h.values, yerr=h.errors, label=label)
            else:
                plot_hist_lines(ax, xbinning=h.xbinning, y=h.values, label=label)
        return

    def _label_axes(self, ax, hists, axislabels, binlabels):
        x, y, z = None, None, None
        if axislabels is not None:
            x = axislabels.x
            y = axislabels.y
            z = axislabels.z
        else:
            l = hists.itervalues().next().label.axislabels
            if len(l) > 0:
                x = l[0]
            if len(l) > 1:
                y = l[1]
            if len(l) > 2:
                z = l[2]
        if x:
            ax.set_xlabel(x)
        if y:
            ax.set_ylabel(y)
        #if z:
        #    ax.set_zlabel(z)
        self._set_bin_labels(ax, hists, binlabels)
        return

    def _set_bin_labels(self, ax, hists, binlabels):
        if binlabels is not None:
            xlabels = binlabels.xlabels
            h = hists.itervalues().next()
            if xlabels == drawoptions.BinLabels.auto:
                #get labels from the first histogram
                xlabels = h.label.axislabels[0]
            if xlabels is not None:
                xbinning = h.xbinning
                xbc = _get_bin_centres(xbinning)
                ax.set_xticks(xbc)
                ax.set_xticklabels(xlabels, rotation=90)
        return


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
    b = get_binning_array(xbinning)
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
    binedges = get_binning_array(xbinning)
    x = _get_bin_centres(binedges)
    weights = y
    return ax.hist(x, bins=binedges, weights=weights,
                   *args, **kwargs)

###############################################################################

def plot_hist_lines(ax, xbinning, y, *args, **kwargs):
    binedges = get_binning_array(xbinning)
    xline = []
    yline = []
    for ii in xrange(len(binedges)):
            eb = binedges[ii]
            xline.append(eb)
            xline.append(eb)
            if ii > 0:
                low = y[ii-1]
            else:
                low = 0.0
            if ii < len(y):
                high = y[ii]
            else:
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
