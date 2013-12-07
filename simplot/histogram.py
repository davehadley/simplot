import numpy
import math

###############################################################################

class Stats:
    poisson = "poisson"
    gaussian = "gaussian" 

###############################################################################

def _get_binning_array(xbinning):
    try:
        b = xbinning.binedges
    except AttributeError:
        #object does not have binedges attribute, 
        #assume that we were given an array-like of bin edges. 
        b = xbinning
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
            yerr = numpy.sqrt(y) 
        #if no yerr is given assume poisson statistics
        elif yerrmode == Stats.poisson:
            #TODO implement Poisson errors
            raise NotImplementedError()
    if not show_empty_bins:
        #filter out bins with zero content
        filtered = filter(lambda v: v[2]!=0.0, zip(x, xerr, y, yerr))
        x, xerr, y, yerr = map(numpy.array, zip(*filtered))
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

def _test_plot_hist():
    import numpy
    from nuOscillation.model import binning
    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    dataset = numpy.random.normal(size=100)
    xbinning = binning.Binning(10, -5.0, 5.0)
    y, binedges = numpy.histogram(dataset, xbinning.binedges)
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
