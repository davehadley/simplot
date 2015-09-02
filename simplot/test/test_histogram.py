import unittest
from simplot.mplot.binning import Binning
from simplot.mplot.histogram import HistogramND, HistogramNDLabel
from simplot.mplot.histpaint import HistogramCollectionPainter
from simplot.rootplot.drawoptions import *
from simplot.mplot.io import FigureWriter
import itertools
import numpy as np

_DEFAULT_BINNING = [Binning(10, -5.0, 5.0),
                    Binning(10, 0.0, 5.0),
                    Binning(10, 5.0, 10.0),
]

class TestHistogram(unittest.TestCase):
    def test_fillhist(self):
        b = _DEFAULT_BINNING
        hist = HistogramND(b)
        #fill histogram
        for binnum in itertools.product(*b):
            coord = [b[i].bincentre(n) for i, n in enumerate(binnum)]
            hist.fill(coord, sum(binnum))
        #check individual bins
        for binnum in itertools.product(*b):
            coord = [b[i].bincentre(n) for i, n in enumerate(binnum)]
            value = hist.eval(coord)
            self.assertEquals(sum(binnum), value)
        return

    def test_fillhist(self):
        b = _DEFAULT_BINNING
        hist = HistogramND(b)
        #fill histogram
        for binnum in itertools.product(*b):
            coord = [b[i].bincentre(n) for i, n in enumerate(binnum)]
            index = hist.binindex(coord)
            hist.setbin(index, sum(binnum))
        #check individual bins
        for binnum in itertools.product(*b):
            coord = [b[i].bincentre(n) for i, n in enumerate(binnum)]
            value = hist.eval(coord)
            self.assertEquals(sum(binnum), value)
        return

class TestPaintHistogram(unittest.TestCase):
    def setUp(self):
        #create some dummy histograms
        hists = {}
        binning = Binning(25, -5.0, 5.0)
        for ii in xrange(5):
            k = "mc_%s" % ii
            h = HistogramND(binning)
            for _ in xrange(10000):
                h.fill(np.random.normal(loc=ii))
            hists[k] = h
        h = HistogramND(binning)
        for _ in xrange(1000):
            h.fill(np.random.normal(loc=np.random.randint(5)))
        hists["data"] = h
        self._hists = hists
        return

    def test_plotnorm(self):
        h = self._hists
        p = HistogramCollectionPainter()
        for mode in [Normalisation.mcToData, Normalisation.noNormalisation, Normalisation.unitArea, Normalisation.totalUnitArea]:
            name = "test_mode_plot_norm_%s" % mode
            fig = p.paint(h, Normalisation(mode), TreatAsData("data"))
            out = FigureWriter()
            out(name, fig)
        return

    def test_plotstack(self):
        h = self._hists
        p = HistogramCollectionPainter()
        for mode in [Stack.overlay, Stack.stackMC]:
            name = "test_mode_plot_stack_%s" % mode
            fig = p.paint(h, Normalisation(Normalisation.mcToData), Stack(mode), TreatAsData("data"))
            out = FigureWriter()
            out(name, fig)
        for mode in [None, Stack.orderByInsertion, Stack.orderByIntegral]:
            name = "test_mode_plot_stackorder_%s" % mode
            fig = p.paint(h, Normalisation(Normalisation.mcToData), Stack(Stack.stackMC, orderMode=mode), TreatAsData("data"))
            out = FigureWriter()
            out(name, fig)
        return

    def test_axislabels(self):
        h = self._hists
        p = HistogramCollectionPainter()
        for name, axislabels in [("nolabels", AxisLabels()), ("withlabels", AxisLabels(x="x label", y="y label", z="z label"))]:
            name = "test_mode_plot_axislabel_%s" % name
            fig = p.paint(h, axislabels, TreatAsData("data"))
            out = FigureWriter()
            out(name, fig)
        return

    def test_canvassize(self):
        h = self._hists
        p = HistogramCollectionPainter()
        for name, axislabels in [("default", CanvasSize()), ("long", CanvasSize(1200, 300))]:
            name = "test_mode_plot_canvassize_%s" % name
            fig = p.paint(h, axislabels, TreatAsData("data"))
            out = FigureWriter()
            out(name, fig)
        return

    def test_binlabels(self):
        h = self._hists
        p = HistogramCollectionPainter()
        labels = ["x=%.2f" % low for low in h.values()[0].xbinning[:-1]]
        for name, binlabels in [("nobinlabels", BinLabels()), ("withbinlabels", BinLabels(labels))]:
            name = "test_mode_plot_binlabels_%s" % name
            fig = p.paint(h, binlabels, TreatAsData("data"))
            out = FigureWriter()
            out(name, fig)
        return

    def test_axisscale(self):
        h = self._hists
        p = HistogramCollectionPainter()
        labels = ["x=%.2f" % low for low in h.values()[0].xbinning[:-1]]
        for name, axisscale in [("logY", AxisScale(logY=True)), ("logXlogY", AxisScale(logX=True, logY=True))]:
            name = "test_mode_plot_axisscale_%s" % name
            fig = p.paint(h, axisscale, TreatAsData("data"))
            out = FigureWriter()
            out(name, fig)
        return

    def test_framerange(self):
        h = self._hists
        p = HistogramCollectionPainter()
        labels = ["x=%.2f" % low for low in h.values()[0].xbinning[:-1]]
        for name, framerange in [("auto", FrameRange()), ("xy", FrameRange(xMin=-3, xMax=3., yMin=-10.0, yMax=1000.0)), ("x", FrameRange(xMin=-3., xMax=3.))]:
            name = "test_mode_plot_framerange_%s" % name
            fig = p.paint(h, framerange, TreatAsData("data"))
            out = FigureWriter()
            out(name, fig)
        return

    def test_legendpos(self):
        h = self._hists
        p = HistogramCollectionPainter()
        labels = ["x=%.2f" % low for low in h.values()[0].xbinning[:-1]]
        for name, legendpos in [("topright", LegendPosition(x=0.6, y=0.6)),
                                ("topleft", LegendPosition(x=0.0, y=1.0)),
                                ("bottomleft", LegendPosition(x=0.0, y=0.0)),
                                ("bottomright", LegendPosition(x=1.0, y=0.0)),
        ]:
            name = "test_mode_plot_legendpos_%s" % name
            fig = p.paint(h, legendpos, TreatAsData("data"))
            out = FigureWriter()
            out(name, fig)
        return

    def test_showoverflow(self):
        h = self._hists
        p = HistogramCollectionPainter()
        fig = p.paint(h, ShowOverflow(True), TreatAsData("data"))
        out = FigureWriter()
        name = "test_mode_plot_showoverflow"
        out(name, fig)
        return

            


def main():
    unittest.main()
    return

if __name__ == "__main__":
    main()
