from argparse import ArgumentParser
from collections import OrderedDict
import glob
import os

import ROOT

from simplot.rootplot.drawtools import MultiTreePainter, TreePainter, HistogramCollectionPainter
from simplot.rootplot.drawoptions import *

def parsecml():
    parser = ArgumentParser()
    parser.add_argument("-p", "--plot", help="Plot command.", type=str, default="1")
    parser.add_argument("-c", "--cut", help="Cut command.", type=str, default=None)
    parser.add_argument("-s", "--split", help="Split command in the form \"label:cut,label:cut\".", type=str, default=None)
    parser.add_argument("-n", "--norm", help="Normalisation mode", choices=Normalisation.ALL_MODES, default=Normalisation.noNormalisation)
    parser.add_argument("-b", "--binning", help="Binning in the form \"Nbins,xlow,xhigh\".", default=None)
    parser.add_argument("-t", "--tree", help="Name of input tree.", type=str, default="ntuple")
    parser.add_argument("-r", "--root-draw-opt", help="ROOT TH1 draw option", type=str, default=None)
    parser.add_argument("-o", "--output", help="Output file.", type=str, default=None)
    parser.add_argument("--labels", help="A comma separated list of dataset labels", type=str, default=None)
    parser.add_argument("input_files", help="A list of input files.", type=str, nargs="+")
    return parser.parse_args()

def loadtree(inputfiles, treename):
    chain = ROOT.TChain(treename)
    for fpattern in inputfiles:
        for f in glob.glob(fpattern):
            chain.AddFile(f)
    return chain

def str_to_dict(string):
    result = OrderedDict()
    for s in string.split(","):
        key, val = s.split(":")
        result[key] = val
    return result

def builddrawoptions(args):
    result = []
    if args.cut:
        print "INFO: cut: ", args.cut
        result.append(Cut(args.cut))
    if args.split:
        split = SplitDataset()
        for name, cut in str_to_dict(args.split).iteritems():
            print "INFO: split:", name, cut
            split.add(name, cut)
        result.append(split)
    if args.norm:
        print "INFO: normalisation:", args.norm
        result.append(Normalisation(args.norm))
    if args.binning:
        nbins, xlow, xhigh = args.binning.split(",")
        nbins, xlow, xhigh = int(nbins), float(xlow), float(xhigh)
        print "INFO: binning", nbins, xlow, xhigh
        result.append(UniformBinning(nbins, xlow, xhigh))
    if args.root_draw_opt:
        print "INFO: ROOT draw option", args.root_draw_opt
        result.append(HistogramDrawOption(str_to_dict(args.root_draw_opt)))
    return result

def plot(datasets, treename, plotcmd, drawoptions, output=None):
    hist = _makehistograms(datasets, treename, plotcmd, drawoptions)
    _makeplot(hist, drawoptions, output)
    return

def _makehistograms(datasets, treename, plotcmd, drawoptions):
    name = "plot"
    painter = MultiTreePainter(datasets, treeName=treename)
    hist = painter.makeHistograms(name, plotcmd, plotcmd, *drawoptions)
    hist._painter = painter # prevent automatic deletion of objects
    return hist

def _makeplot(hist, drawoptions, output):
    histpainter = HistogramCollectionPainter()
    canvas = histpainter.paint(hist, *drawoptions)
    if output:
        canvas.SaveAs(output)
    else:
        raw_input("waiting")
    return

def expandfilelist(pattern):
    return glob.glob(
        os.path.expandvars(
            os.path.expanduser(pattern)
        )
    )

def splitdatasets(args):
    if args.labels:
        datasets = {label : expandfilelist(fpattern) for label, fpattern in zip(args.labels.split(","), args.input_files)}
    else:
        flist = []
        for fname in args.input_files:
            flist.extend(expandfilelist(fname))
        datasets = {"" : flist }
    return datasets

def main():
    args = parsecml()
    if args.output:
        ROOT.gROOT.SetBatch(True)
    datasets = splitdatasets(args)
    drawoptions = builddrawoptions(args)
    plot(datasets, args.tree, args.plot, drawoptions, output=args.output)
    return

if __name__ == "__main__":
    main()
