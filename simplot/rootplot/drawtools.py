'''Provides classes for filling and drawing histograms.

Note classes for specifying draw options are in the drawoptions module.

The main classes are:

* HistogramCollection is a container for all histograms that contribute to a plot
  for example it can contain a histogram for each event type.
* HistogramCollectionPainter formats and draws a HistogramCollection to a canvas.
* TreePainter fills and draws histograms from a TTree using the ROOT.TTree.Draw interface.

The easiest way to make a plot is to use TreePainter to draw directly from an 
ntuple. 

Features still to be implemented:

* POT weighting.
* Event weighting.

'''

import ROOT
import drawoptions
import math
import itertools
import pprint
import StringIO
import uuid
import copy
import sys
import numpy

try:
    from collections import OrderedDict
except:
    #not availble in some old versions of python.
    OrderedDict = dict

###############################################################################

def _findOption(optionType, options, default=None):
    """Given a container of drawoptions, retrieve the first one of/inheriting 
    from a specific class.
    
    :param optionType: the type or class of the requested draw option.
    :type optionType: classobj
    :param options: container of draw options
    :type options: iterable
    :param default: the default value to return if no match is found.
    :returns: first draw option of type optionType in options, or default. 
    """
    result = default
    for opt in options:
        if isinstance(opt,optionType):
            result = opt
            break
    return result

###############################################################################

def _findAllOptions(optionType, options, default=[]):
    """Given a container of drawoptions, retrieve the first one of/inheriting 
    from a specific class.
    
    :param optionType: the type or class of the requested draw option.
    :type optionType: classobj
    :param options: container of draw options
    :type options: iterable
    :param default: the default value to return if no match is found.
    :returns: first draw option of type optionType in options, or default. 
    """
    result = []
    for opt in options:
        if isinstance(opt,optionType):
            result.append(opt)
    if len(result)==0:
        result = default
    return result

###############################################################################

def _checkAllOptionsAreValid(options):
    success = True
    for opt in options:
        match = False
        for optionsType in drawoptions._allValidOptionsTypes:
            if isinstance(opt,optionsType):
                match = True
                break
        success = success and match
    return success

###############################################################################

class TreePainter:
    """Easily makes plots from flat TTrees.
    
    Uses the ROOT.TTree.Draw language to create histograms from an ntuple. 
    These are drawn to a TCanvas using HistogramCollectionPainter.
    
    For example, if you have an ntuple containing the variable "muonMomentum"
    you can plot it easily,
    
    >>> tp = TreePainter(myInputDataset,"MyAnalysisNtuple")
    >>> canvas = tp.paint("muonMomentumPlot","plot of muon momentum", "muonMomentum")
    
    The drawoptions module contains an extensive list of draw options to control
    and format the output plot.  
    
    """
    def __init__(self, inputDataset=None, treeName=None, tree=None):
        """
        :param inputDataset: dataset object containing the input ntuple files.
        :type inputDataset: analysis.datasets.Dataset
        :param treeName: name of input ntuple TTree.
        :type treeName: str 
        """
        if inputDataset is not None and treeName is not None:
            self.treeName = treeName
            self.inputDataset = inputDataset
            self.tree = self._getChain()
        elif tree is not None:
            self.treeName = tree.GetName()
            self.inputDataset = inputDataset
            self.tree = tree
        else:
            raise ValueError("TreePainter initialised with invalid arguments. Either inputDataset and treeName must be set or a tree must be provided.", inputDataset, treeName, tree)
        #get a unique name for the canvas (avoids warnings about duplicate canvas)
        canvName = str(uuid.uuid4()).replace("-","_")
        self.dummyCanv = ROOT.TCanvas(canvName,"",800,600)
           
    def paint(self, name, title, command, *options):
        """Create histograms and paint them.
        
        :param name: the name of the returned canvas.
        :type name: str
        :param title: the title of the returned canvas.
        :type title: str
        :param command: the draw command (in the format of the ROOT.TTree.Draw mini-language).
        :type command: str
        :param options: optional number of draw options (see drawoptions module).
        :returns: the plot (ROOT.TCanvas).
        """
        _checkAllOptionsAreValid(options)
        self.histCol = self._generateHistograms(name, title, command,*options)
        self.histPainter = HistogramCollectionPainter()
        self.canv = self.histPainter.paint(self.histCol, *options)
        return self.canv
    
    def paintLater(self, name, title, command, *options):
        """Returns an iterable over the plotted canvas. 
        paint will only called once you iterate. This is useful if you want to 
        delay a cpu intensive command until the end of the script, or to ensure 
        that memory intensive draw commands only stay in memory for as short a 
        time as possible.
        
        :param name: the name of the returned canvas.
        :type name: str
        :param title: the title of the returned canvas.
        :type title: str
        :param command: the draw command (in the format of the ROOT.TTree.Draw mini-language).
        :type command: str
        :param options: optional number of draw options (see drawoptions module).
        :returns: an iterable over ROOT.TCanvas.
        """
        yield self.paint(name, title, command, *options)
    
    def makeHistograms(self, name, title, command, *options):
        """Create histograms but does not draw them.
        
        :param name: the name of the returned canvas.
        :type name: str
        :param title: the title of the returned canvas.
        :type title: str
        :param command: the draw command (in the format of the ROOT.TTree.Draw mini-language).
        :type command: str
        :param options: optional number of draw options (see drawoptions module).
        :returns: the histogram collection.
        """
        _checkAllOptionsAreValid(options)
        return self._generateHistograms(name, title, command,*options)
        
    def _getChain(self):
        chain = ROOT.TChain(self.treeName,self.treeName)
        for fileName in self.inputDataset:
            try:
                chain.AddFile(fileName)
            except TypeError:
                #could have been given a TFile
                chain.AddFile(fileName.GetName())
        return chain
    
    def _generateHistograms(self, name, title, command, *options):
        nDimensions = _findOption(drawoptions.NDimensions, options, default=drawoptions.NDimensions(1))
        if nDimensions.is1D():
            result = self._generateHistograms1D(name, title, command, *options)
        elif nDimensions.is2D():
            result = self._generateHistograms2D(name, title, command, *options)
        elif nDimensions.is3D():
            result = self._generateHistograms3D(name, title, command, *options)
        else:
            raise Exception("plot dimensions not implemented",nDimensions)
        return result
    
    def _generateHistograms1D(self, name, title, command, *options):
        inputName = name
        inputTitle = title
        histCol = HistogramCollection(name, title)
        
        binning = _findOption(drawoptions.BaseBinning, options, default=None)
        splitDataset = _findOption(drawoptions.SplitDataset, options, default=drawoptions.SplitDataset({"all":drawoptions.Cut("")}))
        #userCut = _findOption(drawoptions.Cut, options, default=drawoptions.Cut(""))
        userCutList = _findAllOptions(drawoptions.Cut, options, default=[drawoptions.Cut("")])
        userCut = drawoptions.Cut.mergeListOfCuts(userCutList)
        eventWeight = _findOption(drawoptions.EventWeight, options, default=drawoptions.EventWeight(""))
        profile = _findOption(drawoptions.Profile, options, default=drawoptions.Profile(False))        
        for name,splitCut in splitDataset:
            histogram = None
            #histName = "hDummy"
            histName = "hist_" + str(uuid.uuid4()).replace("-","_")
            if binning and binning.getBinningArray():
                nBins = binning.getNBins()
                array = binning.getBinningArray()
                if profile.flag:
                    histogram = ROOT.TProfile(histName, name, nBins, array, profile.erroroption)
                else:
                    histogram = ROOT.TH1D(histName, name,nBins,array)
            #get cut string
            fullCut = None
            scs = splitCut.getCutString()
            ucs = userCut.getCutString()
            if len(scs) and len(ucs):
                ucs = "("+ucs+")"
                scs = "("+scs+")"
            if len(scs) and len(ucs):
                fullCut = " && ".join([scs,ucs])
                #print "DEBUG ",fullCut
            elif len(scs):
                fullCut = scs
            elif len(ucs):
                fullCut = ucs
            else:
                fullCut = ""
            #add event weight
            ews = eventWeight.getEventWeightString()
            if len(fullCut) and len(ews):
                fullCut = "("+ews+")*("+fullCut+")"
            elif len(ews):
                fullCut = ews
            #get option string
            fullOption = ""
            if profile.flag:
                fullOption += "prof" + profile.erroroption
            #create histogram by calling the tree Draw command
            fullCommand = command + ">>" + histName
            if isinstance(binning, drawoptions.AutoBinning):
                fullCommand += "({0})".format(binning.getNBins())
            self.dummyCanv.cd()
            self.tree.Draw(fullCommand, fullCut, fullOption)
            if not histogram:
                #retrieve the histogram
                histogram = ROOT.gPad.GetPrimitive(histName)
                if not histogram:
                    raise Exception("TreePainter error : failed to draw histogram. Check that your draw command is valid.", name, title, command, histName, fullCommand, fullCut)
                histogram = histogram.Clone(histName)
            #store the histogram in histogram collection
            histogram.SetDirectory(0)
            histogram.Sumw2()
            histCol.addHistogram(name,histogram)
            if (not binning) or isinstance(binning, drawoptions.AutoBinning):
                #binning not set, try to automatically set it from this
                nBins = histogram.GetNbinsX()
                boundaries = numpy.zeros(nBins+1)
                histogram.GetXaxis().GetLowEdge(boundaries)
                boundaries[nBins] = histogram.GetXaxis().GetXmax()
                rooBinning = ROOT.RooBinning(nBins,boundaries)
                binning = drawoptions.RooBinningWrapper(rooBinning)
        return histCol

    def _generateHistograms2D(self, name, title, command, *options):    
        xBinning = _findOption(drawoptions.BaseBinning, options, default=None)
        yBinningWrapper = _findOption(drawoptions.YBinning, options, default=None)
        histCol = HistogramCollection2D(name, title, xBinning, yBinningWrapper.getBinningObject())
        splitDataset = _findOption(drawoptions.SplitDataset, options, default=drawoptions.SplitDataset({"all":drawoptions.Cut("")}))
        #userCut = _findOption(drawoptions.Cut, options, default=drawoptions.Cut(""))
        userCutList = _findAllOptions(drawoptions.Cut, options, default=[drawoptions.Cut("")])
        userCut = drawoptions.Cut.mergeListOfCuts(userCutList)
        eventWeight = _findOption(drawoptions.EventWeight, options, default=drawoptions.EventWeight(""))
        #sanitise options
        if not (xBinning and yBinningWrapper):
            raise Exception("drawing 2D histograms without specifying x and y binning is not supported.", xBinning, yBinningWrapper)
        #make histograms
        for name,splitCut in splitDataset:
            histogram = None
            #histName = "hDummy"
            histName = "hist_" + str(uuid.uuid4()).replace("-","_") 
            if xBinning and yBinningWrapper:
                yBinning = yBinningWrapper.getBinningObject()
                nBinsX = xBinning.getNBins()
                arrayX = xBinning.getBinningArray()
                nBinsY = yBinning.getNBins()
                arrayY = yBinning.getBinningArray()
                histogram = ROOT.TH2F(histName,title,nBinsX,arrayX,nBinsY,arrayY)
            #get cut string
            fullCut = None
            scs = splitCut.getCutString()
            ucs = userCut.getCutString()
            if len(scs) and len(ucs):
                ucs = "("+ucs+")"
                scs = "("+scs+")"
            if len(scs) and len(ucs):
                fullCut = " && ".join([scs,ucs])
                #print "DEBUG ",fullCut
            elif len(scs):
                fullCut = scs
            elif len(ucs):
                fullCut = ucs
            else:
                fullCut = ""
            #add event weight
            ews = eventWeight.getEventWeightString()
            if len(fullCut) and len(ews):
                fullCut = "("+ews+")*("+fullCut+")"
            elif len(ews):
                fullCut = ews
            #print "DEBUG",fullCut
            #create histogram by calling the tree Draw command
            fullCommand = command + ">>" + histName
            self.dummyCanv.cd()
            self.tree.Draw(fullCommand, fullCut)
            if not histogram:
                #retrieve the histogram
                histogram = ROOT.gPad.GetPrimitive(histName)
                histogram = histogram.Clone(histName)
            #store the histogram in histogram collection
            histogram.SetDirectory(0)
            histogram.Sumw2()
            histCol.addHistogram(name,histogram)
        return histCol

    def _generateHistograms3D(self, name, title, command, *options):    
        xBinning = _findOption(drawoptions.BaseBinning, options, default=None)
        yBinningWrapper = _findOption(drawoptions.YBinning, options, default=None)
        zBinningWrapper = _findOption(drawoptions.ZBinning, options, default=None)
        histCol = HistogramCollection3D(name, title, xBinning, yBinningWrapper.getBinningObject(), zBinningWrapper.getBinningObject())
        splitDataset = _findOption(drawoptions.SplitDataset, options, default=drawoptions.SplitDataset({"all":drawoptions.Cut("")}))
        #userCut = _findOption(drawoptions.Cut, options, default=drawoptions.Cut(""))
        userCutList = _findAllOptions(drawoptions.Cut, options, default=[drawoptions.Cut("")])
        userCut = drawoptions.Cut.mergeListOfCuts(userCutList)
        eventWeight = _findOption(drawoptions.EventWeight, options, default=drawoptions.EventWeight(""))
        #sanitise options
        if not (xBinning and yBinningWrapper and zBinningWrapper):
            raise Exception("drawing 3D histograms without specifying x and y and z binning is not supported.", xBinning, yBinningWrapper, zBinningWrapper)
        #make histograms
        for name,splitCut in splitDataset:
            histogram = None
            #histName = "hDummy"
            histName = "hist_" + str(uuid.uuid4()).replace("-","_") 
            if xBinning and yBinningWrapper and zBinningWrapper:
                yBinning = yBinningWrapper.getBinningObject()
                zBinning = zBinningWrapper.getBinningObject()
                nBinsX = xBinning.getNBins()
                arrayX = xBinning.getBinningArray()
                nBinsY = yBinning.getNBins()
                arrayY = yBinning.getBinningArray()
                nBinsZ = zBinning.getNBins()
                arrayZ = zBinning.getBinningArray()
                histogram = ROOT.TH3F(histName,title,nBinsX,arrayX,nBinsY,arrayY,nBinsZ,arrayZ)
            #get cut string
            fullCut = None
            scs = splitCut.getCutString()
            ucs = userCut.getCutString()
            if len(scs) and len(ucs):
                ucs = "("+ucs+")"
                scs = "("+scs+")"
            if len(scs) and len(ucs):
                fullCut = " && ".join([scs,ucs])
                #print "DEBUG ",fullCut
            elif len(scs):
                fullCut = scs
            elif len(ucs):
                fullCut = ucs
            else:
                fullCut = ""
            #add event weight
            ews = eventWeight.getEventWeightString()
            if len(fullCut) and len(ews):
                fullCut = "("+ews+")*("+fullCut+")"
            elif len(ews):
                fullCut = ews
            #print "DEBUG",fullCut
            #create histogram by calling the tree Draw command
            fullCommand = command + ">>" + histName
            self.dummyCanv.cd()
            self.tree.Draw(fullCommand, fullCut)
            if not histogram:
                #retrieve the histogram
                histogram = ROOT.gPad.GetPrimitive(histName)
                histogram = histogram.Clone(histName)
            #store the histogram in histogram collection
            histogram.SetDirectory(0)
            histogram.Sumw2()
            histCol.addHistogram(name,histogram)
        return histCol

###############################################################################

class MultiTreePainter:
    def __init__(self, datasets, treeName=None):
        """
        :param datasets: A dictionary-like with key=dataset name and value=tree or filename or filelist.
        :type datasets: dict
        """
        self._datasets = datasets
        self._painters = {}
        for label, value in datasets.iteritems():
            if isinstance(value, basestring):
                # single input file
                self._painters[label] = TreePainter([value], treeName=treeName)
            else:
                self._painters[label] = TreePainter(value, treeName=treeName)

    def paint(self, name, title, command, *options):
        """Create histograms and paint them.

        :param name: the name of the returned canvas.
        :type name: str
        :param title: the title of the returned canvas.
        :type title: str
        :param command: the draw command (in the format of the ROOT.TTree.Draw mini-language).
        :type command: str
        :param options: optional number of draw options (see drawoptions module).
        :returns: the plot (ROOT.TCanvas).
        """
        _checkAllOptionsAreValid(options)
        self.histCol = self._generateHistograms(name, title, command, *options)
        self.histPainter = HistogramCollectionPainter()
        self.canv = self.histPainter.paint(self.histCol, *options)
        return self.canv

    def _generateHistograms(self, name, title, command,*options):
        self.histCol = HistogramCollection(name, title)
        for label, p in self._painters.iteritems():
            hc = p._generateHistograms(name, title, command,*options)
            for key, h in hc:
                self.histCol.addHistogram(label + " " + key, h)
        return self.histCol

    def makeHistograms(self, name, title, command, *options):
        """Create histograms but does not draw them.

        :param name: the name of the returned canvas.
        :type name: str
        :param title: the title of the returned canvas.
        :type title: str
        :param command: the draw command (in the format of the ROOT.TTree.Draw mini-language).
        :type command: str
        :param options: optional number of draw options (see drawoptions module).
        :returns: the histogram collection.
        """
        return self._generateHistograms(name, title, command,*options)

    
###############################################################################

class HistogramCollection:
    """Container for multiple histograms that go into a single plot.
    
    Histograms are stored in a dictionary with a string key. The primary use 
    case is histograms split by event type created with TreePainter. 
    
    This class does not implment direct filling. The histograms must be created 
    and filled elsewhere and added to the collection with the addHistogram 
    method.
    
    Note: the child class (eg HistogramCollection1D) implement direct filling.
    
    """
    def __init__(self, name, title):
        self.name = name
        self.title = title
        self.hist = OrderedDict()
        self.nameList = list()
        
    def __iter__(self):
        for n in self.nameList:
            yield n,self.hist[n]
            
    def keys(self):
        return self.hist.keys()
            
    def getFirst(self):
        return iter(self).next()[1]
    
    def addHistogram(self, key, histogram):
        if key in self.hist:
            self.hist[key].Add(histogram)
        else:
            self.hist[key] = histogram
            self.nameList.append(key)
    
    def getTitle(self):
        return self.title
    
    def getName(self):
        return self.name
    
    def getXMin(self):
        xMin = [ h.GetXaxis().GetXmin() for k,h in self ]
        return min(xMin)

    def getXMax(self):
        xMax = [ h.GetXaxis().GetXmax() for k,h in self ]
        return max(xMax)

    def getYMin(self):
        yMin = [ h.GetMinimum() for k,h in self ]
        yMin.append(0.0)
        return min(yMin)

    def getYMax(self):
        yMax = [ h.GetMaximum() for k,h in self ]
        yMax.append(0.0)
        return max(yMax)

    #def getZMin(self):
    #    return self.zMin        

    #def getZMax(self):
    #    return self.zMax

    def __reversed__(self):
        return reversed(self.hist.items())
    
    def __getitem__(self, name):
        return self.hist[name]
    
    def __len__(self):
        return len(self.nameList)
    
    def __str__(self):
        sio = StringIO.StringIO()
        print >>sio,"HistogramCollection(name=\"",self.name,"\", title=\"",self.title,"\", ",len(self),"entries"
        for n,h in self:
            print >>sio,n,"=",str(h),", entries =",h.GetEntries(),", integral =",h.Integral()
        print >>sio,")"
        return sio.getvalue()
    
    def __add__(self, rhs):
        clone = self.clone()
        clone += rhs
        return clone
    
    def __iadd__(self, rhs):
        for key,hist in rhs:
            self.addHistogram(key,hist)
        return self
    
    def clone(self):
        """Make a clone, also cloning the constituent histograms.
        :returns: a deep copy of self.
        """
        result = copy.deepcopy(self)
        for k in result.hist.keys():
            result.hist[k] = result.hist[k].Clone()
        return result

###############################################################################
    
class HistogramCollection1D(HistogramCollection):
    """A fillable container for 1D multiple histograms that go into a single plot. 
    
    Histograms are stored in a dictionary with a string key. It can be filled 
    directly and will create empty histograms on demand when encountering new 
    keys. For example,
    >>> from commonAnalysis.plot import drawoptions
    >>> xBinning = drawoptions.UniformBinning(10,0.0,2000.0)
    >>> hist = HistogramCollection1D("muonMomentumPlot","plot of muon momentum",)
    >>> for event in tree:
    >>>     muonMomentum = event.muonMomentum
    >>>     eventType = "background"
    >>>     if event.isSignal:
    >>>         eventType = "signal"
    >>>     hist.fill(eventType, eventWeight, muonMomentum)
    >>> 
    
    """
    def __init__(self, name, title, xBinning):
        HistogramCollection.__init__(self,name,title)
        self.xBinning = xBinning
    
    def fill(self, key, eventWeight, x):
        """Fill the histogram. 
        Fills the histogram matching this key. If no matching histogram exists 
        it is created.
        Note that the eventWeight argument is before the x,y,z value argument.
        
        :param key: the key used to separate histograms (eg an event type)
        :type key: str
        :param eventWeight: a weight to fill the histogram with. Set to 1.0 if you have no event weighting.
        :type eventWeight: float
        :param x: the x coordinate (determines which bin to fill)
        :type x: float  
        """
        try:
            hist = self.hist[key]
        except:
            hist = ROOT.TH1D(str(key),str(key),self.xBinning.numBins(),self.xBinning.array())
            hist.SetDirectory(0)
            self.hist[key] = hist
            self.nameList.append(key)
        hist.Fill(x,eventWeight)
        return
    
###############################################################################
    
class HistogramCollection2D(HistogramCollection):
    def __init__(self, name, title, xBinning, yBinning):
        HistogramCollection.__init__(self,name,title)
        self.xBinning = xBinning
        self.yBinning = yBinning

    def fill(self, key, eventWeight, x, y):
        try:
            hist = self.hist[key]
        except:
            hist = ROOT.TH2F(str(key),str(key),self.xBinning.numBins(),self.xBinning.array(),self.yBinning.numBins(),self.yBinning.array())
            hist.SetDirectory(0)
            self.hist[key] = hist
            self.nameList.append(key)
        hist.Fill(x,y,eventWeight)
        return
    
    def getXMin(self):
        return self.xBinning.lowBound()

    def getXMax(self):
        return self.xBinning.highBound()

    def getYMin(self):
        return self.yBinning.lowBound()

    def getYMax(self):
        return self.yBinning.highBound()

###############################################################################
    
class HistogramCollection3D(HistogramCollection):
    def __init__(self, name, title, xBinning, yBinning, zBinning):
        HistogramCollection.__init__(self,name,title)
        self.xBinning = xBinning
        self.yBinning = yBinning
        self.zBinning = zBinning

    def fill(self, key, eventWeight, x, y, z):
        try:
            hist = self.hist[key]
        except:
            hist = ROOT.TH2F(str(key),str(key),self.xBinning.numBins(),self.xBinning.array(),self.yBinning.numBins(),self.yBinning.array())
            hist.SetDirectory(0)
            self.hist[key] = hist
            self.nameList.append(key)
        hist.Fill(x,y,z,eventWeight)
        return
    
    def getXMin(self):
        return self.xBinning.lowBound()

    def getXMax(self):
        return self.xBinning.highBound()

    def getYMin(self):
        return self.yBinning.lowBound()

    def getYMax(self):
        return self.yBinning.highBound()
    
    def getZMin(self):
        return self.zBinning.lowBound()

    def getZMax(self):
        return self.zBinning.highBound()    

###############################################################################
    
class HistogramCollectionPainter:
    """Formats and draws a set of histograms to a TCanvas.
     
    """
    def __init__(self):
        #nothing to do here.
        pass
    
    def paint(self, histCollection, *options):
        """Formats and draws a collection of histograms to t a TCanvas. 
        The behaviour is controlled through the drawoptions defined in 
        commonAnalysis.plot.drawoptions.
        :param histCollection: the set of histograms to be drawn.
        :type histCollection: commonAnalysis.plot.drawTools.HistogramCollection
        :options: one of the draw options in commonAnalysis.plot.drawoptions. 
        """
        _checkAllOptionsAreValid(options)
        #make a copy of the histogram so that we don't modify the originals.
        self.histCol = histCollection.clone()
        #store a reference of the options in this class so we don't have to pass 
        #them to every method
        self.options = options
        self._handleOverflow()
        self._normaliseHistograms()
        self._buildHistogramsStackMC()
        self._calculateFrameLimits()
        self._makeCanvas()
        self._drawHistograms()
        self._formatHistograms()
        self.canv.RedrawAxis()
        self._drawLegend()
        self.canv.Update()
        return self.canv
    
    def _makeCanvas(self):
        #find relevent options
        axisLabels = self._findOption(drawoptions.AxisLabels, default=drawoptions.AxisLabels(x=self.histCol.name))
        axisScale = self._findOption(drawoptions.AxisScale, default=drawoptions.AxisScale())
        canvasSize = self._findOption(drawoptions.CanvasSize, default=drawoptions.CanvasSize())
        frameBinLabels = self._findOption(drawoptions.BinLabels, default=None)
        #create the canvas
        canv = ROOT.TCanvas(self.histCol.getName(), self.histCol.getTitle(),canvasSize.getX(),canvasSize.getY())
        frameTitle = axisLabels.getTitle()
        xMin,yMin,xMax,yMax = self.frameLimits
        frame = canv.DrawFrame(xMin,yMin,xMax,yMax,frameTitle)
        if frameBinLabels:
            self._setFrameBinLabels(frameBinLabels, frame)
        self.canv = canv
        if axisScale.isLogY():
            self.canv.SetLogy(1)
        if axisScale.isLogX():
            self.canv.SetLogx(1)
        if axisScale.isLogZ():
            self.canv.SetLogz(1)
        self.frame = frame
        return

    def _getbinedges(self, hist):
        numbins = hist.GetNbinsX()
        binedges = []
        for ii in xrange(numbins + 1):
            binedges.append(hist.GetXaxis().GetBinLowEdge(ii + 1))
        return numpy.array(binedges, dtype=float)

    def _setFrameBinLabels(self, frameBinLabels, frame):
        hist = self.histCol.getFirst()
        #get bin edges
        binedges = self._getbinedges(hist)
        #copy binning
        frame.GetXaxis().Set(len(binedges)-1, binedges)
        #set bin labels
        if frameBinLabels.xlabels is not None:
            xlabels = frameBinLabels.xlabels
            if xlabels == drawoptions.BinLabels.auto:
                #copy labels from the input histogram
                xlabels = []
                for ii in xrange(hist.GetNbinsX()):
                    l = hist.GetXaxis().GetBinLabel(ii + 1)
                    xlabels.append(l)
            #set labels
            for ii, label in enumerate(xlabels):
                frame.GetXaxis().SetBinLabel(ii + 1, label)
        return
    
    def _handleOverflow(self):
        showopt = self._findOption(drawoptions.ShowOverflow, default=drawoptions.ShowOverflow(False))
        if showopt.flag:
            for k,h in self.histCol:
                if h.GetDimension()==1: #only implemented for 1D histograms
                    underflow = h.GetBinContent(0)
                    overflow = h.GetBinContent(h.GetNbinsX()+1)
                    firstbin = h.GetBinContent(1)
                    lastbin = h.GetBinContent(h.GetNbinsX())
                    h.SetBinContent(1, underflow + firstbin)
                    h.SetBinContent(h.GetNbinsX(), overflow + lastbin)
                    h.SetBinContent(0, 0.0)
                    h.SetBinContent(h.GetNbinsX()+1, 0.0)
                    #TODO : implement errors
        return
    
    def _normaliseHistograms(self):
        """Normalise the histograms according to the option provided.
        Behaviour is governed by commonAnalysis.drawoptions.Normalisation.
        Currently implemented:
        * Normalisation.noNormalisation (histograms untouched, default behaviour).
        * Normalisation.mcToData (enforces total MC integral == total data integral)
        * Normalisation.unitArea (each histogram is normalised to unit area)
        Currently not implemented:
        * Normalisation.pot 
        """
        norm = self._findOption(drawoptions.Normalisation, default=drawoptions.Normalisation(drawoptions.Normalisation.noNormalisation))
        treatAsData = self._findOption(drawoptions.TreatAsData, default=drawoptions.TreatAsData())
        mode = norm.getMode()
        if mode == drawoptions.Normalisation.noNormalisation:
            #do nothing
            pass
        elif mode == drawoptions.Normalisation.pot:
            raise NotImplemented
        elif mode == drawoptions.Normalisation.mcToData:
            totalData = 0.0
            totalMc = 0.0
            for k,h in self.histCol:
                n = h.Integral(0,h.GetNbinsX()+1)
                if k in treatAsData:
                    totalData += n
                else:
                    totalMc += n
            scale = 1.0
            if (not totalData == 0.0) and (not totalMc == 0.0):
                scale = totalData / totalMc
            for k,h in self.histCol:
                if k not in treatAsData:
                    h.Scale(scale)
        elif mode == drawoptions.Normalisation.totalUnitArea:
            totalData = 0.0
            totalMc = 0.0
            for k,h in self.histCol:
                n = h.Integral(0,h.GetNbinsX()+1)
                if k in treatAsData:
                    totalData += n
                else:
                    totalMc += n
            scaleData = 1.0
            scaleMC = 1.0
            if (not totalMc == 0.0):
                scaleMC = 1.0 / totalMc
            if (not totalData == 0.0):
                scaleData = 1.0 / totalData
            for k,h in self.histCol:
                if k not in treatAsData:
                    h.Scale(scaleMC)
                else:
                    h.Scale(scaleData)
        elif mode == drawoptions.Normalisation.unitArea:
            for k,h in self.histCol:
                n = h.Integral(0,h.GetNbinsX()+1)
                scale = 0.0
                if n>0:
                    scale = 1./n
                h.Scale(scale)
        else:
            raise Exception("invalid normalisation option",norm)
        return
    
    def _calculateFrameLimits(self):
        """Automatically determines the x and y ranges for the plot.
        Alternatively the drawoptions.FrameRange can be used to override the 
        automatic behaviour.
        """
        stackOption = self._findOption(drawoptions.Stack, default=drawoptions.Stack(drawoptions.Stack.overlay))
        stackMode = stackOption.getStackMode()
        frameRange = self._findOption(drawoptions.FrameRange, default=drawoptions.FrameRange())
        xMin = frameRange.getXMin()
        xMax = frameRange.getXMax()
        yMin = frameRange.getYMin()
        yMax = frameRange.getYMax()
        auto = drawoptions.FrameRange.auto
        if frameRange.getXMin() == auto:
            xMin = self.histCol.getXMin()
        if frameRange.getXMax() == auto:
            xMax = self.histCol.getXMax()
        if stackMode == drawoptions.Stack.stackMC:
            if frameRange.getYMin() == auto:
                yMin = self.stackedCol.getYMin()
            if frameRange.getYMax() == auto:
                yMax = self.stackedCol.getYMax()
                headRoom = frameRange.getHeadroom()
                yMax = yMax + headRoom*math.fabs(yMax)
        else:
            if frameRange.getYMin() == auto:
                yMin = self.histCol.getYMin()
            if frameRange.getYMax() == auto:
                yMax = self.histCol.getYMax()
                headRoom = frameRange.getHeadroom()
                yMax = yMax + headRoom*math.fabs(yMax)
        self.frameLimits = (xMin,yMin,xMax,yMax)
        return
        
    
    def _formatHistograms(self):
        """Format marker, line and fill attributes of the histogram.
        By default uses drawoptions.Format but this can be overided by 
        providing a drawoptions object that inherits from BaseFormat.  
        """
        treatAsData = self._findOption(drawoptions.TreatAsData, default=drawoptions.TreatAsData())
        format = self._findOption(drawoptions.BaseFormat, default=drawoptions.Format())
        for name,hist in self.drawnHistograms.iteritems():
            isData = name in treatAsData
            format.format(name,hist,isData)
        return
    
    def _drawHistograms(self):
        """Draw the histograms on the canvas. 
        
        Histogram drawing depends on the stack options. See 
        _drawHistogramsOverlay and _drawHistogramsStackMC for implementations.  
        """
        self.canv.cd()
        stackOption = self._findOption(drawoptions.Stack, default=drawoptions.Stack(drawoptions.Stack.overlay))
        nDimensions = self._findOption(drawoptions.NDimensions, default=drawoptions.NDimensions(1))
        if nDimensions.is2D():
            self._drawHistograms2D()
        else:
            stackMode = stackOption.getStackMode()
            if stackMode == drawoptions.Stack.overlay:
                self._drawHistogramsOverlay()
            elif stackMode == drawoptions.Stack.stackMC:
                self._drawHistogramsStackMC()   
        return
    
    def _drawHistograms2D(self):
        self.drawnHistograms = OrderedDict()
        for name,hist in self.histCol:
            opts = ["SAME","COLZ"]
            hist.Draw(",".join(opts))
            self.drawnHistograms[name] = hist
        return
        
    def _drawHistogramsOverlay(self):
        """Draw histograms overlayed.
        
        By default, data will be drawn with errors shown.
        MC is drawn as a line.
        """
        treatAsData = self._findOption(drawoptions.TreatAsData, default=drawoptions.TreatAsData())
        histogramDrawOption = self._findOption(drawoptions.HistogramDrawOption, default=drawoptions.HistogramDrawOption())        
        self.drawnHistograms = OrderedDict()
        for name, hist in reversed(self.histCol):
            #print "drawing",key,hist,hist.GetEntries()
            opts = ["SAME"]
            if name in histogramDrawOption:
                opts.append(histogramDrawOption[name])
            elif name in treatAsData:
                opts.append("E1")
            else:
                opts.append("HIST")
            hist.Draw(",".join(opts))
            self.drawnHistograms[name] = hist
        return
    
    def _buildHistogramsStackMC(self):
        """Construct a set of stacked histograms.
        
        In order to draw stacked histograms, each histogram is filled with the 
        events of each histogram beneath it. This method computes these summed 
        histograms.
        """
        treatAsData = self._findOption(drawoptions.TreatAsData, default=drawoptions.TreatAsData())
        stackOption = self._findOption(drawoptions.Stack, default=drawoptions.Stack(drawoptions.Stack.overlay))
        stackMode = stackOption.getStackMode()
        stackOrderMode = stackOption.getOrderMode()
        if not stackMode == drawoptions.Stack.stackMC:
            #don't build the stack if it's not required in the options.
            return
        #first get list of histograms and their order.
        allNames = set( name for name,hist in self.histCol)
        self.orderedNames = []
        #user defined names go first
        userDefinedOrder = stackOption.getUserOrder()
        for name in userDefinedOrder:
            if name in allNames:
                self.orderedNames.append(name)
        if stackOrderMode == drawoptions.Stack.orderByIntegral:
            #next order the remaining histograms by their total integral
            remainingNames = list(allNames.difference(self.orderedNames))
            remainingNames.sort(key=lambda n:self.histCol[n].Integral(), reverse=True)
            self.orderedNames += remainingNames
        else:
            for name,hist in self.histCol:
                if name not in self.orderedNames: 
                    self.orderedNames.append(name)
        #calculate stacked histograms
        self.stackedCol = HistogramCollection(self.histCol.getName(),self.histCol.getTitle())
        for i,name in enumerate(self.orderedNames):
            hist = self.histCol[name]
            if name not in treatAsData:
                hist = hist.Clone()
                #sum up histograms
                for j in xrange(i+1,len(self.orderedNames)):
                    jthName = self.orderedNames[j]
                    jthHist = self.histCol[jthName]
                    if jthName not in treatAsData:
                        hist.Add(jthHist)
            self.stackedCol.addHistogram(name,hist)
        return
        
    
    def _drawHistogramsStackMC(self):
        """Draw stacked histograms.
        
        To be displayed correctly the top histograms have to be drawn first
        (to avoid covering up the smaller histograms beneath them.
        """
        treatAsData = self._findOption(drawoptions.TreatAsData, default=drawoptions.TreatAsData())
        histogramDrawOption = self._findOption(drawoptions.HistogramDrawOption, default=drawoptions.HistogramDrawOption())        
        #now draw the stacked histograms
        self.drawnHistograms = OrderedDict()
        for name in reversed(self.orderedNames):
            hist = self.stackedCol[name]
            opts = ["SAME"]
            if name in histogramDrawOption:
                opts.append(histogramDrawOption[name])
            elif name in treatAsData:
                opts.append("E1")
            else:
                opts.append("HIST")
            hist.Draw(",".join(opts))
            self.drawnHistograms[name] = hist
        return
    
    def _drawLegend(self):
        """Draw the legend.
        
        The legend position is automatically calculated unless the user draw 
        option drawoptions.LegendPosition is set. 
        """
        treatAsData = self._findOption(drawoptions.TreatAsData, default=drawoptions.TreatAsData())
        legendPosition = self._findOption(drawoptions.LegendPosition, default=drawoptions.LegendPosition())
        if legendPosition.doDraw():
            if legendPosition.hasUserLimits():
                xLow,yLow,xHigh,yHigh = legendPosition.calculateLegendLimits()
                opt = "brNDC"
            else:
                xLow,yLow,xHigh,yHigh = self._calculateLegendLimits()
                opt = "br"
            #print xLow,yLow,xHigh,yHigh
            leg = ROOT.TLegend(xLow,yLow,xHigh,yHigh, "", opt)
            leg.SetFillStyle(1001)
            leg.SetFillColor(ROOT.kWhite)
            for name in self.drawnHistograms.iterkeys():
                hist = self.drawnHistograms[name]
                opt = "LF"
                if name in treatAsData:
                    opt = "LP"
                leg.AddEntry(hist,name, opt)
            leg.Draw()
            #store the legend in the histogram collection
            self.histCol.legend = leg
        return
    
    def _calculateLegendLimits(self):
        listOfHistograms = self.drawnHistograms.values()
        return _calculateLegendLimits(self.frameLimits, listOfHistograms=listOfHistograms)
        
    def _findOption(self, optionType, default):
        """A convienient wrapper around the global function _findOption.
        We use a wrapper and global function because multiple classes in this module use this function.
        """
        return _findOption(optionType,self.options,default=default)
    
###############################################################################

class EfficiencyRejectionGraph:
    def __init__(self, 
                 name, 
                 title,
                 inputDataset, 
                 treeName, 
                 command,
                 signalCut, 
                 backgroundCut,
                 binning,
                 isLessThan,
                 treePainterOptions=[]
                  ):
        signalName = "signal"
        backgroundName = "background"
        self.name = name
        self.title = title
        #split data set into signal and background
        sd = drawoptions.SplitDataset()
        sd.add(signalName,signalCut)
        sd.add(backgroundName,backgroundCut)
        #make histograms
        treePainter = TreePainter(inputDataset, treeName)
        histograms = treePainter.makeHistograms(name, title, command,sd,binning,*treePainterOptions)
        #get input histograms
        self.signalHistogram = histograms[signalName]
        self.backgroundHistogram = histograms[backgroundName]
        #create cumulative distribution
        self._convertToCumulativeAndScale(self.signalHistogram)
        self._convertToCumulativeAndScale(self.backgroundHistogram)
        #create the plot
        self._makeGraph(isLessThan=isLessThan)

    def _makeGraph(self, isLessThan):
        sh = self.signalHistogram
        bh = self.backgroundHistogram
        points = []
        for i in xrange(1,sh.GetNbinsX()+1):
            s = sh.GetBinContent(i)
            b = bh.GetBinContent(i)
            sErr = sh.GetBinError(i)
            bErr = bh.GetBinError(i)
            if isLessThan:
               signalEff = s
               bgEff = b
            else:
               signalEff = 1.0-s
               bgEff = 1.0-b
            bgRej = 1.0-bgEff
            x = signalEff
            y = bgRej
            ex = sErr
            ey = bErr
            p = ( x, y, ex, ey)
            points.append( p )
        nPoints = len(points)
        graph = ROOT.TGraphErrors(nPoints)
        graph.SetTitle(";signal eff;background rejection")
        for i,p in enumerate(points):
            x,y,ex,ey = p
            graph.SetPoint(i,x,y)
            graph.SetPointError(i,ex,ey)
        self.graph = graph
        return
    
    def _calculateBinomialError(self,n,p):
        err = 0.0
        if n>0:
            n = float(n)
            p = float(p)
            err = math.sqrt(p*(1.-p)/n)
        return err

    def _convertToCumulativeAndScale(self, histogram):
        nEntries = histogram.GetEntries()
        #scale
        #scale = histogram.Integral(0,histogram.GetNbinsX()+1)
        #if scale!=0.0:
        #    scale = 1.0/scale
        #histogram.Scale(scale)
        #calc cumulative
        histogram.ComputeIntegral()
        integral = histogram.GetIntegral()
        histogram.SetContent(integral)
        #calculate errors
        for i in xrange(1,histogram.GetNbinsX()+1):
            p = histogram.GetBinContent(i)
            e = self._calculateBinomialError(nEntries, p)
            histogram.SetBinError(i,e)
        return
    
    def getGraph(self):
        return self.graph
        

###############################################################################

class EfficiencyRejectionPainter:
    
    def __init__(self):
        self.canvas = None
    
    def paint(self, name, title, effRejGraphs,  *options):
        _checkAllOptionsAreValid(options)
        self.options = options
        self.effRejGraphs = effRejGraphs
        self._createCanvas(name,title)
        self._formatGraphs()
        self._drawPlot(effRejGraphs)
        self.canv.RedrawAxis()
        self._drawLegend()
        self.canv.Update()
        return self.canv
    
    def _createCanvas(self, name, title):
        #determine frame range
        frameRange = self._findOption(drawoptions.FrameRange, default=drawoptions.FrameRange())
        xMin = frameRange.getXMin()
        xMax = frameRange.getXMax()
        yMin = frameRange.getYMin()
        yMax = frameRange.getYMax()
        auto = drawoptions.FrameRange.auto
        if xMin==auto:
            xMin = 0.0
        if yMin==auto:
            yMin = 0.0
        if xMax==auto:
            xMax=1.0
        if yMax==auto:
            yMax=1.0
        self.frameLimits = (xMin,yMin,xMax,yMax)
        #determine canvas size
        canvasSize = self._findOption(drawoptions.CanvasSize, default=drawoptions.CanvasSize())
        xPixel = canvasSize.getX()
        yPixel = canvasSize.getY()
        #determine axis labels
        axisLabels = self._findOption(drawoptions.AxisLabels, default=drawoptions.AxisLabels(x="signal efficiency",y="background rejection"))
        frameTitle = axisLabels.getTitle()
        self.canv = ROOT.TCanvas(name, title, xPixel, yPixel)
        self.canv.DrawFrame(xMin,yMin,xMax,yMax, frameTitle)
        return

    def _drawPlot(self, effRejGraphs):
        self.canv.cd()
        for effRej in effRejGraphs:
            graph = effRej.getGraph()
            #col = ROOT.kRed
            #graph.SetFillColor(col)
            #graph.SetFillStyle(3001)
            #graph.SetLineColor(col)
            #graph.SetMarkerSize(0)
            graph.Draw("L2,SAME")
        return
    
    def _formatGraphs(self):
        """Format marker, line and fill attributes of the histogram.
        By default uses drawoptions.Format but this can be overided by 
        providing a drawoptions object that inherits from BaseFormat.  
        """
        format = self._findOption(drawoptions.BaseFormat, default=drawoptions.Format(fillArea=True))
        for effRej in self.effRejGraphs:
            name = effRej.name
            graph = effRej.getGraph()
            format.format(name,graph,isData=False)
        return
    
    def _findOption(self, optionType, default):
        """A convienient wrapper around the global function _findOption.
        We use a wrapper and global function because multiple classes in this module use this function.
        """
        return _findOption(optionType,self.options,default=default)
    
    def _drawLegend(self):
        """Draw the legend.
        
        The legend position is automatically calculated unless the user draw 
        option drawoptions.LegendPosition is set. 
        """
        treatAsData = self._findOption(drawoptions.TreatAsData, default=drawoptions.TreatAsData())
        legendPosition = self._findOption(drawoptions.LegendPosition, default=drawoptions.LegendPosition(doDraw=True))
        if legendPosition.doDraw():
            if legendPosition.hasUserLimits():
                xLow,yLow,xHigh,yHigh = legendPosition.calculateLegendLimits()
                opt = "brNDC"
            else:
                xLow,yLow,xHigh,yHigh = self._calculateLegendLimits()
                opt = "br"
            #print xLow,yLow,xHigh,yHigh
            leg = ROOT.TLegend(xLow,yLow,xHigh,yHigh, "", opt)
            leg.SetFillStyle(1001)
            leg.SetFillColor(ROOT.kWhite)
            for effRej in self.effRejGraphs:
                graph = effRej.getGraph()
                opt = "LF"
                leg.AddEntry(graph,effRej.title, opt)
            leg.Draw()
            #store the legend in the histogram collection
            self.legend = leg
        return

    def _calculateLegendLimits(self):
        listOfGraphs = [g.getGraph() for g in self.effRejGraphs]
        return _calculateLegendLimits(self.frameLimits, listOfGraphs=listOfGraphs, targetLeft=True, targetBottom=True)

###############################################################################

class EfficiencyGraph:
    def __init__(self, 
                 name, 
                 title,
                 inputDataset=None, 
                 treeName=None, 
                 command=None,
                 cut=None, 
                 treePainterOptions=[],
                 tree=None
                  ):
        self.name = name
        self.title = title
        self.xMin = None
        self.yMin = None
        self.xMax = None
        self.yMax = None
        if command is None or cut is None:
            raise ValueError("At least one of \"command\" and \"cut\" arguments not set.")
        if (bool(tree) == bool(inputDataset and treeName)):
            raise ValueError("You must set \"tree\" or \"inputDataset\" and \"treeName\"")
        self._makeHistograms(inputDataset, treeName, command, cut, treePainterOptions, tree)
        self._makePlots()
    
    @classmethod
    def fromDataset(name, title, inputDataset, treeName, command, cut, 
                    treePainterOptions=[]):
        return EfficiencyGraph(name=name, title=title, 
                               inputDataset=inputDataset, treeName=treeName,
                               command=command, cut=cut, 
                               treePainterOptions=treePainterOptions)

    @classmethod
    def fromTree(name, title, tree, command, cut, 
                    treePainterOptions=[]):
        return EfficiencyGraph(name=name, title=title, 
                               tree=tree,
                               command=command, cut=cut, 
                               treePainterOptions=treePainterOptions)
    
    def _makeHistograms(self, inputDataset, treeName, command, cut, treePainterOptions, tree):
        selectedName = "selectedName"
        totalName = "totalName"
        #split data set into selected and all.
        sd = drawoptions.SplitDataset()
        sd.add(selectedName,cut)
        sd.add(totalName,"1")
        #make histograms
        treePainter = TreePainter(inputDataset, treeName, tree)
        #DEBUG : self.binning,
        histograms = treePainter.makeHistograms(self.name, self.title, command,sd,*treePainterOptions)
        #get input histograms
        self.selectedHistogram = histograms[selectedName]
        self.totalHistogram = histograms[totalName]
        return
    
    def _makePlots(self):
        #create the plot
         
        if hasattr(ROOT,"TEfficiency") and (not ROOT.TEfficiency.CheckConsistency(self.selectedHistogram, self.totalHistogram, "w")):
            raise Exception("EfficiencyGraph failed to create efficiency object due to inconsistent histograms",
                            self.selectedHistogram, 
                            self.totalHistogram
                            )
        #self.efficiencyObj = ROOT.TEfficiency(self.selectedHistogram, self.totalHistogram)
        opt = "cl=0.683 b(1,1) mode"
        self.graph = ROOT.TGraphAsymmErrors(self.selectedHistogram, self.totalHistogram, opt) 
        return
    
    def getGraph(self):
        return self.graph
    
    def getXMin(self):
        if self.xMin is None:
            self._calculateFrameLimits()
        return self.xMin
    
    def getXMax(self):
        if self.xMax is None:
            self._calculateFrameLimits()
        return self.xMax 
    
    def getYMin(self):
        if self.yMin is None:
            self._calculateFrameLimits()
        return self.yMin

    def getYMax(self):
        if self.yMax is None:
            self._calculateFrameLimits() 
        return self.yMax
    
    def _calculateFrameLimits(self):
        graph = self.graph
        nPoints = graph.GetN()
        xLowList = []
        yLowList = []
        xHighList = []
        yHighList = []
        for i in xrange(nPoints):
            x,y = ROOT.Double(0.0),ROOT.Double(0.0)
            graph.GetPoint(i,x,y)
            xHigh = graph.GetErrorXhigh(i)
            xLow = graph.GetErrorXlow(i)
            yHigh = graph.GetErrorYhigh(i)
            yLow = graph.GetErrorYlow(i)
            #print i,x,y,xhigh,xlow,yhigh,ylow
            xLowList.append( x - xLow )
            xHighList.append( x + xHigh )
            yLowList.append( y - yLow )
            yHighList.append( y + yHigh )
        try:
            self.xMin = min(xLowList)
            self.xMax = max(xHighList)
            self.yMin = min(yLowList)
            self.yMax = max(yHighList)
        except ValueError:
            #no points, min and max failed.
            self.xMin = 0.0
            self.yMin = 0.0
            self.xMax = 0.0
            self.yMax = 0.0
        return
    
###############################################################################

class EfficiencyPainter:
    def __init__(self):
        self.canvas = None
    
    def paint(self, name, title, efficiencyGraphs,  *options):
        _checkAllOptionsAreValid(options)
        self.name = name
        self.title = title
        self.options = options
        self.efficiencyGraphs = efficiencyGraphs
        self._calculateFrameLimits()
        self._createCanvas(name,title)
        self._formatGraphs()
        self._drawPlot(efficiencyGraphs)
        self.canv.RedrawAxis()
        self._drawLegend()
        self.canv.Update()
        return self.canv
    
    def _createCanvas(self, name, title):
        xMin,yMin,xMax,yMax = self.frameLimits
        #determine canvas size
        canvasSize = self._findOption(drawoptions.CanvasSize, default=drawoptions.CanvasSize())
        xPixel = canvasSize.getX()
        yPixel = canvasSize.getY()
        #determine axis labels
        axisLabels = self._findOption(drawoptions.AxisLabels, default=drawoptions.AxisLabels(x=self.title,y="efficiency"))
        frameTitle = axisLabels.getTitle()
        self.canv = ROOT.TCanvas(name, title, xPixel, yPixel)
        self.canv.DrawFrame(xMin,yMin,xMax,yMax, frameTitle)
        return
    
    def _calculateFrameLimits(self):
        """Automatically determines the x and y ranges for the plot.
        Alternatively the drawoptions.FrameRange can be used to override the 
        automatic behaviour.
        """
        frameRange = self._findOption(drawoptions.FrameRange, default=drawoptions.FrameRange(headroom=0.05))
        xMin = frameRange.getXMin()
        xMax = frameRange.getXMax()
        yMin = frameRange.getYMin()
        yMax = frameRange.getYMax()
        auto = drawoptions.FrameRange.auto
        headRoom = frameRange.getHeadroom()
        if frameRange.getXMin() == auto:
            xMin = min([g.getXMin() for g in self.efficiencyGraphs])
        if frameRange.getXMax() == auto:
            xMax = max([g.getXMax() for g in self.efficiencyGraphs])
        if frameRange.getYMin() == auto:
            yMin = min([g.getYMin() for g in self.efficiencyGraphs])
            yMin = yMin - headRoom*math.fabs(yMin)
            yMin = max(0.0,yMin)
        if frameRange.getYMax() == auto:
            yMax = max([g.getYMax() for g in self.efficiencyGraphs])
            yMaxPlusHeadRoom = yMax + headRoom*math.fabs(yMax)
            if yMax <= 1.0:
                yMax = min(yMaxPlusHeadRoom,1.0)
            else:
                yMax = yMaxPlusHeadRoom
        self.frameLimits = (xMin,yMin,xMax,yMax)
        return

    def _drawPlot(self, efficiencyGraphs):
        self.canv.cd()
        for effRej in efficiencyGraphs:
            graph = effRej.getGraph()
            #col = ROOT.kRed
            #graph.SetFillColor(col)
            #graph.SetFillStyle(3001)
            #graph.SetLineColor(col)
            #graph.SetMarkerSize(0)
            graph.Draw("P,SAME")
        return
    
    def _formatGraphs(self):
        """Format marker, line and fill attributes of the histogram.
        By default uses drawoptions.Format but this can be overided by 
        providing a drawoptions object that inherits from BaseFormat.  
        """
        format = self._findOption(drawoptions.BaseFormat, default=drawoptions.Format(fillArea=True))
        for effRej in self.efficiencyGraphs:
            name = effRej.name
            graph = effRej.getGraph()
            format.format(name,graph,isData=False)
        return
    
    def _findOption(self, optionType, default):
        """A convienient wrapper around the global function _findOption.
        We use a wrapper and global function because multiple classes in this module use this function.
        """
        return _findOption(optionType,self.options,default=default)
    
    def _drawLegend(self):
        """Draw the legend.
        
        The legend position is automatically calculated unless the user draw 
        option drawoptions.LegendPosition is set. 
        """
        treatAsData = self._findOption(drawoptions.TreatAsData, default=drawoptions.TreatAsData())
        legendPosition = self._findOption(drawoptions.LegendPosition, default=drawoptions.LegendPosition())
        if legendPosition.hasUserLimits():
            xLow,yLow,xHigh,yHigh = legendPosition.calculateLegendLimits()
            opt = "brNDC"
        else:
            xLow,yLow,xHigh,yHigh = self._calculateLegendLimits()
            opt = "br"
        #print xLow,yLow,xHigh,yHigh
        leg = ROOT.TLegend(xLow,yLow,xHigh,yHigh, "", opt)
        leg.SetFillStyle(1001)
        leg.SetFillColor(ROOT.kWhite)
        for effRej in self.efficiencyGraphs:
            graph = effRej.getGraph()
            opt = "LP"
            leg.AddEntry(graph,effRej.title, opt)
        leg.Draw()
        #store the legend in the histogram collection
        self.legend = leg
        return

    def _calculateLegendLimits(self):
        listOfGraphs = [g.getGraph() for g in self.efficiencyGraphs]
        return _calculateLegendLimits(self.frameLimits, listOfGraphs=listOfGraphs, targetLeft=True, targetBottom=True)

###############################################################################

def _calculateLegendLimits(frameLimits, listOfHistograms=None, listOfGraphs=None, useExpandAlgorithm = False, targetBottom=False, targetLeft=False):
        """Implements the algorithm to automatically position the legend.
        
        These algorithms can probably be improved, but I have more important 
        things to do right now.
        """
        xMin,yMin,xMax,yMax = frameLimits
        n = 100
        #create empty grid
        grid = dict()
        for i,j in itertools.product(xrange(n),xrange(n)):
            grid[ (i,j) ] = 0
        #fill in black spaces
        for xi in xrange(0,n):
            xVal = xMin + float(xi)*(xMax-xMin)/float(n)
            if listOfHistograms:
                yValues = [ hist.GetBinContent(hist.FindBin(xVal)) + hist.GetBinError(hist.FindBin(xVal)) for hist in listOfHistograms ]
                yLimit = max(yValues)
                for yi in xrange(0,n):
                    yVal = yMin + float(yi)*(yMax-yMin)/float(n)
                    if yVal <= yLimit:
                        grid[ (xi,yi) ] = 1
            if listOfGraphs:
                for graph in listOfGraphs:
                    gv = graph.Eval(xVal)
                    minError = (yMax - yMin) / 20.
                    #maxError = max( [ graph.GetErrorY(p) for p in xrange(graph.GetN()) ] + [minError] )
                    nearestLow = None
                    nearestHigh = None
                    nearestLowError = 0.0
                    nearestHighError = 0.0
                    for pi in xrange(graph.GetN()):
                        px = ROOT.Double(0.0)
                        py = ROOT.Double(0.0)
                        graph.GetPoint(i,px,py)
                        if nearestLow is None or (xVal - px)<nearestLow:
                            nearestLow = xVal - px
                            nearestLowError = graph.GetErrorY(pi)
                        if nearestHigh is None or (xVal - px)>nearestHigh:
                            nearestHigh = xVal - px
                            nearestHighError = graph.GetErrorY(pi)
                    maxError = max([nearestLowError,nearestHighError,minError]) 
                    for yi in xrange(0,n):
                        yVal = yMin + float(yi)*(yMax-yMin)/float(n)
                        if abs(yVal-gv) <= maxError:
                            grid[ (xi,yi) ] = 1
                            
#                for yi in xrange(0,n):
#                    yVal = xMin + float(xi)*(xMax-xMin)/float(n)
#                    covered = False
#                    for graph in listOfGraphs:
#                        for ip in xrange(graph.GetN()):
#                            px = ROOT.Double(0.0)
#                            py = ROOT.Double(0.0)
#                            ex = graph.GetErrorX(ip)
#                            ey = graph.GetErrorY(ip)
#                            graph.GetPoint(ip,px,py)
#                            if (px-ex) < xVal < (px+ex) and (py-ey) < yVal < (py+ey):
#                                covered = True
#                                break
#                        if covered:
#                            grid[ (xi,yi) ] = 1
#                            break
        initialGrid = dict(grid)
        if useExpandAlgorithm:
            while len(grid)>0:
                #remove squares with value >0.
                #tag neighbours for removal next iteration.
                keys = grid.keys()
                keys.sort(key=lambda x:(x[1],x[0]))
                for (i,j) in keys:
                    if grid[(i,j)]>0:
                        del grid[(i,j)]
                    if (i+1,j) in grid:
                        grid[(i+1,j)] = 1
                    if (i-1,j) in grid:
                        grid[(i-1,j)] = 1
                    if (i,j+1) in grid:
                        grid[(i,j+1)] = 1
                    if (i,j-1) in grid:
                        grid[(i,j-1)] = 1
    
            seedPoint = (i,j)
            #expand around the seed point until hit black 
            limits = [max(0,i-1),max(0,j-1), min(i+1,n-1), min(j+1,n-1) ]
            status = [True,True,True,True]
            while any(status):
                for k in xrange(0,4):
                    if status[k]:
                        if k>1:
                            limits[k] += 1
                        else:
                            limits[k] -= 1
                        if limits[k] >= n:
                            limits[k] = n-1
                            status[k] = False
                        if limits[k] < 0:
                            limits[k] = 0
                            status[k] = False
#                        #check every square for black spots
#                        failed = False
#                        for i,j in itertools.product(xrange(limits[0],limits[2]),xrange(limits[1],limits[3])):
#                            if initialGrid[(i,j)]>0:
#                                failed = True
#                                break
                        #check border for black spots
                        failed = False
                        for i,j in itertools.product(xrange(limits[0],limits[2]),[limits[1],limits[3]]):
                            if initialGrid[(i,j)] > 0:
                                failed = True
                                break
                        if not failed:
                            for i,j in itertools.product([limits[0],limits[2]],xrange(limits[1],limits[3])):
                                if initialGrid[(i,j)] > 0:
                                    failed = True
                                    break
                        if failed:
                            status[k] = False
                            if k>1:
                                limits[k] -= 1
                            else:
                                limits[k] += 1
        else:
            #brute force algorithm
            limits = [0.7,0.7,0.92,0.92]
            maxFigOfMerit = -sys.maxint
            stepSize = 3
            for xHigh in reversed(xrange(1,n,stepSize)):
                for xLow in xrange(1,xHigh,stepSize):
                    for yHigh in reversed(xrange(0,n,stepSize)):
                        for yLow in xrange(0,yHigh,stepSize):
                            xSize = (xHigh-xLow)
                            ySize = (yHigh-yLow)
                            targetSize = 2500
                            #print xSize,ySize,xSize*ySize
                            x = xHigh
                            y = yHigh
                            if targetLeft:
                                x = n-xLow
                            if targetBottom:
                                y = n-yLow
                            figOfMerit = -abs(targetSize - xSize*ySize) - (xSize-ySize)**2 + x + y
                            if figOfMerit<maxFigOfMerit or (figOfMerit==maxFigOfMerit and yHigh<limits[3]):
                                continue
                            #check corners for black spots
                            isValid = True
                            for i,j in itertools.product([xLow,xHigh],[yLow,yHigh]):
                                if grid[(i,j)] > 0:
                                    isValid = False
                                    break                                
                            if not isValid:
                                continue
                            #check border for black spots
                            for i,j in itertools.product(xrange(xLow,xHigh),[yLow,yHigh]):
                                if grid[(i,j)] > 0:
                                    isValid = False
                                    break
                            if not isValid:
                                continue
                            for i,j in itertools.product([xLow,xHigh],xrange(yLow,yHigh)):
                                if grid[(i,j)] > 0:
                                    isValid = False
                                    break
                            if not isValid:
                                continue
                            #satified all conditions, make this the current best
                            maxFigOfMerit = figOfMerit 
                            limits = [xLow,yLow,xHigh,yHigh]
        #add in a small border
        border = 3
        if limits[2] - limits[0] > border*4:
            limits[0] += border
            limits[2] -= border
        if limits[3] - limits[1] > border*4:
            limits[1] += border
            #limits[3] -= border
        #calculate final legend limits 
        xLow = xMin + float(limits[0])*(xMax-xMin)/float(n)
        yLow = yMin + float(limits[1])*(yMax-yMin)/float(n)
        xHigh = xMin + float(limits[2])*(xMax-xMin)/float(n)
        yHigh = yMin + float(limits[3])*(yMax-yMin)/float(n)
        return xLow,yLow,xHigh,yHigh
     
###############################################################################

_unitTestInputFileName = "unitTest_drawTools_inputFile.root"
_unitTestInputTreeName = "unitTest_drawTools_tree"

###############################################################################

def _makeUnitTestFile():
    import os
    from commonAnalysis.tools import ntuple
    if not os.path.exists(_unitTestInputFileName):
        tfile = ROOT.TFile(_unitTestInputFileName,"recreate")
        tree = ROOT.TTree(_unitTestInputTreeName,_unitTestInputTreeName)
        branch_eventNum = ntuple.BranchPrimative("eventNum",tree)
        branch_eventType = ntuple.BranchPrimative("eventType",tree)
        branch_likelihood = ntuple.BranchPrimative("likelihood",tree)
        branch_ann = ntuple.BranchPrimative("ann",tree)
        branch_momentum = ntuple.BranchPrimative("momentum",tree)
        rand = ROOT.TRandom3(1928)
        for i in xrange(100000):
            branch_eventNum.setValue(i)
            #fil signal
            sigLik = rand.Gaus(-1.0,1.0)
            branch_eventType.setValue( 1.0 )
            branch_likelihood.setValue( sigLik )
            branch_ann.setValue( rand.Gaus(1.0,1.0) )
            branch_momentum.setValue( rand.Gaus(2.0+sigLik,1.0) )
            tree.Fill()
            #fill background
            bgLik = rand.Gaus(1.0,1.0)
            branch_eventType.setValue( 0.0 )
            branch_likelihood.setValue( bgLik )
            branch_ann.setValue( rand.Gaus(-1.0,2.0) )
            branch_momentum.setValue( rand.Gaus(2.0-bgLik,1.0) )
            tree.Fill()
        tree.Write()
        tfile.Close()
    return

###############################################################################

class GraphCollectionPainter:

    def __init__(self):
        self.canv = None

    def paint(self, name, title, graphs, *options):
        _checkAllOptionsAreValid(options)
        #make a copy of the graphs so that we don't modify the originals
        graphs = [gr.Clone() for gr in graphs]
        #store a reference of the options in this class so we don't have to pass
        #them to every method
        framelimits = self._calculateFrameLimits(graphs, options)
        canv, frame = self._makeCanvas(name, title, framelimits, options)
        self._formatGraphs(graphs, options)
        self._drawGraphs(graphs, canv, options)
        canv.RedrawAxis()
        leg = self._drawLegend(graphs, framelimits, options)
        canv.Update()
        #store drawn objects in the canvas to prevent automatic deletion
        canv._graphcolpainter_store = [leg, graphs, frame]
        return canv

    def _formatGraphs(self, graphs, options):
        treatAsData = _findOption(drawoptions.TreatAsData, options, default=drawoptions.TreatAsData())
        format = _findOption(drawoptions.BaseFormat, options, default=drawoptions.Format())
        for gr in graphs:
            name = gr.GetName()
            isData = name in treatAsData
            format.format(name, gr, isData)
        return

    def _drawGraphs(self, graphs, canv, options):
        canv.cd()
        for gr in graphs:
            gr.Draw()
        return

    def _drawLegend(self, graphs, frameLimits, options):
        """Draw the legend.

        The legend position is automatically calculated unless the user draw
        option drawoptions.LegendPosition is set.
        """
        leg = None
        #treatAsData = _findOption(drawoptions.TreatAsData, options, default=drawoptions.TreatAsData())
        legendPosition = _findOption(drawoptions.LegendPosition, options, default=drawoptions.LegendPosition())
        if legendPosition.doDraw():
            if legendPosition.hasUserLimits():
                xLow,yLow,xHigh,yHigh = legendPosition.calculateLegendLimits()
                opt = "brNDC"
            else:
                xLow,yLow,xHigh,yHigh = _calculateLegendLimits(frameLimits, listOfGraphs=graphs, targetLeft=True, targetBottom=True)
                opt = "br"
            #print xLow,yLow,xHigh,yHigh
            leg = ROOT.TLegend(xLow,yLow,xHigh,yHigh, "", opt)
            leg.SetFillStyle(1001)
            leg.SetFillColor(ROOT.kWhite)
            for gr in graphs:
                name = gr.GetTitle()
                opt = "LP"
                leg.AddEntry(gr, name, opt)
            leg.Draw()
        return leg



    def _calculateFrameLimits(self, graphs, options):
        """Automatically determines the x and y ranges for the plot.
        Alternatively the drawoptions.FrameRange can be used to override the
        automatic behaviour.
        """
        frameRange = _findOption(drawoptions.FrameRange, options, default=drawoptions.FrameRange())
        xMin = frameRange.getXMin()
        xMax = frameRange.getXMax()
        yMin = frameRange.getYMin()
        yMax = frameRange.getYMax()
        auto = drawoptions.FrameRange.auto
        if frameRange.getXMin() == auto:
            xMin = min(gr.GetXaxis().GetXmin() for gr in graphs)
        if frameRange.getXMax() == auto:
            xMax = min(gr.GetXaxis().GetXmax() for gr in graphs)
        if frameRange.getYMin() == auto:
            yMin = min(gr.GetHistogram().GetMinimum() for gr in graphs)
        if frameRange.getYMax() == auto:
            yMax = max(gr.GetHistogram().GetMaximum() for gr in graphs)
            headRoom = frameRange.getHeadroom()
            yMax = yMax + headRoom*math.fabs(yMax)
        frameLimits = (xMin,yMin,xMax,yMax)
        return frameLimits

    def _makeCanvas(self, name, title, frameLimits, options):
        #find relevent options
        axisLabels = _findOption(drawoptions.AxisLabels, options, default=drawoptions.AxisLabels(x=name))
        axisScale = _findOption(drawoptions.AxisScale, options, default=drawoptions.AxisScale())
        canvasSize = _findOption(drawoptions.CanvasSize, options, default=drawoptions.CanvasSize())
        #create the canvas
        canv = ROOT.TCanvas(name, title, canvasSize.getX(), canvasSize.getY())
        frameTitle = axisLabels.getTitle()
        xMin,yMin,xMax,yMax = frameLimits
        frame = canv.DrawFrame(xMin,yMin,xMax,yMax,frameTitle)
        self.canv = canv
        if axisScale.isLogY():
            canv.SetLogy(1)
        if axisScale.isLogX():
            canv.SetLogx(1)
        if axisScale.isLogZ():
            canv.SetLogz(1)
        return canv, frame

###############################################################################

def _unitTest_EfficiencyPainter():
    from commonAnalysis.analysis import datasets
    inputDataset = datasets.Dataset("unitTestDataset",_unitTestInputFileName)
    effPlotContainer = []
    for (name,isLessThan) in [ ("likelihood",True), ("ann", False) ]:
        effPlot = EfficiencyGraph(
                     name=name, 
                     title=name,
                     inputDataset=inputDataset, 
                     treeName=_unitTestInputTreeName, 
                     command="momentum",
                     cut=name+">0.0",
                     binning = drawoptions.UniformBinning(100,0.0,200.0),
                     treePainterOptions=[drawoptions.Cut("eventType==1")],
                     )
        effPlotContainer.append(effPlot)
    painter = EfficiencyPainter()
    canvas = painter.paint("unitTest_drawTools_efficiencyPlot","",effPlotContainer)
    canvas.SaveAs(canvas.GetName()+".eps")
    return    

###############################################################################

def _unitTest_EfficiencyRejectionPainter():
    from commonAnalysis.analysis import datasets
    inputDataset = datasets.Dataset("unitTestDataset",_unitTestInputFileName)
    effRejContainer = []
    for (name,isLessThan) in [ ("likelihood",True), ("ann", False) ]:
        effRej = EfficiencyRejectionGraph(
                     name=name, 
                     title=name,
                     inputDataset=inputDataset, 
                     treeName=_unitTestInputTreeName, 
                     command=name,
                     signalCut="eventType==1", 
                     backgroundCut="eventType!=1",
                     binning=drawoptions.UniformBinning(1000,-3.0,3.0),
                     isLessThan=isLessThan,
                     treePainterOptions=[drawoptions.Cut("eventNum<1000")],
                     )
        effRejContainer.append(effRej)
    painter = EfficiencyRejectionPainter()
    canvas = painter.paint("unitTest_drawTools_effRejPlot","",effRejContainer)
    canvas.SaveAs(canvas.GetName()+".eps")
    return

###############################################################################

def main():
    ROOT.gROOT.SetBatch(True)
    _makeUnitTestFile()
    _unitTest_EfficiencyPainter()
    #_unitTest_EfficiencyRejectionPainter()
    return

if __name__ == "__main__":
    main()       
