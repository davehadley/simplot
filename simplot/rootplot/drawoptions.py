"""Provides plot drawing options to be supplied to the tools in the drawTools module. 
"""
import ROOT
import StringIO
import itertools
import numpy
from . import palette


###############################################################################

class InvalidDrawOption(Exception):
    """ An exception to be raised if case of bad user input to the draw options classes. 
    """
    pass

###############################################################################

class EventWeight:
    """Define the variable to be used as an event weight.
    
    Fill a histogram with an event weight. For example, if you want to draw 
    reconstructed muon momentum with flux and pot weighting and  your ntuple 
    contains variables called "pMu", "FluxWeight" and "PotWeight" your code 
    would look similar to:
    
    >>> from commonAnalysis.plot import drawTools
    >>> from commonAnalysis.plot import drawOptions
    >>> 
    >>> tp = drawTools.TreePainter(myInputDataset,"MyTree")
    >>> canvas = tp.paint("muonMomentum",
    >>>                   "Muon momementum",
    >>>                   "pMu",
    >>>                   drawOptions.EventWeight("FluxWeight*PotWeight")
    >>>                  )
    >>>

      
    """
    def __init__(self,eventWeight):
        """
        :param eventWeight: TTree::Draw compatible string corresponding to the event weight.
        :type eventWeight: str
        """
        self.eventWeight = eventWeight

    def __str__(self):
        return self.getEventWeightString()
    
    def getEventWeightString(self):
        return self.eventWeight

###############################################################################


class Cut:
    """Define a cut to select events to plot.
    
    Fill a histogram with only events passing this cut. For example, if you want to draw 
    reconstructed muon momentum ("pMu") for events passing a PID cut ("muonLikelihood > 0.5")
    your code would look similar to,
    
    >>> from commonAnalysis.plot import drawTools
    >>> from commonAnalysis.plot import drawOptions
    >>> 
    >>> tp = drawTools.TreePainter(myInputDataset,"MyTree")
    >>> canvas = tp.paint("muonMomentum",
    >>>                   "Muon momementum",
    >>>                   "pMu",
    >>>                   drawOptions.Cut("muonLikelihood > 0.5")
    >>>                  )
    >>>
    
    """
    def __init__(self, cut):
        """
        :param cut: TTree::Draw compatible string corresponding to selection cut to apply.
        :type cut: str
        """
        self.cut = cut
        
    def __str__(self):
        return str(self.cut)
        
    def getCutString(self):
        return self.cut

    @staticmethod
    def convertBinningToCuts(variableName, binning):
        """Convert a binning object to a set of cuts.
        
        Returns a set of cuts from a binning object. For example, if you want to 
        split a data set up into bins of true neutrino energy ("eNu") your code
        would look like,
        
        >>> from commonAnalysis.plot import drawTools
        >>> from commonAnalysis.plot import drawOptions
        >>> 
        >>> binning = ROOT.RooBinning(5,0.0,5.0)
        >>> cuts = drawOptions.Cut.convertBinningToCuts("eNu", binning)
        >>> tp = drawTools.TreePainter(myInputDataset,"MyTree")
        >>> canvas = tp.paint("muonMomentum",
        >>>                   "Muon momementum",
        >>>                   "pMu",
        >>>                   drawOptions.SplitDataset(cuts)
        >>>                  )
        >>>
        
        """
        cuts = dict()
        for ebin in xrange(binning.numBins()):
            eMin = binning.binLow(ebin)
            eMax = binning.binHigh(ebin)
            if ebin == 0:
                #first bin
                cut = variableName+" < "+str(eMax)
            elif ebin == (energyBinning.numBins()-1):
                #last bin
                cut = variableName+" >= "+str(eMin)
            else:
                #all other bins
                cut = str(eMin)+" <= "+variableName+" && "+variableName+" < "+str(eMax)
            cuts[ebin] = Cut(cut)
        return cuts
    
    @staticmethod
    def mergeListOfCuts(listOfCuts):
        cutStr = ""
        if len(listOfCuts)>1:
            cutStr = " && ".join(["("+cut.getCutString()+")" for cut in listOfCuts])
        elif len(listOfCuts)==1:
            cutStr = listOfCuts[0].getCutString()
        result = Cut(cutStr)
        return result
    
    def add(self, other):
        self.cut = Cut.mergeListOfCuts([self,other]).getCutString()
        return
    
###############################################################################

class SplitDataset:
    """Split the total dataset up into subsets matching certain cuts.
    
    The main use for this method is for plotting MC histograms broken up into 
    event types. For example, if you wanted to plot signal and background 
    separately, your code would look like,
    
    >>> from commonAnalysis.plot import drawTools
    >>> from commonAnalysis.plot import drawOptions
    >>> 
    >>> eventTypes = { "signal" : drawOptions.Cut("isSignal==1"),
    >>>                "background" : drawOptions.Cut("isSignal==0"),
    >>>              }
    >>> splitDataset = drawOptions.SplitDataset(eventTypes)
    >>> tp = drawTools.TreePainter(myInputDataset,"MyTree")
    >>> canvas = tp.paint("muonMomentum",
    >>>                   "Muon momementum",
    >>>                   "pMu",
    >>>                   splitDataset
    >>>                  )
    >>>
    
    """
    def __init__(self, cuts=None):
        """Construct a split dataset.
        
        You may provide the a dictionary defining the cuts here, or add them later by calling SplitDataset.add.
        
        :param cuts: an optional dictionary with key=name of sub-dataset and values=a drawOptions.Cut object that defined the selection for this subset of data.
        :type cuts: dict
        """
        self.cuts = dict()
        self.nameList = []
        if cuts:
            keys = cuts.keys()
            keys.sort()
            for name in keys:
                cut = cuts[name]
                self.add(name,cut)
                
    def add(self, name, cut):
        """Define another sub-dataset to be plotted.
        
        :param name: the name of the sub-dataset. This should be unique.
        :type name: str
        :param cut: defines the selection for this sub dataset.
        :type cut: drawOptions.Cut  
        """
        if isinstance(cut,str):
            cut = Cut(cut)
        self.nameList.append(name)
        self.cuts[name] = cut
    
    @classmethod
    def from_integer_map(cls, varname, map):
        '''Requires input which dictionary of category names mapped to integer values.'''
        cuts = dict()
        for key, val in map.iteritems():
            cuts[key] = "{0}=={1}".format(varname, val)
        return SplitDataset(cuts)
        
    def getNames(self):
        """
        :returns: a list of the names of the defined sub-datasets.
        """
        return self.nameList
    
    def __iter__(self):
        """Iterate over name,Cut pairs.
        :returns: iterator over name,cut pairs.
        """
        for n in self.nameList:
            yield (n,self.cuts[n])
            

###############################################################################

class CanvasSize:
    """Set the canvas x and/or y size in pixels.  
    """
    def __init__(self, x = 800, y = 600):
        """
        :param x: x size in pixels.
        :type x: int
        :param y: y size in pixels.
        :type y: int  
        """
        self.x = x
        self.y = y
        
    def getX(self):
        """
        :returns: an integer corresponding to x size in pixels 
        """
        return self.x

    def getY(self):
        """
        :returns: an integer corresponding to y size in pixels 
        """
        return self.y    

###############################################################################

class AxisLabels:
    """Set the plot axis labels.
    """
    def __init__(self, x="", y="", z=""):
        """
        :param x: x-axis label
        :type x: str
        :param y: y-axis label
        :type y: str
        :param z: z-axis label
        :type z: str
        """
        self.x = x
        self.y = y
        self.z = z
    
    def getX(self):
        """
        :returns: the x-axis label (of type str). 
        """
        return self.x

    def getY(self):
        """
        :returns: the y-axis label (of type str). 
        """
        return self.y
    
    def getZ(self):
        """
        :returns: the z-axis label (of type str). 
        """
        return self.z
    
    def getTitle(self):
        """Construct a title string in the format understood by ROOT objects.
        
        In ROOT the axis labels are set by parsing a semi-colon separated list in the title string set when constructing the object. 
        
        For example, when constructing a TH1F, the axis labels are set in the title string.
        
        >>> import ROOT
        >>> hist = ROOT.TH1F("histName","hist title;x axis label;y axis label;z-axis label",10,-3,3)
        >>> hist.Draw()
        
        The corresponding title string return by this method is ";x axis label;y axis label;z-axis label"
        
        :returns: the title string. 
        """
        title = ";"+self.getX()+";"+self.getY()
        if self.getZ():
            title += ";"+self.getZ()
        return title
    
    def __str__(self):
        return "AxisLabels(" + ",".join([self.x,self.y,self.z]) + ")"

###############################################################################
    
class AxisScale:
    
    def __init__(self, logX=False, logY=False, logZ=False):
        self.logX = logX
        self.logY = logY
        self.logZ = logZ
        
    def isLogY(self):
        return self.logY
    
    def isLogX(self):
        return self.logX
    
    def isLogZ(self):
        return self.logZ

###############################################################################

class FrameRange:
    """Force the x,y or z range of the plot axes.
    
    The axis ranges are by default automatically calculated. 
    With this draw option you can force them to take specific values.
    You can also specify the "headroom" (the size of the gap between the 
    maximum data point and the top of the plot, primarily to provide space for 
    the legend). 
    """
    auto = "auto"
    def __init__(self, xMin=auto, xMax=auto, yMin=auto, yMax=auto, zMin=auto, zMax=auto, headroom=0.15):
        """
        :param xMin: set the minimum edge of the x-axis (in units of the axis).
        :type xMin: float (or FrameRange.auto)
        :param xMax: set the maximum edge of the x-axis (in units of the axis).
        :type xMax: float (or FrameRange.auto)
        :param yMin: set the minimum edge of the y-axis (in units of the axis).
        :type yMin: float (or FrameRange.auto)
        :param yMax: set the maximum edge of the y-axis (in units of the axis).
        :type yMax: float (or FrameRange.auto)
        :param headroom: set the size of the gap between the maximum data point and the top edge of the plot (as a fraction of the y-value of the max data point). Primarily to give space for a legend.
        :type headroom: float

        """
        self.xMin = xMin
        self.xMax = xMax
        self.yMin = yMin
        self.yMax = yMax
        self.zMin = zMin
        self.zMax = zMax
        self.headroom = headroom
        
    def getXMin(self):
        return self.xMin        

    def getXMax(self):
        return self.xMax

    def getYMin(self):
        return self.yMin        

    def getYMax(self):
        return self.yMax

    def getZMin(self):
        return self.zMin        

    def getZMax(self):
        return self.zMax
    
    def getHeadroom(self):
        return self.headroom

###############################################################################

class NDimensions:
    """Specify the dimensionality of the variable being plotted.
    
    The drawTools assume that the plotted data is 1D. If you are plotting 2D or 
    ND data then you must manually specify with this draw option.
    """
    def __init__(self, nDimensions):
        """
        :param nDimensions: number of dimensions of the plot (eg 1 or 2).
        :type nDimensions: int
        """
        self.nDimensions = nDimensions
        return
    
    def is1D(self):
        """
        :returns: True if nDimensions is 1, False otherwise.
        """
        return self.nDimensions == 1
    
    def is2D(self):
        """
        :returns: True if nDimensions is 2, False otherwise.
        """
        return self.nDimensions == 2
    
    def is3D(self):
        """
        :returns: True if nDimensions is 3, False otherwise.
        """
        return self.nDimensions == 3
    
    def isHyperD(self):
        """
        :returns: True if nDimensions is greater than 3, False otherwise.
        """
        return self.nDimensions > 3
    
    def getNDimensions(self):
        """
        :returns: an integer nDimensions.
        """
        return self.nDimensions

###############################################################################

class BaseBinning:
    """A base class for the various binning options.
    
    Should not be directly used, use one of the child classes (eg 
    UniformBinning or VariableBinning).
    """
    def getBinningArray(self):
        """
        :returns: a c-style array containing the bin edges.
        """
        return self.rooBinning.array()
        
    def getNBins(self):
        """
        :returns: integer number of bins
        """
        return self.rooBinning.numBins()
    
    def lowBound(self):
        return self.rooBinning.lowBound()
    
    def highBound(self):
        return self.rooBinning.highBound()

class AutoBinning(BaseBinning):
    
    def __init__(self, nbins):
        self._nbins = nbins
    
    def getBinningArray(self):
        """
        :returns: a c-style array containing the bin edges.
        """
        return None
        
    def getNBins(self):
        """
        :returns: integer number of bins
        """
        return self._nbins
    
    def lowBound(self):
        return None
    
    def highBound(self):
        return None

class UniformBinning(BaseBinning):
    """Draw option to set uniform histogram binning.
    """
    def __init__(self, nBins, xMin, xMax):
        """
        :param nBins: number of bins (uniformly distributed within the range)
        :type nBins: int
        :param xMax: high end of the axis range.
        :type xMax: int
        :param xMin: low end of the axis range.
        :type xMin: int                
        """
        self.nBins = nBins
        self.xMin = xMin
        self.xMax = xMax
        self.rooBinning = ROOT.RooBinning(nBins, xMin, xMax)

class VariableBinning(BaseBinning):
    """Draw option to set variable width histogram bins.
    """ 
    def __init__(self, binEdges):
        """
        :param binEdges: an ordered list of bin edges
        :type binEdges: list or tuple of floats.
        """
        nBins = len(binEdges)-1
        self.nBins = nBins
        self.xMin = min(binEdges)
        self.xMin = max(binEdges)
        self.binEdges = binEdges
        self.arrayEdges = numpy.array(binEdges)
        self.rooBinning = ROOT.RooBinning(nBins, self.arrayEdges)
        
class RooBinningWrapper(BaseBinning):
    """Specify binning directly with a ROOT.RooBinning object.
    
    A wrapper around a RooBinning object so that it can be used as a BaseBinning
    type draw option.
    """
    def __init__(self, rooBinning):
        """
        :type rooBinning: ROOT.RooBinning
        """
        self.rooBinning = rooBinning

###############################################################################

class YBinning:
    """Specify the y-binning of a 2D histogram.
    
    By default the binning options apply to x-axis only. To set the y-axis 
    binning you should wrap your binning option in a YBinning option. For example,
    
    
    
    """
    def __init__(self, binningObject):
        self.binningObject = binningObject
    
    def getBinningObject(self):
        return self.binningObject

###############################################################################

class ZBinning:
    """Specify the z-binning of a 3D histogram.
    
    By default the binning options apply to x-axis only. To set the z-axis 
    binning you should wrap your binning option in a ZBinning option.
    """
    def __init__(self, binningObject):
        self.binningObject = binningObject
    
    def getBinningObject(self):
        return self.binningObject

###############################################################################

class Stack:
    """Draw options to set stacking mode of multiple datasets.
    
    Set how to control how plots with multiple datasets are simultaneously plotted.
    The supported modes are overlay (datasets are drawn on top of each other, no stacking) 
    and stackMC (stack the MC datasets, data is never stacked).
    
    When plotting stacked histograms optionally the order in which plots are 
    stacked can be automatically determined depending on the orderMode. 
    "orderbyInsertion" orders by the order histograms were added to the 
    HistogramCollection.
    "orderByIntegral" places the histogram with the largest area on top.
    Finally, the user can set the stacking order by hand by providing a list of 
    histogram/sub-dataset names. 
     
    For example, if I wanted to draw a stacked histogram with my signal on top, 
    the combinatorial-background next and the order the remaining backgrounds 
    by integrated number of events the code would look like,
    
    >>> from commonAnalysis.plot import drawTools
    >>> from commonAnalysis.plot import drawOptions
    >>> cuts = { "signal" : drawOptions.Cut("eventType==1"),
    >>>          "combinatorialBg" : drawOptions.Cut("eventType==2"),
    >>>          "backgroundType1" : drawOptions.Cut("eventType==3"),
    >>>          "backgroundType2" : drawOptions.Cut("eventType==4"),
    >>>        }
    >>> split = drawOptions.SplitDataset(cuts)
    >>> userOrder = ["signal", "combinatorialBg"]
    >>> stack = drawOptions.Stack(drawOptions.Stack.orderByIntegral,userOrder)
    >>> tp = drawTools.TreePainter(myInputDataset,"MyTree")
    >>> canvas = tp.paint("muonMomentum",
    >>>                   "Muon momementum",
    >>>                   "pMu",
    >>>                   split,
    >>>                   stack,
    >>>                  )
    >>>
    
    """
    overlay = "overlay"
    stackMC = "stackMC"
    
    orderByInsertion = "orderByInsertion"
    orderByIntegral = "orderByIntegral"
    
    _validModes = set([overlay,stackMC])
    
    _validOrderingModes = set([orderByInsertion,orderByIntegral]) 
    
    def __init__(self, mode, orderMode=orderByInsertion, userOrder=[]):
        """
        :param mode: Set the stacking mode. This should be one of the modes defined within this class eg Stack.overlay, Stack.stackMC.
        :type mode: str
        :param orderMode: Set automatic ordering mode when stacking. This should be one of the modes defined within this class eg Stack.orderByInsertion, Stack.orderByIntegral.
        :type orderMode: str
        :param userOrder: manually set the stacking order by providing an ordered list of histogram/sub-dataset names. The order is provided from top to bottom.
        :type userOrder: list or tuple
        """
        self.mode = mode
        self.setUserOrder(userOrder)
        self.setOrderMode(orderMode)
        #check inputs
        if mode not in Stack._validModes:
            raise InvalidDrawOption("unknown Stack option set",mode)
        
    def getStackMode(self):
        return self.mode
    
    def setUserOrder(self, listOfNames):
        self.userOrder = listOfNames
        
    def getUserOrder(self):
        return self.userOrder
    
    def setOrderMode(self, orderMode):
        self.orderMode = orderMode
    
    def getOrderMode(self):
        return self.orderMode
    
    def __str__(self):
        sio = StringIO.StringIO()
        print >>sio,"StackOption(",
        print >>sio,"mode =",self.getMode(),
        print >>sio,", order =",str(self.getOrder()),
        print >>sio,")",
        return sio.getvalue()

###############################################################################

class TreatAsData:
    """Draw option to define which datasets are data.
    
    The drawing tools sometimes need to treat data and MC differently. 
    With this drawing option you can inform the drawing tools which datasets are
    data. For example,
    
    >>> from commonAnalysis.plot import drawTools
    >>> from commonAnalysis.plot import drawOptions
    >>> cuts = { "data" : drawOptions.Cut("isData==1"),
    >>>            "mc" : drawOptions.Cut("isData==0"),
    >>>        }
    >>> split = drawOptions.SplitDataset(cuts)
    >>> treatAsData = drawOptions.TreatAsData("data")
    >>> tp = drawTools.TreePainter(myInputDataset,"MyTree")
    >>> canvas = tp.paint("muonMomentum",
    >>>                   "Muon momementum",
    >>>                   "pMu",
    >>>                   split,
    >>>                   treatAsData,
    >>>                  )
    >>>
    
    """
    def __init__(self, *names):
        """
        :param names: the names of the histograms/sub-datasets that are to be plotted as data.
        :type names: str
        """
        self.names = set(names)
    
    def __contains__(self, name):
        return name in self.names
    
    def __iter__(self):
        return iter(self.names)

###############################################################################

class LegendPosition:
    """Manually place the legend.
    
    The drawing tools witll automatically be placed. Use this drawing option 
    to overide the automatic placement and manually set the legend position 
    within the canvas. 
    """
    def __init__(self, x=None, y=None, width=0.2, height=0.2, doDraw=True):
        """
        
        All numbers are given as a fraction of the pad size in the relevant 
        dimension. If the legend will not fit on the canvas it will 
        automatically be shifted to fit on the canvas. 
        
        :param x: x coordinate of bottom left corner of the legend.
        :type x: float
        :param y: y coordinate of bottom left corner of the legend.
        :type y: float
        :param width: x-size of the legend.
        :type width: float
        :param height: y-size of the legend.
        :type height: float
        """
        self.x = x
        self.y = y
        self.width = width
        self.height = height
        self._doDrawFlag = doDraw
        
    def hasUserLimits(self):
        return not(self.x == None or self.y == None)
        
    def calculateLegendLimits(self):
        """
        :returns: legend position values in the format expected by the ROOT.TLegend constructor.
        """
        x = self.x
        y = self.y
        width = self.width
        height = self.height
        #check it will fit inside the area
        if x+width>1.0:
            shift = 1.0-(x+width)
            x += shift
        if y+height>1.0:
            shift = 1.0-(y+height)
            y += shift
        #calculate limits needed by TLegend constuctor
        xLow = x
        xHigh = x + width
        yLow = y
        yHigh = y + height
        return xLow,yLow,xHigh,yHigh
    
    def doDraw(self):
        return self._doDrawFlag

###############################################################################

class BaseFormat:
    """An abstract interface for draw options for formatting histograms.
    
    This is in place for future proofing in case alternative formatting classes 
    are implemented.
    """
    def __init__(self):
        pass
    
    def format(self, name, histogram, isData):
        """Format a histogram.
        
        All formating classes should provide an implementation for this method.
        
        :param name: name of the histogram/data set
        :type name: str
        :param histogram: the histogram to be formatted.
        :type histogram: ROOT.TH1
        :param isData: True if the histogram is to be treated as data. False if it is to be treated as MC.
        :type isData: bool
        """
        return

###############################################################################

class Format(BaseFormat):
    """The nominal formatting class.
    
    Automatically colors and styles lines, markers and fill areas for MC and 
    data.
    """
    def __init__(self, fillArea=False, dashedLines=False):
        """
        :param fillArea: switch on/off setting of histogram fill areas.
        :type fillArea: bool
        :param dashedLines: switch on/off setting of line style.
        :type dashedLines: bool
        """
        BaseFormat.__init__(self)
        self.fillArea = fillArea
        self.colors = palette.ColorPalette.createPaletteCycle()
        self.dataColors = itertools.chain( [1],self.colors )
        self.marker = itertools.cycle(xrange(22,35))
        self.dataMarker = itertools.cycle(xrange(20,35))
        self.fillStyle = itertools.cycle([3145, 3154, 3105, 3195])
        if dashedLines:
            self.lineStyle = itertools.cycle([1,3,4,6,2])
        else:
            self.lineStyle = itertools.cycle([1])
        self.userDefined = dict()
        
    def set(self, name, color=None, markerStyle=None, lineStyle=None):
        """Override default behaviour for a particular histogram/sub-dataset.
        :param name: name of the histogram.
        :type name: str
        :param color: set ROOT color code for fill/line/marker color.
        :type color: int
        :param markerStyle: set ROOT marker style. See ROOT.TAttMarker.
        :type markerStyle: int
        :param lineStyle: set ROOT line style. See ROOT.TAttLine.
        :type lineStyle: int
        """
        self.userDefined[name] = (color,markerStyle,lineStyle)
    
    def format(self, name, histogram, isData):
        """Implementation of auto-formatting algorithm.
        
        See BaseFormat for more details on arguments.
        """
        #get colors options
        if isData:
            col = self.dataColors.next()
            mkr = self.dataMarker.next()
            line = 1
        else:
            col = self.colors.next()
            mkr = self.marker.next()
            line = self.lineStyle.next()
        #get user overrides
        if name in self.userDefined:
            userCol,userMkr,userLine = self.userDefined[name]
            if userCol is not None:
                col = userCol
            if userMkr is not None:
                mkr = userMkr
            if userLine is not None:
                line = userLine                
        #print "DEBUG : format",name,col,mkr,line
        #set histogram attributes
        histogram.SetLineColor(col)
        if not isData and self.fillArea:
            style = self.fillStyle.next()
            histogram.SetFillColor(col)
            histogram.SetFillStyle(style)
            histogram.SetMarkerColor(col)
            histogram.SetMarkerStyle(0)
        elif isData:
            histogram.SetMarkerColor(col)
            histogram.SetMarkerStyle(mkr)
            histogram.SetLineStyle(line)
        else:
            histogram.SetMarkerColor(col)
            histogram.SetMarkerStyle(0)
            histogram.SetLineStyle(line)            
        return

###############################################################################

class Normalisation:
    """TODO : still needs to be implemented.
    """
    noNormalisation = "noNormalisation"
    pot = "pot"
    mcToData = "mcToData"
    unitArea = "unitArea"
    def __init__(self, mode):
        self.mode = mode
        return
    
    def getMode(self):
        return self.mode

###############################################################################

class ShowOverflow:
    def __init__(self, flag):
        self.flag = flag

###############################################################################    
    
_allValidOptionsTypes = set([EventWeight, 
                             Cut, 
                             SplitDataset, 
                             CanvasSize, 
                             AxisLabels, 
                             FrameRange, 
                             NDimensions, 
                             BaseBinning, 
                             UniformBinning,
                             VariableBinning,
                             RooBinningWrapper,
                             YBinning,
                             ZBinning,
                             Stack,
                             TreatAsData,
                             LegendPosition,
                             BaseFormat,
                             Format,
                             Normalisation,
                             AxisScale,
                             ]
                            )
