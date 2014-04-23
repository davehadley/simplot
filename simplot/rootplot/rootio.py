import os
import ROOT # @UnresolvedImport
import itertools
import simplot.table

class CanvasWriter:

    def __init__(self, path="./", extensions=["png","eps", "C"], verbose=True, debug=False):
        """
        :param _path: output directory to store plots (created if it does not exist). Plots are stored in sub directories by file extension eg _path/png
        :type _path: str
        :param extensions: a list of strings eg ["png","eps"]
        :type extensions: list
        :param debug: switch on/off debug messages.
        :type debug: bool 
        """
        self._setpath(path)
        self._setext(extensions)
        self._normal_error_ignore_level = ROOT.kInfo
        self._debug = debug
        self._verbose = verbose
        
    def _setext(self, ext):
        self._extensions = ext
        
    def _setpath(self, path):
        self._path = os.path.expanduser(path)
        
    def save(self, *canvases):
        for it in canvases:
            try:
                iterable = iter(it)
            except TypeError:
                #input is not interable, assume that it is just a single canvas
                iterable = [it]
            for canv in iterable:
                if self._debug:
                    print "CanvasWriter(canv:", canv, ", _path:", self._path, ", extensions:", self._extensions, ")"
                if self._verbose:
                    print "saving plot", canv.GetName()
                self._silence_root()
                try: 
                    iter(canv);
                except TypeError: 
                    canv = [canv] #put inside list
                for c in canv:
                    for ext in self._extensions:
                        if c:
                            fname = self._path+"/"+ext+"/"+c.GetName() + "." + ext
                            try:
                                os.makedirs(os.path.dirname(fname))
                            except OSError:
                                pass #error thrown if directory already exists, but this is fine. 
                            if self._debug:
                                print "saving ",c,"to",fname
                            if ext.upper() == "C":
                                #special case, rename the object so that it is valid c++
                                name = c.GetName()
                                name = name.replace("/", "_")
                                c.SetName(name)
                            c.SaveAs(fname)
                self._unsilence_root()
            
    def _silence_root(self):
        """Prevent ROOT from printing a message everytime a plot is saved."""
        if not self._debug:
            self._normal_error_ignore_level = ROOT.gErrorIgnoreLevel
            ROOT.gErrorIgnoreLevel = ROOT.kError
            
    def _unsilence_root(self):
        """Undo _silenceRoot."""
        if not self._debug:
            ROOT.gErrorIgnoreLevel = self._normal_error_ignore_level
            

###############################################################################

class TableHistogramOutput(simplot.table.TableOutputBase):
    def __init__(self, formatter=None):
        super(TableHistogramOutput, self).__init__(formatter)
    
    def write(self, table, filename):
        try:
            os.makedirs(os.path.dirname(filename))
        except os.error:
            pass
        self._savecanvas(table, filename)
        return
    
    def _savecanvas(self, table, filename, textformat=None):
        canv = ROOT.TCanvas(table.get_name(), table.get_name(), 800, 600)
        if textformat:
            startingTextFormat = ROOT.gStyle.GetPaintTextFormat()
            ROOT.gStyle.SetPaintTextFormat(textformat)
        hist = self._colz_plot(table)
        hist.Draw("COLZ,TEXT")
        canv.SaveAs(filename)
        if textformat:
            ROOT.gStyle.SetPaintTextFormat(startingTextFormat)
        return

    def _hasrowlabels(self, table):
        return any( (len(r)>0 and type(r[0]) is str for r in table.get_rows()) )

    def _colz_plot(self, table, minZ=None, maxZ=None):
        name = table.get_name()
        nRows = table.get_nrows()
        nCols = table.get_ncols()
        hasRowLabels = self._hasrowlabels(table)
        nY = nRows
        nX = nCols
        if hasRowLabels:
            nX = nCols - 1
        hist = ROOT.TH2F("hist"+name,"",nX,0,nX,nY,0,nY)
        hist.SetDirectory(0)
        autoMax = minZ
        rows = table.get_rows()
        for iX, iY in itertools.product(range(nX),range(nY)):
            if hasRowLabels:
                indexX = iX + 1
            else:
                indexX = iX
            indexY = iY
            binX = iX + 1
            binY = nRows - indexY
            entry = rows[indexY][indexX]
            #convert to floating point if we can
            value,error = self._get_value_and_error(entry)
            autoMax = max(value,autoMax)
            hist.SetBinContent(binX,binY,value)
            if error is not None:
                hist.SetBinError(binX,binY,value)
        #set axis labels
        if table.get_headrow():
            for i,xLabel in enumerate(table.get_headrow()):
                hist.GetXaxis().SetBinLabel(i+1,xLabel)
        if hasRowLabels:
            for iY in xrange(nY):
                indexY = iY
                yLabel = rows[indexY][0]
                binY = nRows - indexY
                hist.GetYaxis().SetBinLabel(binY,yLabel)
        if minZ is not None:
            hist.SetMaximum(minZ)
            maxZ = autoMax
        if maxZ is not None:
            hist.SetMaximum(maxZ)
        return hist

###############################################################################