import os
import ROOT

class CanvasWriter:
    '''Manages writing of ROOT TCanvas to files.
    
    For example,
    
    >>> from commonAnalysis.plot import io
    >>> 
    >>> #make a very simple plot
    >>> canvas = ROOT.TCanvas()
    >>> hist = ROOT.TH1F("hist","hist",100,-3,3)
    >>> hist.FillRandom("gaus",10000)
    >>> hist.Draw()
    >>>
    >>> #save the plot to a directory "myPlots"
    >>> output = io.CanvasWriter(outputPath="myPlots", extensions=["png","eps"])
    >>> output.save(canvas)
    '''
    def __init__(self,outputPath="./plots",extensions=["png","eps", "C"], debug=False):
        """
        :param outputPath: output directory to store plots (created if it does not exist). Plots are stored in sub directories by file extension eg outputPath/png
        :type outputPath: str
        :param extensions: a list of strings eg ["png","eps"]
        :type extensions: list
        :param debug: switch on/off debug messages.
        :type debug: bool 
        """
        self.setOutputPath(outputPath)
        self.setExtensions(extensions)
        self.normalErrorIgnoreLevel = ROOT.kInfo
        self.debug = debug
    def setExtensions(self,ext):
        """Set file types to save to.
        :param ext: a list of strings eg ["png","eps"]
        :type ext: list
        """
        self.extensions = ext
    def setOutputPath(self,path):
        """Set folder name to write plots to.
        :param path: output directory to store plots (created if it does not exist). Plots are stored in sub directories by file extension eg outputPath/png
        :type path: str
        """
        self.outputPath = os.path.expanduser(path)
    def save(self,*canvases):
        """Save the canvas. Multiple copies are created, 1 for each file extension.
        :param canvases: a plot or list of plots to write out.
        :type: ROOT.TCanvas or list of ROOT.TCanvas
        """
        for canv in canvases:
            if self.debug:
                print "CanvasWriter(canv:",canv,", path:",self.outputPath,", extensions:", self.extensions,")"
            self._silenceRoot()
            try: iter(canv);
            except TypeError: canv = [canv] #put inside list
                #if ext == "root":
                #    pass # TODO
                #else:
            for c in canv:
                for ext in self.extensions:
                    if c:
                        fname = self.outputPath+"/"+ext+"/"+c.GetName() + "." + ext
                        try:
                            os.makedirs(os.path.dirname(fname))
                        except OSError:
                            pass #error thrown if directory already exists, but this is fine. 
                        if self.debug:
                            print "saving ",c,"to",fname                       
                        c.SaveAs(fname)
            self._unsilenceRoot()
    def _silenceRoot(self):
        """Prevent ROOT from printing a message everytime a plot is saved."""
        if not self.debug:
            self.normalErrorIgnoreLevel = ROOT.gErrorIgnoreLevel
            ROOT.gErrorIgnoreLevel = ROOT.kError
    def _unsilenceRoot(self):
        """Undo _silenceRoot."""
        if not self.debug:
            ROOT.gErrorIgnoreLevel = self.normalErrorIgnoreLevel
