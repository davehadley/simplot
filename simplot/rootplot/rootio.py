import os
import ROOT # @UnresolvedImport

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
        for canv in canvases:
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
