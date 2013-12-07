import os

import matplotlib.pyplot as plt

class FigureWriter:
    def __init__(self, path="./plots", extensions=["png", "eps"], show=False):
        self._ext = extensions
        self._path = self._expand(path)
        self._show = show
        
    def __call__(self, *args):
        """Valid call signatures:
         a single 2-element tuple
         iterable that yields 2 element tuples
         
         The tuples must contain:
         name (str) and figure (matplotlib.figure.Figure) in that order."""
        if len(args) == 1:
            iterable = args[0]
            for name, fig in iterable:
                if self._show:
                    plt.show()
                    raw_input("waiting on plot " + name)
                self._save(name, fig)
        elif len(args) == 2:
            name, fig = args
            self._save(name, fig)
        else:
            raise TypeError("FigureWriter requires 1 (iterable over tuples) or 2 arguments (tuple of name,figure)")
    
    def _save(self, name, fig):
        for ext in self._ext:
            fname = "".join((self._path,
                             os.path.sep,
                             ext,
                             os.path.sep,
                             name,
                             ".",
                             ext,
                            )) 
            fname = self._expand(fname)
            #ensure that the parent directory 
            self._make_parent_dir(fname)
            print "Saving",fname
            fig.savefig(fname)
        return
    
    def _make_parent_dir(self, fname):
        dirname = os.path.dirname(fname)
        try:
            os.makedirs(dirname)
        except os.error:
            pass # if path already exists os.error raised, ignore this
        return
    
    def _expand(self, path):
        return os.path.expandvars(os.path.expanduser(path))
