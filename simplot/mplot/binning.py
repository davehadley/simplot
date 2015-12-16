import numpy

class Binning:
    def __init__(self, *args):
        binedges = None
        if len(args)==1:
            #list of bin edges
            binedges = numpy.array(args[0], dtype=float)
        elif len(args)==3:
            nbins, xmin, xmax = args
            binedges = numpy.linspace(xmin, xmax, num=nbins+1)
        else:
            raise TypeError("Binning initialised with wrong number of arguments.", args)
        self._binedges = binedges
        self._numbins = len(binedges) - 1
        self._lowedges = binedges[:-1]
        self._highedges = binedges[1:]
        self._verify()
        
    def _verify(self):
        if not numpy.all(self._binedges == numpy.sort(self._binedges)):
            raise ValueError("Binning binedges constructed with a non-sorted list") 
    
    def __iter__(self):
        return iter(xrange(len(self)))
    
    def __len__(self):
        return self._numbins
    
    def binnumber(self, x):
        #Brute force method
        #this is inefficient as it always tests every bin.
        #binary tree method may be faster.
        #return ((x >= self._lowedges) & (x < self._highedges)).nonzero()[0]
        
        #Binary Search Method
        # returns bin number (0, nbins-1). Ouside range returns -1 for below, and nbins for high. 
        r = self._numbins
        if x < self.max():
            r = numpy.searchsorted(self._binedges, x, side="right")-1
        return r
    
    def min(self):
        return self._binedges[0]
    
    def max(self):
        return self._binedges[-1]
    
    def binlow(self, binnum):
        if binnum < 0 or binnum > self._numbins:
            raise IndexError("out of range", binnum)
        return self._binedges[binnum]
    
    def binhigh(self, binnum):
        if binnum < -1 or binnum >= self._numbins:
            raise IndexError("out of range", binnum)
        return self._binedges[binnum+1]
    
    def bincentre(self, binnum):
        if binnum < 0 or binnum > self._numbins:
            raise IndexError("out of range", binnum)
        return numpy.average(self._binedges[binnum:(binnum+2)])

    def __eq__(self, rhs):
        return numpy.all(self._binedges == rhs._binedges)
    
    def __str__(self):
        edgesstr = ", edges=" + ",".join(("%.2e"%e for e in self._binedges))
        return "".join(("Binning(n=",
                        str(len(self)),
                        edgesstr,
                        ")",
                       ))

    @property
    def binedges(self):
        return self._binedges

