# cython: profile=True

import itertools
import collections
from operator import attrgetter

import numpy as np
import ROOT

from simplot.progress import printprogress

###############################################################################

def calculate_statistics_from_toymc(toymc, statistics, npe, transform=None, name=None):
    if transform is None:
        transform = attrgetter("vec")
    return calculate_statistics(toymc, statistics, npe, transform=transform, name=name)

def calculate_statistics(generator, statistics, npe, transform=None, name=None):
    try:
        iter(statistics)
    except TypeError:
        statistics = [statistics]
    iterable = xrange(npe)
    if name is not None:
        iterable = printprogress(name, npe, iterable, update=True)
    for _ in iterable:
        exp = generator()
        if transform:
            exp = transform(exp)
        for s in statistics:
            s.add(exp)
    return statistics

###############################################################################

#def _set_naninf(arr, val=0.0):
#    #remove NaN and inf values
#    arr[~np.isfinite(arr)] = val
#    return arr

###############################################################################

def safedivide(lhs, rhs):
    '''Divide by zero is set to 0.0.'''
    l = np.array(lhs)
    r = np.array(rhs)
    l[r==0] = 0.0
    r[r==0] = 1.0
    return l / r

###############################################################################

class Mean(object):
    def __init__(self):
        self._rms = StandardDeviation()
    
    def add(self, vec):
        self._rms.add(vec) 
    
    def eval(self):
        return self._rms._mean()
    
    def err(self):
        rms = self._rms.eval()
        return np.divide(rms, np.sqrt(self._rms._count))

###############################################################################

class StandardDeviation(object):
    def __init__(self, ndbinning=None, projection=None):
        self._ndbinning = ndbinning
        self._projection = projection
        self._sumw = None
        self._sumw2 = None
        self._count = 0
    
    def add(self, vec):
        self._count += 1
        try:
            self._sumw += vec
            self._sumw2 += np.power(vec, 2)
        except TypeError:
            self._sumw = np.copy(vec)
            self._sumw2 = np.power(vec, 2)
        return
    
    def _mean(self):
        return np.divide(self._sumw, float(self._count))
        
    def eval(self):
        mean = self._mean()
        mean2 = np.power(mean, 2)
        avgsum2 = self._sumw2 / float(self._count)
        rms = np.sqrt(np.abs(avgsum2 - mean2))
        return rms

    def err(self):
        rms = self.eval()
        return np.divide(rms, np.sqrt(2.0*self._count))
    
###############################################################################

class Covariance(object):
    def __init__(self, fractional=False):
        self._sumw = None
        self._sumw2 = None
        self._count = 0
        self._fractional = fractional
    
    def add(self, vec):
        vec = np.array(vec, copy=False) # convert to vec
        self._count += 1
        try:
            self._sumw += vec
            self._sumw2 += self._product(vec, vec)
        except TypeError:
            self._sumw = np.copy(vec)
            self._sumw2 = self._product(vec, vec)
        return
    
    def _product(self, a, b):
        #return np.array([bi * a for bi in b])
        ax = a.reshape((a.shape[0], 1))
        bx = b.reshape((b.shape[0], 1))
        return np.dot(ax, bx.T)
    
    def _mean(self):
        return np.divide(self._sumw, float(self._count))
        
    def eval(self):
        mean = self._mean()
        mean2 = self._product(mean, mean)
        avgsum2 = self._sumw2 / float(self._count)
        rms = avgsum2 - mean2
        if self._fractional:
            rms /= mean2
        return rms

    def err(self):
        cov = self.eval()
        result = np.copy(cov)
        for ii, jj in itertools.product(xrange(cov.shape[0]), xrange(cov.shape[1])):
            result[ii,jj] = self._calcerr(cov, ii, jj, self._count)
        return result

    def _calcerr(self, cov, ii, jj, N):
        if ii == jj:
            return np.multiply(np.sqrt(2), np.divide(cov[ii][jj], np.sqrt(N)))
        else:
            #no correlations
            #return np.divide(np.sqrt(cov[ii][ii]*cov[jj][jj]), np.sqrt(N))
            #with correlations
            uncorrterm = np.divide(np.array(cov[ii][ii]*cov[jj][jj]), N)
            corr = cov[ii][jj] / np.sqrt(cov[ii][ii]*cov[jj][jj])
            errsq = (1.0 + abs(corr)) * uncorrterm
            ret = np.sqrt(errsq)
            return ret

###############################################################################

class FractionalStandardDeviation(StandardDeviation):
    def __init__(self, ndbinning=None, projection=None):
        super(FractionalStandardDeviation, self).__init__(ndbinning, projection)
    
    def eval(self):
        rms = super(FractionalStandardDeviation, self).eval()
        return safedivide(rms, self._mean())

###############################################################################

class RootHistogram(object):
    def __init__(self, name, title, nbins=100, range=None, allowrebin=True):
        self._name = name
        self._title = title
        self._nbins = nbins
        self._range = range
        self._allowrebin = allowrebin
        self._hist = collections.OrderedDict()
        
    def _make_hist(self, i):
        name = "_".join([self._name, str(i)])
        title = "_".join([self._title, str(i)])
        nbins = self._nbins
        if self._range is not None:
            try:
                xmin, xmax = self._range[i]
            except TypeError:
                xmin, xmax = self._range
            hist = ROOT.TH1F(name, title, nbins, xmin, xmax)
        else:
            hist = ROOT.TH1F(name, title, nbins, 0.0, 0.0)
        if self._allowrebin:
            hist.SetBit(ROOT.TH1.kCanRebin)
        hist.SetDirectory(0)
        hist.Sumw2()
        return hist

    def add(self, vec):
        for i, v in enumerate(vec):
            try:
                self._hist[i].Fill(v)
            except KeyError:
                h = self._make_hist(i)
                h.Fill(v)
                self._hist[i] = h
        return
    
    def eval(self):
        return self._hist.values()

