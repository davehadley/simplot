import ROOT
import StringIO
import numpy
import time
import datetime
from mcmcfitter.lib import cmcmc

from simplot.rootplot import ntuple
from simplot import progress
import itertools
import collections

##############################################################################

_gaus = ROOT.TMath.Gaus
_rand = ROOT.TRandom3(14142)


##############################################################################

class McMcSetupError(Exception):
    pass

##############################################################################

class Parameter(object):
    def __init__(self, name, start, low, high):
        self.name = name
        self.start = start
        self.low = low
        self.high = high

    def infoString(self):
        return "Parameter(\"{}\", start={:.2g}, low={:.2g}, high={:.2g})".format(self.name, self.start, self.low, self.high)

    def __str__(self):
        return self.infoString()

##############################################################################

class NamedParameters(object):
    def __init__(self, pars=None):
        self._pars = collections.OrderedDict()
        if pars:
            for n in sorted(pars.iterkeys()):
                s, l, h = pars[n]
                self.add(n, s, l, h)
        return

    def __len__(self):
        return len(self._pars)
                
    def start(self):
        return [v.start for v in self._pars.itervalues()]
    
    def add(self, name, start, low, high):
        self._pars[name] = Parameter(name, start, low, high)
        return
    
    def parameter_ranges(self):
        return [(v.low, v.high) for v in self._pars.itervalues()]
    
    def __iter__(self):
        return self._pars.itervalues()

    def infoString(self):
        sio = StringIO.StringIO()
        print >>sio, "NamedParameters"
        for i, (n, s, r) in enumerate(zip(self._pars.keys(), self.start(), self.parameter_ranges())):
            print >>sio, "    ", i, ":", (n, s, r)
        return sio.getvalue()

    def __str__(self):
        return self.infoString()

##############################################################################

class ListDataset(object):
    def __init__(self, burnIn=1000):
        self.burnIn = burnIn
        self.dataset = []
        self.likelihood = []
        
    def add(self, parameters, likelihood):
        self.dataset.append(parameters)
        self.likelihood.append(likelihood)
        
    def get(self, axis=None):
        if axis is None:
            r = self.dataset[self.burnIn:]
        else:
            r = []
            for i in xrange(self.burnIn,len(self.dataset)):
                x = self.dataset[i]
                y = [x[a] for a in axis]
                r.append(y)
        return r
    
    def getVariable(self, index):
        r = []
        for i in xrange(self.burnIn, len(self.dataset)):
            x = self.dataset[i]
            y = x[index]
            r.append(y)
        return r
    
    def finalize(self):
        return

##############################################################################

class RootDataSet(object):

    def __init__(self, file_name, parameters, tree_name="mcmc", prescale=None):
        self._tfile = self._open_file(file_name)
        self._branches = []
        self._tree = self._create_tree(self._tfile, tree_name, parameters)
        self._prescale = prescale
        self._counter = 0

    def _open_file(self, file_name):
        tfile = ROOT.TFile(file_name, "RECREATE")
        if not tfile.IsOpen():
            raise IOError("cannot open ROOT file", file_name)
        return tfile

    def add(self, parameters, likelihood):
        if (not self._prescale) or (self._counter%self._prescale == 0):
            self._likelihood.setvalue(likelihood)
            for v, b in itertools.izip(parameters, self._branches):
                b.setvalue(v)
            self._filltree()
            for b in self._branches:
                b.reset()
        self._counter += 1
        return

    def _filltree(self):
        self._tree.Fill()

    def _create_tree(self, tfile, tree_name, parameters):
        tfile.cd()
        tree = ROOT.TTree(tree_name, tree_name)
        self._likelihood = ntuple.BranchPrimitive("likelihood", tree, 0.0)
        for p in parameters:
            b = ntuple.BranchPrimitive(p.name, tree, p.start)
            self._branches.append(b)
        return tree

    def finalize(self):
        self._tree.Write()
        self._tfile.Close()
        return

##############################################################################

class CRootDataSet(object):
    '''A faster (but harder to debug) C++ based implementation of RootDataSet.'''
    _default_tree_name = "mcmc"
    def __init__(self, file_name, parameters, tree_name=None, prescale=None, maxpars=None):
        if tree_name is None:
            tree_name = CRootDataSet._default_tree_name
        pvec = ROOT.std.vector(ROOT.std.string)()
        for p in parameters:
            s = ROOT.std.string(p.name)
            pvec.push_back(s)
        if prescale is None:
            prescale = 1
        if maxpars is None:
            maxpars = len(parameters)
        self._maxpars = maxpars
        self._cobj = cmcmc.RootDataSet(file_name, pvec, tree_name, prescale)

    def add(self, parameters, likelihood):
        self._cobj.add(max(len(parameters),self._maxpars), parameters, likelihood)
        return

    def finalize(self):
        self._cobj.finalize()
        return

    class open:
        def __init__(self, file_name, parameters, tree_name=None, prescale=None, maxpars=None):
            self._file_name = file_name
            self._parameters = parameters
            self._tree_name = tree_name
            self._prescale = prescale
            self._cobj = None
            self._maxpars = maxpars
        def __enter__(self):
            self._cobj = CRootDataSet(self._file_name, self._parameters, self._tree_name, prescale=self._prescale, maxpars=self._maxpars)
            return self._cobj
        def __exit__(self, exc_type, exc_val, exc_tb):
            if self._cobj:
                self._cobj.finalize()


##############################################################################

class MetropolisHastingsMC(object):

    def __init__(self, alg, dataset, quiet=False):
        self.alg = alg
        self.burnIn = 0
        self.dataset = dataset
        self._quiet = quiet

    def setBurnIn(self, n):
        self.burnIn = n

    def generate(self, n, maxseconds=None):
        if maxseconds is None:
            self._generatefixed(n)
        else:
            self._generatetimed(maxseconds=maxseconds, maxevents=n)

    def _generatetimed(self, maxseconds, maxevents=None):
        start = time.time()
        end = start + maxseconds
        count = 0
        if not self._quiet:
            fmt = "%Y-%m-%d %H:%M:%S"
            startstr = datetime.datetime.fromtimestamp(start).strftime(fmt)
            endstr = datetime.datetime.fromtimestamp(end).strftime(fmt)
            print "MetropolisHastingsMC start =", startstr, " expected end =", endstr
        while time.time() < end:
            parameters, likelihood = self.alg.generate()
            self.dataset.add(parameters, likelihood)
            if maxevents is not None and count > maxevents:
                break
        if not self._quiet:
            print "efficiency =",(self.alg.efficiency() * 100.), "%"
        self.dataset.finalize()
        return

    def _generatefixed(self, n):
        if self._quiet:
            iterN = xrange(n)
        else:
            iterN = progress.printprogress("MetropolisHastingsMC", n, xrange(n))
        for i in iterN:
            parameters, likelihood = self.alg.generate()
            self.dataset.add(parameters, likelihood)
        if not self._quiet:
            print "efficiency =",(self.alg.efficiency() * 100.), "%"
        self.dataset.finalize()
        return

##############################################################################

class MetropolisHastingsAlgorithm(object):
    def __init__(self, function, proposal, start, parameter_range):
        self._start = numpy.array(start)
        self._low = numpy.array([x[0] for x in parameter_range])
        self._high = numpy.array([x[1] for x in parameter_range])
        self._theta = self._start
        #if not self._function:

        self._function = function
        self._proposal = proposal
        self._total = 0
        self._success = 0
        #compute the intial likelihood values
        self._likelihood = self._function(self._theta) 
        return
    
    # def generate(self):
    #     theta0 = self._theta
    #     thetaProposal = self._proposal.generate(theta0)
    #     #self._print_vec(theta0, "theta0")
    #     #self._print_vec(thetaProposal, "thetaProposal"),
    #     #self._print_vec(self._low, "low")
    #     #self._print_vec(self._high, "high")
    #     #raw_input("wait")
    #     if cmcmc.check_range(len(thetaProposal), thetaProposal, len(self._low), self._low, len(self._high), self._high):
    #         pt = self._function(thetaProposal)
    #         jt = self._proposal.logDensity(thetaProposal, theta0)
    #         pt0 = self._function(theta0)
    #         jt0 = self._proposal.logDensity(theta0, thetaProposal)
    #         r = pt - jt - pt0 + jt0
    #         prob = min(numpy.exp(r), 1.0)
    #         if numpy.random.uniform() < prob:
    #             thetaNext = numpy.copy(thetaProposal)
    #             likelihood = pt
    #             self._success += 1
    #         else:
    #             thetaNext = theta0
    #             likelihood = pt0
    #         self._theta = thetaNext
    #         self._likelihood = likelihood
    #     self._total += 1
    #     return self._theta, self._likelihood

    def generate(self):
        theta0 = self._theta
        while True:
            thetaProposal = self._proposal.generate(theta0)
            #self._print_vec(theta0, "theta0")
            #self._print_vec(thetaProposal, "thetaProposal"),
            #self._print_vec(self._low, "low")
            #self._print_vec(self._high, "high")
            #raw_input("wait")
            if cmcmc.check_range(len(thetaProposal), thetaProposal, len(self._low), self._low, len(self._high), self._high):
                pt = self._function(thetaProposal)
                jt = self._proposal.logDensity(thetaProposal, theta0)
                pt0 = self._function(theta0)
                jt0 = self._proposal.logDensity(theta0, thetaProposal)
                r = pt - jt - pt0 + jt0
                prob = min(numpy.exp(r), 1.0)
                if numpy.random.uniform() < prob:
                    thetaNext = numpy.copy(thetaProposal)
                    likelihood = pt
                    self._success += 1
                else:
                    thetaNext = theta0
                    likelihood = pt0
                    #if self._success == 0:
                    #    #kept trying until we find successful proposal if this is the first time
                    #    continue
                self._theta = thetaNext
                self._likelihood = likelihood
                self._total += 1
                break
        return self._theta, self._likelihood

    def efficiency(self):
        eff = 0.0
        if self._total > 0:
            eff = float(self._success) / float(self._total)
        return eff
    
    def _print_vec(self, vec, name):
        sep = "----"
        print "DEBUG", sep, name, "start", sep
        for i, val in enumerate(vec):
            print i, val
        print "DEBUG", sep, name, "stop", sep
        return

##############################################################################

class AdaptiveMetropolisHastingsAlgorithm(MetropolisHastingsAlgorithm):
    def __init__(self, *args, **kwargs):
        self._period = 10000
        self._converged = False
        self._save_all = False
        self._data = numpy.zeros(shape=(self._period, len(args[2])))
        super(AdaptiveMetropolisHastingsAlgorithm, self).__init__(*args, **kwargs)

    def generate(self):
        if self._save_all:
            return self._generate_all()
        else:
            return self._generate_only_converged()

    def _generate_only_converged(self):
        while not self._converged:
            self._adaptivegenerate()
        return MetropolisHastingsAlgorithm.generate(self)

    def _generate_all(self):
        if not self._converged:
            return self._adaptivegenerate()
        else:
            #converged, behave as normal
            return MetropolisHastingsAlgorithm.generate(self)

    def _adaptivegenerate(self):
        theta, likelihood = MetropolisHastingsAlgorithm.generate(self)
        self._data[self._total-1] = theta # total - 1 as total is incremented in generate
        #not converged, try to adapt
        if self._total % self._period == 0:
            if self._total > 0:
                self._converged = self._proposal.adapt(self.efficiency(), self._data)
                #reset efficiency counters
                self._success = 0
                self._total = 0
        return theta, likelihood
