import StringIO

import numpy as np

################################################################################

class MonteCarloException(Exception):
    pass

class MonteCarloParameterMismatch(MonteCarloException):
    @classmethod
    def compare_parameters_message(cls, lhs, rhs):
        result = ""
        if set(lhs) == set(rhs) and not (lhs == rhs):
            result = "Wrong order:" + ",\n".join("%s: %s != %s)"%(i, l, r) for i,(l,r) in enumerate(zip(lhs,rhs)) if not l==r)
        else:
            slhs = set(lhs)
            srhs = set(rhs)
            result = "Different parameter names: in left=" + str(slhs - srhs) + "\n in right=" + str(srhs - slhs) 
        return result

################################################################################

class ToyMCExperiment:
    def __init__(self, pars, vec):
        self.pars = pars
        self.vec = vec

    def clone(self):
        return ToyMCExperiment(np.copy(self.pars), np.copy(self.vec))

################################################################################

class ToyMC:
    def __init__(self, ratevector, generator):
        self.ratevector = ratevector
        self.generator = generator
        self._verify(ratevector, generator)
        
    def _verify(self, ratevector, generator):
        #check that the inputs match
        if not ratevector.parameter_names == generator.parameter_names:
            raise MonteCarloParameterMismatch("ToyMC given ratevector and generator with mis-matched names.", MonteCarloParameterMismatch.compare_parameters_message(ratevector.parameter_names, generator.parameter_names))
        return
    
    def __call__(self):
        return self._generate()
    
    def asimov(self):
        '''Return the prediction with all parameters at their starting values and no statistical variations.'''
        pars = np.copy(self.generator.start_values)
        vec = np.copy(self.ratevector(pars))
        return ToyMCExperiment(pars, vec)

    def _generate(self):
        pars = self.generator()
#         if __debug__:
#             logging.debug("ToyMC generating parameters")
#             for n, v in itertools.izip_longest(self.generator.parameter_names, pars):
#                 logging.debug("parameter value %s = %s", n, v)
        vec = self.ratevector(pars)
        return ToyMCExperiment(pars, vec)
    
    def _infostring(self):
        sio = StringIO.StringIO()
        asimov = self.asimov()
        npars = len(asimov.pars)
        nbins = len(asimov.vec)
        total = sum(asimov.vec)
        print >>sio, "ToyMC(npars={}, nbins={}, sumbins={})".format(npars, nbins, total) 
        print >>sio, str(self.generator)
        print >>sio, str(self.ratevector)
        return sio.getvalue() 
    
    def __str__(self):
        return self._infostring()

################################################################################

def generate_timed(self, toymc, maxseconds, maxevents=None, name=None):
        start = time.time()
        end = start + maxseconds
        count = 0
        if name is not None:
            fmt = "%Y-%m-%d %H:%M:%S"
            startstr = datetime.datetime.fromtimestamp(start).strftime(fmt)
            endstr = datetime.datetime.fromtimestamp(end).strftime(fmt)
            print name + " start =", startstr, " expected end =", endstr
        while time.time() < end:
            yield toymc()
            if maxevents is not None and count > maxevents:
                break
        return

################################################################################

def generate_events(self, toymc, n, name=None):
        if name is None:
            iterN = xrange(n)
        else:
            iterN = progress.printprogress(name, n, xrange(n))
        for i in iterN:
            yield toymc()
        return
