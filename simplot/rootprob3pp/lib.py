import ROOT  # @UnresolvedImport
import os
import pkg_resources

_hasloaded = False

def load():
    global _hasloaded
    if (not _hasloaded) and (not hasattr(ROOT, "crootprob3pp")):
        libname = "libProb3ppShared"
        if not _checklibexists(libname):
            raise Exception("Cannot find library in LD_LIBRARY_PATH. Did you run the install script?", libname)
        ROOT.gSystem.Load(libname)
        #compile prob3++
        probdir = pkg_resources.resource_filename("simplot.rootprob3pp", "Prob3++.20121225")
        ROOT.gSystem.SetIncludePath(ROOT.gSystem.GetIncludePath()+ " -I{} ".format(probdir))
        for fname in [#"mosc.c", "mosc3.c",
                      "EarthDensity.cc", "NeutrinoPropagator.h", "BargerPropagator.cc"]:
            fullfname = os.sep.join([probdir, fname])
            ROOT.gROOT.ProcessLine(".L {}+".format(fullfname))
        fname = pkg_resources.resource_filename("simplot.rootprob3pp", "Prob3ppProbability.cxx")  # @UndefinedVariable
        ROOT.gROOT.ProcessLine(".L {}+".format(fname))
        _hasloaded = True
    return

def _checklibexists(libname):
    ldpath = os.environ["LD_LIBRARY_PATH"]
    for dirname in ldpath.split(":"):
        for ext in [".so"]:
            fname = os.sep.join([dirname, libname + ext])
            if os.path.exists(fname):
                return True
    return False

def _testprob():
    class Flavour:
        NU_E = 1
        NU_MU = 2
        NU_TAU = 3
    prob = ROOT.crootprob3pp.Probability()
    for i in xrange(200, 1000, 100):
        e = float(i) / 1000.0
        print "numu -> nue, e={}, P={}".format(e, prob.prob(Flavour.NU_MU, Flavour.NU_E, e))
    for i in xrange(200, 1000, 100):
        e = float(i) / 1000.0
        print "numu -> numu, e={}, P={}".format(e, prob.prob(Flavour.NU_MU, Flavour.NU_MU, e))
    for i in xrange(200, 1000, 100):
        e = float(i) / 1000.0
        print "numu -> nutau, e={}, P={}".format(e, prob.prob(Flavour.NU_MU, Flavour.NU_TAU, e))
    return

def main():
    _testprob()
    print "done."
    return

if __name__ == "__main__":
    main()
