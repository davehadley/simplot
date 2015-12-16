import lib
import ROOT

lib.load()
crootprob3pp = ROOT.crootprob3pp
Probability = crootprob3pp.Probability

class Flavour:
    NU_E = 1;
    NU_MU = 2;
    NU_TAU = 3;

class CP:
    MATTER = 1;
    ANTI_MATTER = -1;
