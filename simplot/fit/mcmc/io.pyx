from simplot.mc.montecarlo import generate_timed, generate_events

import numpy as np
cimport numpy as np
import ROOT

def write_mc_to_root_file(mc, filename, treename, nevents, seconds=None, burnin=None):
    out = _RootOutput(mc, filename, treename)
    if burnin is None:
        burnin  = 0
    out.write(nevents, seconds=seconds, burnin=burnin)
    return

cdef class _RootOutput:
    cdef object _mc
    cdef str _filename
    cdef str _treename
    cdef list _arrays
    cdef np.ndarray _likelihood

    def __init__(self, mc, filename, treename="mc"):
        self._mc = mc
        self._arrays = []
        self._filename = filename
        self._treename = treename

    def _setup_output_file(self):
        tfile = ROOT.TFile(self._filename, "RECREATE")
        tree = ROOT.TTree(self._treename, self._treename)
        self._setup_branches(tree)
        return tree, tfile

    def _setup_branches(self, tree):
        arr = np.zeros(dtype=float, shape=(1,))
        tree.Branch("likelihood", arr, "likelihood/D")
        self._likelihood = arr
        for parname in self._mc.parameter_names:
            arr = np.zeros(dtype=float, shape=(1,))
            tree.Branch(parname, arr, parname + "/D")
            self._arrays.append(arr)
        return

    def _iterable(self, nevents, seconds):
        mc = self._mc
        name = "writing mc to " + self._filename
        if seconds is None:
            iterable = generate_events(mc, nevents, name=name)
        else:
            iterable = generate_timed(mc, maxseconds=seconds, name=name)
        return iterable

    def _set_values(self, np.ndarray[dtype=double, ndim=1] parameters, double likelihood):
        for ii in xrange(parameters.shape[0]):
            self._arrays[ii][0] = parameters[ii]
        self._likelihood[0] = likelihood
        return

    def write(self, nevents, seconds=None, int burnin=0):
        tree, tfile = self._setup_output_file()
        iterable = self._iterable(nevents, seconds)
        cdef int count = 0
        for parameters, likelihood in iterable:
            if count >= burnin:
                self._set_values(parameters, likelihood)
                tree.Fill()
            count += 1
        tree.Write()
        tfile.Close()
        return
