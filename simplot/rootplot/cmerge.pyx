import ROOT
import glob
from simplot.progress import printprogress
import numpy as np
cimport numpy as np

def cmerge(outfilename, treename, flist=None):
    flist = [_openrootfile(f) for f in flist]
    tlist = [_gettree(f, treename) for f in flist]
    itree = tlist[0]
    for t in tlist[1:]:
        itree.AddFriend(t)
    ofile = ROOT.TFile(outfilename, "recreate")
    otree = ROOT.TTree(treename, treename)
    #setup input and output branches
    cdef list outputcontainers = []
    branches = [br for br in itree.GetListOfBranches()]
    branches += [br for fr in itree.GetListOfFriends() for br in fr.GetTree().GetListOfBranches()]
    for br in branches:
        arr = np.zeros(shape=(1,), dtype=float)
        itree.SetBranchAddress(br.GetName(), arr)
        otree.Branch(br.GetName(), arr, br.GetName() + "/D")
        outputcontainers.append(arr)
    #loop over entries
    for entrynum in printprogress("merge", itree.GetEntries(), xrange(itree.GetEntries())):
        itree.GetEntry(entrynum)
        otree.Fill()
    otree.Write()
    ofile.Close()
    return

#merge("test.root", "titus_weights", glob.glob("weights_*.root"))

def _openrootfile(fname):
    # check we were successful
    f = ROOT.TFile(fname, "READ")
    if not (f and f.IsOpen()):
        raise ValueError("unable to open file", fname)
    return f

def _gettree(tfile, treename):
    tree = tfile.Get(treename)
    if not tree:
        raise ValueError("unable to find tree in file", treename, tfile.GetName())
    return tree
