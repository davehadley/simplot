from argparse import ArgumentParser
import glob

import ROOT
from simplot.rootplot.ntuple import BranchPrimitive
from simplot.progress import printprogress

from simplot.rootplot.cmerge import cmerge

def merge(outfilename, treename, flist=None):
    flist = [ROOT.TFile(f) for f in flist]
    tlist = [f.Get(treename) for f in flist]
    itree = tlist[0]
    for t in tlist[1:]:
        itree.AddFriend(t)
    ofile = ROOT.TFile(outfilename, "recreate")
    otree = ROOT.TTree(treename, treename)
    branches = [BranchPrimitive(br.GetName(), otree) for br in itree.GetListOfBranches()]
    branches += [BranchPrimitive(br.GetName(), otree) for fr in itree.GetListOfFriends() for br in fr.GetTree().GetListOfBranches()]
    for entry in printprogress("merge", itree.GetEntries(), itree):
        for br in branches:
            br.setvalue(getattr(entry, br.name))
        otree.Fill()
    otree.Write()
    ofile.Close()
    return

def parsecml():
    parser = ArgumentParser()
    parser.add_argument("-o", "--output-file", type=str, help="Output filename", default="merged.root")
    parser.add_argument("-t", "--tree", type=str, help="Input/output treename", default="titus_weights")
    parser.add_argument("input_files", nargs="+", help="Input file list")
    return parser.parse_args()

def main():
    args = parsecml()
    infilelist = []
    for pat in args.input_files:
        match = glob.glob(pat)
        if len(match) == 0:
            raise Exception("ERROR: no matching files, check command line arguments", pat)
        infilelist.extend(match)
    cmerge(args.output_file, args.tree, infilelist)
    return

if __name__ == "__main__":
    main()
