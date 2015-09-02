import unittest
import tempfile

import ROOT
import numpy as np

from simplot.rootplot.merge import merge, cmerge

class TestMerge(unittest.TestCase):
    def setUp(self):
        #create some temp files
        self._treename = "TestTree"
        self._branchnames = []
        N = 2
        tempfiles = [tempfile.NamedTemporaryFile(suffix=".root") for _ in xrange(N)]
        for ii, f in enumerate(tempfiles):
            f = ROOT.TFile(f.name, "recreate")
            tree = ROOT.TTree(self._treename, self._treename)
            name = chr(ord("a") + ii)
            arr = np.zeros(shape=(1,), dtype=float)
            tree.Branch(name, arr, name + "/D")
            for jj in xrange(1000):
                arr[0] = ii + jj
                tree.Fill()
            tree.Write()
            f.Close()
        # hold on until the end
        self._tempfiles = tempfiles
        self._tempfilenames = [t.name for t in tempfiles]
        return


    def test_merge(self):
        outfile = tempfile.NamedTemporaryFile(suffix=".root")
        merge(outfile.name, self._treename, self._tempfilenames)
        self._checkmerged(outfile.name)
        return

    def test_cmerge(self):
        outfile = tempfile.NamedTemporaryFile(suffix=".root")
        cmerge(outfile.name, self._treename, self._tempfilenames)
        self._checkmerged(outfile.name)
        return

    def _checkmerged(self, fname):
        tfile = ROOT.TFile(fname)
        tree = tfile.Get(self._treename)
        for ii, entry in enumerate(tree):
            for jj, br in enumerate(self._branchnames):
                x = getattr(entry, br)
                self.assertEquals(ii + jj, x)
        return

def main():
    unittest.main()

if __name__ == "__main__":
    main()
