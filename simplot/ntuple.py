import itertools
import ROOT
import numpy
import operator
import os
import glob

###############################################################################

class BranchPrimitive:
    def __init__(self, name, tree, start=None):
        """
        :param name: name of the output variable.
        :type name: str
        :param tree: the output TTree.
        """
        self.name = name
        if start is None:
            start = 0.0
        self.value = numpy.array([start], dtype=float)
        self.start = start
        self.tree = tree
        self._setbranch()
        
    def _setbranch(self):
        branch = self.tree.FindBranch(self.name)
        if not branch:
            branch = self.tree.Branch(self.name, self.value, self.name+"/D" )
        else:
            #self.tree.SetBranchAddress(self.name, self.value, self.name+"/D" )
            self.tree.SetBranchAddress(self.name, self.value)
        self.branch = branch
    
    def _clearbranch(self):
        if self.branch:
            self.branch.SetAddress(0)
    
    def setvalue(self, value):
        """Set the value. This should be done prior to calling TTree.Fill().
        
        :param value: value to be set for current entry in TTree.
        :type value: float 
        """
        self.value[0] = value
    
    def getvalue(self):
        """
        :returns: current value.
        """
        return self.value[0]
    
    def reset(self):
        self.setvalue(self.start)
    
    def __str__(self):
        return "(",self.name,self.value[0],")"
    
###############################################################################
    
class BranchObject:
    """Writes C++ objects to a ROOT TTree. 
    
    The object must have a default constructor that takes no arguments.
    The dictionaries for this object must be loaded.
    """
    def __init__(self,ObjectType, name, tree, start=None):
        """
        :param ObjectType: the class of this type eg ROOT.TLorentzVector
        :type ObjectType: class
        :param name: name of the output branch.
        :type name: str
        :param tree: the output TTree.
        :type tree: ROOT.TTree
        """
        self.name = name
        self.ObjectType = ObjectType
        if start is None:
            start = ObjectType() 
        self.value = start
        self.start = start
        self.tree = tree
        self._setbranch()
        
    def _setbranch(self):
        self.branch = self.tree.FindBranch(self.name)
        if not self.branch:
            self.branch = self.tree.Branch(self.name, self.value)
        else:
            #tree.SetBranchAddress(self.name, self.value)
            self.branch.SetAddress(ROOT.AddressOf(self.value))
    
    def _clearbranch(self):
        if self.branch:
            self.branch.SetAddress(0)

    def setvalue(self,value):
        '''Set reference to the output object.
        
        Takes reference of the object to be writen. 
        Preferred method provided value is guaranteed to exist at time of ROOT.TTree.Fill.
        
        :param value: reference to the output object to be filled.
        '''
        self.value = value
        #update branch address
        self.branch.SetAddress(ROOT.AddressOf(self.value))
    
    def copyvalue(self,value):
        '''Tries to copy the value. 
        
        Use in case value is not guaranteed to exist at the time of ROOT.TTree.Fill.
        This is potentially a slow method so set should be used if possible.
        
        :param value: reference to the output object to be copied then filled.
        '''
        self.value = self.ObjectType(value)
        #update branch address
        self.branch.SetAddress(ROOT.AddressOf(self.value))
        
    def reset(self):
        self.setvalue(self.start)
    
    def getvalue(self):
        """
        :returns: reference to the current object.
        """
        return self.value
    
    def __str__(self):
        return "BranchObject("+self.name+")"

###############################################################################

class TreeCopier:
    
    def __init__(self,treeName,inputDataset,outputFileName):
        self.treeName = treeName
        self.inputDataset = inputDataset
        self.outputFileName = outputFileName
        #create input chain
        self.inputChain = ROOT.TChain(self.treeName,self.treeName)
        for fileName in self.inputDataset:
            self.inputChain.AddFile(fileName)
        #create output file
        self.outputFile = ROOT.TFile(self.outputFileName,"recreate")
        assert self.outputFile.IsOpen()
        self.outputTree = self.inputChain.CloneTree(0)
        self.outputFile.cd()
        
    def __iter__(self):
        for event in self.inputChain:
            yield event
            
    def getOutputTree(self):
        return self.outputTree

    def fill(self):
        self.outputTree.Fill()
    
    def close(self):
        self.outputTree.AutoSave()
        self.outputFile.Close()
        
    def numEntries(self):
        return self.inputChain.GetEntries()

###############################################################################

class Algorithm(object):
    '''Define the interface for algorithms applied to ROOT tree by the ProcessTree runnable.'''
    
    def begin(self):
        pass
    
    def file(self, tfile):
        pass
    
    def event(self, event):
        pass
    
    def end(self, end):
        pass

###############################################################################

class AlgorithmList(Algorithm):
    
    def __init__(self, algs):
        self._algs = algs
    
    def begin(self):
        for a in self._algs:
            a.begin()
    
    def file(self, tfile):
        for a in self._algs:
            a.file(tfile)
    
    def event(self, event):
        for a in self._algs:
            a.event(event)
    
    def end(self):
        for a in self._algs:
            a.end()

###############################################################################

class TreeFillerAlgorithm(Algorithm):
    
    def __init__(self, outfilename, treename, branches):
        self._outfilename = outfilename
        self._outtreename = treename
        self._branches = branches
        self._outfile = None
        self._outtree = None

    def begin(self):
        self._outfile = ROOT.TFile(self._outfilename, "RECREATE")
        self._outtree = ROOT.TTree(self._outtreename, self._outtreename)
        for b in self._branches:
            b.createbranch(self._outtree)
        return

    def file(self, tfile):
        return

    def event(self, event):
        for br in self._branches:
            br.event(event)
        self._outtree.Fill()
        return

    def end(self):
        self._outtree.AutoSave()
        self._outfile.Close()
        return

###############################################################################

class FileObjectGetter:
    def __init__(self, name):
        self._name = name
    def __call__(self, tfile):
        return tfile.Get(self._name)

###############################################################################

class ProcessTree:
    def __init__(self, infilelist, treename, alg):
        '''Iterates over the input tree and applies the provided algorithm.
        '''
        self._filelist = _expand_file_patterns(infilelist)
        self._alg = alg
        self._tree_getter = FileObjectGetter(treename)
        return
    
    def run(self):
        self._alg.begin()
        for tfile in _iter_root_files(self._filelist):
            self._alg.file(tfile)
            tree = self._tree_getter(tfile)
            for event in tree:
                self._alg.event(event)
        self._alg.end()
        return

###############################################################################

class BranchFiller:
    def __init__(self, name, function, start_value=0.0, ):
        '''BranchFiller handles setting branch values on each event and is designed to be provided to a TreeFillerAlgorithm.
        The required arguments are the name of the output branch and a callable object 
        or function that is applied to the event that returns a double.
        Alternatively, function if is a string, a simple attribute getter is constructed with it. 
        Optionally the start_value can be set.
        '''
        self._name = name
        self._start_value = start_value
        self._branch = None
        if isinstance(function, str):
            function = operator.attrgetter(function)
        self._function = function
        if not callable(self._function):
            raise Exception("BranchFiller given a non-callable function object.")
    
    def createbranch(self, tree):
        self._branch = BranchPrimitive(name=self._name, tree=tree, start=self._start_value)
        return
    
    def event(self, event):
        self._branch.setvalue(self._function(event))

###############################################################################

def _find_file_matches(pattern):
    pattern = os.path.expanduser(os.path.expandvars(pattern))
    matches = glob.glob(pattern)
    return matches

###############################################################################

def _expand_file_patterns(pattern_list):
    result = []
    for p in pattern_list:
        result.extend(_find_file_matches(p))
    return result

###############################################################################

def _iter_root_files(filelist):
    for fname in filelist:
        yield ROOT.TFile(fname)

###############################################################################

_testfilename = "unitTestBranchObj.root"
_testtreename = "unitTestBranchObj"

def _unitTestBranchObjWrite():
    print "--- _unitTestBranchObjWrite"
    #generate dictionaries
    ROOT.gInterpreter.GenerateDictionary("std::vector<TLorentzVector>","TLorentzVector.h")
    #create output tree
    tfile = ROOT.TFile(_testfilename,"recreate")
    tree = ROOT.TTree("unitTestBranchObj","unitTestBranchObj")
    #define branches
    x = BranchPrimitive("x", tree)
    y = BranchPrimitive("y", tree)
    z = BranchPrimitive("z", tree)
    testLorentzVec = BranchObject(ROOT.TLorentzVector,"testLorentzVec", tree)
    testStdVec = BranchObject(ROOT.std.vector("double"),"testStdVec", tree)
    #fill tree with some random data
    StdVecHlv = ROOT.std.vector("TLorentzVector")
    testStdVecLorentz = BranchObject(StdVecHlv,"testStdVecLorentz", tree)
    rand = ROOT.TRandom3()
    print "unitTestBranchObj begining event loop"
    for i in xrange(1000):
        x.setvalue(rand.Gaus(1,1))
        y.setvalue(rand.Gaus(10,1))
        z.setvalue(rand.Gaus(20,1))
        e = rand.Gaus(1000,100)
        hlv = ROOT.TLorentzVector(1.0,2.0,3.0,4.0)
        testLorentzVec.setvalue(hlv)
        vec = ROOT.std.vector("double")()
        for i in xrange(10):
            vec.push_back(rand.Gaus(-1,1))
        testStdVec.setvalue(vec)
        vecHlv = StdVecHlv()
        hlv1 = ROOT.TLorentzVector(1.0,2.0,3.0,4.0)
        hlv2 = ROOT.TLorentzVector(5.0,6.0,7.0,8.0)
        vecHlv.push_back(hlv1)
        vecHlv.push_back(hlv2)
        vecHlv.push_back(ROOT.TLorentzVector(9.0,10.0,11.0,12.0))
        testStdVecLorentz.setvalue(vecHlv)
        tree.Fill()
    tree.Write()
    tfile.Close()
    print "--- _unitTestBranchObjWrite... done." 
    return

##################################################

def _unitTestBranchObjRead():
    print "--- _unitTestBranchObjRead"
    tfile = ROOT.TFile(_testfilename,"read")
    tree = tfile.Get(_testtreename)
    print tree
    for i,event in enumerate(tree):
        print "TEST double = ",event.x
        event.testLorentzVec.Print()
        print "TEST testStdVec = ",event.testStdVec.size(),[ x for x in event.testStdVec ]
        print "TEST vec Lorentz Vec:",event.testStdVecLorentz.size(),"entries"
        for hlv in event.testStdVecLorentz:
            hlv.Print()
        if i>=2:
            break
    print "--- _unitTestBranchObjRead ..done"
    return

##################################################

def _unitTestAlgorithm():
    #create and output file
    _unitTestBranchObjWrite()
    branches = [BranchFiller("x", "x", start_value=0.0),
                BranchFiller("x2", lambda x: x.x*2, start_value=0.0),
                ]
    treefiller = TreeFillerAlgorithm(outfilename="test_alg_ntuple.root", treename="test_alg", branches=branches)
    process = ProcessTree(infilelist=[_testfilename], treename=_testtreename, alg=treefiller)
    process.run()
    return

##################################################

def main():
    """Run some unit tests that check reading and writing of a ttree still works."""
    #_unitTestBranchObjWrite()
    #_unitTestBranchObjRead()
    _unitTestAlgorithm()
    return

if __name__ == "__main__":
    main()
