import os
import simplot.filelock as filelock

import ROOT

_hasloadedlib = False

def makeproject(filename, name="", oaanalysis=False):
    global _hasloadedlib
    global _libname
    #do sanity checks
    namesMatch = (not _hasloadedlib) or _libname == name
    if not namesMatch:
        raise Exception("already loaded an oaAnalaysis library! Loading multiple libraries within the same session is not supported", name, _libname)
    #open input file
    inputFileName = filename
    ROOT.TFile.SetReadStreamerInfo(False) #suppress missing dictionary warnings by not reading in streamer info.
    tfile = ROOT.TFile(inputFileName,"read")
    ROOT.TFile.SetReadStreamerInfo(True)
    if oaanalysis:
        #Automatically determine the software version from header
        basicHeader = tfile.Get("HeaderDir/BasicHeader")
        if basicHeader and basicHeader.GetEntries() > 0:
            basicHeader.GetEntry(0)
            name = basicHeader.SoftwareVersion
            #remove null characters from SoftwareVersion
            name = "".join(name.split("\x00"))
            if len(name) == 0:
                name = "unknownversion"
            if int(str(basicHeader.IsMC)): # string conversion neccessary as in some versions IsMC is of char type
                name += "_mc"
            else:
                name += "_data"
    outputPath = "libMakeProject"
    if len(name)>0:
        outputPath += "_"+name
    libraryFileName = outputPath+"/"+outputPath+".so"
#    libraryFileName = "./"+libraryFileName
    workingDir = os.getcwd()
    libraryFileName = workingDir+"/"+libraryFileName
    # get lock to prevent multiple jobs simultaneously trying to make library
    lock = filelock.FileLock(".lockFile_libReadoaAnalysis_"+name)    
    lock.lock()
    if not os.path.exists(libraryFileName):
        if not tfile.IsOpen():
            raise Exception("loadLibReadOaAnalysis could not open input file",inputFileName, name)
        tfile.ReadStreamerInfo() #manually read streamer info since it was disabled above
        tfile.MakeProject(outputPath,"*","update+")
    ROOT.gSystem.Load("libCint.so")
    ROOT.gSystem.Load(libraryFileName)
#    ROOT.gSystem.Load("libReadoaAnalysis/libReadoaAnalysis.so")
    _hasloadedlib = True
    _libname = name
    #load hack scripts
    lock.release()
    return outputPath

def makeproject_new_process(filename, name="", oaanalysis=False):
    '''As makeproject above but does the compilation in a separate process.'''
    import subprocess
    cmd = ["python", "-m", "simplot.rootplot.makeproject", filename, "--name", name]
    if oaanalysis:
        cmd.append("--oaanalysis")
    subprocess.check_call(cmd)
    return makeproject(filename, name, oaanalysis)

def _parsecml():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("fname", help="Input file name")
    parser.add_argument("--name", help="Name of library file to generate.", type=str, default="")
    parser.add_argument("--oaanalysis", help="oaAnalysis mode.", action="store_true")
    return parser.parse_args() 

def main():
    args = _parsecml()
    makeproject(args.fname, args.name, args.oaanalysis)
    return

if __name__ == "__main__":
    main()
