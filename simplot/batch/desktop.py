import os
import time
import stat
import itertools
import subprocess
import threading
from time import strftime
from Queue import Queue, Empty

###############################################################################

def _addwarwickdomain(l):
    return [n + ".epp.warwick.ac.uk" for n in l]

class DesktopDefs:
    #see http://www2.warwick.ac.uk/fac/sci/physics/research/epp/internal/computing/hardware/cluster/batch/
    desktops4Core = [
                 "angliru", "gavia", "mortirolo", "pordoi", "zoncolan", "magpie", 
                 "bosberg", "kronplatz", "redroute", "soulor", "stockeu", "tourmalet",
                 "galibier", "kapelmuur", "hautacam", "iseran", "izoard", "rosier", "stelvio",
                 "aspin",
                 "aubisque",
                ]
    
    t2kMachines = ["bosberg", "magpie", "tourmalet", "zoncolan", "oropa", "redroute", "soulor", "aubisque"] # "aubisque","simplon" are out of action

    #t2kMachines = _addwarwickdomain(t2kMachines)
    #desktops4Core = _addwarwickdomain(desktops4Core)

###############################################################################

class DesktopJobException(Exception):
    pass

###############################################################################

class DesktopJob:
    
    #a globally accessible time string that all jobs created this session may use to label themselves. 
    initialiseTime = time.strftime("%Y-%m-%d.%H%M")
    
    def __init__(self, jobName, cmd, workingDirectory="./", memLimit=1024, vmemLimit=1024, isTest=False, isVerbose=False, isQuiet=False):
        self.cmd = cmd
        self.jobName = jobName
        self.workingDirectory = self.getAbsolutePath(workingDirectory)
        self.isTest = isTest
        self.isQuiet = isQuiet
        self.isVerbose = isVerbose
        self.scriptPath = None
        self.timeStr = DesktopJob.initialiseTime
        #set memory limits
        self.memoryLimit = memLimit
        self.swapLimit = vmemLimit
        self.sanitiseInputs()
        
    def __str__(self):
        optStr = ",".join([self.jobName,self.queue])
        return "DesktopJob("+optStr+")"
    
    def getAbsolutePath(self,path):
        return os.path.abspath(os.path.expanduser(path))
    
    def sanitiseInputs(self):
        #ensure that working directory exists
        if not ( os.path.exists(self.workingDirectory) and os.path.isdir(self.workingDirectory) ):
            raise DesktopJobException("working directory does not exist",self.workingDirectory)
        #ensure that job name is not null
        if not self.jobName:
            raise DesktopJobException("job name not correctly set: it should be a none null string",self.jobName)
        
    def generateScript(self):
        #create a unique file to write the script to.
        #fileHandle = tempfile.NamedTemporaryFile(prefix=self.jobName, suffix=".sh", dir=self.workingDirectory, delete=False)
        #fileName = fileHandle.name
        fileHandle = None
        counter = 0
        while not fileHandle:
            fileName = "%(pwd)s/%(name)s.%(time)s.%(count)03d.sh" % {"pwd":self.workingDirectory,
                                                                  "name":self.jobName,
                                                                  "time":self.timeStr,
                                                                  "count":counter,
                                                                  }
            counter += 1
            if not os.path.exists(fileName):
                fileHandle = open(fileName,"w")
        #make sure we are using bash
        print >>fileHandle,"#!/bin/bash"
        self.stdoutFile = os.path.splitext(fileName)[0]+".out"
        #copy environment into script
        for k,v in os.environ.iteritems():
            if "BASH_FUNC" in k:
                #environment variables containing BASH_FUNC cause problems on CSC.
                continue
            print >>fileHandle,"export {k}=\"{v}\"".format(k=k,v=v)
        print >>fileHandle,"ulimit -m "+str(self.memoryLimit*1000)
        #change working directory to current working directory
        if self.workingDirectory is None:
            self.workingDirectory = os.getcwd()
        if not self.workingDirectory == "":
            print >>fileHandle,"cd "+self.workingDirectory
        #finally, add the command to be run
        print >>fileHandle,self.cmd
        #clean up resources and record the script filename
        fileHandle.close()
        #make file executable
        os.chmod(fileName, stat.S_IXUSR | stat.S_IRUSR | stat.S_IWUSR )
        #return the absolute path of the file we created
        self.scriptPath = fileName
        return fileName
        
###############################################################################

class DesktopJobRunner:
    
    def __init__(self, jobs, machineList=DesktopDefs.t2kMachines, beNice=True):
        '''Run several shell commands in on multiple desktop machines. 
        '''
        self._jobs = jobs
        self._beNice = beNice
        self._machineList = machineList

    def submit(self):
        #Fill jobs into thread safe queue.
        jobqueue = Queue()
        for subjobnum, j in enumerate(self._jobs):
            script = j.generateScript()
            logfile = j.stdoutFile
            jobqueue.put((subjobnum, script, logfile))
        #Create _NodeJobRunner for each node.
        nodes = [_NodeJobRunner(m, jobqueue) for m in self._machineList]
        #Run each NodeJobRunner until the Queue is empty.
        threads = []
        for n in nodes:
            t = threading.Thread(target=n)
            threads.append(t)
            t.start()
        #Wait until all threads are finished
        for t in threads:
            t.join()
        return        

class _NodeJobRunner:
    def __init__(self, node, jobqueue):
        self.node = node
        self.joblist = jobqueue
        self._totaljobs = jobqueue.qsize()

    def __call__(self):
        return self.run()

    def run(self):
        node = self.node
        joblist = self.joblist
        njobs = self._totaljobs
        try:
            while True:
                subjobnum, scriptpath, logfile = joblist.get(block=False)
                cmd = "nice -n 15 "+ scriptpath + " >> " + logfile + " 2>&1 "
                cmd = "ssh " + node + " -C \"" + cmd + "\""
                print "DesktopJobRunner(%s, %s,  starting, %s / %s, %s)" % (strftime("%Y-%m-%d %H:%M:%S"), node, subjobnum, njobs, logfile)
                try:
                    subprocess.check_call(cmd, shell=True)
                except Exception as ex:
                    print "ERROR ", node, "jobs failed."
                    print ex
                else:
                    print "DesktopJobRunner(%s, %s, completed, %s / %s, %s)" % (strftime("%Y-%m-%d %H:%M:%S"), node, subjobnum, njobs, logfile)
        except Empty:
            # queue is empty, nothing left to do
            pass
        return

###############################################################################

def _unitTestDesktopJob():
    jobs = []
    machines = DesktopDefs.t2kMachines
    for i in xrange(3*len(machines)):
        cmd = "echo Test desktop job submission "+str(i)
        cmd += "\nsleep 5;\necho Job "+str(i)+" done."
        j = DesktopJob("job"+str(i), cmd, workingDirectory="./")
        jobs.append(j)
    job = DesktopJobRunner(jobs, machineList=machines)
    job.submit()

###############################################################################

def main():
    _unitTestDesktopJob()
    return    

if __name__ == "__main__":
    main()
