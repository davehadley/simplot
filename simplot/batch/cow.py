import os
import tempfile
import time
import socket

from ..tools import hosts

###############################################################################

class CowJobException(Exception):
    pass

###############################################################################

class CowJob:
    
    #a globally accessible time string that all jobs created this session may use to label themselves. 
    initialiseTime = time.strftime("%Y-%m-%d.%H%M")
    
    def __init__(self, jobName, queue, cmd, workingDirectory="./", walltime=2, memoryLimit=1024, setupND280=False, isTest=False, isVerbose=False, isQuiet=False):
        self.cmd = cmd
        self.jobName = jobName
        self.queue = queue
        self.workingDirectory = self.getAbsolutePath(workingDirectory)
        self.doSetupND280 = setupND280
        self.isTest = isTest
        self.isQuiet = isQuiet
        self.isVerbose = isVerbose
        self.scriptPath = None
        self.timeStr = CowJob.initialiseTime
        #set memory limits
        self.memoryLimit = memoryLimit # in MB
        self.swapLimit = max(512,memoryLimit) # in MB
        self.queue = "serial"
        if self.memoryLimit>256:
            self.queue = "taskfarm"
        self.walltime = str(walltime) #walltime in hours
        self.sanitiseInputs()
        
    def __str__(self):
        optStr = ",".join([self.jobName,self.queue])
        return "ClusterJob("+optStr+")"
    
    def getAbsolutePath(self,path):
        return os.path.abspath(os.path.expanduser(path))
    
    def sanitiseInputs(self):
        #ensure that working directory exists
        if not ( os.path.exists(self.workingDirectory) and os.path.isdir(self.workingDirectory) ):
            raise CowJobException("working directory does not exist",self.workingDirectory)
        #ensure that job name is not null
        if not self.jobName:
            raise CowJobException("job name not correctly set: it should be a none null string",self.jobName)
        
    def submit(self):
        self.generateScript()
        cmd = "qsub < "+self.scriptPath
        if self.isTest or not self.isQuiet:
            print str(self),"executing:",cmd
        os.system(cmd)
            
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
        #set batch system options
        self.stdoutFile = os.path.splitext(fileName)[0]+".stdout.out"
        self.stderrFile = os.path.splitext(fileName)[0]+".stderr.out"
        
        options = ["-V", #export the parent shells environment
                   "-N "+self.jobName,
                   "-q "+self.queue,
                   " -o "+self.stdoutFile,
                   " -e "+self.stderrFile,
                   "-l nodes=1:ppn=1,mem={memoryLimit}mb,vmem={swapLimit}mb,walltime={walltime}:00:00".format(memoryLimit=self.memoryLimit, swapLimit=self.swapLimit, walltime=self.walltime)
                   ]
        for opt in options:
            print >>fileHandle,"#PBS",opt
        #setup nd280 software if required
        if self.doSetupND280:
            try:
                cmthome = os.environ["MYCMTHOME"]
                #pianalysisroot = os.environ["PIANALYSISROOT"]
            except KeyError:
                raise CowJobException("error cannot setup ND280, MYCMTHOME environment variable is not set!")
            setupCmd = '''#setting up nd280 software
cd %(cmthome)s && source setup.sh;
''' % { "cmthome":cmthome }
            print >>fileHandle,setupCmd
        #change working directory to current working directory
        if self.workingDirectory is None:
            self.workingDirectory = os.getcwd()
        if not self.workingDirectory == "":
            print >>fileHandle,"cd "+self.workingDirectory
        #finally, add the command to be run
        print >>fileHandle,self.cmd
        #clean up resources and record the script filename
        fileHandle.close()
        #return the absolute path of the file we created
        self.scriptPath = fileName
        return fileName
    
###############################################################################

def _unitTestCowJob():
    cmd = "echo Hello World && env && sleep 5 && echo Done"
    job = CowJob("unitTestCowJob", "express",cmd)
    job.submit()

###############################################################################

def main():
    _unitTestCowJob()
    return    

if __name__ == "__main__":
    main()

