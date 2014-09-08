import os
import time
import stat
import itertools

###############################################################################

class DesktopDefs:
    #see http://www2.warwick.ac.uk/fac/sci/physics/research/epp/internal/computing/hardware/cluster/batch/
    desktops4Core = [
                 "angliru", "gavia", "mortirolo", "pordoi", "zoncolan", "magpie", 
                 "bosberg", "kronplatz", "redroute", "soulor", "stockeu", "tourmalet",
                 "galibier", "kapelmuur", "hautacam", "iseran", "izoard", "rosier", "stelvio",
                 "aspin",
                 "aubisque",
                ]
    
    t2kMachines = ["bosberg", "magpie", "tourmalet", "zoncolan", "grimsel", "redroute", "soulor", ] # "aubisque","simplon" are out of action

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
        #create sub-job scripts
        listOfSubJobScripts = []
        for j in self._jobs:
            scriptPath = j.generateScript()
            logFile = j.stdoutFile
            listOfSubJobScripts.append((scriptPath,logFile))
        #Assign jobs equally between input machines
        jobAssignments = dict([(m,list()) for m in self._machineList])
        nodeCycle = itertools.cycle(jobAssignments.itervalues())
        while len(listOfSubJobScripts)>0:
            subJob = listOfSubJobScripts.pop()
            jobList = nodeCycle.next()
            jobList.append(subJob,)
        #Create command for each machine
        commands = []
        for workerNode,jobList in jobAssignments.iteritems():
            if len(jobList)>0:
                cmd = "ssh "+workerNode+" -f 'screen -d -m  /bin/bash -c \""
                for i,(scriptPath,logFile) in enumerate(jobList):
                    if i>0:
                        cmd += " && "
                    cmd += "nice -n 15 "+scriptPath + " >> "+logFile+ " 2>&1 "
                cmd += " \" ' "
                print "Submitting ",len(jobList),"jobs to",workerNode
                commands.append(cmd)
        #Finally, submit commands
        for cmd in commands:
            print cmd
            os.system(cmd)
        return

###############################################################################

def _unitTestDesktopJob():
    jobs = []
    for i in xrange(10):
        cmd = "echo Test desktop job submission "+str(i)
        j = DesktopJob("job"+str(i), cmd, workingDirectory="./")
        jobs.append(j)
    cmd = "echo Hello World && sleep 5"
    job = DesktopJobRunner(jobs, machineList=DesktopDefs.t2kMachines)
    job.submit()

###############################################################################

def main():
    _unitTestDesktopJob()
    return    

if __name__ == "__main__":
    main()
