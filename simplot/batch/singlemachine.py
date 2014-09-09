import subprocess
import time
import re
import multiprocessing

###############################################################################

class Job:
    def __init__(self, name, cmd, logFileName="", memLimit=1024, vmemLimit=4*1024):
        '''Define a single job.
        
        Set logFileName to None to switch off log files.
        If logFileName is an empty string it will be automatically assigned to
        "log."+name+".txt".
        
        :param name: name of job
        :type name: str
        :param command: shell command to run job
        :type command: str
        :param logFileName: log file name.
        :type logFileName: str
        :param memLimit: physical memory limit in units of mega-bytes.
        :type memLimit: int
        :param vmemLimit: virtual memory limit in units of mega-bytes.
        :type vmemLimit: int
        '''
        self.name = name
        self.command = cmd
        if logFileName is not None and logFileName == "":
            logFileName = "log.{0}.txt".format(name)
        self.logFileName = logFileName
        self.memLimit = memLimit
        self.vmemLimit = vmemLimit

###############################################################################

class ParallelJobRunner:
    def __init__(self, jobs, ncores=None, beNice=True):
        '''Run several shell commands in parallel on a single machine. 
        
        If provided, the stdout and stderr of the jobs will be written to the 
        job log file name.
        If nCores is None then it is automatically runs on all cores on this machine.
        
        :param jobs: list of jobs.
        :type commands: list of :class:`Job`
        :param nCores: maximum number of simultaneous jobs to run.
        :type nCores: int
        '''
        self._jobs = jobs
        self._beNice = beNice
        maxCores = multiprocessing.cpu_count()
        if ncores is None:
            #Automatically run all cores
            ncores = max(1,maxCores)
        self._nCores = ncores
        #Check inputs are valid
        if ncores <= 0:
            raise ValueError("number of cores must be greater than 0. ({0} given)".format(ncores))
        if  ncores > maxCores:
            raise ValueError("number of cores must be less than or equal to the number of cores available on this machine ({0} given, {1} available on this host).".format(ncores,maxCores))

    def run(self):
        self._runCommandsInParallel()
        return
    

    def _runCommandsInParallel(self):
        beNice = self._beNice
        msgPrefix = "ParallelJobRunner :"
        jobs = list(self._jobs)
        print msgPrefix,len(jobs),"jobs to run."
        jobs = jobs[::-1] # reverse the list
        allNames = [j.name for j in jobs]
        allJobs = [j for j in jobs]
        commands = [j.command for j in jobs]
        nCores = self._nCores
        logFileNames = [j.logFileName for j in jobs]
        processes = [None]*nCores
        names = [None]*nCores
        waitTime = 0
        timeSinceUpdate = 0
        updateMilestone = 0
        while True:
            waitFlag = True
            #Start jobs
            for i in xrange(nCores):
                #Check if job is complete
                if processes[i] and processes[i].poll() is not None:
                    print "Job complete:",names[i]
                    processes[i] = None
                    names[i] = None
                #Start new job if there is a free core
                if processes[i] is None and len(commands)>0:
                    names[i] = allNames.pop()
                    job = allJobs.pop()
                    cmd = commands.pop()
                    if beNice and "nice -n" not in cmd:
                        cmd = "nice -n 15 " + cmd
                    #Add memory limits
                    bytes = 1000
                    cmd = "ulimit -m {mem} && ulimit -v {vmem} && {cmd}".format(mem=job.memLimit*bytes,vmem=job.vmemLimit*bytes,cmd=cmd)
                    #Create log file
                    outFile = None
                    if logFileNames:
                        fname = logFileNames.pop()
                        outFile = open(fname,"w")
                    print >>outFile,"Running command:"
                    print >>outFile,cmd
                    processes[i] = subprocess.Popen(cmd,shell=True, stdout=outFile, stderr=outFile)
                    waitFlag = False
            #Check whether all jobs have finished
            if all( (p is None for p in processes)  ):
                break
            #Wait before trying again
            if waitFlag:
                if timeSinceUpdate==0:
                    running = [n for n in names if n is not None]
                    print msgPrefix,len(running),"running jobs ",running,",",len(commands),"queued, (",time.asctime(),")"
                    if waitTime<10:
                        #slowly increase wait time as the job goes on
                        waitTime += 1
                    if updateMilestone < 300:
                        #slowly increase period of updates to prevent too many messages on long jobs
                        # stop at 5*60=300s between messages.
                        updateMilestone += 10
                timeSinceUpdate += waitTime
                if timeSinceUpdate >= updateMilestone:
                    timeSinceUpdate = 0
                time.sleep(waitTime)
        print msgPrefix,"complete."
        return

###############################################################################

def _example():
    commands = ["echo Hello World",
                "sleep 1",
                "sleep 2",
                "sleep 3",
                "sleep 4",
                "sleep 5",
                "sleep 6",
                "sleep 7",
                "sleep 8",
                "sleep 9",
                "sleep 10",
                ]
    #commands = ["sleep 600"]*3
    jobs = [Job(i,c) for i,c in enumerate(commands)]
    cores = None
    job = ParallelJobRunner(jobs,nCores=cores)
    job.run()
    return

###############################################################################

def _memoryLeakExample():
    njobs = 4
    commands = ['python -c "import time; l=[range(10**6) for _ in xrange(10)]; print range(10); time.sleep(100);" 2>&1 '] * njobs
    jobs = [Job(i,c,logFileName="log.memoryLeakTest.txt", memLimit=128, vmemLimit=128) for i, c in enumerate(commands)]
    cores = None
    job = ParallelJobRunner(jobs,nCores=cores)
    job.run()
    return

###############################################################################

def main():
    _example()
    #_memoryLeakExample()

###############################################################################

if __name__ == "__main__":
    main()