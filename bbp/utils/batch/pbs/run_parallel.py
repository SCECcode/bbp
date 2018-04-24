#!/usr/bin/env python
"""
Copyright 2010-2017 University Of Southern California

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""
import sys
import os
import time
import socket
import subprocess

class RunParallel:

    def __init__(self, envscript):
        self.envscript = envscript
        return


    def runMultiSSH(self, nodefile, numcores, cmdlist):
        hostname = socket.gethostname()

        # Read in nodefile
        nodelist = []
        if (nodefile == 'localhost'):
            for i in xrange(0, numcores):
                nodelist.append(hostname)
        else:
            ip = open(nodefile)
            lines = ip.read().splitlines()
            ip.close()
            for line in lines:
                if ((line != hostname) and (line != '')):
                    for i in xrange(0, numcores):
                        nodelist.append(line)
                elif (line == hostname):
                    for i in xrange(0, numcores):
                        nodelist.append(hostname)

        if (len(nodelist) == 0):
            print "No compute nodes available"
            return(1)
        else:
            print "Running on %s cores" % (len(nodelist))

        # Execute each station's xml file
        proclist = []
        for c in cmdlist:
            while (len(nodelist) == 0):
                # Free up some nodes
                newproc = []
                for proc in proclist:
                    if (proc[0].poll() != None):
                        if (proc[0].wait() != 0):
                            print "Process on node %s failed" % (proc[1])
                            return(1)
                        else:
                            nodelist.append(proc[1])
                    else:
                        newproc.append(proc)

                time.sleep(5)
                proclist = newproc

            # Allocate next available node
            node = nodelist.pop()
            # Make sure we set TMPDIR and PBS_JOBID
            if not "TMPDIR" in os.environ:
                os.environ["TMPDIR"] = ("/tmp/%s" %
                                        (os.environ["PBS_JOBID"]))
            cmd = "/usr/bin/ssh %s \"/bin/sh -c \'TMPDIR=%s;PBS_JOBID=%s;source %s;%s\'\"" % (node, os.environ["TMPDIR"], os.environ["PBS_JOBID"], self.envscript, c)
            print "Running on %s: %s" % (node, cmd)
            proclist.append([subprocess.Popen(cmd,shell=True), node])
            # Ensure unique simids
            time.sleep(5)

        # Wait for all child processes to finish
        if (len(proclist) > 0):
            for proc in proclist:
                proc[0].wait()

        return(0)


if __name__ == '__main__':

    envscript = sys.argv[1]
    cmdfile = sys.argv[2]
    nodefile = sys.argv[3]
    numcores = int(sys.argv[4])

    # Read in command list
    ip = open(cmdfile)
    cmdlist = ip.read().splitlines()
    ip.close()

    # Run the commands
    runobj = RunParallel(envscript)
    runobj.runMultiSSH(nodefile, numcores, cmdlist)

    sys.exit(0)
