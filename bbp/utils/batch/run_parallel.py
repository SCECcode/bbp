#!/usr/bin/env python

import sys
import os
import time
import socket
import subprocess

class RunParallel:

    def __init__(self, envscript):
        self.envscript = envscript
        return

        
    def runMultiSSH(self, nodes, numcores, cmdlist):
        hostname = socket.gethostname()

        # Process list of nodes
        nodelist = []
        nodes = nodes.split(',')

        for node in nodes:
            node = node.strip()
            if node != '':
                for i in xrange(0, numcores):
                    nodelist.append(node)

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
            # Make sure we set TMPDIR and SLURM_JOB_ID
            if not "TMPDIR" in os.environ:
                os.environ["TMPDIR"] = ("/tmp/%s" %
                                        (os.environ["SLURM_JOB_ID"]))
            cmd = "/usr/bin/ssh %s \"/bin/sh -c \'TMPDIR=%s;SLURM_JOB_ID=%s;source %s;%s\'\"" % (node, os.environ["TMPDIR"], os.environ["SLURM_JOB_ID"], self.envscript, c)
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
    nodelist = sys.argv[3]
    numcores = int(sys.argv[4])

    # Read in command list
    ip = open(cmdfile)
    cmdlist = ip.read().splitlines()
    ip.close()

    # Run the commands
    runobj = RunParallel(envscript)
    runobj.runMultiSSH(nodelist, numcores, cmdlist)
    
    sys.exit(0)
