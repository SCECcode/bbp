#!/usr/bin/env python

from __future__ import division, print_function

import sys
import os
import time
import socket
import subprocess

class CommandParallel:

    def __init__(self, envscript):
        self.envscript = envscript
        return

    def runMultiSSH(self, remotecmd, nodes):
        hostname = socket.gethostname()

        # Read in nodefile
        nodelist = []
        nodes = nodes.split(',')
        for node in nodes:
            node = node.strip(' \n\t')
            if ((node != hostname) and (node != '')):
                if node in nodelist:
                    continue
                nodelist.append(node)
            elif node == hostname:
                if 'localhost' in nodelist:
                    continue
                nodelist.append('localhost')

        print("nodelist = ", nodelist)

        if (len(nodelist) == 0):
            print("No compute nodes available")
            return(1)
        else:
            print("Running on %s nodes" % (len(nodelist)))

        # Execute the copy command on each remote node
        proclist = []
        while (len(nodelist) > 0):
            # Use next node
            node = nodelist.pop()
            # Make sure we set TMPDIR and SLURM_JOB_ID
            if not "TMPDIR" in os.environ:
                os.environ["TMPDIR"] = ("/tmp/%s" %
                                        (os.environ["SLURM_JOB_ID"]))
            if (node == 'localhost'):
                cmd = ("TMPDIR=%s;SLURM_JOB_ID=%s;source %s;%s" %
                       (os.environ["TMPDIR"], os.environ["SLURM_JOB_ID"],
                        self.envscript, remotecmd))
            else:
                cmd = "/usr/bin/ssh -o \"ServerAliveInterval 60\" %s \"/bin/sh -c \'TMPDIR=%s;SLURM_JOB_ID=%s;source %s;%s\'\"" % (node, os.environ["TMPDIR"], os.environ["SLURM_JOB_ID"], self.envscript, remotecmd)
            
            print("Running on %s: %s" % (node, cmd))
            proclist.append([subprocess.Popen(cmd,shell=True), node])
            
            # Sleep a bit
            time.sleep(5)
            
        # Wait for all child processes to finish
        if (len(proclist) > 0):
            for proc in proclist:
                proc[0].wait()
            
        return(0)

if __name__ == '__main__':
    envscript = sys.argv[1]
    remotecmd = sys.argv[2]
    nodelist = sys.argv[3]

    # Run the commands
    runobj = CommandParallel(envscript)
    runobj.runMultiSSH(remotecmd, nodelist)
    
    sys.exit(0)
