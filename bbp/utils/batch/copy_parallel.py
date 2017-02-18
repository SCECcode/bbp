#!/usr/bin/env python

import sys
import os
import time
import socket
import subprocess

class CopyParallel:

    def __init__(self, envscript):
        self.envscript = envscript
        return

        
    def runMultiSSH(self, remotedir, localdir, nodefile):
        hostname = socket.gethostname()

        # Read in nodefile
        nodelist = []
        if (nodefile == 'localhost'):
            nodelist.append('localhost')
        else:
            nodes = open(nodefile)
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
            nodes.close()

        print "nodelist = ", nodelist

        if (len(nodelist) == 0):
            print "No compute nodes available"
            return(1)
        else:
            print "Running on %s nodes" % (len(nodelist))

        # Command to execute remotely
        c = ("cp -frp %s/* %s/." % (remotedir, localdir))

        # Execute the copy command on each remote node
        proclist = []
        while (len(nodelist) > 0):
            # Use next node
            node = nodelist.pop()
     
            # Make sure we set TMPDIR and PBS_JOBID
            if not "TMPDIR" in os.environ:
                os.environ["TMPDIR"] = ("/tmp/%s" %
                                        (os.environ["PBS_JOBID"]))
            if (node == 'localhost'):
                cmd = ("TMPDIR=%s;PBS_JOBID=%s;source %s;%s" %
                       (os.environ["TMPDIR"], os.environ["PBS_JOBID"],
                        self.envscript, c))
#                cmd = "%s" % (c)
            else:
#                cmd = "/usr/bin/ssh %s \"/bin/sh -c \'source %s;%s\'\"" % (node, self.envscript, c)
                cmd = "/usr/bin/ssh %s \"/bin/sh -c \'TMPDIR=%s;PBS_JOBID=%s;source %s;%s\'\"" % (node, os.environ["TMPDIR"], os.environ["PBS_JOBID"], self.envscript, c)
            
            print "Running on %s: %s" % (node, cmd)   
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
    remotedir = sys.argv[2]
    localdir = sys.argv[3]
    nodefile = sys.argv[4]

    # Run the commands
    runobj = CopyParallel(envscript)
    runobj.runMultiSSH(remotedir, localdir, nodefile)
    
    sys.exit(0)
