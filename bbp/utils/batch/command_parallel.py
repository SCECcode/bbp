#!/usr/bin/env python
"""
BSD 3-Clause License

Copyright (c) 2021, University of Southern California
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

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
