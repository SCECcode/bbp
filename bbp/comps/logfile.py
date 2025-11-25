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

Class to log file using a sequence number and return a pathname to the log file.
"""
from __future__ import division, print_function

import os
import sys
import time
import datetime
import pwd
import seqnum
from install_cfg import InstallCfg

class LogFile(object):
    """
    You can create a log file by constructing with no sequence number.
    This file will assign a new random number (based on epoch time), or
    You can pass in a sequence number and this will create a file with
    with that sequence number. You might want to do this if you want to
    use the sequence number in other places in your code.
    """

    def __init__(self, snum=None):
        if snum == None:
            self.sim_id = seqnum.get_seq_num()
        else:
            self.sim_id = int(snum)
        install = InstallCfg.getInstance()
        logdir = install.A_OUT_LOG_DIR
        self.outlogfile = logdir + "/%d/%d.txt" % (self.sim_id, self.sim_id)

    def compose_header(self):
        lines = []
        l1 = ("# SCEC BB Platform Log File for Simulation ID: %d\n" %
              (self.sim_id))
        lines.append(l1)
        t = datetime.datetime.utcnow()
        now = datetime.datetime.fromtimestamp(time.mktime(t.timetuple()))
        dat = now.ctime()
        l2 = "# Start time: %s\n" % (dat)
        lines.append(l2)
        uid = os.getuid()
        uname = pwd.getpwuid(uid)
        l3 = "# UserID: %s\n" % (uname[0])
        lines.append(l3)
        return lines

    def getlogfile(self):
        print("Opening logfile:%s" % (self.outlogfile))
        if os.path.exists(self.outlogfile) == 1:
            fin = open(self.outlogfile, "w+")
        else:
            fin = open(self.outlogfile, "w")
        header = self.compose_header()
        if fin:
            for line in header:
                fin.write(line)
            fin.close()
        else:
            print("Error opening log file %s:" % (self.outlogfile))
            sys.exit(-1)
        return self.outlogfile

if __name__ == "__main__":
    logg = LogFile()
    lf = logg.getlogfile()
    mlf = open(lf, "r")
    mylines = mlf.readlines()
    mlf.close()
    oline = logg.compose_header()
    linnum = 0
    for lin in mylines:
        if lin != oline[linnum]:
            print("Error matching header")
            print("this: %s" % (lin))
            print("that: %s" % (oline[linnum]))
        linnum += 1

    sys.exit(0)
