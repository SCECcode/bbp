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

Compares two sets of reference data for equality within
a tolerance.
"""
from __future__ import division, print_function

import os
import sys
import glob
import cmp_bbp

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("Usage: %s <dir1> <dir2>" % (os.path.basename(sys.argv[0])))
        sys.exit(1)

    dir1 = sys.argv[1]
    dir2 = sys.argv[2]

    if not os.path.exists(dir1):
        print("Dir %s does not exist" % (dir1))
        sys.exit(1)

    if not os.path.exists(dir2):
        print("Dir %s does not exist" % (dir2))
        sys.exit(1)

    # Check for missing, mismatching files in dir2
    testdirs = glob.glob("%s/*" % (dir1))
    for test in testdirs:
        testbase = os.path.basename(test)

        # Compare bbp files
        bbpfiles = glob.glob("%s/*.bbp" % (test))
        for bbpfile in bbpfiles:
            filebase = os.path.basename(bbpfile)
            otherfile = "%s/%s/%s" % (dir2, testbase, filebase)
            print("Comparing:")
            print("\t%s" % (bbpfile))
            print("\t%s" % (otherfile))
            if not os.path.exists(otherfile):
                print("File %s not found" % (otherfile))
            if cmp_bbp.cmp_bbp(bbpfile, otherfile, 0.03) != 0:
                print("Compare failed for %s and %s" %
                      (bbpfile, otherfile))

        # Compare rsp files
        rspfiles = glob.glob("%s/*.rsp" % (test))
        for rspfile in rspfiles:
            filebase = os.path.basename(rspfile)
            otherfile = "%s/%s/%s" % (dir2, testbase, filebase)
            print("Comparing:")
            print("\t%s" % (rspfile))
            print("\t%s" % (otherfile))
            if not os.path.exists(otherfile):
                print("File %s not found" % (otherfile))
            if cmp_bbp.cmp_bbp(rspfile, otherfile, 0.03) != 0:
                print("Compare failed for %s and %s" %
                      (rspfile, otherfile))

    # Check for missing files in dir1
    testdirs = glob.glob("%s/*" % (dir2))
    for test in testdirs:
        testbase = os.path.basename(test)
        bbpfiles = glob.glob("%s/*.bbp" % (test))
        for bbpfile in bbpfiles:
            filebase = os.path.basename(bbpfile)
            otherfile = "%s/%s/%s" % (dir1, testbase, filebase)
            if not os.path.exists(otherfile):
                print("File %s not found" % (otherfile))

        rspfiles = glob.glob("%s/*.rsp" % (test))
        for rspfile in rspfiles:
            filebase = os.path.basename(rspfile)
            otherfile = "%s/%s/%s" % (dir1, testbase, filebase)
            if not os.path.exists(otherfile):
                print("File %s not found" % (otherfile))

    sys.exit(0)
