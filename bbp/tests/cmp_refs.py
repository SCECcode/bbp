#!/usr/bin/env python
"""
SCopyright 2010-2018 University Of Southern California

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

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
