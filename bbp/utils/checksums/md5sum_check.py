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

This script calculates md5sums for all the Greens functions and
compares against the expected values.  The path to the greens
functions is set up by install_cfg.py.
"""
# Import Python modules
import os
import sys
import subprocess

# Import Broadband modules
from install_cfg import InstallCfg
import bband_utils

def check_md5s_in_dir(md5_dir, data_dir):
    """
    This function checks the stored md5 sums stored in md5_dir and
    compares them with the files present in data_dir
    """
    return_code = 0
    print "Entering directory %s" % (data_dir)
    for entry in os.listdir(md5_dir):
        # Make sure we have an absolute path
        a_entry = os.path.join(md5_dir, entry)
        if os.path.isdir(a_entry) and entry[0] != '.':
            # Skip hidden dirs
            return_code += check_md5s_in_dir(a_entry,
                                             os.path.join(data_dir, entry))
        elif a_entry.endswith(".md5"):
            # It's a valid md5 file, cut off .md5 extension
            r_gf_filename = entry[0:len(entry) - 4]
            a_gf_filename = os.path.join(data_dir, r_gf_filename)
            try:
                if sys.platform == "linux2":
                    # We are running on linux
                    proc = subprocess.Popen(["md5sum", a_gf_filename],
                                            stdout=subprocess.PIPE)
                    calculated_md5 = proc.communicate()[0].split()[0]
                elif sys.platform == "darwin":
                    # This is a Mac OS X computer
                    proc = subprocess.Popen(["md5", a_gf_filename],
                                            stdout=subprocess.PIPE)
                    calculated_md5 = proc.communicate()[0].split()[3]
                else:
                    print "Unknow platform %s!" % (sys.platform)
                    sys.exit(1)
            except IndexError:
                print "Error calculating md5sum of file %s" % (a_gf_filename)
                sys.exit(1)
            saved_fp = open(a_entry, 'r')
            stored_md5 = saved_fp.readline().split()[0]
            saved_fp.close()
            if calculated_md5 != stored_md5:
                print ("Error: calculated md5sum %s " % (calculated_md5) +
                       "doesn't agree with expected md5sum %s for file %s.\n" %
                       (stored_md5, a_gf_filename))
                return_code += 1
    return return_code

def generate_md5s(data_dir, md5_dir, top_level=False):
    """
    This function generates md5 checksums for the files stored in
    data_dir, and stores the results in md5_dir.
    """

    print "Entering directory %s\n" % data_dir
    # Make sure md5_dir exists
    bband_utils.mkdirs([md5_dir], print_cmd=False)

    for entry in os.listdir(data_dir):
        # Skip top-level checksums directory
        if top_level and entry == "checksums":
            continue
        a_entry = os.path.join(data_dir, entry)
        if os.path.isdir(a_entry) and entry[0] != '.': #skip hidden dirs
            generate_md5s(a_entry, os.path.join(md5_dir, entry))
        else:
            md5_filename = os.path.join(md5_dir, "%s.md5" % (entry))
            filename = os.path.basename(entry)
            if sys.platform == "linux2":
                # We are running on linux
                md5sum = subprocess.Popen(["md5sum", a_entry],
                                          stdout=subprocess.PIPE)
                computed_sum = md5sum.communicate()[0].split()[0]
            elif sys.platform == "darwin":
                # This is a Mac OS X computer
                md5sum = subprocess.Popen(["md5", a_entry],
                                          stdout=subprocess.PIPE)
                computed_sum = md5sum.communicate()[0].split()[3]
            else:
                print "Unknow platform %s!" % (sys.platform)
                sys.exit(1)
            # Now, write md5sum to file
            md5_out = open(md5_filename, 'w')
            md5_out.write("%s %s\n" % (computed_sum, filename))
            md5_out.close()

#-------------------------------------------------------------------------------
# Main program starts here
#-------------------------------------------------------------------------------
CFG = InstallCfg()

if len(sys.argv) == 2:
    if sys.argv[1] == '-g':
        # Generate md5sums instead
        for CHECK_DIR in [CFG.A_GF_DIR, CFG.A_VAL_DIR]:
            SUB_DIRS = bband_utils.list_subdirs(CHECK_DIR)
            for DIRECTORY in SUB_DIRS:
                BASE_DIR = os.path.join(CHECK_DIR, DIRECTORY)
                MD5_DIR = os.path.join(CHECK_DIR, DIRECTORY, "checksums")
                # First, remove checksums directory if it already exists
                if os.path.exists(MD5_DIR):
                    os.system("rm -rf %s" % (MD5_DIR))
                generate_md5s(BASE_DIR, MD5_DIR, True)
        # All done!
        sys.exit(0)

    if sys.argv[1] == '-d':
        # Delete md5sums
        for CHECK_DIR in [CFG.A_GF_DIR, CFG.A_VAL_DIR]:
            SUB_DIRS = bband_utils.list_subdirs(CHECK_DIR)
            for DIRECTORY in SUB_DIRS:
                BASE_DIR = os.path.join(CHECK_DIR, DIRECTORY)
                MD5_DIR = os.path.join(CHECK_DIR, DIRECTORY, "checksums")
                # First, remove checksums directory if it already exists
                if os.path.exists(MD5_DIR):
                    os.system("rm -rf %s" % (MD5_DIR))
        # All done!
        sys.exit(0)

print "Starting verifying checksums..."
STATUS = 0
for CHECK_DIR in [CFG.A_GF_DIR, CFG.A_VAL_DIR]:
    SUB_DIRS = bband_utils.list_subdirs(CHECK_DIR)
    for DIRECTORY in SUB_DIRS:
        BASE_DIR = os.path.join(CHECK_DIR, DIRECTORY)
        MD5_DIR = os.path.join(CHECK_DIR, DIRECTORY, "checksums")
        if os.path.exists(MD5_DIR):
            STATUS += check_md5s_in_dir(MD5_DIR, BASE_DIR)

if STATUS == 0:
    print "All checksums agree!"
sys.exit(STATUS)
