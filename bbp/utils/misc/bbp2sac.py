#!/usr/bin/env python
"""
Copyright 2010-2019 University Of Southern California

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

Utility to convert BBP time history files to SAC format
"""

# Import python modules
import os
import sys

# Import Broadband Modules
import bband_utils
from install_cfg import InstallCfg

if len(sys.argv) < 2:
    print "Usage: bbp2sac input_bbp_file"
    sys.exit(1)

input_file = sys.argv[1]

# First pass, figure DT and number of samples
samples = 0
first_time = 0.0
second_time = 0.0
dt = 0.0
ifile = open(input_file)
for line in ifile:
    # Skip comments
    if line.startswith("#") or line.startswith("%"):
        continue
    pieces = line.split()
    samples = samples + 1
    if samples == 1:
        first_time = float(pieces[0])
    if samples == 2:
        second_time = float(pieces[0])
ifile.close()

# Calculate DT
dt = second_time - first_time

# Create file names
base_file = os.path.splitext(input_file)[0]
ns_file = "%s.000" % (base_file)
ew_file = "%s.090" % (base_file)
ud_file = "%s.ver" % (base_file)

# Second pass, create the output files
ifile = open(input_file)
nsfile = open(ns_file, 'w')
ewfile = open(ew_file, 'w')
udfile = open(ud_file, 'w')

# Write headers
nsfile.write("\t%d   %1.9E\n" % (samples, dt))
ewfile.write("\t%d   %1.9E\n" % (samples, dt))
udfile.write("\t%d   %1.9E\n" % (samples, dt))

for line in ifile:
    # Skip comments
    if line.startswith("#") or line.startswith("%"):
        continue
    pieces = line.split()
    pieces = [float(piece) for piece in pieces]
    nsfile.write(" %1.9E\n" % (pieces[1]))
    ewfile.write(" %1.9E\n" % (pieces[2]))
    udfile.write(" %1.9E\n" % (pieces[3]))

# All done, close everything
nsfile.close()
ewfile.close()
udfile.close()
ifile.close()

# Get pointers
install = InstallCfg.getInstance()

# Now convert them into sac format
for comp in ["000", "090", "ver"]:
    file_in = "%s.%s" % (base_file, comp)
    file_out = "%s.sac" % (file_in)
    cmd = "echo '%s' > tmp" % (file_in)
    bband_utils.runprog(cmd)
    cmd = ("%s < tmp >> /dev/null 2>&1" %
           (os.path.join(install.A_UCSB_BIN_DIR,
                         "BBPtoSAC")))
    bband_utils.runprog(cmd)
    os.rename("output.sac", file_out)
    os.unlink(file_in)

# Remove tmp file
os.unlink("tmp")
