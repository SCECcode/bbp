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
