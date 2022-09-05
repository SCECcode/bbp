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

This converts a RWG format velocity model file to a SDSU version.
Note there are hard coded values in this conversion including the
depth to the moho and the vs values at the moho
"""
from __future__ import division, print_function

# Import Python modules
import sys

def gpvm2sdsuvm(gpfilename, sdsufilename):
    """
    Converts a Graves&Pitarka velocity model into a SDSU velocity
    model (for the srcV1.5 release)
    """
    depth = 0.0
    infile = open(gpfilename, "r")
    outfile = open(sdsufilename, "w")
    # Write header
    outfile.write("% sdsu format velocity model for BB Platform\n")
    outfile.write("% Depth(km)  Vp(km/s)  Vs(km/s)  rho(g/cm3)   Qp       Qs\n")
    # Get the first layer and write it to the output file
    for line in infile:
        elems = line.split()
        # Skip lines without the correct number of items
        if len(elems) == 6:
            break
    outstr = ("%7.4f    %7.4f    %7.4f   %7.4f    %7.4f   %7.4f\n" %
              (depth, float(elems[1]), float(elems[2]),
               float(elems[3]), float(elems[4]), float(elems[5])))
    outfile.write(outstr)
    depth = depth + float(elems[0])
    # Continue processing the file
    for line in infile:
        # Skip lines without the correct number of items
        if len(line.split()) != 6:
            continue
        outstr = ("%7.4f    %7.4f    %7.4f   %7.4f    %7.4f   %7.4f\n" %
                  ((depth - 0.0001), float(elems[1]), float(elems[2]),
                   float(elems[3]), float(elems[4]), float(elems[5])))
        outfile.write(outstr)
        elems = line.split()
        outstr = ("%7.4f    %7.4f    %7.4f   %7.4f    %7.4f   %7.4f\n" %
                  ((depth + 0.0001), float(elems[1]), float(elems[2]),
                   float(elems[3]), float(elems[4]), float(elems[5])))
        outfile.write(outstr)
        depth = depth + float(elems[0])
    # Now, add the last line
    outstr = ("%7.4f    %7.4f    %7.4f   %7.4f    %7.4f   %7.4f\n" %
              (depth, float(elems[1]), float(elems[2]),
               float(elems[3]), float(elems[4]), float(elems[5])))
    outfile.write(outstr)

    # All done, close files
    infile.close()
    outfile.close()

def gpvm2sdsuvm_old(gpfilename, sdsufilename):
    """
    Converts a Graves&Pitarka velocity model into a SDSU velocity
    model (for the srcV1.4 release, but with 6 columns)
    """
    depth = 0.0
    infile = open(gpfilename, "r")
    lines = infile.readlines()
    infile.close()
    outfile = open(sdsufilename, "w")
    outfile.write("% sdsu format velocity model for BB Platform\n")
    #outfile.write("%   z     vp      vs     rho       Qp       Qs\n")
    outfile.write("% Depth(km)  Vp(km/s)  Vs(km/s)  rho(g/cm3)   Qp       Qs\n")
    # Skip first line
    for line in lines[1:-1]:
        elems = line.split()
        #str = ("%7.3f %7.3f %7.3f %7.3f %7.0f %7.0f\n" %
        #       (dd, float(elems[1]), float(elems[2]),
        #        float(elems[3]), float(elems[4]), float(elems[5])))
        outstr = ("%7.3f    %7.3f    %7.3f   %7.3f    %7.3f   %7.3f\n" %
                  (depth, float(elems[1]), float(elems[2]),
                   float(elems[3]), float(elems[4]), float(elems[5])))
        depth = float(elems[0]) + depth
        outfile.write(outstr)
    # Now add the last line
    line = lines[-1]
    elems = line.split()
    outstr = ("%7.2f    %7.3f    %7.3f   %7.3f    %7.3f   %7.3f\n" %
              (1024.00, float(elems[1]),
               float(elems[2]), float(elems[3]),
               float(elems[4]), float(elems[5])))
    outfile.write(outstr)
    outfile.close()

def velmodellayers_sdsu(ifile):
    """
    Input a SDSU type velocity model file and return number of layers
    as defined by SDSU. This is designed for use in the SDSU code
    which required the number of layers in a file.
    """
    lincount = 0
    infile = open(ifile, "r")
    for _ in infile:
        lincount = lincount + 1
    infile.close()
    if lincount < 3:
        return 0
    else:
        return lincount-2

def gpvm2ucsbvm(gpfilename, ucsbfilename):
    """
    Converts a Graves&Pitarka velocity model into a UCSB velocity model
    """
    # dd = 0.0
    infile = open(gpfilename, "r")
    lines = infile.readlines()
    infile.close()
    ofile = ucsbfilename
    outfile = open(ofile, "w")
    # lent = len(lines)
    for i, line in enumerate(lines):
        if i == 0:
            elem = line.split()
            outstr = "%d 1.0\n" % (int(elem[0]))
            outfile.write(outstr)
        else:
            elems = []
            elems = line.split()
            outstr = ("%.2f %.2f %.2f %.2f %f %f\n" %
                      (float(elems[1]), float(elems[2]), float(elems[3]),
                       float(elems[0]), float(elems[4]), float(elems[5])))
            outfile.write(outstr)
    # Close file
    outfile.close()

if __name__ == "__main__":
    print("Test Config Class: %s" % sys.argv[0])
    gpvm2ucsbvm(sys.argv[1], sys.argv[2])
