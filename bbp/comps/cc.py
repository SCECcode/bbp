#!/usr/bin/env python3
"""
BSD 3-Clause License

Copyright (c) 2023, University of Southern California
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

Coordinate Conversion methods for BB platform
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys
import math

# Import Broadband modules
from install_cfg import InstallCfg

def oll2xy(slo, ofile, mlon, mlat, x_azim):
    """
    invoke Rob's cc program to convert from lat/lon to xy
    This inputs a station list object.
    This will work as long at the station list has both lat/lon and XY in it.
    """
    # install = InstallCfg()
    site_list = slo.get_station_list()
    out_file = open(ofile, "w")
    out_file.write("% station coordinates in x,y system, "
                   "origin at fault center\n")
    out_file.write("%    X (north)     Y (east)\n")
    for site in site_list:
        x_north = round(float(site.north))
        y_east = round(float(site.east))
        out_str = "%10d  %10d\n" % (x_north, y_east)
        out_file.write(out_str)
    out_file.close()
    return ofile

def ll2xy(ifile, ofile, mlon, mlat, x_azim):
    """
    Invoke Rob's cc program to convert from lat/lon to xy
    This inputs a station list file.
    Currently, the slo contains lat/lon and x/y. But since the user may
    not be able to do this conversion, we will convert it here again for them.
    """
    install = InstallCfg()
    ll2xy_bin = os.path.join(install.A_GP_BIN_DIR, "ll2xy")
    cmd = ("%s mlon=%f mlat=%f xazim=%f < %s" %
           (ll2xy_bin, mlon, mlat, x_azim, ifile))
    print(cmd)
    oput = os.popen(cmd)
    out_file = open(ofile, "w")
    out_file.write("% station coordinates in x,y system, "
                   "origin at fault center\n")
    out_file.write("%    X (north)     Y (east)\n")
    for line in oput:
        elems = []
        elems = line.split()
        x_north = round(float(elems[0]))
        y_east = round(float(elems[1]))
        out_str = "%10d  %10d\n" % (x_north, y_east)
        out_file.write(out_str)

    out_file.close()
    return ofile

def xy2ll(ifile, ofile, mlon, mlat, x_azim):
    """
    invoke Rob's cc program to convert from xy to ,lat/lon
    Input file name will be returned with a changed extension.
    Error returned if file does not end with xy
    """
    install = InstallCfg()
    ll2xy_bin = os.path.join(install.A_GP_BIN_DIR, "ll2xy")
    cmd = ("%s mlon=%f mlat=%f xazim=%f < %s" %
           (ll2xy_bin, mlon, mlat, x_azim, ifile))
    print(cmd)
    oput = os.popen(cmd)

    out_file = open(ofile, "w")
    for line in oput:
        out_file.write(line)
    out_file.close()
    return ofile

def find_fx_fy_fz(hypalongstrike, faultlen, dip, hypdowndip, depthtotop):
    """
    This function finds the hypocenter's x/y/z
    """
    mcfg = []
    fault_x = hypalongstrike - (0.5*faultlen)
    mcfg.append(fault_x)
    fault_y = math.cos(dip*math.pi/180.0) * hypdowndip
    mcfg.append(fault_y)
    fault_z = depthtotop + (math.sin(dip*math.pi/180.0) * hypdowndip)
    mcfg.append(fault_z)
    return mcfg

if __name__ == "__main__":
    ll2xy(sys.argv[1], sys.argv[2], float(sys.argv[3]),
          float(sys.argv[4]), float(sys.argv[5]))
