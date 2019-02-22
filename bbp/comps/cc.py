#!/usr/bin/env python
"""
Copyright 2010-2018 University Of Southern California

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

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
    site_list = slo.getStationList()
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
