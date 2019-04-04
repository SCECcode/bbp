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

UCSB makeFault Global (mfg) file
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys
import math
import numpy
import random

# Import Broadband modules
from ucrmg_cfg import UCrmgCfg
from install_cfg import InstallCfg

def uc_create_ffsp_inp(a_ffsp_inp, a_srcfile, i_vmodel_name):
    """
    This function creates the ffsp.inp file with the parameters
    specified in the SRC file
    """
    # Read source file
    cfg = UCrmgCfg(i_vmodel_name, a_srcfile)

    # Calculate some parameters
    deg2rad = numpy.pi / 180.
    cxp = cfg.CFGDICT["fault_length"] / 2. + cfg.CFGDICT["hypo_along_stk"]
    cyp = cfg.CFGDICT["hypo_down_dip"]
    depth_hypc = (cfg.CFGDICT["depth_to_top"] +
                  cfg.CFGDICT["hypo_down_dip"] *
                  numpy.sin(cfg.CFGDICT["dip"] * deg2rad))

    # Top center of fault is in the origin
    xref_hypc = (numpy.cos(cfg.CFGDICT["strike"] * deg2rad) *
                 cfg.CFGDICT["hypo_along_stk"] -
                 numpy.cos(cfg.CFGDICT["dip"] * deg2rad) *
                 numpy.sin(cfg.CFGDICT["strike"] * deg2rad)
                 * cfg.CFGDICT["hypo_down_dip"])
    yref_hypc = (numpy.sin(cfg.CFGDICT["strike"] * deg2rad) *
                 cfg.CFGDICT["hypo_along_stk"] +
                 numpy.cos(cfg.CFGDICT["dip"] * deg2rad) *
                 numpy.cos(cfg.CFGDICT["strike"] * deg2rad)
                 * cfg.CFGDICT["hypo_down_dip"])

    moment = 10.**(1.5 * cfg.CFGDICT['magnitude'] + 9.05)

    nsubx = int(2**math.ceil(math.log(cfg.CFGDICT["fault_length"] /
                                      cfg.CFGDICT["dlen"] + 1, 2)))
    nsuby = int(2**math.ceil(math.log(cfg.CFGDICT["fault_width"] /
                                      cfg.CFGDICT["dwid"] + 1, 2)))

    # Calculate random numbers
    random.seed(cfg.CFGDICT['seed'])
    rdm1 = 100000 * random.random()
    rdm2 = 100000 * random.random()
    rdm3 = 100000 * random.random()

    # Now, write the ffsp_inp file
    ffsp_inp_file = open(a_ffsp_inp, "w")
    ffsp_inp_file.write("8 %.2f\n" % (cfg.FMAX / 2.0))
    ffsp_inp_file.write("%.2f %.2f\n" % (cfg.CFGDICT["fault_length"],
                                         cfg.CFGDICT["fault_width"]))
    ffsp_inp_file.write("%.2f %.2f %.2f\n" % (cxp, cyp, depth_hypc))
    ffsp_inp_file.write("%.2f, %.2f\n" % (xref_hypc, yref_hypc))
    ffsp_inp_file.write("%1.2e %.2f %.2f\n" %
                        (moment, cfg.CFGDICT['corner_freq'], cfg.RV_AVG))
    ffsp_inp_file.write("%.2f\n" % (cfg.TP_TR))
    ffsp_inp_file.write("%.f. %.f. %.f.\n" % (cfg.CFGDICT["strike"],
                                              cfg.CFGDICT["dip"],
                                              cfg.CFGDICT["rake"]))
    ffsp_inp_file.write("20.0   30.0\n")
    ffsp_inp_file.write("%d   %d\n" % (nsubx, nsuby))
    ffsp_inp_file.write("5  5  5  5\n")
    ffsp_inp_file.write("%d %d %d\n" % (-rdm1, -rdm2, -rdm3))
    ffsp_inp_file.write("1 1\n")
    ffsp_inp_file.write("%s\n" % (os.path.basename(cfg.A_UC_LF_VELMODEL)))
    ffsp_inp_file.write("0.0\n")
    ffsp_inp_file.write("1\n")
    ffsp_inp_file.write("%s\n" % (cfg.FFSP_OUTPUT_PREFIX))
    ffsp_inp_file.close()

def uc_create_fault_global(a_faultfile, sim_id, r_srcfile,
                           i_vmodel_name, r_velmodel, r_srffile):
    """
    This fuction creates the faultfile with the parameters specified
    in the src_file and velocity model
    """
    install = InstallCfg.getInstance()
    a_srcfile = os.path.join(install.A_IN_DATA_DIR, str(sim_id), r_srcfile)
    # since KinModel appends .srf automatically, make sure it's not
    # already there
    if r_srffile.endswith(".srf"):
        r_srffile = r_srffile[0:len(r_srffile) - 4]

    # Read source file
    cfg = UCrmgCfg(i_vmodel_name, a_srcfile)

    # Write configuration file
    print("Creating %s" % (a_faultfile))

    fault_file = open(a_faultfile, "w")
    fault_file.write("%.3f %.3f %.1f\n" % (cfg.CFGDICT["lon_top_center"],
                                           cfg.CFGDICT["lat_top_center"],
                                           cfg.CFGDICT["depth_to_top"]))
    fault_file.write("%.2f %.2f\n" % (cfg.CFGDICT["fault_length"],
                                      cfg.CFGDICT["fault_width"]))
    fault_file.write("%.f. %.f. %.f.\n" % (cfg.CFGDICT["strike"],
                                           cfg.CFGDICT["dip"],
                                           cfg.CFGDICT["rake"]))
    fault_file.write("%.1f %.1f\n" % (cfg.CFGDICT["hypo_along_stk"],
                                      cfg.CFGDICT["hypo_down_dip"]))
    fault_file.write("%4.2f\n" % (cfg.CFGDICT['magnitude']))
    fault_file.write("%.3f %.3f\n" % (cfg.CFGDICT["dlen"],
                                      cfg.CFGDICT["dwid"]))
    fault_file.write("%d\n" % (cfg.CFGDICT['seed']))
    fault_file.write("%.2f\n" % (cfg.DT))
    fault_file.write("%.2f\n" % (cfg.CFGDICT['corner_freq']))
    fault_file.write("%s\n" % (r_velmodel))
    fault_file.write("%s\n" % (r_srffile))
    fault_file.close()

def uc_create_kin_model(a_fault_file, a_kin_model_file,
                        corner_wave_number,
                        max_rand_dist):
    """
    This function reads the faultGlobal file (a_fault_file) and writes
    the KinModel.inp file (a_kin_model_file). It replaces the call to
    the old convertFaultInputKM program.
    """
    # Reads faultGlobal file
    print("Reading input file: %s" % (a_fault_file))
    # Open input file, read it line by line
    fault_f = open(a_fault_file)
    [top_cen_lon,
     top_cen_lat,
     top_cen_z] = [float(val) for val in fault_f.readline().split()]
    top_cen_x = 0.0
    top_cen_y = 0.0
    [rup_len, dd_width] = [float(val) for val in fault_f.readline().split()]
    [strike,
     dip,
     rake] = [float(val) for val in fault_f.readline().split()]
    [hypo_strike,
     hypo_dip] = [float(val) for val in fault_f.readline().split()]
    [mw] = [float(val) for val in fault_f.readline().split()]
    [dx, dy] = [float(val) for val in fault_f.readline().split()]
    [random_seed] = [int(val) for val in fault_f.readline().split()]
    [dt] = [float(val) for val in fault_f.readline().split()]
    [fc] = [float(val) for val in fault_f.readline().split()]
    velmod = fault_f.readline().split()[0]
    name = fault_f.readline().split()[0]
    # Close input file
    fault_f.close()

    # Perform some convertions
    rad = math.pi / 180.0
    hypo_x = top_cen_x + math.cos(strike * rad) * hypo_strike
    hypo_x = hypo_x - math.cos(dip * rad) * math.sin(strike * rad) * hypo_dip
    hypo_y = top_cen_y + math.sin(strike * rad) * hypo_strike
    hypo_y = hypo_y + math.cos(dip * rad) * math.cos(strike * rad) * hypo_dip
    hypo_z = top_cen_z + math.sin(dip * rad) * hypo_dip
    m0 = math.pow(10, 1.5 * mw + 9.05)
    if random_seed < 0:
        random_seed = abs(random_seed)
    dx2 = dx * 1000.0
    # nxx = math.floor(rup_len / dx)
    # nyy = math.floor(dd_width / dy)
    # nx = 2
    # while (nx < nxx):
    #     nx = nx * 2.0
    # ny = 2
    # while (ny < nyy):
    #     ny = ny * 2.0
    # dx = 1.0 * rup_len / 1.0 * nx
    # dy = 1.0 * dd_width / 1.0 * ny
    hypo_strike = hypo_strike + rup_len / 2.0

    # Writes KinModel.inp file
    print("Writing output file: %s" % (a_kin_model_file))
    kin_f = open(a_kin_model_file, 'w')
    kin_f.write("%lf %lf\n" % (rup_len * 1000.0, dd_width * 1000.0))
    kin_f.write("%lf %lf\n" % (hypo_strike * 1000.0, hypo_dip * 1000.0))
    kin_f.write("%lf %lf %lf\n" % (hypo_x * 1000.0,
                                   hypo_y * 1000.0,
                                   hypo_z * 1000.0))
    kin_f.write("%g %lf\n" % (m0, fc))
    kin_f.write("%lf %lf %lf\n" % (strike, dip, rake))
    kin_f.write("%lf %lf\n" % (dx2, dt))
    kin_f.write("1\n")
    kin_f.write("%d %d %d\n" % (random_seed + 23123,
                                random_seed + 123213,
                                random_seed + 454532))
    kin_f.write("%s %.2f %.2f\n" % (velmod, corner_wave_number, max_rand_dist))
    # Close output file
    kin_f.close()

if __name__ == "__main__":
    print("Test Config Class: %s" % (os.path.basename(sys.argv[0])))
