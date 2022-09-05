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

Broadband Platform Version of UCSB site response
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys
import shutil

# Import Broadband modules
import bband_utils
import stas2files
import fault_utils
from install_cfg import InstallCfg
from uc_site_cfg import UCSiteCfg
from station_list import StationList

class UCSite(object):
    """
    Implements UCSB site response as a Broadband Component
    """

    def __init__(self, i_r_velocity, i_r_srcfile,
                 i_r_stations, sim_id=0):
        self.sim_id = sim_id
        self.r_srcfile = i_r_srcfile
        self.r_velocity = i_r_velocity
        self.r_stations = i_r_stations
        self.cfg = None

    def create_fault_global_in(self, a_fault_file):
        """
        This function creates the faultglobal.in file from data
        obtained in the src file
        """
        out_fp = open(a_fault_file, "w")
        out_fp.write("%.3f %.3f %2.3f\n" %
                     (self.cfg.CFGDICT['lon_top_center'],
                      self.cfg.CFGDICT['lat_top_center'],
                      self.cfg.CFGDICT["depth_to_top"]))
        out_fp.write("%.2f %.2f\n" %
                     (self.cfg.CFGDICT["fault_length"],
                      self.cfg.CFGDICT["fault_width"]))
        out_fp.write("%.1f %.1f %.1f\n" %
                     (self.cfg.CFGDICT['strike'],
                      self.cfg.CFGDICT['dip'],
                      self.cfg.CFGDICT['rake']))
        out_fp.write("%.2f %.2f\n" %
                     (self.cfg.CFGDICT["hypo_along_stk"],
                      self.cfg.CFGDICT["hypo_down_dip"]))
        out_fp.write("%.2f\n" % (self.cfg.CFGDICT['magnitude']))
        out_fp.write("-1 -1\n")
        out_fp.write("%d\n" % (self.cfg.CFGDICT['seed']))
        out_fp.write("-1\n")
        out_fp.write("-1\n")
        out_fp.write("%s\n" % (self.r_velocity))
        out_fp.write("FFSP_OUTPUT\n")
        out_fp.close()

    def separate_stats(self, install, a_uc_vs30, input_dt):
        """
        need to construct modified station list with only
        non-class A stations, since class A stations
        shouldn't be stitched.
        """
        stations_to_stitch = []
        vs30_fp = open(a_uc_vs30, 'r')
        data = vs30_fp.readlines()
        vs30_fp.close()
        comps = ['000', '090', 'ver']
        for line in data:
            pieces = line.split()
            name = pieces[0]
            vs30 = float(pieces[1])
            if vs30 >= 1524:
                my_class = 'A'
            elif vs30 >= 762:
                my_class = 'B'
                my_type = 'rock'
            elif vs30 >= 555:
                my_class = 'BC'
                my_type = 'rock'
            elif vs30 >= 366:
                my_class = 'C'
                my_type = 'sand'
            elif vs30 >= 270:
                my_class = 'CD'
                my_type = 'sand'
            elif vs30 >= 183:
                my_class = 'D'
                my_type = 'sand'
            else:
                my_class = 'E'
                my_type = 'sand'

            r_toSAC = 'toSAC'
            a_toSAC = os.path.join(install.A_UCSB_BIN_DIR, r_toSAC)

            r_noah = 'noah_w'
            a_noah = os.path.join(install.A_UCSB_BIN_DIR, r_noah)

            if not my_class == 'A':
                stations_to_stitch.append(name)

                cmd = "cp %s/input%s.iwan input.iwan" % (install.A_UCSB_DATA_DIR, my_class)
                bband_utils.runprog(cmd)

                cmd = "cp %s/Q%s.par Q.par" % (install.A_UCSB_DATA_DIR, my_class)
                bband_utils.runprog(cmd)
                cmd = "cp %s/tremorka.%s ." % (install.A_UCSB_DATA_DIR, my_type)
                bband_utils.runprog(cmd)
                cmd = "cp %s/intermed.par ." % (install.A_UCSB_DATA_DIR)
                bband_utils.runprog(cmd)

                # dt * iskip = input_dt in intermed.par
                # Lines 1 and 8 of input?.iwan should be adjusted
                #
                intermed_in = open("intermed.par", "r")
                intermed_data = intermed_in.readlines()
                intermed_in.close()
                intermed_out = open("intermed.par.tmp", "w")
                for line in intermed_data:
                    if line.split()[1] == 'dt':
                        dt = float(line.split()[0])
                        intermed_out.write(line)
                    elif line.split()[1] == 'iskip':
                        iskip = int(input_dt / dt)
                        if not dt * iskip == input_dt:
                            print(dt, iskip, input_dt)
                        intermed_out.write(" %d\t\tiskip\n" % iskip)
                    else:
                        intermed_out.write(line)
                intermed_out.flush()
                intermed_out.close()
                shutil.move("intermed.par.tmp", "intermed.par")

                input_in = open("input.iwan", 'r')
                input_out = open("input.iwan.tmp", 'w')
                input_data = input_in.readlines()
                input_in.close()

                for i in range(0, len(input_data)):
                    if i == 0:
                        val = min([10, 0.8 * 1.0 / (2 * input_dt)])
                        input_out.write("     %f\n" % val)
                    elif i == 7:
                        val = min([15, 0.8 * 1.0 / (2 * input_dt)])
                        input_out.write("     %f\n" % val)
                    else:
                        input_out.write(input_data[i])
                input_out.flush()
                input_out.close()

                shutil.move("input.iwan.tmp", "input.iwan")

            for comp in comps:
                cmd = "echo '%s.%s.de1D.001' > tmp" % (name, comp)
                bband_utils.runprog(cmd)
                cmd = "%s < tmp >> %s 2>&1" % (a_toSAC, self.log)
                bband_utils.runprog(cmd)

                if my_class == 'A':
                    continue

                cmd = ("%s tremorka.%s %s.%s.de1D.001.sac %s.%s.nl1D.001 >> %s 2>&1" %
                       (a_noah, my_type, name, comp, name, comp, self.log))
                bband_utils.runprog(cmd)

        return stations_to_stitch

    def run(self):
        """
        Runs the UCSB site response program
        """
        print("UCSB Site".center(80, '-'))

        #
        # Global installation parameters
        #
        install = InstallCfg.getInstance()
        #
        # Required Inputs are sim_id, SRC file, and station list
        #

        sim_id = self.sim_id
        sta_base = os.path.basename(os.path.splitext(self.r_stations)[0])
        self.log = os.path.join(install.A_OUT_LOG_DIR,
                                str(sim_id),
                                "%d.uc_site_%s.log" % (sim_id, sta_base))
        a_indir = os.path.join(install.A_IN_DATA_DIR, str(sim_id))
        a_tmpdir = os.path.join(install.A_TMP_DATA_DIR, str(sim_id))
        a_tmpdir_mod = os.path.join(install.A_TMP_DATA_DIR, str(sim_id),
                                    "uc_site_%s" % (sta_base))
        a_outdir = os.path.join(install.A_OUT_DATA_DIR, str(sim_id))

        a_velocity = os.path.join(a_indir, self.r_velocity)

        #
        # Make sure the output and tmp directories exist
        #
        bband_utils.mkdirs([a_tmpdir, a_tmpdir_mod, a_outdir],
                           print_cmd=False)

        # Parse SRC file
        if self.r_srcfile is None or self.r_srcfile == "":
            raise bband_utils.ParameterError("SRC file not defined!")

        a_srcfile = os.path.join(a_indir, self.r_srcfile)
        self.cfg = UCSiteCfg(a_srcfile)
        cfg = self.cfg

        # Store cwd and change over to tmpdir so the executable can
        # find the files
        old_cwd = os.getcwd()
        os.chdir(a_tmpdir_mod)

        # Copy velocity file to tmpdir_mod
        shutil.copy2(a_velocity, os.path.join(a_tmpdir_mod, self.r_velocity))

        # Read station list
        a_stations = os.path.join(a_indir, self.r_stations)
        print(a_stations)
        slo = StationList(a_stations)
        site_list = slo.getStationList()

        # This is not a UCSB format station list, convert station
        # list to UCSB format, generating the station file and the
        # vs30 file
        a_uc_stations = os.path.join(a_tmpdir_mod, cfg.R_UC_STATION_FILE)
        a_uc_vs30 = os.path.join(a_tmpdir_mod, cfg.R_UC_VS30_FILE)
        stas2files.gp2uc_stalist(slo, a_uc_stations, a_uc_vs30)

        #
        # The UCSB codes require fixed input names.  So here, we copy
        # the UCSB file over to the expected name "stations.ll"
        #
        cmd = ("cp %s %s" % (a_uc_stations,
                             os.path.join(a_tmpdir_mod, "stations.ll")))
        bband_utils.runprog(cmd)

        # Copy .bbp files over to .3comp
        # If we have anything but just a hybrid file, combine them first
        # Use site 0 as the dummy

        for site in site_list:
            if os.path.exists("%s/%d.%s.bbp" % (a_tmpdir, sim_id, site.scode)):
                shutil.copy2("%s/%d.%s.bbp" % (a_tmpdir, sim_id, site.scode),
                             "%s/%s.3comp" % (a_tmpdir_mod, site.scode))
            elif os.path.exists("%s/%s.3comp" % (a_tmpdir, site.scode)):
                shutil.copy2("%s/%s.3comp" % (a_tmpdir, site.scode),
                             "%s/%s.3comp" % (a_tmpdir_mod, site.scode))

        # determine dt for input seismogram
        bbp_fp = open("%s/%s.3comp" % (a_tmpdir_mod, site.scode), 'r')
        bbp_data = bbp_fp.readlines()
        bbp_fp.close()
        i = 0
        while bbp_data[i][0] == '%' or bbp_data[i][0] == '#':
            i += 1

        t0 = float(bbp_data[i].split()[0])
        t1 = float(bbp_data[i+1].split()[0])
        input_dt = t1-t0
        print("input_dt: %f\n" % (input_dt))

        #
        # Create deconvBBP.inp, stitchBBP.inp, VMname.list
        #

        dBBP_in = open("deconvBBP.inp", "w")
        dBBP_in.write("%s\n" % self.r_velocity)
        dBBP_in.write("0.29\n")
        dBBP_in.write("1\n")
        dBBP_in.write("%d\n" % len(site_list))
        for site in site_list:
            dBBP_in.write("%s\n" % site.scode)
        dBBP_in.flush()
        dBBP_in.close()

        sBBP_in = open("stitchBBP.inp", "w")
        sBBP_in.write("stations.xy\n")
        sBBP_in.write("VMname.list\n")
        sBBP_in.write("./\n")
        sBBP_in.write("./\n")
        sBBP_in.write("1.0, 15.0\n")

        # depth of hypocenter
        hypo_dep = fault_utils.calculate_hypo_depth(a_srcfile)

        sBBP_in.write("%s\n" % hypo_dep)
        sBBP_in.write("1\n")
        sBBP_in.write("2\n")
        sBBP_in.flush()
        sBBP_in.close()

        vMname_in = open("VMname.list", "w")
        for site in site_list:
            vMname_in.write("%s\n" % self.r_velocity)
        vMname_in.flush()
        vMname_in.close()

        #
        # Create stations.xy if it doesn't exist yet
        #
        if not os.path.exists("stations.xy"):
            #
            # Create faultGlobal.in
            #
            r_faultfile = "faultGlobal.in"
            a_faultfile = os.path.join(a_tmpdir_mod, r_faultfile)
            self.create_fault_global_in(a_faultfile)

            #
            # Convert stations to xy
            #
            cmd = "%s >> %s 2>&1" % (self.cfg.A_SLL2XY, self.log)
            bband_utils.runprog(cmd)

        #
        # Deconvolve
        #
        cmd = "%s >> %s 2>&1" % (self.cfg.A_UC_DECON_EXE, self.log)
        bband_utils.runprog(cmd)

        #
        # Logic of separateStats.csh pulled out into function
        #
        stations_to_stitch = self.separate_stats(install, a_uc_vs30, input_dt)

        #
        # Stitch
        #
        # Update station files to only stitch non class A stations
        # (Class A stations don't have a non-linear component)
        # Must use 'stations.xy' because it's in stitchBBP.inp
        #
        shutil.copy2("stations.xy", "stations.xy.orig")
        shutil.copy2("stations.ll", "stations.ll.orig")
        station_in = open("stations.xy.orig", 'r')
        station_ll_in = open("stations.ll.orig", "r")
        station_ll_data = station_ll_in.readlines()
        station_data = station_in.readlines()
        station_in.close()
        station_ll_in.close()
        station_out = open("stations.xy", "w")
        station_ll_out = open("stations.ll", "w")
        pieces = station_data[0].split()
        station_out.write("%d %f %f %f\n" % (len(stations_to_stitch),
                                             float(pieces[1]),
                                             float(pieces[2]),
                                             float(pieces[3])))
        station_ll_out.write("%d\n" % len(stations_to_stitch))
        i = 1
        while i < len(station_data):
            inList = False
            stat_data_name = station_data[i].strip()
            for site in stations_to_stitch:
                if stat_data_name == site:
                    inList = True
                    break
            if inList:
                station_out.write("%s\n" % stat_data_name)
                station_out.write("%s" % station_data[i+1])
                station_ll_out.write("%s" % station_ll_data[(i+1)//2])
            i += 2
        station_out.flush()
        station_ll_out.flush()
        station_out.close()
        station_ll_out.close()

        cmd = "%s >> %s 2>&1" % (self.cfg.A_STITCH, self.log)
        bband_utils.runprog(cmd)

        #
        # Copy original stations file back in
        #
        shutil.copy2("stations.xy", "stations.xy.stitch")
        shutil.copy2("stations.xy.orig", "stations.xy")

        # Convert to 3-component seismograms
        #
        cmd = "%s/conv3CompBB >> %s 2>&1" % (install.A_UCSB_BIN_DIR, self.log)
        bband_utils.runprog(cmd)

        shutil.copy2("stations.ll", "stations.ll.stitch")
        shutil.copy2("stations.ll.orig", "station.ll")

        # Move the results to the output directory, as bbp format
        for result_file in os.listdir(a_tmpdir_mod):
            dot_index = result_file.rfind('.3comp')
            if dot_index > -1:
                basename = result_file[0:dot_index]
                shutil.copy2(result_file, "%s/%d.%s.vel.bbp" % (a_outdir,
                                                                sim_id,
                                                                basename))
                shutil.copy2(result_file, "%s/%d.%s.vel.bbp" % (a_tmpdir,
                                                                sim_id,
                                                                basename))
                shutil.copy2(result_file, "%s/%s.3comp" % (a_tmpdir, basename))

                # Create acceleration seismogram

                # Create path names and check if their sizes are
                # within bounds
                nsfile = os.path.join(a_tmpdir,
                                      "%d.%s.000" % (sim_id, basename))
                ewfile = os.path.join(a_tmpdir,
                                      "%d.%s.090" % (sim_id, basename))
                udfile = os.path.join(a_tmpdir,
                                      "%d.%s.ver" % (sim_id, basename))
                bbpfile = os.path.join(a_tmpdir,
                                       "%d.%s.vel.bbp" % (sim_id, basename))

                bband_utils.check_path_lengths([nsfile, ewfile, udfile],
                                               bband_utils.GP_MAX_FILENAME)

                cmd = ("%s/wcc2bbp " % (install.A_GP_BIN_DIR) +
                       "nsfile=%s ewfile=%s udfile=%s " %
                       (nsfile, ewfile, udfile) +
                       "wcc2bbp=0 < %s >> %s 2>&1" %
                       (bbpfile, self.log))
                bband_utils.runprog(cmd, abort_on_error=True)

                for comp in cfg.COMPS:
                    # Differentiate each component
                    filein = os.path.join(a_tmpdir,
                                          "%d.%s.%s" %
                                          (sim_id, basename, comp))
                    fileout = os.path.join(a_tmpdir,
                                           "%d.%s.acc.%s" %
                                           (sim_id, basename, comp))

                    bband_utils.check_path_lengths([filein, fileout],
                                                   bband_utils.GP_MAX_FILENAME)

                    cmd = ("%s/integ_diff diff=1 filein=%s fileout=%s" %
                           (install.A_GP_BIN_DIR, filein, fileout))
                    bband_utils.runprog(cmd, abort_on_error=True)

                # Create path names and check if their sizes are
                # within bounds
                nsfile = os.path.join(a_tmpdir,
                                      "%d.%s.acc.000" % (sim_id, basename))
                ewfile = os.path.join(a_tmpdir,
                                      "%d.%s.acc.090" % (sim_id, basename))
                udfile = os.path.join(a_tmpdir,
                                      "%d.%s.acc.ver" % (sim_id, basename))
                bbpfile = os.path.join(a_tmpdir,
                                       "%d.%s.acc.bbp" % (sim_id, basename))

                bband_utils.check_path_lengths([nsfile, ewfile, udfile],
                                               bband_utils.GP_MAX_FILENAME)

                cmd = ("%s/wcc2bbp " % (install.A_GP_BIN_DIR) +
                       "nsfile=%s ewfile=%s udfile=%s " %
                       (nsfile, ewfile, udfile) +
                       "units=cm/s/s wcc2bbp=1 > %s 2>> %s" %
                       (bbpfile, self.log))
                bband_utils.runprog(cmd, abort_on_error=True)

                # Copy acceleration bbp file to outdir
                shutil.copy2(os.path.join(a_tmpdir, "%d.%s.acc.bbp" %
                                          (sim_id, basename)),
                             os.path.join(a_outdir, "%d.%s.acc.bbp" %
                                          (sim_id, basename)))

        os.chdir(old_cwd)

        print("UCSB Site Completed".center(80, '-'))

if __name__ == "__main__":
    print("Testing Module: %s" % (os.path.basename(sys.argv[0])))
    ME = UCSite(sys.argv[1], sys.argv[2], sys.argv[3],
                sim_id=int(sys.argv[4]))
    ME.run()
