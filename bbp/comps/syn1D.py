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

Broadband Platform Module for the UCSB wave propagation codes
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
import uc_fault_utils
from install_cfg import InstallCfg
from syn1D_cfg import Syn1DCfg
from station_list import StationList

class Syn1D(object):
    """
    Implement UCSB syn1D as a Broadband Component
    """

    def __init__(self, i_r_velmodel, i_r_srcfile, i_r_srffile,
                 i_r_stations, i_vmodel_name, sim_id=0):
        """
        Initialize basic class parameters
        """
        self.sim_id = sim_id
        self.r_velmodel = i_r_velmodel
        self.r_srcfile = i_r_srcfile
        self.r_srffile = i_r_srffile
        self.r_stations = i_r_stations
        self.vmodel_name = i_vmodel_name
        self.cfg = None
        self.slo = None
        self.a_indir = None
        self.a_tmpdir = None

    def run_syn1d(self, a_tmpdir_mod, a_velmodel,
                  a_greenfile, a_greensoil, a_lahfile):
        """
        Run the Syn1D simulator with the parameters provided
        """
        # Store cwd and change over to tmpdir so the executable can
        # find the files
        old_cwd = os.getcwd()
        os.chdir(a_tmpdir_mod)

        # Copy velocity model
        r_velmodel = os.path.basename(a_velmodel)
        shutil.copy2(a_velmodel, os.path.join(a_tmpdir_mod, r_velmodel))

        # This is not a UCSB format station list, convert station
        # list to UCSB format, generating the station file and the
        # vs30 file
        a_uc_stations = os.path.join(a_tmpdir_mod, self.cfg.R_UC_STATION_FILE)
        a_uc_vs30 = os.path.join(a_tmpdir_mod, self.cfg.R_UC_VS30_FILE)
        stas2files.gp2uc_stalist(self.slo, a_uc_stations, a_uc_vs30)

        #
        # The UCSB codes require fixed input names.  So here, we copy
        # the UCSB file over to the expected name "stations.ll"
        #
        shutil.copy2(a_uc_stations,
                     os.path.join(a_tmpdir_mod, "stations.ll"))

        # Now copy source_model.list
        shutil.copy2(os.path.join(self.a_tmpdir, self.cfg.R_UC_SOURCE_MODEL),
                     os.path.join(a_tmpdir_mod, self.cfg.R_UC_SOURCE_MODEL))

        #
        # Define SRF file as files in the working directory
        #
        a_ffspfile = os.path.join(self.a_indir, self.cfg.R_FFSP_FILE)
        shutil.copy2(a_ffspfile, a_tmpdir_mod)
        a_ffspfile = os.path.join(a_tmpdir_mod, self.cfg.R_FFSP_FILE)

        #
        # Move a Green_Bank.inf and associated velocity model
        # into the working directory
        #
        shutil.copy2(a_greenfile, os.path.join(a_tmpdir_mod,
                                               os.path.basename(a_greenfile)))

        # Symlink green soil (too big to copy)
        if (not os.path.exists(os.path.join(a_tmpdir_mod,
                                            os.path.basename(a_greensoil)))):
            os.symlink(a_greensoil, os.path.join(a_tmpdir_mod,
                                                 os.path.basename(a_greensoil)))

        #
        # Move syn_1d.inp into working directory
        #
        shutil.copy2(a_lahfile, os.path.join(a_tmpdir_mod,
                                             os.path.basename(a_lahfile)))

        #
        # Create faultGlobal.in
        #
        r_faultfile = "faultGlobal.in"
        a_faultfile = os.path.join(a_tmpdir_mod, r_faultfile)
        uc_fault_utils.uc_create_fault_global(a_faultfile, self.sim_id,
                                              self.r_srcfile, self.vmodel_name,
                                              self.r_velmodel,
                                              self.cfg.R_FFSP_FILE)

        #
        # Convert stations to xy
        #
        cmd = "%s >> %s 2>&1" % (self.cfg.A_SLL2XY, self.log)
        bband_utils.runprog(cmd)

        #
        # Run the bb codes
        #
        cmd = "%s >> %s 2>&1" % (self.cfg.A_SYN1D, self.log)
        bband_utils.runprog(cmd, abort_on_error=True)

        # Restore previous directory
        os.chdir(old_cwd)

    def run_stitch(self, a_tmpdir_stitch, a_tmpdir_lf,
                   a_tmpdir_hf, a_velmodel):
        """
        Run the UCSB Stitch code to merge LF and HF seismograms,
        creating BB seismograms
        """
        # Copy velocity model
        r_velmodel = os.path.basename(a_velmodel)
        shutil.copy2(a_velmodel, os.path.join(a_tmpdir_stitch, r_velmodel))

        # Write VMname.list file
        vmname_file = open(os.path.join(a_tmpdir_stitch, "VMname.list"), 'w')
        site_list = self.slo.getStationList()
        for _ in site_list:
            vmname_file.write("%s\n" % (r_velmodel))
        vmname_file.close()

        #
        # Create faultGlobal.in
        #
        r_faultfile = "faultGlobal.in"
        a_faultfile = os.path.join(a_tmpdir_stitch, r_faultfile)
        uc_fault_utils.uc_create_fault_global(a_faultfile, self.sim_id,
                                              self.r_srcfile, self.vmodel_name,
                                              self.r_velmodel,
                                              self.cfg.R_FFSP_FILE)

        # Copy station list xy to the stitch directory
        r_station_file = "stations.xy"
        shutil.copy2(os.path.join(a_tmpdir_lf, r_station_file),
                     os.path.join(a_tmpdir_stitch, r_station_file))

        # Figure out hypocenter depth
        if self.r_srcfile is not None and self.r_srcfile != "":
            a_srcfile = os.path.join(self.a_indir, self.r_srcfile)
            hypo_dep = fault_utils.calculate_hypo_depth(a_srcfile)
        elif self.r_srffile is not None and self.r_srffile != "":
            sta_base = os.path.basename(os.path.splitext(self.r_stations)[0])
            a_srffile = os.path.join(self.a_indir, self.r_srffile)
            hypo_dep = fault_utils.get_hypocenter(a_srffile, sta_base)[2]
        else:
            # No way to determine hypocenter depth, existing
            print("No way to determine hypocenter depth, exiting!")
            sys.exit(1)
        # Write stitch.inp file
        stitch_inp_file = open(os.path.join(a_tmpdir_stitch, "stitch.inp"), 'w')
        stitch_inp_file.write("%s\n" % (r_station_file))
        stitch_inp_file.write("VMname.list\n")
        stitch_inp_file.write("%s%s\n" % (a_tmpdir_hf, os.sep))
        stitch_inp_file.write("%s%s\n" % (a_tmpdir_lf, os.sep))
        # Merging frequency
        stitch_inp_file.write("1.0,  45.0	!merging frequency, fmax\n")
        stitch_inp_file.write("%.7f	!Hypocenter depth\n" % (hypo_dep))
        stitch_inp_file.write("1	!number of sources\n")
        stitch_inp_file.write("2	!displacement (1), velocity (2),"
                              " acceleration (3)\n")
        stitch_inp_file.close()

        # Save old directory
        old_cwd = os.getcwd()
        os.chdir(a_tmpdir_stitch)

        # All good, run the code!
        cmd = "%s >> %s 2>&1" % (self.cfg.A_STITCH, self.log)
        bband_utils.runprog(cmd, abort_on_error=True)

        # Restore old directory
        os.chdir(old_cwd)

    def run(self):
        """
        Runs the UCSB Syn1D simulator
        """
        print("UCSB Syn1D".center(80, '-'))

        #
        # Global installation parameters
        #
        install = InstallCfg.getInstance()
        #
        # Required inputs are sim_id, the src file, the FFSP output
        # and station list
        #
        sim_id = self.sim_id
        sta_base = os.path.basename(os.path.splitext(self.r_stations)[0])
        self.log = os.path.join(install.A_OUT_LOG_DIR,
                                str(sim_id),
                                "%d.syn1d_%s.log" % (sim_id, sta_base))

        self.a_indir = os.path.join(install.A_IN_DATA_DIR, str(sim_id))
        self.a_tmpdir = os.path.join(install.A_TMP_DATA_DIR, str(sim_id))
        a_tmpdir_lf = os.path.join(install.A_TMP_DATA_DIR, str(sim_id),
                                   "syn1D_lf_%s" % (sta_base))
        a_tmpdir_hf = os.path.join(install.A_TMP_DATA_DIR, str(sim_id),
                                   "syn1D_hf_%s" % (sta_base))
        a_tmpdir_stitch = os.path.join(install.A_TMP_DATA_DIR, str(sim_id),
                                       "stitch_%s" % (sta_base))
        a_outdir = os.path.join(install.A_OUT_DATA_DIR, str(sim_id))

        #
        # Make sure the output and tmp directories exist
        #
        bband_utils.mkdirs([self.a_tmpdir, a_tmpdir_lf, a_tmpdir_hf,
                            a_tmpdir_stitch, a_outdir],
                           print_cmd=False)

        # Parse SRC file
        a_srcfile = os.path.join(self.a_indir, self.r_srcfile)
        self.cfg = Syn1DCfg(self.vmodel_name, a_srcfile)

        # Read station list
        a_stations = os.path.join(self.a_indir, self.r_stations)
        print(a_stations)
        self.slo = StationList(a_stations)
        site_list = self.slo.getStationList()

        # Make sure syn1D can handle our station list
        if len(site_list) > self.cfg.MAX_STATIONS:
            raise bband_utils.ParameterError("Too many stations in "
                                             "the station list: %d. " %
                                             (len(site_list)) +
                                             "Maximum limit is %d." %
                                             (self.cfg.MAX_STATIONS))

        # Run Syn1D for LF
        print("Running Syn1D for LF...")
        self.run_syn1d(a_tmpdir_lf, self.cfg.A_UC_LF_VELMODEL,
                       self.cfg.A_UC_GREENBANK, self.cfg.A_UC_GREEN_SOIL,
                       self.cfg.A_UC_SYN1D_INP_FILE)

        # Run Syn1D for HF
        print("Running Syn1D for HF...")
        self.run_syn1d(a_tmpdir_hf, self.cfg.A_UC_HF_VELMODEL,
                       self.cfg.A_UC_HF_GREENBANK, self.cfg.A_UC_HF_GREEN_SOIL,
                       self.cfg.A_UC_SYN1D_INP_FILE)

        # Run Stitch to combine LF and HF
        print("Running Stitch...")
        self.run_stitch(a_tmpdir_stitch, a_tmpdir_lf, a_tmpdir_hf,
                        self.cfg.A_UC_LF_VELMODEL)

        #
        # Convert the outputs to BB format
        #

        # Copy station list ll to the stitch directory
        r_station_file = "stations.ll"
        shutil.copy2(os.path.join(a_tmpdir_lf, r_station_file),
                     os.path.join(a_tmpdir_stitch, r_station_file))

        # Save old directory
        old_cwd = os.getcwd()
        os.chdir(a_tmpdir_stitch)

        cmd = "%s >> %s 2>&1" % (self.cfg.A_CONV, self.log)
        bband_utils.runprog(cmd)

        # Restore old directory
        os.chdir(old_cwd)

        #
        # Move the results to the tmpdir directory. Use the stations
        # list to determine the names of the output file the system
        # should have produced.  Define an output name for each
        # station BB file.  Read each line in the file as a station.
        #
        for stat in site_list:
            a_tmpfile = os.path.join(a_tmpdir_stitch, "%s.3comp" % (stat.scode))
            expected_file = os.path.join(self.a_tmpdir,
                                         "%d.%s.bbp" % (sim_id, stat.scode))
            shutil.copy2(a_tmpfile, expected_file)

        if self.r_srcfile == "":
            # calculate magnitude and write to file
            mag = fault_utils.get_magnitude(os.path.join(self.a_indir,
                                                         self.r_velmodel),
                                            os.path.join(self.a_indir,
                                                         self.r_srffile),
                                            sta_base)
            mag_file = open(os.path.join(self.a_indir,
                                         "magnitude_%s" % (sta_base)), 'w')
            mag_file.write("%.2f" % mag)
            mag_file.flush()
            mag_file.close()

        print("UCSB Syn1D Completed".center(80, '-'))

if __name__ == "__main__":
    print("Testing Module: %s" % (os.path.basename(sys.argv[0])))
    INSTALL = InstallCfg()
    ME = Syn1D(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4],
               sys.argv[5], sim_id=int(sys.argv[6]))
    ME.run()
