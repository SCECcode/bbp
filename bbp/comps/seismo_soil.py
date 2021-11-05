#!/usr/bin/env python
"""
Copyright 2010-2021 University Of Southern California

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys
import tempfile

import numpy as np

import matplotlib as mpl
mpl.use("AGG")

# Import Broadband modules
import bband_utils
from station_list import StationList
from install_cfg import InstallCfg
from seismo_soil_cfg import SeismoSoilCfg
from PySeismoSoil.class_ground_motion import Ground_Motion
from PySeismoSoil.class_Vs_profile import Vs_Profile
from PySeismoSoil.class_site_effect_adjustment import Site_Effect_Adjustment

class SeismoSoil(object):
    """
    This class contains the glue code needed to interface the
    Broadband Platform and the PySeismoSoil module
    """
    def __init__(self, i_r_srcfile, i_r_velmodel, method,
                 i_r_stations, vmodel_name, sim_id=0, debug=False):
        """
        This function initializes basic class objects
        """
        self.sim_id = sim_id
        self.r_srcfile = i_r_srcfile
        self.r_velmodel = i_r_velmodel
        self.method = method
        self.r_stations = i_r_stations
        self.vmodel_name = vmodel_name
        self.stat_list = None
        self.install = None
        self.config = None
        self.log = None
        self.vs_profile = None
        self.debug = debug

    def read_velocity_model(self, a_velmodel):
        """
        This function parses the velocity model file and
        converts h and vs to meters.
        """
        # Initialize velocity model structure
        vmodel = {'h': [],
                  'vs': []}

        vel_file = open(a_velmodel, 'r')
        for line in vel_file:
            line = line.strip()
            pieces = line.split()
            # Skip lines without the 6 values
            if len(pieces) != 6:
                continue
            pieces = [float(piece) for piece in pieces]
            # Convert to meters
            thick = pieces[0] * 1000
            vs = pieces[2] * 1000
            if vs > 1000:
                # Stop when vs > 1000
                break
            vmodel['h'].append(thick)
            vmodel['vs'].append(vs)
        vel_file.close()

        vmodel['h'].append(0.0)
        vmodel['vs'].append(1000.0)

        vs_profile = np.zeros((len(vmodel['h']), 2))
        for index, thick, vs in zip(range(len(vmodel['h'])),
                                          vmodel['h'], vmodel['vs']):
            vs_profile[index][0] = thick
            vs_profile[index][1] = vs

        self.vs_profile = Vs_Profile(vs_profile)

        if self.debug:
            a_param_outdir = os.path.join(self.install.A_OUT_DATA_DIR,
                                          str(self.sim_id),
                                          "param_files")
            vs_profile_fn = os.path.join(a_param_outdir,
                                         "seismosoil.%s.vs_profile.txt" %
                                         (str(self.sim_id)))
            print(self.vs_profile, file=open(vs_profile_fn, 'w'))

    def write_single_component_seismograms(self, input_bbp_file):
        """
        This functions reads the input_bbp_file and writes the
        individual components separately, it returns three filenames
        for the first, second, and third component.
        """
        f_input = open(input_bbp_file, 'r')

        # Create temporary files for the three components
        f_comp1 = tempfile.NamedTemporaryFile(mode='w', delete=False)
        f_comp2 = tempfile.NamedTemporaryFile(mode='w', delete=False)
        f_comp3 = tempfile.NamedTemporaryFile(mode='w', delete=False)

        for line in f_input:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#") or line.startswith("%"):
                continue
            pieces = line.split()
            f_comp1.write("%s\t%s\n" % (pieces[0], pieces[1]))
            f_comp2.write("%s\t%s\n" % (pieces[0], pieces[2]))
            f_comp3.write("%s\t%s\n" % (pieces[0], pieces[3]))

        # All done, save filenames and close files
        f1_name = f_comp1.name
        f2_name = f_comp2.name
        f3_name = f_comp3.name

        f_input.close()
        f_comp1.close()
        f_comp2.close()
        f_comp3.close()

        # Return filenames
        return f1_name, f2_name, f3_name

    def convert_to_bbp(self, f_comp1, f_comp2, f_comp3, out_acc_file):
        """
        This function merged the three input files into a single
        BBP acceleration file
        """
        comp1_file = open(f_comp1, 'r')
        comp2_file = open(f_comp2, 'r')
        comp3_file = open(f_comp3, 'r')
        bbp_file = open(out_acc_file, 'w')
        # Add header to BBP file
        bbp_file.write("#    time(sec)      N-S(cm/s/s)      "
                       "E-W(cm/s/s)      U-D(cm/s/s)\n")
        for line1, line2, line3 in zip(comp1_file, comp2_file, comp3_file):
            line1 = line1.strip()
            line2 = line2.strip()
            line3 = line3.strip()
            pieces1 = line1.split()
            pieces1 = [float(val) for val in pieces1]
            pieces2 = line2.split()
            pieces2 = [float(val) for val in pieces2]
            pieces3 = line3.split()
            pieces3 = [float(val) for val in pieces3]
            bbp_file.write("%15.6e%15.6e%15.6e%15.6e\n" %
                           (pieces1[0], pieces1[1], pieces2[1], pieces3[1]))
        bbp_file.close()
        comp1_file.close()
        comp2_file.close()
        comp3_file.close()

    def create_vel_bbp(self, stat):
        """
        This function derives a velocity bbp file from an acceleration
        file
        """
        install = self.install
        sim_id = self.sim_id
        a_tmpdir = os.path.join(install.A_TMP_DATA_DIR, str(sim_id))
        a_outdir = os.path.join(install.A_OUT_DATA_DIR, str(sim_id))
        in_acc_file = os.path.join(a_outdir, "%d.%s.acc.bbp" %
                                   (sim_id, stat))
        out_vel_file = os.path.join(a_outdir, "%d.%s.vel.bbp" %
                                    (sim_id, stat))

        # Since we have acceleration files, we need to integrate to
        # get to velocity

        # Create path names and check if their sizes are within bounds
        nsfile = os.path.join(a_tmpdir,
                              "%d.%s.acc.000" % (sim_id, stat))
        ewfile = os.path.join(a_tmpdir,
                              "%d.%s.acc.090" % (sim_id, stat))
        udfile = os.path.join(a_tmpdir,
                              "%d.%s.acc.ver" % (sim_id, stat))
        bband_utils.check_path_lengths([nsfile, ewfile, udfile],
                                       bband_utils.GP_MAX_FILENAME)

        cmd = ("%s" % (os.path.join(install.A_GP_BIN_DIR, "wcc2bbp")) +
               " wcc2bbp=0 nsfile=%s ewfile=%s udfile=%s < %s >> %s 2>&1" %
               (nsfile, ewfile, udfile, in_acc_file, self.log))
        bband_utils.runprog(cmd, abort_on_error=True, print_cmd=False)
        # Now we need to integrate to get to velocity
        for comp in ['000', '090', 'ver']:
            file_in = os.path.join(a_tmpdir,
                                   "%d.%s.acc.%s" % (sim_id, stat, comp))
            file_out = os.path.join(a_tmpdir,
                                    "%s.%s.vel.%s" % (sim_id, stat, comp))
            bband_utils.check_path_lengths([file_in, file_out],
                                           bband_utils.GP_MAX_FILENAME)
            cmd = ("%s" % (os.path.join(install.A_GP_BIN_DIR, "integ_diff")) +
                   " integ=1 filein=%s fileout=%s" % (file_in, file_out))
            bband_utils.runprog(cmd, abort_on_error=True, print_cmd=False)
        # Now we put together the file again as a velocity bbp file
        nsfile = os.path.join(a_tmpdir,
                              "%d.%s.vel.000" % (sim_id, stat))
        ewfile = os.path.join(a_tmpdir,
                              "%d.%s.vel.090" % (sim_id, stat))
        udfile = os.path.join(a_tmpdir,
                              "%d.%s.vel.ver" % (sim_id, stat))
        bband_utils.check_path_lengths([nsfile, ewfile, udfile],
                                       bband_utils.GP_MAX_FILENAME)

        cmd = ("%s" % (os.path.join(install.A_GP_BIN_DIR, "wcc2bbp")) +
               " wcc2bbp=1 nsfile=%s ewfile=%s udfile=%s > %s 2>> %s" %
               (nsfile, ewfile, udfile, out_vel_file, self.log))
        bband_utils.runprog(cmd, abort_on_error=True, print_cmd=False)

    def run(self):
        """
        This function processes a station list and applies the
        PySeismoSoil site amplification module to the two
        horizontal components of each station.
        """
        print("PySeismoSoil".center(80, '-'))

        self.install = InstallCfg.getInstance()
        install = self.install
        self.config = SeismoSoilCfg(self.vmodel_name, self.method)
        config = self.config
        sim_id = self.sim_id

        sta_base = os.path.basename(os.path.splitext(self.r_stations)[0])
        self.log = os.path.join(install.A_OUT_LOG_DIR,
                                str(sim_id),
                                "%d.pyseismosoil_%s.log" % (sim_id, sta_base))
        a_indir = os.path.join(install.A_IN_DATA_DIR, str(sim_id))
        a_tmpdir = os.path.join(install.A_TMP_DATA_DIR, str(sim_id))
        a_outdir = os.path.join(install.A_OUT_DATA_DIR, str(sim_id))
        a_param_outdir = os.path.join(a_outdir, "param_files")
        #
        # Make sure the output and tmp directories exist
        #
        bband_utils.mkdirs([a_indir, a_tmpdir, a_outdir, a_param_outdir])

        a_velmodel = os.path.join(a_indir, self.r_velmodel)
        a_stations = os.path.join(a_indir, self.r_stations)
        self.stat_list = StationList(a_stations)

        # Read velocity model
        self.read_velocity_model(a_velmodel)
        vref = config.LF_VREF

        # In Debug mode we output the velocity profile
        if self.debug:
            print(self.vs_profile)

        # In Debug mode we save a plot of the velocity profile
        if self.debug:
            import pylab
            fig = pylab.plt.figure(figsize=(3.0, 4.0), dpi=100)
            self.vs_profile.plot(fig=fig) #fig=fig, figsize=(2.6, 3.2), dpi=100)
            fig.savefig(os.path.join(a_param_outdir,
                                     "seismosoil.%d.vs_profile.png" % (sim_id)),
                        format="png")
            mpl.pyplot.close('all')

        #
        # Read and parse the station list with this call
        #
        slo = StationList(a_stations)
        site_list = slo.getStationList()

        for site in site_list:
            station_name = site.scode
            station_vs30 = site.vs30

            print("Station: %s - Vs30: %f ..."  % (station_name, station_vs30))
            bbp_file = os.path.join(a_outdir,
                                    "%d.%s.acc.bbp" %
                                    (sim_id, station_name))

            f_comp1, f_comp2, f_comp3 = self.write_single_component_seismograms(bbp_file)

            # Read input timeseries
            gm1 = Ground_Motion(f_comp1, unit='cm/s/s', motion_type='accel')
            gm2 = Ground_Motion(f_comp2, unit='cm/s/s', motion_type='accel')
            gm3 = Ground_Motion(f_comp3, unit='cm/s/s', motion_type='accel')

            if vref < 1000:
                # Input motion may already contain some soil response, need
                # to cancel out this soil response by performing linear
                # deconvolution (as per PySeismoSoil examples/Pipeline_01_Site
                # Effects_Adjustments.ipynb)
                gm1 = gm1.deconvolve(self.vs_profile, boundary='elastic',
                                     show_fig=False)
                gm2 = gm2.deconvolve(self.vs_profile, boundary='elastic',
                                     show_fig=False)
            else:
                # On the other hand, if the motion was simulated with a
                # VREF > 1000, we may be missing some site amplification
                # from propagating throught the rock layers
                gm1 = gm1.amplify(self.vs_profile, show_fig=False)
                gm2 = gm2.amplify(self.vs_profile, show_fig=False)

            # Now apply the non linear site effects
            sea_1 = Site_Effect_Adjustment(gm1, station_vs30, lenient=True)
            sea_2 = Site_Effect_Adjustment(gm2, station_vs30, lenient=True)
            # REMOVE
            # gm1.save_accel(os.path.join(a_tmpdir, "test_acc.txt"), unit='cm/s/s')
            output_motion_1, fig1, _ = sea_1.run(show_fig=True, return_fig_obj=True)
            output_motion_2, fig2, _ = sea_2.run(show_fig=True, return_fig_obj=True)

            # In Debug mode, we save these plots
            if self.debug:
                fig1.savefig(os.path.join(a_param_outdir,
                                          "seismosoil.%d.%s-comp1.png" %
                                          (sim_id, station_name)),
                             format="png")
                fig2.savefig(os.path.join(a_param_outdir,
                                          "seismosoil.%d.%s-comp2.png" %
                                          (sim_id, station_name)),
                             format="png")
                mpl.pyplot.close('all')

            # Write output BBP file
            output_motion_1.save_accel(f_comp1, unit='cm/s/s')
            output_motion_2.save_accel(f_comp2, unit='cm/s/s')
            gm3.save_accel(f_comp3, unit='cm/s/s')

            # Need to recreate acceleration BBP file
            self.convert_to_bbp(f_comp1, f_comp2, f_comp3, bbp_file)
            self.create_vel_bbp(station_name)

            # Delete temp files
            os.unlink(f_comp1)
            os.unlink(f_comp2)
            os.unlink(f_comp3)

        print("PySeismoSoil Completed".center(80, '-'))

if __name__ == "__main__":
    print("Testing Module: %s" % (sys.argv[0]))
    if len(sys.argv) < 7:
        print("Usage: %s src_file vel_model_file method "
              "station_list vmodel_name sim_id [debug_flag]" % (sys.argv[0]))
        sys.exit(-1)
    DEBUG = False
    if len(sys.argv) > 7:
        DEBUG = int(sys.argv[7]) > 0
    ME = SeismoSoil(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4],
                    sys.argv[5], sim_id=int(sys.argv[6]), debug=DEBUG)
    ME.run()
