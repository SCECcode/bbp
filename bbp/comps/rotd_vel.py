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
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys

# Import Broadband modules
import bband_utils
import install_cfg
import bbp_formatter
from station_list import StationList

class RotDVel(object):
    """
    BBP module implementation of rotd50 for velocity provided by UCB.
    Rotd50 inputs seismograms and outputs response spectra
    """
    @staticmethod
    def do_rotd_vel(workdir, peer_input_e_file, peer_input_n_file,
                  peer_input_z_file, output_rotd_file,
                  logfile):
        """
        This function runs the rotd50 command inside workdir, using
        the inputs and outputs specified
        """
        install = install_cfg.InstallCfg.getInstance()

        # Make sure we don't have absolute path names
        peer_input_e_file = os.path.basename(peer_input_e_file)
        peer_input_n_file = os.path.basename(peer_input_n_file)
        peer_input_z_file = os.path.basename(peer_input_z_file)
        output_rotd_file = os.path.basename(output_rotd_file)

        # Save cwd, change back to it at the end
        old_cwd = os.getcwd()
        os.chdir(workdir)

        # Make sure we remove the output files first or Fortran will
        # complain if they already exist
        try:
            os.unlink(output_rotd_file)
        except OSError:
            pass

        #
        # write config file for rotd100 program
        rd50_conf = open("rotd100_inp.cfg", 'w')
        # This flag indicates inputs acceleration
        rd50_conf.write("2 interp flag\n")
        # This flag indicate we are processing two input files (horizontals)
        rd50_conf.write("1 Npairs\n")
        # Number of headers in the file
        rd50_conf.write("6 Nhead\n")
        rd50_conf.write("%s\n" % peer_input_e_file)
        rd50_conf.write("%s\n" % peer_input_n_file)
        rd50_conf.write("%s\n" % output_rotd_file)
        # Close file
        rd50_conf.close()

        progstring = ("%s/rotd100 >> %s 2>&1" % (install.A_UCB_BIN_DIR, logfile))
        bband_utils.runprog(progstring, abort_on_error=True, print_cmd=False)

        # Restore working directory
        os.chdir(old_cwd)

    def __init__(self, i_r_stations, sim_id=0):
        """
        Initializes class variables
        """
        self.sim_id = sim_id
        self.r_stations = i_r_stations

    def run(self):
        print("RotDVel".center(80, '-'))
        #
        # convert input bbp acc files to peer format acc files
        #

        install = install_cfg.InstallCfg.getInstance()
        sim_id = self.sim_id
        sta_base = os.path.basename(os.path.splitext(self.r_stations)[0])
        self.log = os.path.join(install.A_OUT_LOG_DIR, str(sim_id),
                                "%d.rotdvel_%s.log" % (sim_id, sta_base))
        a_statfile = os.path.join(install.A_IN_DATA_DIR,
                                  str(sim_id),
                                  self.r_stations)
        a_tmpdir = os.path.join(install.A_TMP_DATA_DIR, str(sim_id))
        a_outdir = os.path.join(install.A_OUT_DATA_DIR, str(sim_id))

        #
        # Make sure the tmp and out directories exist
        #
        bband_utils.mkdirs([a_tmpdir, a_outdir], print_cmd=False)

        slo = StationList(a_statfile)
        site_list = slo.getStationList()

        for site in site_list:
            stat = site.scode
            print("==> Processing station: %s" % (stat))

            # Create path names and check if their sizes are within bounds
            nsfile = os.path.join(a_tmpdir,
                                  "%d.%s.000" % (sim_id, stat))
            ewfile = os.path.join(a_tmpdir,
                                  "%d.%s.090" % (sim_id, stat))
            udfile = os.path.join(a_tmpdir,
                                  "%d.%s.ver" % (sim_id, stat))
            bbpfile = os.path.join(a_outdir,
                                   "%d.%s.vel.bbp" % (sim_id, stat))

            bband_utils.check_path_lengths([nsfile, ewfile, udfile],
                                           bband_utils.GP_MAX_FILENAME)

            cmd = ("%s/wcc2bbp " % (install.A_GP_BIN_DIR) +
                   "wcc2bbp=0 nsfile=%s ewfile=%s udfile=%s < %s >> %s 2>&1" %
                   (nsfile, ewfile, udfile, bbpfile, self.log))
            bband_utils.runprog(cmd, abort_on_error=True, print_cmd=False)

            # Now we need to convert to peer format
            out_n_vel = os.path.join(a_tmpdir,
                                     "%d.%s.peer_n.vel" % (sim_id, stat))
            out_e_vel = os.path.join(a_tmpdir,
                                     "%d.%s.peer_e.vel" % (sim_id, stat))
            out_z_vel = os.path.join(a_tmpdir,
                                     "%d.%s.peer_z.vel" % (sim_id, stat))
            bbp_formatter.bbp2peer(bbpfile, out_n_vel, out_e_vel, out_z_vel, accel=False)

            # Let's have rotD50 create these output files
            out_rotd_base = "%d.%s.rdvel" % (sim_id, stat)
            tmp_rotd = os.path.join(a_tmpdir, out_rotd_base)
            out_rotd = os.path.join(a_outdir, out_rotd_base)

            # Run the rotD50 program
            self.do_rotd_vel(a_tmpdir, out_e_vel, out_n_vel, out_z_vel,
                           out_rotd, self.log)

            cmd = "cp %s %s" % (tmp_rotd, out_rotd)
            bband_utils.runprog(cmd, abort_on_error=True, print_cmd=False)

        # All done!
        print("RotDVel Completed".center(80, '-'))

if __name__ == '__main__':
    print("Testing Module: %s" % (os.path.basename(sys.argv[0])))
    ME = RotDVel(sys.argv[1], sim_id=int(sys.argv[2]))
    ME.run()
    sys.exit(0)
