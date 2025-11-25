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

class RotD50(object):
    """
    BBP module implementation of rotd50 provided by UCB.
    Rotd50 inputs seismograms and outputs response spectra
    """
    @staticmethod
    def do_rotd50(workdir, peer_input_e_file, peer_input_n_file,
                  peer_input_z_file, output_rotd50_file,
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
        output_rotd50_file = os.path.basename(output_rotd50_file)

        # Save cwd, change back to it at the end
        old_cwd = os.getcwd()
        os.chdir(workdir)

        # Make sure we remove the output files first or Fortran will
        # complain if they already exist
        try:
            os.unlink(output_rotd50_file)
        except OSError:
            pass

        #
        # write config file for rotd50 program
        rd50_conf = open("rotd50_inp.cfg", 'w')
        # This flag indicates inputs acceleration
        rd50_conf.write("2 interp flag\n")
        # This flag indicate we are processing two input files (horizontals)
        rd50_conf.write("1 Npairs\n")
        # Number of headers in the file
        rd50_conf.write("6 Nhead\n")
        rd50_conf.write("%s\n" % peer_input_e_file)
        rd50_conf.write("%s\n" % peer_input_n_file)
        rd50_conf.write("%s\n" % output_rotd50_file)
        # Close file
        rd50_conf.close()

        progstring = ("%s/rotd50 >> %s 2>&1" % (install.A_UCB_BIN_DIR, logfile))
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
        print("RotD50".center(80, '-'))
        #
        # convert input bbp acc files to peer format acc files
        #

        install = install_cfg.InstallCfg.getInstance()
        sim_id = self.sim_id
        sta_base = os.path.basename(os.path.splitext(self.r_stations)[0])
        self.log = os.path.join(install.A_OUT_LOG_DIR, str(sim_id),
                                "%d.rotd50_%s.log" % (sim_id, sta_base))
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
        site_list = slo.get_station_list()

        for site in site_list:
            stat = site.scode
            print("==> Processing station: %s" % (stat))

            # Since we have velocity files, we need to differentiate
            # to get to acceleration

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

            for c in ["090", "000", "ver"]:
                # Differentiate to get from velocity to accl needed by rotd50
                # Create path names and check if their sizes are within bounds
                filein = os.path.join(a_tmpdir,
                                      "%d.%s.%s" %
                                      (sim_id, stat, c))
                fileout = os.path.join(a_tmpdir,
                                       "%d.%s.acc.%s" %
                                       (sim_id, stat, c))

                bband_utils.check_path_lengths([filein, fileout],
                                               bband_utils.GP_MAX_FILENAME)

                cmd = ("%s/integ_diff diff=1 " % (install.A_GP_BIN_DIR) +
                       "filein=%s fileout=%s >> %s 2>&1" %
                       (filein, fileout, self.log))
                bband_utils.runprog(cmd, abort_on_error=True, print_cmd=False)

                # Check file length
                bband_utils.check_path_lengths(["%s/%d.%s.acc.%s" %
                                                (a_tmpdir, sim_id, stat, c)],
                                               bband_utils.GP_MAX_FILENAME)

            # Now we need to convert them back to bbp
            # Create path names and check if their sizes are
            # within bounds
            nsfile = os.path.join(a_tmpdir,
                                  "%d.%s.acc.000" % (sim_id, stat))
            ewfile = os.path.join(a_tmpdir,
                                  "%d.%s.acc.090" % (sim_id, stat))
            udfile = os.path.join(a_tmpdir,
                                  "%d.%s.acc.ver" % (sim_id, stat))
            bbpfile = os.path.join(a_tmpdir,
                                   "%d.%s.acc.bbp" % (sim_id, stat))

            bband_utils.check_path_lengths([nsfile, ewfile, udfile],
                                           bband_utils.GP_MAX_FILENAME)

            cmd = ("%s/wcc2bbp " % (install.A_GP_BIN_DIR) +
                   "nsfile=%s ewfile=%s udfile=%s " %
                   (nsfile, ewfile, udfile) +
                   "units=cm/s/s wcc2bbp=1 > %s 2>> %s" %
                   (bbpfile, self.log))
            bband_utils.runprog(cmd, abort_on_error=True, print_cmd=False)

            # Now we need to convert to peer format
            out_n_acc = os.path.join(a_tmpdir,
                                     "%d.%s.peer_n.acc" % (sim_id, stat))
            out_e_acc = os.path.join(a_tmpdir,
                                     "%d.%s.peer_e.acc" % (sim_id, stat))
            out_z_acc = os.path.join(a_tmpdir,
                                     "%d.%s.peer_z.acc" % (sim_id, stat))
            bbp_formatter.bbp2peer(bbpfile, out_n_acc, out_e_acc, out_z_acc)

            # Let's have rotD50 create these output files
            out_rotd50_base = "%d.%s.rd50" % (sim_id, stat)
            tmp_rotd50 = os.path.join(a_tmpdir, out_rotd50_base)
            out_rotd50 = os.path.join(a_outdir, out_rotd50_base)

            # Run the rotD50 program
            self.do_rotd50(a_tmpdir, out_e_acc, out_n_acc, out_z_acc,
                           out_rotd50, self.log)

            cmd = "cp %s %s" % (tmp_rotd50, out_rotd50)
            bband_utils.runprog(cmd, abort_on_error=True, print_cmd=False)

        # All done!
        print("RotD50 Completed".center(80, '-'))

if __name__ == '__main__':
    print("Testing Module: %s" % (os.path.basename(sys.argv[0])))
    ME = RotD50(sys.argv[1], sim_id=int(sys.argv[2]))
    ME.run()
    sys.exit(0)
