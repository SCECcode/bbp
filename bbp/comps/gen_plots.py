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

Broadband module that generates rotd50 and seismogram overlay plots
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys

# Import Broadband modules
import bband_utils
import bbp_formatter
import arias_duration
import plot_seismograms
import plot_rotd50
from install_cfg import InstallCfg
from station_list import StationList

class GenPlots(object):

    def __init__(self, i_r_stations, i_format,
                 i_comparison_label, sim_id=0):
        self.sim_id = sim_id
        self.r_stations = i_r_stations
        self.format = i_format
        self.comp_label = i_comparison_label

    def run(self):
        print("Generating Plots".center(80, '-'))

        # Initialize basic variables
        install = InstallCfg.getInstance()
        sim_id = self.sim_id
        sta_base = os.path.basename(os.path.splitext(self.r_stations)[0])

        self.log = os.path.join(install.A_OUT_LOG_DIR,
                                str(sim_id),
                                "%d.gen_plots.log" % (sim_id))


        # Input, tmp, and output directories
        a_tmpdir = os.path.join(install.A_TMP_DATA_DIR, str(sim_id))
        a_outdir = os.path.join(install.A_OUT_DATA_DIR, str(sim_id))
        a_tmpdir_seis = os.path.join(install.A_TMP_DATA_DIR, str(sim_id),
                                     "obs_seis_%s" % (sta_base))

        # Station file
        a_statfile = os.path.join(install.A_IN_DATA_DIR,
                                  str(sim_id),
                                  self.r_stations)
        # List of observed seismogram files
        filelist = os.listdir(a_tmpdir_seis)

        slo = StationList(a_statfile)
        site_list = slo.getStationList()

        for site in site_list:
            stat = site.scode

            # Look for the files we need
            bbpfile = os.path.join(a_tmpdir_seis, "%s.bbp" % stat)
            expected_file = os.path.join(a_outdir, "%d.%s.vel.bbp" %
                                         (sim_id, stat))
            if (not os.path.exists(expected_file) or
                not os.path.exists(bbpfile)):
                # just skip this station
                continue

            print("==> Plotting seismogram comparison for station: %s" % (stat))
            if self.format == 'vel':
                # We have velocity, nothing we need to do
                filename1 = bbpfile
            elif self.format == 'acc':
                # We have acceleration, must integrate first
                # Create path names and check if their sizes are within bounds
                nsfile = os.path.join(a_tmpdir, "temp.acc.000")
                ewfile = os.path.join(a_tmpdir, "temp.acc.090")
                udfile = os.path.join(a_tmpdir, "temp.acc.ver")
                bband_utils.check_path_lengths([nsfile, ewfile, udfile],
                                               bband_utils.GP_MAX_FILENAME)

                cmd = ("%s/wcc2bbp " % (install.A_GP_BIN_DIR) +
                       "nsfile=%s ewfile=%s udfile=%s " %
                       (nsfile, ewfile, udfile) +
                       "wcc2bbp=0 < %s >> %s 2>&1" %
                       (bbpfile, self.log))
                bband_utils.runprog(cmd, abort_on_error=True, print_cmd=False)

                for comp in ['000', '090', 'ver']:
                    # Create path names and check if their sizes are
                    # within bounds
                    filein = os.path.join(a_tmpdir,
                                          "temp.acc.%s" % (comp))
                    fileout = os.path.join(a_tmpdir,
                                           "temp.vel.%s" % (comp))

                    bband_utils.check_path_lengths([filein, fileout],
                                                   bband_utils.GP_MAX_FILENAME)

                    cmd = ("%s/integ_diff integ=1 " % (install.A_GP_BIN_DIR) +
                           "filein=%s fileout=%s >> %s 2>&1" %
                           (filein, fileout, self.log))
                    bband_utils.runprog(cmd, abort_on_error=True,
                                        print_cmd=False)

                # Create path names and check if their sizes are within bounds
                nsfile = os.path.join(a_tmpdir, "temp.vel.000")
                ewfile = os.path.join(a_tmpdir, "temp.vel.090")
                udfile = os.path.join(a_tmpdir, "temp.vel.ver")
                vel_bbp_file = os.path.join(a_tmpdir, "temp.%s.vel" % stat)

                bband_utils.check_path_lengths([nsfile, ewfile, udfile],
                                               bband_utils.GP_MAX_FILENAME)

                cmd = ("%s/wcc2bbp wcc2bbp=1 " % install.A_GP_BIN_DIR +
                       "nsfile=%s ewfile=%s udfile=%s > %s 2>> %s" %
                       (nsfile, ewfile, udfile, vel_bbp_file, self.log))
                bband_utils.runprog(cmd, abort_on_error=True, print_cmd=False)
                filename1 = vel_bbp_file

            # Generate arias duration files for calculated data
            calc_acc = os.path.join(a_outdir, "%d.%s.acc.bbp" %
                                    (sim_id, stat))
            calc_peer_n = os.path.join(a_tmpdir, "%d.%s_N.acc" %
                                       (sim_id, stat))
            calc_peer_e = os.path.join(a_tmpdir, "%d.%s_E.acc" %
                                       (sim_id, stat))
            calc_peer_z = os.path.join(a_tmpdir, "%d.%s_Z.acc" %
                                       (sim_id, stat))
            # Convert calculated acc seismogram into peer format
            bbp_formatter.bbp2peer(calc_acc, calc_peer_n,
                                   calc_peer_e, calc_peer_z)

            # Now calculate arias duration for each component
            for comp in ["N", "E", "Z"]:
                file_in = os.path.join(a_tmpdir, "%d.%s_%s.acc" %
                                       (sim_id, stat, comp))
                file_out = os.path.join(a_tmpdir, "%d.%s_%s.arias" %
                                        (sim_id, stat, comp))
                arias_duration.ad_from_acc(file_in, file_out)

            # Generate arias duration files for observed data
            obs_acc = os.path.join(a_tmpdir_seis, "%s.bbp" % stat)
            obs_peer_n = os.path.join(a_tmpdir, "obs.%s_N.acc" %
                                      (stat))
            obs_peer_e = os.path.join(a_tmpdir, "obs.%s_E.acc" %
                                      (stat))
            obs_peer_z = os.path.join(a_tmpdir, "obs.%s_Z.acc" %
                                      (stat))
            # Convert observed acc seismogram into peer format
            bbp_formatter.bbp2peer(obs_acc, obs_peer_n,
                                   obs_peer_e, obs_peer_z)

            # Now calculate arias duration for each component
            for comp in ["N", "E", "Z"]:
                file_in = os.path.join(a_tmpdir, "obs.%s_%s.acc" %
                                       (stat, comp))
                file_out = os.path.join(a_tmpdir, "obs.%s_%s.arias" %
                                        (stat, comp))
                arias_duration.ad_from_acc(file_in, file_out)

            # Plot seismograms with arias duration
            filename2 = os.path.join(a_outdir, "%d.%s.vel.bbp" %
                                     (sim_id, stat))
            outfile = os.path.join(a_outdir, "%s_%d_%s_overlay.png" %
                                   (self.comp_label,
                                    sim_id, stat))
            obs_arias_n = os.path.join(a_tmpdir, "obs.%s_N.arias" %
                                       (stat))
            obs_arias_e = os.path.join(a_tmpdir, "obs.%s_E.arias" %
                                       (stat))
            obs_arias_z = os.path.join(a_tmpdir, "obs.%s_Z.arias" %
                                       (stat))
            calc_arias_n = os.path.join(a_tmpdir, "%d.%s_N.arias" %
                                        (sim_id, stat))
            calc_arias_e = os.path.join(a_tmpdir, "%d.%s_E.arias" %
                                        (sim_id, stat))
            calc_arias_z = os.path.join(a_tmpdir, "%d.%s_Z.arias" %
                                        (sim_id, stat))

            plot_seismograms.plot_overlay_with_arias(stat, filename1,
                                                     filename2,
                                                     obs_arias_n,
                                                     obs_arias_e,
                                                     obs_arias_z,
                                                     calc_arias_n,
                                                     calc_arias_e,
                                                     calc_arias_z,
                                                     self.comp_label,
                                                     "run %d" % sim_id,
                                                     outfile)

        # Now create rd50 comparison plots
        for site in site_list:
            stat = site.scode
            print("==> Plotting RotD50 comparison for station: %s" % (stat))

            # Now process rd50 files
            expected_rd50_file = os.path.join(a_outdir, "%d.%s.rd50" %
                                              (sim_id, stat))
            if not os.path.exists(expected_rd50_file):
                # just skip it
                print("Skipping rotd50/psa5 for station %s..." % (stat))
                continue

            # See if .rd50 file exists for comparison. If it don't
            # exist, skip it
            rd50_file = None
            if ("%s.rd50" % (stat)) in filelist:
                rd50_file = "%s.rd50" % (stat)
            else:
                # Skip this station
                continue

            # Plot rotd50 results
            rd50_filename1 = os.path.join(a_tmpdir_seis, rd50_file)
            rd50_filename2 = os.path.join(a_outdir, "%d.%s.rd50" %
                                          (sim_id, stat))
            outfile = os.path.join(a_outdir, "%s_%d_%s_rotd50.png" %
                                   (self.comp_label, sim_id, stat))

            plot_rotd50.plot_rd50(stat, rd50_filename1, rd50_filename2,
                                  self.comp_label, sim_id, outfile,
                                  site.low_freq_corner,
                                  site.high_freq_corner,
                                  quiet=True)

        print("Generating Plots Completed".center(80, '-'))

if __name__ == "__main__":
    PROG_BASE = os.path.basename(sys.argv[0])
    if len(sys.argv) != 6:
        print("Usage: %s " % (PROG_BASE) +
              "station_list obs_dir obs_format comp_label sim_id")
        sys.exit(1)
    print("Testing Module: %s" % (PROG_BASE))
    ME = GenPlots(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4],
                  sim_id=int(sys.argv[5]))
    ME.run()
    sys.exit(0)
