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

Broadband Platform script to generate GMPE station comparisons and GOF plots
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys
# import shutil

# Import Broadband modules
import plot_gmpe
import bband_utils
import gmpe_config
from install_cfg import InstallCfg
from station_list import StationList

class GMPEPlot(object):

    def __init__(self, i_r_stations, i_format,
                 i_comparison_label, i_gmpe_group_name,
                 sim_id=0):
        self.sim_id = sim_id
        self.r_stations = i_r_stations
        self.format = i_format
        self.comp_label = i_comparison_label
        self.gmpe_group_name = i_gmpe_group_name
        # Make sure gmpes are in the right format
        if i_format != "gmpe":
            raise bband_utils.ParameterError("Format %s for " %
                                             (self.format) +
                                             "gmpe results "
                                             "not supported")

    def run(self):
        """
        This function creates GMPE plots for all stations
        """
        print("GMPE Plot".center(80, '-'))

        # Initialize basic variables
        install = InstallCfg.getInstance()
        sim_id = self.sim_id
        sta_base = os.path.basename(os.path.splitext(self.r_stations)[0])

        self.log = os.path.join(install.A_OUT_LOG_DIR, str(sim_id),
                                "%d.gmpe_gof.log" % (sim_id))

        # Input, tmp, and output directories
        a_tmpdir = os.path.join(install.A_TMP_DATA_DIR, str(sim_id))
        a_outdir = os.path.join(install.A_OUT_DATA_DIR, str(sim_id))
        a_outdir_gmpe = os.path.join(install.A_OUT_DATA_DIR, str(sim_id),
                                     "gmpe_data_%s" % (sta_base))

        #
        # Make sure the output and tmp directories exist
        #
        dirs = [a_tmpdir, a_outdir, a_outdir_gmpe]
        bband_utils.mkdirs(dirs, print_cmd=False)

        # Figure out gmpe labels
        gmpe_group = gmpe_config.GMPES[self.gmpe_group_name]
        gmpe_labels = gmpe_group["labels"]

        # Station file
        a_statfile = os.path.join(install.A_IN_DATA_DIR,
                                  str(sim_id),
                                  self.r_stations)
        # List of gmpe files
        filelist = os.listdir(a_outdir_gmpe)

        slo = StationList(a_statfile)
        site_list = slo.get_station_list()

        # Go through each station
        for site in site_list:
            stat = site.scode
            print("==> Generating plot for station: %s" % (stat))
            # Since we're using the GP station list, make sure the
            # .rd50 for the station exists.  It might not if we ran the
            # validation with a shorter station list
            sim_file = os.path.join(a_outdir, "%d.%s.rd50" % (sim_id, stat))
            if not os.path.exists(sim_file):
                # just skip it
                print("Couldn't find file %s. " % (sim_file) +
                      "This is not necessarily an error, as you may have " +
                      "run with a subset of a stations. Continuing " +
                      "with available stations.")
                continue

            # Ok, we have the calculated rd50 for this station
            # Look for the gmpe file
            r_gmpe_file = "%s-gmpe.ri50" % (stat)
            if r_gmpe_file not in filelist:
                # No gmpe file for this station
                continue
            a_gmpe_file = os.path.join(a_outdir_gmpe, r_gmpe_file)

            # Plot GMPE rotd50 results
            outfile = os.path.join(a_outdir, "%s_%d_%s_gmpe.png" %
                                   (self.comp_label, sim_id, stat))

            plot_gmpe.plot_gmpe(stat, sim_file, a_gmpe_file, gmpe_labels,
                                sim_id, self.comp_label, outfile)

        print("GMPE Plot Completed".center(80, '-'))

if __name__ == "__main__":
    PROG_BASE = os.path.basename(sys.argv[0])
    if len(sys.argv) != 7:
        print("Usage: %s " % (PROG_BASE) +
              "station_list gmpe_format magnitude "
              "comp_label gmpe_group_name sim_id")
        sys.exit(1)
    print("Testing Module: %s" % (PROG_BASE))
    ME = GMPEPlot(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4],
                  sys.argv[5], sim_id=int(sys.argv[6]))
    ME.run()
    sys.exit(0)
