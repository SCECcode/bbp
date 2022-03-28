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

This serves as a dummy module - it copies in the SDSU seismograms,
and converts them to .bbp format
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys
import shutil

# Import Broadband modules
import bband_utils
from station_list import StationList
from install_cfg import InstallCfg

class SDSUSeismograms(object):
    """
    This module brings SDSU low-frequency seismograms from the bin
    file into the right place so that a high-frequency module can
    start
    """

    def __init__(self, i_r_binfile, i_r_stations, i_r_full_stations, sim_id=0):
        """
        Initialize SDSUSeismograms class

        stations = stations you want seismograms for
        full_stations = stations in the binary file
        This is needed to support subsets for validation runs
        """
        self.r_binfile = i_r_binfile
        self.r_stations = i_r_stations
        if i_r_full_stations is None or i_r_full_stations == "":
            self.r_full_stations = self.r_stations
        else:
            self.r_full_stations = i_r_full_stations
        self.sim_id = sim_id

    def run(self):
        """
        Extracts needed seismograms from the bin file
        """
        print("SDSU Seismograms".center(80, '-'))

        install = InstallCfg.getInstance()
        sim_id = self.sim_id
        sta_base = os.path.basename(os.path.splitext(self.r_stations)[0])
        a_indir = os.path.join(install.A_IN_DATA_DIR, str(sim_id))
        a_outdir = os.path.join(install.A_OUT_DATA_DIR, str(sim_id))
        a_tmpdir = os.path.join(install.A_TMP_DATA_DIR, str(sim_id))
        a_tmpdir_mod = os.path.join(install.A_TMP_DATA_DIR, str(sim_id),
                                    "sdsu_seismograms_%s" % (sta_base))

        binfile = os.path.join(a_indir, self.r_binfile)

        #
        # Make sure the output and tmp directories exist
        #
        bband_utils.mkdirs([a_tmpdir, a_tmpdir_mod, a_outdir],
                           print_cmd=False)

        a_full_stations = os.path.join(a_indir, self.r_full_stations)
        a_stations = os.path.join(a_indir, self.r_stations)

        # Copy station files to the tmpdir_mod directory
        cmd = "cp %s %s" % (a_full_stations,
                            os.path.join(a_tmpdir_mod, self.r_full_stations))
        bband_utils.runprog(cmd)
        cmd = "cp %s %s" % (a_stations,
                            os.path.join(a_tmpdir_mod, self.r_stations))
        bband_utils.runprog(cmd)

        #
        # Make sure path names are within the limits accepted by the
        # Fortran code
        #
        if len(binfile) >= bband_utils.SDSU_MAX_FILENAME:
            raise ValueError("binfile is %d characters long, maximum is %d" %
                             (len(binfile), bband_utils.SDSU_MAX_FILENAME))

        old_cwd = os.getcwd()
        os.chdir(a_tmpdir_mod)

        # Get number of stations in seismogram file, this is a
        # variation of the code in station_list.py
        stat_names = {}
        num_stations = 0
        stat_fp = open(a_full_stations, 'r')
        for line in stat_fp:
            if line.startswith('#'):
                continue
            sta = line.split()
            if len(sta) >= 3:
                scode = sta[2]
                num_stations = num_stations + 1
                stat_names[scode] = num_stations
        stat_fp.close()

        # Create list of stations to save
        slo = StationList(a_stations)
        site_list = slo.getStationList()
        save_stat_names = []
        for stat in site_list:
            save_stat_names.append(stat.scode)

        # Convert to bbp format
        cmd = "%s/bin2bbp %s %d" % (install.A_SDSU_BIN_DIR,
                                    binfile, len(stat_names))
        bband_utils.runprog(cmd)

        # Copy over the names
        for stat in save_stat_names:
            if not stat in stat_names:
                continue
            sta_id = stat_names[stat]
            shutil.copy2("%s/%d.bbp" % (a_tmpdir_mod, sta_id),
                         "%s/%d.%s-lf.bbp" %
                         (a_tmpdir, sim_id, stat))
            del stat_names[stat]

        # Delete the ones you don't need
        for stat in stat_names.keys():
            os.remove("%s/%d.bbp" %
                      (a_tmpdir_mod, stat_names[stat]))

        os.chdir(old_cwd)

        print("SDSU Seismograms Completed".center(80, '-'))

if __name__ == "__main__":
    print("Testing Module: %s" % (os.path.basename(sys.argv[0])))
    ME = SDSUSeismograms(sys.argv[1], sys.argv[2], sys.argv[3],
                         sim_id=int(sys.argv[4]))
    ME.run()
