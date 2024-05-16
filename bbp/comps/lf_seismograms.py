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
import glob
import shutil

# Import Broadband modules
from install_cfg import InstallCfg
from station_list import StationList

def find_bbp_file(input_dir, station_name, label=None):
    """
    Looks into input_dir for a seismogram for station station_name
    """
    # Find input file, try with label if specified
    if label is not None:
        input_list = glob.glob("%s%s*%s.%s*.bbp" %
                               (input_dir, os.sep, label, station_name))
    else:
        input_list = []

    if not len(input_list):
        # Try to match filename without the label
        input_list = glob.glob("%s%s*%s*.bbp" %
                               (input_dir, os.sep, station_name))
        if not len(input_list):
            print("[ERROR]: Can't find input file for station %s" % (station_name))
            sys.exit(1)
    if len(input_list) > 1:
        # Found more than one file
        print("[ERROR]: Found multiple input files for station %s" % (station_name))
        sys.exit(1)

    input_file = os.path.basename(input_list[0])

    return input_file

class LFSeismograms(object):
    """
    This module copies pre-computed low-frequency seismograms to the
    tmpdir directory
    """

    def __init__(self, i_seis_dir, i_r_stations,
                 allow_negative_timestamps=False, sim_id=0):
        """
        Initialize class variables
        """
        self.seis_dir = i_seis_dir
        self.r_stations = i_r_stations
        self.allow_negative_timestamps = allow_negative_timestamps
        self.sim_id = sim_id

    def run(self):
        """
        Goes through the station list and copy each low-frequency
        seismogram from the seis_dir to the simulation's tmpdir
        """
        install = InstallCfg.getInstance()
        sim_id = self.sim_id

        a_tmpdir = os.path.join(install.A_TMP_DATA_DIR, str(sim_id))
        a_indir = os.path.join(install.A_IN_DATA_DIR, str(sim_id))
        a_stations = os.path.join(a_indir, self.r_stations)

        print("Looking for seismograms in: %s" % (self.seis_dir))

        slo = StationList(a_stations)
        station_list = slo.get_station_list()
        for station in station_list:
            station_name = station.scode
            # Look for bbp seismogram, make a copy in the tmpdata folder
            input_file = find_bbp_file(self.seis_dir, station_name)
            a_input_file = os.path.join(self.seis_dir, input_file)
            a_output_file = os.path.join(a_tmpdir,
                                         "%d.%s-lf.bbp" % (int(self.sim_id),
                                                               station_name))
            print("Copying LF file for station %s..." % (station_name))

            if self.allow_negative_timestamps:
                shutil.copy2(a_input_file, a_output_file)
            else:
                # Copy LF seismograms but remove points with t<0
                fp_in = open(a_input_file, 'r')
                fp_out = open(a_output_file, 'w')
                for line in fp_in:
                    pieces = line.split()
                    try:
                        if pieces[0] == '#' or pieces[0] == '%':
                            fp_out.write(line)
                        elif float(pieces[0]) < -0.0001:
                            continue
                        elif float(pieces[0]) < 0.0001:
                            fp_out.write("0.0\t%s\t%s\t%s\n" % (pieces[1],
                                                                pieces[2],
                                                                pieces[3]))
                        else:
                            fp_out.write(line)
                    except ValueError:
                        fp_out.write(line)
                fp_in.close()
                fp_out.flush()
                fp_out.close()
