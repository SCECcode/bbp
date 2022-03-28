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
"""
from __future__ import division, print_function

# Import Python modules
import os

# Import Broadband modules
from install_cfg import InstallCfg
from station_list import StationList

class LFSeismograms(object):
    """
    This module copies pre-computed low-frequency seismograms to the
    tmpdir directory
    """

    def __init__(self, i_seis_dir, i_r_stations, sim_id=0):
        """
        Initialize class variables
        """
        self.seis_dir = i_seis_dir
        self.r_stations = i_r_stations
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

        print(self.seis_dir)

        slo = StationList(a_stations)
        stat_list = slo.getStationList()
        for stat in stat_list:
            # Look for bbp seismogram, copy in
            print("%s/%s-lf.bbp" % (self.seis_dir, stat.scode))
            if os.path.exists("%s/%s-lf.bbp" % (self.seis_dir, stat.scode)):
                print("Copying for site %s" % (stat.scode))
                # Need to eliminate negative times
                fp_in = open("%s/%s-lf.bbp" % (self.seis_dir, stat.scode), 'r')
                fp_out = open("%s/%d.%s-lf.bbp" %
                              (a_tmpdir, sim_id, stat.scode), 'w')
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
            else:
                print("Could not find LF seismogram for station %s!" %
                      (stat.scode))
