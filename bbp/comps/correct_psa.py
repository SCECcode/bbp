#!/usr/bin/env python
"""
Copyright 2010-2018 University Of Southern California

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

This module applies the user-supplied correction factors to the
outputs of RotD50/RotD100
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys

# Import Broadband modules
import seqnum
import bband_utils
from install_cfg import InstallCfg
from station_list import StationList

class CorrectPSA(object):
    """
    This class corrects the output of RotD50/RotD100 using
    user-supplied correction factors
    """

    def __init__(self, i_r_stations, i_extension, i_corr_file,
                 i_a_proc_dir, sim_id=0):
        """
        Initialize basic class parameters
        """
        self.sim_id = sim_id
        self.r_stations = i_r_stations
        self.extension = i_extension
        self.corr_file = i_corr_file
        self.proc_dir = i_a_proc_dir
        # Start empty for now
        self.periods = []
        self.factors = {}
        # Load the correction factors from the file
        self.load_correction_factors()

    def load_correction_factors(self):
        """
        This function loads the correction factors from the corr_file
        provided
        """
        try:
            cfile = open(self.corr_file, 'r')
        except IOError:
            raise bband_utils.ParameterError("Cannot read correction file %s" %
                                             (self.corr_file))

        # We are looking for the header first
        headers = None
        # Loop through the lines
        for line in cfile:
            if line.startswith("#StaName"):
                headers = line.split()
                break

        # Make sure we got the header line
        if headers is None:
            cfile.close()
            raise bband_utils.ProcessingError("Cannot find header line in "
                                              "correction file %s" %
                                              (self.corr_file))

        skip_headers = 0
        # Pick up the periods
        while len(headers) > 0:
            try:
                tmp = float(headers[0])
            except:
                # Skip this one, and remove from list
                skip_headers = skip_headers + 1
                headers.pop(0)
            else:
                # Found first period, get out
                break
        # Make sure we have at least 1 period
        if not headers:
            cfile.close()
            raise bband_utils.ProcessingError("Cannot find any periods in "
                                              "correction file %s" %
                                              (self.corr_file))

        # Convert periods to floats
        self.periods = [float(value) for value in headers]

        # Now read the rest of the correction file
        for line in cfile:
            if line.startswith("#"):
                continue
            factors = line.split()
            station = factors[0]
            to_skip = skip_headers
            # Remove everything other than the correction factors
            while to_skip > 0:
                factors.pop(0)
                to_skip = to_skip - 1

            factors = [float(value) for value in factors]
            # Make sure we have the proper number of correction factors
            if len(factors) != len(self.periods):
                cfile.close()
                raise bband_utils.ProcessingError("Station %s has %d periods" %
                                                  (station, len(factors)) +
                                                  ", expecting %s periods" %
                                                  (len(self.periods)))

            self.factors[station] = factors

        # All done
        cfile.close()

    def correct_file(self, factors, input_file, output_file):
        """
        This function reads input_file and writes output_file after
        applying the correction factors
        """

        # Open files
        ifile = open(input_file, 'r')
        ofile = open(output_file, 'w')

        for line in ifile:
            if line.startswith("#"):
                # Comment line, just write to output
                ofile.write(line)
            else:
                # This line needs correction
                pieces = line.split()
                period = float(pieces[0])
                pieces.pop(0)
                values = [float(value) for value in pieces]
                try:
                    index = self.periods.index(period)
                except ValueError:
                    raise bband_utils.ProcessingError("Period %f not on the " %
                                                      (period) +
                                                      "correction list!")
                # Apply correction value
                values = [value * factors[index] for value in values]
                ostr = " %10.5E" % (period)
                for value in values:
                    ostr = "%s %10.5E" % (ostr, value)
                # Write to output file
                ofile.write("%s\n" % (ostr))

        # All done!
        ifile.close()
        ofile.close()

    def correct_station(self, station, extension):
        """
        This function applies the user-provided correction factors to
        the station's amplitudes, and outputs the corrected file
        """
        if not station in self.factors:
            raise bband_utils.ParameterError("Unknown station %s!" % (station))

        orig_file = os.path.join(self.proc_dir, "%s-orig.%s" % (station, extension))
        corr_file = os.path.join(self.proc_dir, "%s.%s" % (station, extension))

        # Make sure input files exist
        if not os.path.exists(orig_file):
            raise bband_utils.ProcessingError("File %s not found!" %
                                              (orig_file))

        # Pick set of correction factors
        factors = self.factors[station]

        # Correct rd50 file
        self.correct_file(factors, orig_file, corr_file)

    def run(self):
        """
        Corrects the amplitudes from all stations found in the station
        list according to the correction coefficients provided by the user
        """
        print("Correct PSA".center(80, '-'))

        # Initialize basic variables
        install = InstallCfg.getInstance()
        sim_id = self.sim_id
        sta_base = os.path.basename(os.path.splitext(self.r_stations)[0])

        self.log = os.path.join(install.A_OUT_LOG_DIR,
                                str(sim_id),
                                "%d.obs_seis.log" %
                                (sim_id))

        # Input, tmp, and output directories
        a_indir = os.path.join(install.A_IN_DATA_DIR, str(sim_id))

        # Station file
        a_statfile = os.path.join(a_indir, self.r_stations)

        slo = StationList(a_statfile)
        site_list = slo.getStationList()

        # Go through each station
        # print "Working dir: %s" % (self.proc_dir)
        for site in site_list:
            stat = site.scode
            print("==> Correcting amplitudes for station: %s" % (stat))
            self.correct_station(stat, self.extension)

        print("Correct PSA Completed".center(80, '-'))

if __name__ == "__main__":
    print("Testing module: %s" % os.path.basename((sys.argv[0])))
    if len(sys.argv) < 5:
        print("usage: %s " % (sys.argv[0]) +
              "station_file corr_file extension"
              " proc_dir [sim_id]")
        sys.exit(1)
    if len(sys.argv) == 5:
        SIM_ID = int(seqnum.get_seq_num())
    else:
        SIM_ID = int(sys.argv[5])
    CORR_PSA = CorrectPSA(sys.argv[1], sys.argv[2], sys.argv[3],
                          sys.argv[4], SIM_ID)
    CORR_PSA.run()
