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
import sys
import glob
import numpy
import shutil
import optparse
import tempfile

# Import Broadband modules
import bband_utils
from station_list import StationList

# Import Pynga and its utilities
import pynga.utils as putils

# --------------------------------------------------------------------------
# Functions
# --------------------------------------------------------------------------

def read_resid_file(resid_file, resid_data):
    """
    This function reads the resid file and copies all
    data we will need for averaging to resid_data
    """
    input_file = open(resid_file, 'r')
    header = input_file.readline()
    header = header.strip()
    pieces = header.split()
    periods = pieces[13:]
    for line in input_file:
        line = line.strip()
        pieces = line.split()
        station = pieces[2]
        comp = pieces[12]

        #print("%s %s %s" % (station, comp, " ".join(pieces[13:])))

        if station not in resid_data:
            resid_data[station] = {}
        if comp not in resid_data[station]:
            resid_data[station][comp] = {}
            resid_data[station][comp]['header'] = pieces[0:12]
            resid_data[station][comp]['periods'] = {}
        for period, new_data in zip(periods, pieces[13:]):
            if period not in resid_data[station][comp]['periods']:
                resid_data[station][comp]['periods'][period] = []
            resid_data[station][comp]['periods'][period].append(float(new_data))
    input_file.close()

    return header, periods

def combine_realizations_data(input_dir, temp_dir,
                              input_indir, combined_file):
    """
    This function averages the residuals across all realizations
    """
    # Get realizations
    realizations = sorted(os.listdir(input_dir))
    one_realization = realizations[0]
    basedir = os.path.join(input_dir, one_realization)
    # Figure out what our stations are
    rd50_files = glob.glob("%s%s%s.*.rd50" % (basedir,
                                              os.sep,
                                              one_realization))
    rd50_files = [os.path.basename(each_file) for each_file in rd50_files]
    stations = [station.split(".")[1] for station in rd50_files]
    # Capture event_label
    bias_file = glob.glob("%s%s*.bias" % (basedir, os.sep))
    if len(bias_file) < 1:
        raise bband_utils.ProcessingError("Cannot find event label!")
    bias_file = bias_file[0]
    # Let's capture the event label
    event_label = os.path.basename(bias_file).split("-")[0]

    # Get source file
    basedir =  os.path.join(input_indir, one_realization)
    a_srcfile = glob.glob("%s%s*.src" % (basedir, os.sep))
    if len(a_srcfile) != 1:
        a_srcfile = select_src_file(a_srcfile)
    else:
        a_srcfile = a_srcfile[0]

    # Parse it!
    src_keys = bband_utils.parse_src_file(a_srcfile)
    
    resid_header = None
    resid_periods = None
    resid_data = {}
    # Now walk through each realization and collect what we need
    for cur_realization in realizations:
        cur_dir = os.path.join(input_dir, cur_realization)
        resid_file = glob.glob("%s%s*rd50-resid.txt" % (cur_dir, os.sep))
        if len(resid_file) < 1:
            raise bband_utils.ProcessingError("Cannot find residuals file!")
        resid_file = resid_file[0]
        header, periods = read_resid_file(resid_file, resid_data)
        if resid_header is None:
            resid_header = header
        if resid_periods is None:
            resid_periods = periods

    # Got all the data we need, now write the combined residuals file

    # Copy header for first file, set logfile
    if os.path.isfile(combined_file):
        # But not, if file already exists
        copy_header = 0
    else:
        copy_header = 1

    output_file = open(combined_file, 'a')
    if copy_header:
        output_file.write("%s\n" % (header))
    for station in resid_data:
        for comp in resid_data[station]:
            comp_header = resid_data[station][comp]["header"]
            comp_periods = resid_data[station][comp]["periods"]
            slon = float(comp_header[3])
            slat = float(comp_header[4])

            # Calculate Rrup
            origin = (src_keys['lon_top_center'],
                      src_keys['lat_top_center'])
            dims = (src_keys['fault_length'], src_keys['dlen'],
                    src_keys['fault_width'], src_keys['dwid'],
                    src_keys['depth_to_top'])
            mech = (src_keys['strike'], src_keys['dip'],
                    src_keys['rake'])

            site_geom = [float(slon), float(slat), 0.0]
            (fault_trace1, up_seis_depth,
             low_seis_depth, ave_dip,
             dummy1, dummy2) = putils.FaultTraceGen(origin, dims, mech)
            _, rrup, _ = putils.DistanceToSimpleFaultSurface(site_geom,
                                                             fault_trace1,
                                                             up_seis_depth,
                                                             low_seis_depth,
                                                             ave_dip)
            comp_header[7] = "%.2f" % (rrup)
            output_file.write("\t".join(comp_header))
            output_file.write("\t%s" % (comp))
            for period in resid_periods:
                per_data = comp_periods[period]
                output_file.write("\t%.5e" % (numpy.mean(per_data)))
            output_file.write("\n")
    output_file.close()

    return event_label, len(realizations), len(stations)

def select_src_file(a_srcfiles):
    """
    Returns the src file that contains the hypocenter
    """
    for a_srcfile in a_srcfiles:
        src_keys = bband_utils.parse_src_file(a_srcfile)
        if 'true_hypo' in src_keys:
            # If event has true_hypo key and it is set to 1, use this segment
            if int(float(src_keys['true_hypo'])) == 1:
                return a_srcfile
        fault_len = float(src_keys['fault_length']) / 2.0
        hypo_along_strike = abs(float(src_keys['hypo_along_stk']))
        if hypo_along_strike <= fault_len:
            # Use this segment
            return a_srcfile

    raise bband_utils.ProcessingError("Cannot parse src files!")

#
# Main
# --------------------------------------------------------------------------
def main():
    """
    Get output file and list of directories from where we should read
    data
    """
    if len(sys.argv) < 3:
        print("Usage: %s output_file dir1 [dir2 dir3... dirn]" % (sys.argv[0]))
        sys.exit(1)

    output_file = sys.argv[1]
    dirs = sys.argv[2:]
    for input_dir in dirs:
        if not os.path.isdir(input_dir):
            print("Skipping directory: %s!" % (input_dir))
            continue
        if not "Sims" in os.listdir(input_dir):
            print("Skipping invalid directory: %s!" % (input_dir))
            continue
        print("Processing %s..." % (input_dir))
        input_outdir = os.path.join(input_dir, "Sims" , "outdata")
        input_tmpdir = os.path.join(input_dir, "Sims" , "tmpdata")
        input_indir = os.path.join(input_dir, "Sims" , "indata")

        # Create temp dir
        tmpdir = tempfile.mkdtemp(prefix="bbp-")

        # Combine realizations' data
        comp_label, _, _ = combine_realizations_data(input_outdir, tmpdir,
                                                     input_indir, output_file)

        # Clean-up, all done!
        shutil.rmtree(tmpdir)

if __name__ == "__main__":
    main()
