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

def combine_station_data(station, input_dir, temp_dir):
    """
    This function combines data for a given station across multiple
    realizations, writting a single output file in temp_dir
    """
    data = {}
    # Get realizations
    realizations = sorted(os.listdir(input_dir))
    for realization in realizations:
        basedir = os.path.join(input_dir, realization)
        data_file = glob.glob("%s%s%s.%s.rd50" % (basedir, os.sep,
                                                  realization,
                                                  station))
        if len(data_file) != 1:
            raise bband_utils.ProcessingError("Data for station %s " %
                                              (station) +
                                              "not found for "
                                              "realization %s!" %
                                              (realization))
        data_file = data_file[0]
        in_data = open(data_file, 'r')
        for line in in_data:
            line = line.strip()
            # Skip comments
            if line.startswith("#"):
                continue
            pieces = line.split()
            pieces = [float(piece) for piece in pieces]
            key = pieces[0]
            pieces = pieces[1:]
            if not key in data:
                # Key is new to dictionary
                empty_set = [[] for _ in pieces]
                data[key] = empty_set
            for idx, value in enumerate(pieces):
                data[key][idx].append(value)
        in_data.close()

    # Now, write the output file
    out_file = open((os.path.join(temp_dir, "%s.rd50" % (station))), 'w')
    keys = sorted(data.keys())
    for key in keys:
        out_file.write("%10.4f" % (key))
        for comp in data[key]:
            out_file.write(" %10.5e" % (numpy.mean(comp)))
        out_file.write("\n")

def combine_realizations_data(input_dir, temp_dir):
    """
    This function creates a single file averaging the rd50 files for
    each of the stations across all realizations
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

    # Now walk through all realizations and combine stations data
    for station in stations:
        print("working on station: %s" % (station))
        combine_station_data(station, input_dir, temp_dir)

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

def create_resid_data_file(comp_label, input_indir, input_obsdir,
                           combined_file, temp_dir):
    """
    This function creates a file containing the combined residuals
    from the simulation data from all stations
    """
    # Copy header for first file, set logfile
    if os.path.isfile(combined_file):
        # But not, if file already exists
        copy_header = 0
    else:
        copy_header = 1
    logfile = os.path.join(temp_dir, "log.txt")

    # Figure out where out binaries are
    if "BBP_DIR" in os.environ:
        install_root = os.path.normpath(os.environ["BBP_DIR"])
    else:
        raise bband_utils.ProcessingError("BBP_DIR is not set!")
    gp_bin_dir = os.path.join(install_root, "src", "gp", "bin")

    # Get realizations
    realizations = sorted(os.listdir(input_indir))
    one_realization = realizations[0]
    basedir = os.path.join(input_indir, one_realization)

    # Get the station list
    a_statfile = glob.glob("%s%s*.stl" % (basedir, os.sep))
    if len(a_statfile) != 1:
        raise bband_utils.ProcessingError("Cannot get station list!")
    a_statfile = a_statfile[0]
    slo = StationList(a_statfile)
    site_list = slo.get_station_list()

    # Get source file
    a_srcfile = glob.glob("%s%s*.src" % (basedir, os.sep))
    if len(a_srcfile) != 1:
        a_srcfile = select_src_file(a_srcfile)
    else:
        a_srcfile = a_srcfile[0]

    # Parse it!
    src_keys = bband_utils.parse_src_file(a_srcfile)

    # Get the obsdir
    print(input_obsdir)
    realizations = sorted(os.listdir(input_obsdir))
    one_realization = realizations[0]
    basedir = os.path.join(input_obsdir, one_realization)
    obs_dir = glob.glob("%s%sobs_seis*" % (basedir, os.sep))
    if len(obs_dir) != 1:
        raise bband_utils.ProcessingError("Cannot get observation dir!")
    obs_dir = obs_dir[0]

    # Go through all stations
    for site in site_list:
        slon = float(site.lon)
        slat = float(site.lat)
        stat = site.scode

        # Calculate Rrup
        origin = (src_keys['lon_top_center'],
                  src_keys['lat_top_center'])
        dims = (src_keys['fault_length'], src_keys['dlen'],
                src_keys['fault_width'], src_keys['dwid'],
                src_keys['depth_to_top'])
        mech = (src_keys['strike'], src_keys['dip'],
                src_keys['rake'])

        site_geom = [float(site.lon), float(site.lat), 0.0]
        (fault_trace1, up_seis_depth,
         low_seis_depth, ave_dip,
         dummy1, dummy2) = putils.FaultTraceGen(origin, dims, mech)
        _, rrup, _ = putils.DistanceToSimpleFaultSurface(site_geom,
                                                         fault_trace1,
                                                         up_seis_depth,
                                                         low_seis_depth,
                                                         ave_dip)

        simfile1 = os.path.join(temp_dir, "%s.rd50" % (stat))
        datafile1 = os.path.join(obs_dir, "%s.rd50" % (stat))

        cmd = ("%s/gen_resid_tbl_3comp bbp_format=1 " % (gp_bin_dir) +
               "datafile1=%s simfile1=%s " % (datafile1, simfile1) +
               "comp1=psa5n comp2=psa5e comp3=rotd50 " +
               "eqname=%s mag=0.0 stat=%s lon=%.4f lat=%.4f " %
               (comp_label, stat, slon, slat) +
               "vs30=%d cd=%.2f " % (site.vs30, rrup) +
               "flo=%f fhi=%f " % (site.low_freq_corner,
                                   site.high_freq_corner) +
               "print_header=%d >> %s 2>> %s" %
               (copy_header, combined_file, logfile))
        bband_utils.runprog(cmd, abort_on_error=True)

        if copy_header == 1:
            copy_header = 0

# --------------------------------------------------------------------------
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
        comp_label, _, _ = combine_realizations_data(input_outdir, tmpdir)

        # Create data files with both gmpe and simulation data
        create_resid_data_file(comp_label, input_indir, input_outdir,
                               output_file, tmpdir)

        # Clean-up, all done!
        shutil.rmtree(tmpdir)

if __name__ == "__main__":
    main()
