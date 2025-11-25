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

Python version of Ronnie Kamai Matlab scripts to generate a combined
bias plot. It collects information from the rd50 files for each
realization, groups the results by station (averaging) and then
generates a residuals file using the rd50 files from the recorded
data. This single residuals file uses the same resid2uncer_varN
program used in single bias plots to generate data for the combined
plot.
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
import plot_config
from PlotGOF import PlotGoF
from gp_gof_cfg import GPGofCfg
from station_list import StationList

# Import Pynga and its utilities
import pynga.utils as putils

# --------------------------------------------------------------------------
# Functions
# --------------------------------------------------------------------------

def summarize_rotd50(tmpdir, outdir, combined_resid_file,
                     comp_label, num_stations, num_realization,
                     codebase, output_mode):
    """
    This function summarizes all rotd50 data and creates the combined
    rotd50 GOF plot
    """
    config = GPGofCfg()
    method = codebase.lower()
    freq_ranges = plot_config.PSA_GOF_FREQ
    lfreq=freq_ranges[method]['freq_low']
    hfreq=freq_ranges[method]['freq_high']

    # Set output mode
    if output_mode == "periods":
        mode = "rd50-single"
    elif output_mode == "freq":
        mode = "rd50-single-freq"

    # Figure out where out binaries are
    if "BBP_DIR" in os.environ:
        install_root = os.path.normpath(os.environ["BBP_DIR"])
    else:
        raise bband_utils.ProcessingError("BBP_DIR is not set!")
    gp_bin_dir = os.path.join(install_root, "src", "gp", "bin")

    logfile = os.path.join(tmpdir, "log.txt")

    for comp in config.COMPS_PSA5:
        # Build paths and check lengths
        fileroot = os.path.join(tmpdir, "%s-%s-combined-rd50-%s" %
                                (codebase, comp_label, comp))
        bband_utils.check_path_lengths([combined_resid_file, fileroot],
                                       bband_utils.GP_MAX_FILENAME)

        cmd = ("%s " % (os.path.join(gp_bin_dir, "resid2uncer_varN")) +
               "residfile=%s fileroot=%s " % (combined_resid_file, fileroot) +
               "comp=%s nstat=%d nper=63 " % (comp, num_stations) +
               " >> %s 2>&1" % (logfile))
        bband_utils.runprog(cmd, abort_on_error=True)

    plottitle = ("Combined GOF Plot for %s\n%d Realizations - %s Method" %
                 (comp_label, num_realization, codebase.upper()))
    fileroot = "%s-%s-combined-rd50" % (codebase, comp_label)
    plotter = PlotGoF()
    plotter.plot(plottitle, fileroot, tmpdir, outdir,
                 cutoff=0, mode=mode, lfreq=lfreq, hfreq=hfreq,
                 colorset="combined")

    print("Stations used: %s" % (num_stations))

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

def combine_resid_data(input_dir, temp_dir, combined_file):
    """
    This function combines the residuals from multiple
    realizations by averaging the individual residuals
    from each realization period by period
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
    output_file = open(combined_file, 'w')
    output_file.write("%s\n" % (header))
    for station in resid_data:
        for comp in resid_data[station]:
            comp_header = resid_data[station][comp]["header"]
            comp_periods = resid_data[station][comp]["periods"]
            output_file.write("\t".join(comp_header))
            output_file.write("\t%s" % (comp))
            for period in resid_periods:
                per_data = comp_periods[period]
                output_file.write("\t%.5e" % (numpy.mean(per_data)))
            output_file.write("\n")
    output_file.close()

    return event_label, len(realizations), len(stations)

# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

PARSER = optparse.OptionParser()
PARSER.add_option("-d", "--dir", dest="input_dir",
                  help="Input directory containing simulation results")
PARSER.add_option("-o", "--output_dir", dest="output_dir",
                  help="Output file")
PARSER.add_option("-c", "--codebase", dest="codebase",
                  help="Method used for the simulation")
PARSER.add_option("-f", "--freq", action="store_true", dest="freq",
                  help="Use frequencies on X axis instead of periods")
(OPTIONS, ARGS) = PARSER.parse_args()


if OPTIONS.input_dir is None:
    PARSER.error("Please specify the input directory!")
TOP_INPUT_DIR = OPTIONS.input_dir
if not os.path.isdir(TOP_INPUT_DIR):
    PARSER.error("Invalid input directory!")
if not "Sims" in os.listdir(TOP_INPUT_DIR):
    PARSER.error("Please provide the top-level simulation directory!\n"
                 "This is the directory given to the cluster script")

# Select output format
OUTPUT_MODE = "periods"
if OPTIONS.freq:
    OUTPUT_MODE = "freq"

INPUT_OUTDIR = os.path.join(TOP_INPUT_DIR, "Sims" , "outdata")
INPUT_TMPDIR = os.path.join(TOP_INPUT_DIR, "Sims" , "tmpdata")
INPUT_INDIR = os.path.join(TOP_INPUT_DIR, "Sims" , "indata")

if OPTIONS.output_dir is None:
    PARSER.error("error specify output directory!")
else:
    OUTPUT_DIR = OPTIONS.output_dir
    if not os.path.isdir(OUTPUT_DIR):
        PARSER.error("Invalid output directory!")

if OPTIONS.codebase is None:
    PARSER.error("Please specify codebase!")

# Create temp dir
TMPDIR = tempfile.mkdtemp(prefix="bbp-")
COMBINED_FILE = os.path.join(TMPDIR,
                             "bbp-rd50-resid-combined.txt")

# Combine data from multiple realizations
(COMP_LABEL,
 NUM_REALIZATIONS,
 NUM_STAT) = combine_resid_data(INPUT_OUTDIR, TMPDIR,
                                COMBINED_FILE)

# Generate plots
summarize_rotd50(TMPDIR, OUTPUT_DIR,
                 COMBINED_FILE,
                 COMP_LABEL,
                 NUM_STAT,
                 NUM_REALIZATIONS,
                 OPTIONS.codebase,
                 OUTPUT_MODE)

print("All Done!")
# Clean-up, all done!
shutil.rmtree(TMPDIR)
