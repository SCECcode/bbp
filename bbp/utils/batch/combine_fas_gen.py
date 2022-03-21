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

Script to combine FAS data from multiple realizations and create
a combined FAS GOF plot.
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
from fas_gof_cfg import FASGofCfg, resid2uncer_py
from station_list import StationList

# Import Pynga and its utilities
import pynga.utils as putils

# --------------------------------------------------------------------------
# Functions
# --------------------------------------------------------------------------

def summarize_fas(tmpdir, outdir, fas_residfile,
                  comp_label, num_stations, num_realization,
                  codebase):
    """
    This function summarizes all fas data and creates the combined fas
    GOF plot

    """
    config = FASGofCfg()
    min_cdst = 0
    max_cdst = 1e+15

    method = codebase.lower()
    freq_ranges = plot_config.FAS_GOF_FREQ
    lfreq=freq_ranges[method]['freq_low']
    hfreq=freq_ranges[method]['freq_high']

    logfile = os.path.join(tmpdir, "log.txt")

    for comp in config.COMPS_FAS:
        # Build paths and check lengths
        fileroot = os.path.join(tmpdir, "%s-%s-combined-fas-%s" %
                                (codebase, comp_label, comp))
        bband_utils.check_path_lengths([fas_residfile, fileroot],
                                       bband_utils.GP_MAX_FILENAME)
        resid2uncer_py(fas_residfile, fileroot, comp,
                       num_stations, min_cdst, max_cdst)

    plottitle = ("Combined GOF Plot for %s\n%d Realizations - %s Method" %
                 (comp_label, num_realization, codebase.upper()))
    fileroot = "%s-%s-combined-fas" % (codebase, comp_label)
    plotter = PlotGoF()
    plotter.plot_fas_gof(plottitle, fileroot, tmpdir, outdir,
                         cutoff=0, lfreq=lfreq, hfreq=hfreq,
                         colorset="combined")

    print("Stations used: %s" % (num_stations))

def read_fas_resid_file(resid_file, resid_data):
    """
    This function reads the resid file and copies all
    data we will need for averaging to resid_data
    """
    input_file = open(resid_file, 'r')
    header = input_file.readline()
    header = header.strip()
    pieces = header.split()
    frequencies = pieces[13:]
    for line in input_file:
        line = line.strip()
        pieces = line.split()
        station = pieces[2]
        comp = pieces[12]

        if station not in resid_data:
            resid_data[station] = {}
        if comp not in resid_data[station]:
            resid_data[station][comp] = {}
            resid_data[station][comp]['header'] = pieces[0:12]
            resid_data[station][comp]['frequencies'] = {}
        for freq, new_data in zip(frequencies, pieces[13:]):
            if freq not in resid_data[station][comp]['frequencies']:
                resid_data[station][comp]['frequencies'][freq] = []
            resid_data[station][comp]['frequencies'][freq].append(float(new_data))
    input_file.close()

    return header, frequencies

def combine_resid_data(outdir, combined_file):
    """
    This function combines the residuals from multiple
    realizations by averaging the individual residuals
    from each realization period by frequency
    """
    # Get realizations
    realizations = sorted(os.listdir(outdir))
    one_realization = realizations[0]
    basedir = os.path.join(outdir, one_realization)
    # Figure out what our stations are
    fas_files = glob.glob("%s%sFAS%sobs*.col" % (basedir,
                                                 os.sep,
                                                 os.sep))
    fas_files = [os.path.basename(each_file) for each_file in fas_files]
    stations = [station.split(".")[1] for station in fas_files]
    # Capture event_label
    bias_file = glob.glob("%s%s*.bias" % (basedir, os.sep))
    if len(bias_file) < 1:
        raise bband_utils.ProcessingError("Cannot find event label!")
    bias_file = bias_file[0]
    # Let's capture the event label
    event_label = os.path.basename(bias_file).split("-")[0]

    resid_header = None
    resid_frequencies = None
    resid_data = {}
    # Now walk through each realizations and collect what we need
    for cur_realization in realizations:
        cur_dir = os.path.join(outdir, cur_realization)
        resid_file = glob.glob("%s%s*fas-resid.txt" % (cur_dir, os.sep))
        if len(resid_file) < 1:
            raise bband_utils.ProcessingError("Cannot find residuals file!")
        resid_file = resid_file[0]
        header, frequencies = read_fas_resid_file(resid_file, resid_data)
        if resid_header is None:
            resid_header = header
        if resid_frequencies is None:
            resid_frequencies = frequencies

    # Got all the data we need, now write the combined residuals file
    output_file = open(combined_file, 'w')
    output_file.write("%s\n" % (header))
    for station in resid_data:
        for comp in resid_data[station]:
            comp_header = resid_data[station][comp]["header"]
            comp_frequencies = resid_data[station][comp]["frequencies"]
            output_file.write("\t".join(comp_header))
            output_file.write("\t%s" % (comp))
            for freq in resid_frequencies:
                per_data = comp_frequencies[freq]
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
(OPTIONS, ARGS) = PARSER.parse_args()

if OPTIONS.input_dir is None:
    PARSER.error("Please specify the input directory!")
TOP_INPUT_DIR = OPTIONS.input_dir
if not os.path.isdir(TOP_INPUT_DIR):
    PARSER.error("Invalid input directory!")
if not "Sims" in os.listdir(TOP_INPUT_DIR):
    PARSER.error("Please provide the top-level simulation directory!\n"
                 "This is the directory given to the cluster script")
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
                             "bbp-fas-resid-combined.txt")

# Combine data from multiple realizations
(COMP_LABEL,
 NUM_REALIZATIONS,
 NUM_STAT) = combine_resid_data(INPUT_OUTDIR,
                                COMBINED_FILE)

summarize_fas(TMPDIR, OUTPUT_DIR,
              COMBINED_FILE,
              COMP_LABEL,
              NUM_STAT,
              NUM_REALIZATIONS,
              OPTIONS.codebase)

print("All Done!")
# Clean-up, all done!
shutil.rmtree(TMPDIR)
