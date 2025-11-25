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

This program creates a Vs30-based GOF, combining information from
all realizations into a single plot.
"""
from __future__ import division, print_function

# Import Python modules
import os
import glob
import optparse
import matplotlib
if (matplotlib.get_backend() != 'agg'):
    matplotlib.use('Agg') # Disables use of Tk/X11
import pylab
import numpy

# Import Broadband modules
import plot_config
import bband_utils

# Constants
MIN_Y_AXIS = -1.75
MAX_Y_AXIS = 1.75
VS30_PERIODS = [0.01, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0]

# --------------------------------------------------------------------------
# Functions
# --------------------------------------------------------------------------

def combine_realization_data(tmpdir, period):
    """
    This function reads the resid files from all realizations and
    returns the combined data set
    """
    data = {}

    realizations = sorted(os.listdir(tmpdir))
    for realization in realizations:
        basedir = os.path.join(tmpdir, realization)
        resid_file = glob.glob("%s%s*-resid-vs30-%.3f-rotd50.txt" %
                               (basedir, os.sep, period))
        if len(resid_file) != 1:
            raise bband_utils.ProcessingError("Residuals file not found for "
                                              "realization %s!" % (realization))
        resid_file = resid_file[0]
        input_file = open(resid_file, 'r')
        for line in input_file:
            line = line.strip()
            # Skip comments and empty lines
            if line.startswith("#") or line.startswith("%") or not line:
                continue
            pieces = line.split()
            # Make sure we have enough tokens
            if len(pieces) != 2:
                continue
            # Convert to floats
            pieces = [float(piece) for piece in pieces]
            vs30 = pieces[0]
            val = pieces[1]
            if vs30 in data:
                data[vs30].append(val)
            else:
                data[vs30] = [val]
        input_file.close()
    # Ok, processed all realizations, now combine the data
    sta_vs30_data = []
    sta_min_data = []
    sta_max_data = []
    sta_resid_data = []

    for item in data:
        sta_vs30_data.append(item)
        sta_min_data.append(numpy.std(data[item]))
        sta_max_data.append(numpy.std(data[item]))
        sta_resid_data.append(numpy.mean(data[item]))
    # Return the data we found
    return sta_vs30_data, sta_min_data, sta_max_data, sta_resid_data

def plot_combined_vs30_gof(indir, outdir, codebase):
    """
    This function reads data from the residuals files from multiple
    realizations and plots a vs30 gof plot with a number of periods
    """
    # Capture number of realizations and event label
    num_realizations = len(os.listdir(indir))
    basedir = os.path.join(indir, os.listdir(indir)[0])
    resid_file = glob.glob("%s%s*-resid-vs30-%.3f-rotd50.txt" %
                           (basedir, os.sep, VS30_PERIODS[0]))[0]
    event_label = os.path.basename(resid_file).split("-")[0]

    # Collect all the data from the residuals file
    all_sta_vs30_data = []
    all_sta_min_data = []
    all_sta_max_data = []
    all_sta_resid_data = []
    for period in VS30_PERIODS:
        (sta_vs30_data,
         sta_min_data,
         sta_max_data,
         sta_resid_data) = combine_realization_data(indir, period)
        all_sta_vs30_data.append(sta_vs30_data)
        all_sta_min_data.append(sta_min_data)
        all_sta_max_data.append(sta_max_data)
        all_sta_resid_data.append(sta_resid_data)

    # Now create the vs30 GOF
    outfile = os.path.join(outdir, "gof-vs30-combined-%s-%s-rotd50.png" %
                           (codebase, event_label))
    create_combined_vs30_gof(all_sta_resid_data, all_sta_vs30_data,
                             all_sta_min_data, all_sta_max_data,
                             event_label, num_realizations, codebase, outfile)

def create_combined_vs30_gof(all_sta_resid_data, all_sta_vs30_data,
                             all_sta_min_data, all_sta_max_data,
                             event_label, num_realizations, codebase, outfile):
    """
    Creates a combined Vs30 gof plot for all the data provided
    """
    plottitle = ("GOF Comparison for %s\n%d Realizations\n%s Method" %
                 (event_label, num_realizations, codebase.upper()))

    # Set up ticks to match matplotlib 1.x style
    matplotlib.rcParams['xtick.direction'] = 'in'
    matplotlib.rcParams['ytick.direction'] = 'in'
    matplotlib.rcParams['xtick.top'] = True
    matplotlib.rcParams['ytick.right'] = True
    # And errorbars too
    matplotlib.rcParams['errorbar.capsize'] = 3

    # Create figure
    num_plots = len(VS30_PERIODS)
    if len(VS30_PERIODS) % 2:
        num_plots = num_plots + 1
    num_columns = num_plots // 2
    fig, axs = pylab.plt.subplots(2, num_columns)
    fig.set_size_inches(18, 8)
    fig.subplots_adjust(left = 0.05, right = 0.95, top = 0.86, bottom = 0.06)
    
    # Find max, min values for x_axis
    min_x = 0
    max_x = 0
    for vs30s in all_sta_vs30_data:
        # Check if not empty
        if len(vs30s):
            max_x = max(max_x, max(vs30s))
    # If no data, set it to 900 (will get rounded to 1000)
    if max_x == 0:
        max_x = 900
    # Round to the next 10'
    max_x = max_x + (10 - (max_x % 10))
    # y-axis is fixed
    min_y = MIN_Y_AXIS
    max_y = MAX_Y_AXIS

    # Convert to list
    subfigs = []
    for y_subplot in range(0, 2):
        for x_subplot in range(0, num_columns):
            subfigs.append(axs[y_subplot, x_subplot])

    # Good, now walk through each subfig
    for (subfig, sta_vs30_data,
         sta_min_data, sta_max_data,
         sta_resid_data, period) in zip(subfigs, all_sta_vs30_data,
                                        all_sta_min_data, all_sta_max_data,
                                        all_sta_resid_data, VS30_PERIODS):
        subfig.set_xlim(min_x, max_x)
        subfig.set_ylim(min_y, max_y)
        subfig.set_title("Period = %.3f s" % (period), size=8)
        if VS30_PERIODS.index(period) % num_columns == 0:
            subfig.set_ylabel("ln (data/model)", size=8)
        subfig.tick_params(labelsize=8)
        subfig.errorbar(sta_vs30_data, sta_resid_data,
                        yerr=[sta_min_data, sta_max_data],
                        fmt='o', color='black', ecolor='grey',
                        label='_nolegend_')
        subfig.plot([min_x, max_x], [0.0, 0.0],
                    color='grey', label='_nolegend_')
        # Only add label to last row
        if VS30_PERIODS.index(period) >= (2 * num_columns) / 2:
            subfig.set_xlabel("Vs30 (m/s)", size=8)

    fig.suptitle('%s' % (plottitle), size=12)
    print("Saving Vs30 GoF plot to %s" % (outfile))
    fig.savefig(outfile, format="png", transparent=False, dpi=plot_config.dpi)

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

if OPTIONS.output_dir is None:
    PARSER.error("error specify output directory!")
else:
    OUTPUT_DIR = OPTIONS.output_dir
    if not os.path.isdir(OUTPUT_DIR):
        PARSER.error("Invalid output directory!")

if OPTIONS.codebase is None:
    PARSER.error("Please specify codebase!")

# Creates combined Vs30 gof plot
plot_combined_vs30_gof(INPUT_OUTDIR, OUTPUT_DIR, OPTIONS.codebase)

# All done!
print("All Done!")
