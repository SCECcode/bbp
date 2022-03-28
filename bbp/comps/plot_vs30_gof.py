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

Creates a Vs30 gof plot for a list of periods
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys
import matplotlib as mpl
if (mpl.get_backend() != 'agg'):
    mpl.use('Agg') # Disables use of Tk/X11
import pylab

# Import Broadband modules
from install_cfg import InstallCfg
import plot_config

# Constants
MIN_Y_AXIS = -1.75
MAX_Y_AXIS = 1.75
COMP_EXT_RD50 = 'rotd50'
COMP_TITLE_RD50 = 'RotD50'
PLOT_PERIODS = [0.01, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0]

def read_resid(resid_file, period, summary_output):
    """
    Reads the residual file resid_file and returns all data for the
    requested period
    """
    # Start empty
    data = []
    vs30s = []

    # Read residuals file and get information we need
    input_file = open(resid_file, 'r')
    # Look over header and figure out which column contains the period
    # we need to plot
    header = input_file.readline()
    header = header.strip()
    items = header.split()
    index = -1
    for idx, item in enumerate(items):
        try:
            val = float(item)
            if val == period:
                # Found period, save index
                index = idx
                break
        except:
            pass

    if index < 0:
        # If we don't have this period, nothing to do
        print("Residuals file %s does not have data for period %f" %
              (resid_file, period))
        # Close input file
        input_file.close()
        # Return empty sets
        return data, vs30s

    # Index #6 has vs30
    # Index #12 has component
    # Indexes #10 and #11 have period range for valid data

    # Read the rest of the file
    for line in input_file:
        items = line.split()
        comp = items[12]
        vs30 = items[6]
        tmin = items[10]
        tmax = items[11]
        value = items[index]
        # Skip components we don't know
        if comp != COMP_EXT_RD50:
            continue
        if period >= float(tmin) and period <= float(tmax):
            # Data within range, take it
            data.append(float(value))
            vs30s.append(float(vs30))

    # Done reading the file
    input_file.close()

    # Write summary output for later processing
    output_file = open(summary_output, 'w')
    for vs30, val in zip(vs30s, data):
        output_file.write("%f %f\n" % (vs30, val))
    output_file.close()

    # Return the data we found
    return data, vs30s

def plot_vs30_gof(resid_file, comp_label, sim_id):
    """
    Reads data from resid_file and creates a GoF Vs30 plot
    for all periods
    """
    # Get directory names
    install = InstallCfg.getInstance()
    a_outdir = os.path.join(install.A_OUT_DATA_DIR, str(sim_id))

    # Collect all the data
    all_data = []
    all_vs30s = []
    # Read the residuals data
    for period in PLOT_PERIODS:
        summary_output = os.path.join(a_outdir, "%s-%d-resid-vs30-%.3f-%s.txt" %
                                      (comp_label, sim_id,
                                       period, COMP_EXT_RD50))
        data, vs30s = read_resid(resid_file, period, summary_output)
        all_data.append(data)
        all_vs30s.append(vs30s)

    # Now create the GoF plot
    outfile = os.path.join(a_outdir, "gof-vs30-%s-%d-rotd50.png" %
                           (comp_label, sim_id))
    create_vs30_gof(all_data, all_vs30s,
                    comp_label, sim_id, outfile)

def create_vs30_gof(all_data, all_vs30s,
                    comp_label, sim_id, outfile):
    """
    Creates a Vs30 GoF plot
    """

    plottitle = ("GOF Comparison between %s and simulation %d" %
                 (comp_label, sim_id))

    # Create figure
    num_plots = len(PLOT_PERIODS)
    if len(PLOT_PERIODS) % 2:
        num_plots = num_plots + 1
    num_columns = num_plots // 2
    fig, axs = pylab.plt.subplots(2, num_columns)
    fig.set_size_inches(18, 8)
    fig.subplots_adjust(left=0.05)
    fig.subplots_adjust(right=0.95)
    fig.subplots_adjust(hspace=0.25)

    # Find max, min values for x_axis
    min_x = 0
    max_x = 0
    for vs30s in all_vs30s:
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
    for subfig, data, vs30s, period in zip(subfigs,
                                           all_data,
                                           all_vs30s,
                                           PLOT_PERIODS):
        subfig.set_xlim(min_x, max_x)
        subfig.set_ylim(min_y, max_y)
        subfig.set_title("Period = %.3f s" % (period), size=8)
        subfig.set_ylabel("ln (data/model)", size=8)
        subfig.tick_params(labelsize=7)
        subfig.plot(vs30s, data, 'o', color='black',
                    label='_nolegend_')
        subfig.plot([min_x, max_x], [0.0, 0.0],
                    color='grey', label='_nolegend_')
        subfig.set_xlabel("Vs30 (m/s)", size=8)

    fig.suptitle('%s' % (plottitle), size=12)
    print("==> Created Vs30 GoF plot: %s" % (outfile))
    fig.savefig(outfile, format="png", transparent=False, dpi=plot_config.dpi)
    pylab.close()

def usage():
    """
    Prints usage information
    """
    print("usage: %s <resid_file> <label> <sim_id>" % (sys.argv[0]))

if __name__ == '__main__':
    if len(sys.argv) < 4:
        usage()
        sys.exit(-1)

    plot_vs30_gof(sys.argv[1], sys.argv[2], int(sys.argv[3]))
