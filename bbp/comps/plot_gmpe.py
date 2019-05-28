#!/usr/bin/python
"""
Copyright 2010-2019 University Of Southern California

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

This module contains functions to plot GMPE results
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys
import matplotlib as mpl
mpl.use('AGG', warn=False)
import pylab

# Import plot config file
import plot_config

def plot_gmpe(stat, sim_file, gmpe_file, gmpe_labels, label1, label2, outfile):
    """
    This function generates comparison plots between the simulated
    results and the gmpe data
    """
    periods1 = []
    rd50_aa1 = []

    periods2 = []
    gmpe_ri50 = []

    # Read simulated rd50 data file
    simfile = open(sim_file, 'r')
    for line in simfile:
        if line.startswith("#") or line.startswith("%"):
            continue
        pieces = line.split()
        periods1.append(float(pieces[0]))
        rd50_aa1.append(float(pieces[3]))
    simfile.close()

    if periods1 == []:
        print("Input file %s is missing data! Aborting..." % (sim_file))
        sys.exit(1)

    # Read gmpe data file
    subplot_titles = []
    gmpefile = open(gmpe_file, 'r')
    for line in gmpefile:
        if line.startswith("#") or line.startswith("%"):
            continue
        pieces = line.split()
        periods2.append(float(pieces[0]))
        # Figure out if we need to create list
        if not gmpe_ri50:
            for item in pieces[1:]:
                gmpe_ri50.append([])
        for item, dst in zip(pieces[1:], gmpe_ri50):
            dst.append(float(item))
    gmpefile.close()

    if not periods2:
        print("Input file %s is missing data! Aborting..." % (gmpe_file))
        sys.exit(1)

    subplot_titles = gmpe_labels

    # Start plot
    num_plots = len(gmpe_ri50)
    if len(gmpe_ri50) % 2:
        num_plots = num_plots + 1
    num_columns = num_plots // 2
    fig, axs = pylab.plt.subplots(2, num_columns)
    fig.set_size_inches(8, 7)
    pylab.subplots_adjust(left=0.075)
    pylab.subplots_adjust(right=0.975)
    pylab.subplots_adjust(hspace=0.3)

    # Figure out min and max values
    min_x = min([min((periods1)), min(periods2)])
    max_x = max([max((periods1)), max(periods2)])
    min_y = min([min(gmpe_values) for gmpe_values in gmpe_ri50]) / 1.1
    max_y_gmpes = max([max(gmpe_values) for gmpe_values in gmpe_ri50])
    max_y = 1.1 * max(max_y_gmpes, max(rd50_aa1))

    # Convert to list
    subfigs = []
    for y_subplot in range(0, 2):
        for x_subplot in range(0, num_columns):
            subfigs.append(axs[y_subplot, x_subplot])

    # Now walk through each subfig, if we have more subfigs that we
    # have data for, the for loop will only create enough subfigs to
    # cover the data, ignoring the extra one
    for subfig, subplot_title, gmpe_values in zip(subfigs,
                                                  subplot_titles,
                                                  gmpe_ri50):
        subfig.set_xlim(min_x, max_x)
        subfig.set_ylim(min_y, max_y)
        subfig.set_title("%s" % subplot_title, fontsize='small')
        subfig.plot(periods1, rd50_aa1, label=str(label1))
        subfig.plot(periods2, gmpe_values, label=str(label2))
        subfig.set_ylabel("PSA (g)")
        subfig.set_xlabel("Period (s)")
        subfig.legend(prop=mpl.font_manager.FontProperties(size=8))
        subfig.set_xscale('log')

    # All done, label and save it!
    fig.suptitle('GMPE Comparison for station %s, sim %s' %
                 (stat, label1), size=14)
    fig.savefig(outfile, format="png", dpi=plot_config.dpi)
    pylab.close()

if __name__ == '__main__':
    PROG_BASE = os.path.split(sys.argv[0])[1]
    plot_gmpe(sys.argv[1], sys.argv[2], sys.argv[3],
              sys.argv[4], sys.argv[5], sys.argv[6])
