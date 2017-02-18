#!/usr/bin/python
"""
Southern California Earthquake Center Broadband Platform
Copyright 2010-2016 Southern California Earthquake Center

This module contains functions to plot GMPE results
$Id: plot_gmpe.py 1719 2016-08-18 21:44:13Z fsilva $
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
    #print("Plotting %s vs %s" % (sim_file, gmpe_file))

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
        # # Catch plot titles
        # if line.startswith("#period") or line.startswith("%period"):
        #    subplot_titles = line.split()[1:]
        #    continue
        # Other comments
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

    #if not subplot_titles:
    #    print "Couldn't find header line with GMPE labels! Aborting..."
    #    sys.exit(1)

    # Start plot
    num_plots = len(gmpe_ri50)
    if len(gmpe_ri50) % 2:
        num_plots = num_plots + 1
    num_columns = num_plots // 2
    fig, axs = pylab.plt.subplots(2, num_columns)
    fig.set_size_inches(8, 7)
#    fig.subplots_adjust(left = 0.05, right = 0.95)

    # Figure out min and max values
    min_x = min([min((periods1)), min(periods2)])
    max_x = max([max((periods1)), max(periods2)])
    min_y = min([min(gmpe_values) for gmpe_values in gmpe_ri50]) / 1.1
    max_y = 1.1 * max([max(gmpe_values) for gmpe_values in gmpe_ri50])

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
        subfig.set_title("%s" % subplot_title)
        subfig.plot(periods1, rd50_aa1, label=str(label1))
        subfig.plot(periods2, gmpe_values, label=str(label2))
        # Add ylabel only to leftmost plots
        if subplot_titles.index(subplot_title) % num_columns == 0:
            subfig.set_ylabel("PSA (g)")
        # Add xlabel only to bottom plots
        if subplot_titles.index(subplot_title) >= (2 * num_columns) // 2:
            subfig.set_xlabel("Period (s)")
        subfig.legend(prop=mpl.font_manager.FontProperties(size=10))
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
