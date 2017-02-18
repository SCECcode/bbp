#!/usr/bin/env python
"""
Southern California Earthquake Center Broadband Platform
Copyright 2010-2016 Southern California Earthquake Center

Creates a distance gof plot for a list of periods
$Id: plot_dist_gof.py 1718 2016-08-16 21:55:36Z fsilva $
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys
import matplotlib as mpl
if (mpl.get_backend() != 'agg'):
    mpl.use('Agg') # Disables use of Tk/X11
import numpy
import pylab

# Import Broadband modules
#import bband_utils
from install_cfg import InstallCfg
import plot_config

# Constants
MIN_Y_AXIS = -1.75
MAX_Y_AXIS = 1.75
COMP_EXT_RD50 = 'rotd50'
COMP_TITLE_RD50 = 'RotD50'
DIST_PERIODS = [0.01, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0]
NGA_MODELS = ['AS', 'BA', 'CB', 'CY']

def read_gmpe_file(resid_file, period):
    """
    Reads the gmpe residuals file and returns all the data
    """
    gmpe_data = []

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
        return gmpe_data

    # Read the rest of the file
    # Index #2 has station name
    # Index #7 has distance
    for line in input_file:
        items = line.split()
        stat = items[2]
        dist = items[7]
        value = items[index]
        gmpe_data.append((stat, dist, value))

    # Done reading the file
    input_file.close()

    return gmpe_data

def read_resid(resid_file, period, summary_output):
    """
    Reads the residual file resid_file and returns all data for the
    requested period
    """
    # Start empty
    data = []
    distance = []

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
        return data, distance

    # Index #7 has distance
    # Index #12 has component
    # Indexes #10 and #11 have period range for valid data

    # Read the rest of the file
    for line in input_file:
        items = line.split()
        comp = items[12]
        dist = items[7]
        tmin = items[10]
        tmax = items[11]
        value = items[index]
        # Skip components we don't know
        if comp != COMP_EXT_RD50:
            continue
        if period >= float(tmin) and period <= float(tmax):
            # Data within range, take it
            data.append(float(value))
            distance.append(float(dist))

    # Done reading the file
    input_file.close()

    # Write summary output for later processing
    output_file = open(summary_output, 'w')
    for dist, val in zip(distance, data):
        output_file.write("%f %f\n" % (dist, val))
    output_file.close()

    # Return the data we found
    return data, distance

def read_gmpe_resid(a_outdir_gmpe, sim_id, period, summary_output):
    """
    This function reads the 4 GMPE residuals files, averages the
    results for each station using the data for the desired period
    """
    # Start empty
    data = []
    distance = []

    all_gmpe_data = {}
    merged_gmpe_data = []
    for gmpe in NGA_MODELS:
        resid_file = os.path.join(a_outdir_gmpe, "%s-%d.resid.txt" %
                                  (gmpe.lower(), sim_id))
        gmpe_data = read_gmpe_file(resid_file, period)
        for item in gmpe_data:
            # item[0] is station name
            # item[1] is station distance
            # item[2] is residual value
            if item[0] not in all_gmpe_data:
                all_gmpe_data[item[0]] = ((float(item[1]),
                                           [float(item[2])]))
            else:
                all_gmpe_data[item[0]][1].append(float(item[2]))

    # Now combine the data
    for station in all_gmpe_data:
        merged_gmpe_data.append((all_gmpe_data[station][0],
                                numpy.mean(all_gmpe_data[station][1])))

    # Almost there, sort it and split in the data and distance variables!
    for item in sorted(merged_gmpe_data):
        distance.append(item[0])
        data.append(item[1])

    # Write summary output for later processing
    output_file = open(summary_output, 'w')
    for dist, val in zip(distance, data):
        output_file.write("%f %f\n" % (dist, val))
    output_file.close()

    # Return the data we found
    return data, distance

def plot_dist_gof(resid_file, comp_label, a_outdir_gmpe, sim_id):
    """
    Reads data from resid_file and plots a gof distance plot all
    periods
    """
    # Get directory names
    install = InstallCfg.getInstance()
    a_outdir = os.path.join(install.A_OUT_DATA_DIR, str(sim_id))

    # Collect all the data
    all_data = []
    all_distances = []
    all_gmpe_data = []
    all_gmpe_distances = []
    # Read the residuals data
    for period in DIST_PERIODS:
        summary_output = os.path.join(a_outdir, "%s-%d-resid-%.3f-%s.txt" %
                                      (comp_label, sim_id,
                                       period, COMP_EXT_RD50))
        data, distance = read_resid(resid_file, period, summary_output)
        all_data.append(data)
        all_distances.append(distance)
    # Now do the same for the GMPE data
#    for period in DIST_PERIODS:
#        if os.path.isdir(a_outdir_gmpe):
#            summary_output = os.path.join(a_outdir,
#                                          "%s-%d-resid-gmpe-%.3f-%s.txt" %
#                                          (comp_label, sim_id,
#                                           period, COMP_EXT_RD50))
#            data, distance = read_gmpe_resid(a_outdir_gmpe, sim_id,
#                                             period, summary_output)
#            all_gmpe_data.append(data)
#            all_gmpe_distances.append(distance)
#        else:
#            all_gmpe_data.append([])
#            all_gmpe_distances.append([])

    # Now create the 2 plots, 1 linear and 1 log
    outfile = os.path.join(a_outdir, "gof-dist-linear-%s-%d-rotd50.png" %
                           (comp_label, sim_id))
    create_dist_gof(all_data, all_distances,
                    all_gmpe_data, all_gmpe_distances,
                    comp_label, sim_id, outfile)

    outfile = os.path.join(a_outdir, "gof-dist-log-%s-%d-rotd50.png" %
                           (comp_label, sim_id))
    create_dist_gof(all_data, all_distances,
                    all_gmpe_data, all_gmpe_distances,
                    comp_label, sim_id, outfile, log_scale=True)

def create_dist_gof(all_data, all_distances,
                    all_gmpe_data, all_gmpe_dist,
                    comp_label, sim_id, outfile, log_scale=False):
    """
    Creates a gof distance plots for all the data and distances
    provided
    """

    plottitle = ("GOF Comparison between %s and simulation %d" %
                 (comp_label, sim_id))

    # Create figure
    num_plots = len(DIST_PERIODS)
    if len(DIST_PERIODS) % 2:
        num_plots = num_plots + 1
    num_columns = num_plots // 2
    fig, axs = pylab.plt.subplots(2, num_columns)
    fig.set_size_inches(18, 8)
    fig.subplots_adjust(left=0.05, right=0.95)

    # Find max, min values for x_axis
    if log_scale:
        min_x = 1
    else:
        min_x = 0
    max_x = 0
    for dist in all_distances:
        # Check if not empty
        if len(dist):
            max_x = max(max_x, max(dist))
    # If no data, set it to 90 (will get rounded to 100)
    if max_x == 0:
        max_x = 90
    # Round to the next 10'
    max_x = max_x + (10 - (max_x % 10))
    if log_scale and max_x > 100:
        # Round to the next 100'
        max_x = max_x + (100 - (max_x % 100))
    # y-axis is fixed
    min_y = MIN_Y_AXIS
    max_y = MAX_Y_AXIS

    # Convert to list
    subfigs = []
    for y_subplot in range(0, 2):
        for x_subplot in range(0, num_columns):
            subfigs.append(axs[y_subplot, x_subplot])

    # Good, now walk through each subfig
    for subfig, data, dist, period in zip(subfigs,
                                          all_data,
                                          all_distances,
                                          DIST_PERIODS):
        subfig.set_xlim(min_x, max_x)
        subfig.set_ylim(min_y, max_y)
        subfig.set_title("Period = %.3f s" % (period), size=8)
        if DIST_PERIODS.index(period) % num_columns == 0:
            subfig.set_ylabel("ln (data/model)", size=8)
        subfig.tick_params(labelsize=8)
        subfig.plot(dist, data, 'o', color='black',
                    label='_nolegend_')
        subfig.plot([min_x, max_x], [0.0, 0.0],
                    color='grey', label='_nolegend_')
        if log_scale:
            subfig.set_xscale('log')
        # Only add label to last row
        if DIST_PERIODS.index(period) >= (2 * num_columns) // 2:
            subfig.set_xlabel("Distance (km)", size=8)

    fig.suptitle('%s' % (plottitle), size=12)
    print("==> Created Distance GoF plot: %s" % (outfile))
    fig.savefig(outfile, format="png", transparent=False, dpi=plot_config.dpi)
    pylab.close()

def usage():
    """
    Prints usage information
    """
    print("usage: %s <resid_file> <label> <gmpe_dir> <sim_id>" % (sys.argv[0]))

if __name__ == '__main__':
    if len(sys.argv) < 5:
        usage()
        sys.exit(1)

    plot_dist_gof(sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]))
    sys.exit(0)
