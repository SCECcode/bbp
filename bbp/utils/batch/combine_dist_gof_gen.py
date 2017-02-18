#!/usr/bin/env python
"""
This program created a distance-based GOF, combining information from
all realizations into a single plot.

$Id: combine_dist_gof_gen.py 1333 2014-03-27 22:51:03Z fsilva $
"""

# Import Python modules
import os
import glob
import optparse
import matplotlib
if (matplotlib.get_backend() != 'agg'):
    matplotlib.use('Agg') # Disables use of Tk/X11
#import matplotlib.colors as mcolors
#import matplotlib.cm as cm
#from matplotlib.ticker import FormatStrFormatter
import pylab
import numpy

# Import Broadband modules
import plot_config
import bband_utils

# Constants
MIN_Y_AXIS = -1.75
MAX_Y_AXIS = 1.75
DIST_PERIODS = [0.01, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0]

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
        resid_file = glob.glob("%s%s*-resid-%.3f-rotd50.txt" %
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
            dist = pieces[0]
            val = pieces[1]
            if dist in data:
                data[dist].append(val)
            else:
                data[dist] = [val]
        input_file.close()
    # Ok, processed all realizations, now combine the data
    sta_dist_data = []
    sta_min_data = []
    sta_max_data = []
    sta_resid_data = []

    for item in data:
        sta_dist_data.append(item)
        sta_min_data.append(numpy.std(data[item]))
        sta_max_data.append(numpy.std(data[item]))
        sta_resid_data.append(numpy.mean(data[item]))
    # Return the data we found
    return sta_dist_data, sta_min_data, sta_max_data, sta_resid_data
        
def plot_combined_dist_gof(indir, outdir, codebase):
    """
    This function reads data from the residuals files from multiple
    realizations and plots a dist gof plot with a number of periods
    """
    # Capture number of realizations and event label
    num_realizations = len(os.listdir(indir))
    basedir = os.path.join(indir, os.listdir(indir)[0])
    resid_file = glob.glob("%s%s*-resid-%.3f-rotd50.txt" %
                           (basedir, os.sep, DIST_PERIODS[0]))[0]
    event_label = os.path.basename(resid_file).split("-")[0]

    # Collect all the data from the residuals file
    all_sta_dist_data = []
    all_sta_min_data = []
    all_sta_max_data = []
    all_sta_resid_data = []
    for period in DIST_PERIODS:
        (sta_dist_data,
         sta_min_data,
         sta_max_data,
         sta_resid_data) = combine_realization_data(indir, period)
        all_sta_dist_data.append(sta_dist_data)
        all_sta_min_data.append(sta_min_data)
        all_sta_max_data.append(sta_max_data)
        all_sta_resid_data.append(sta_resid_data)

    # Now create the dist GOFs
    outfile = os.path.join(outdir, "gof-dist-linear-combined-%s-%s-rotd50.png" %
                           (codebase, event_label))
    create_combined_dist_gof(all_sta_resid_data, all_sta_dist_data,
                             all_sta_min_data, all_sta_max_data,
                             event_label, num_realizations, codebase, outfile)
    outfile = os.path.join(outdir, "gof-dist-log-combined-%s-%s-rotd50.png" %
                           (codebase, event_label))
    create_combined_dist_gof(all_sta_resid_data, all_sta_dist_data,
                             all_sta_min_data, all_sta_max_data,
                             event_label, num_realizations, codebase, outfile,
                             log_scale=True)

def create_combined_dist_gof(all_sta_resid_data, all_sta_dist_data,
                             all_sta_min_data, all_sta_max_data,
                             event_label, num_realizations, codebase, outfile,
                             log_scale=False):
    """
    Creates a combined gof dist plot for all the data and distances
    provided
    """
    plottitle = ("GOF Comparison for %s\n%d Realizations\n%s Method" %
                 (event_label, num_realizations, codebase.upper()))

    # Create figure
    num_plots = len(DIST_PERIODS)
    if len(DIST_PERIODS) % 2:
        num_plots = num_plots + 1
    num_columns = num_plots / 2
    fig, axs = pylab.plt.subplots(2, num_columns)
    fig.set_size_inches(18, 8)
    fig.subplots_adjust(left = 0.05, right = 0.95, top = 0.86, bottom = 0.06)

    # Find max, min values for x_axis
    if log_scale:
        min_x = 1
    else:
        min_x = 0
    max_x = 0
    for dist in all_sta_dist_data:
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
    for (subfig, sta_dist_data,
         sta_min_data, sta_max_data,
         sta_resid_data, period) in zip(subfigs, all_sta_dist_data,
                                        all_sta_min_data, all_sta_max_data,
                                        all_sta_resid_data, DIST_PERIODS):

#        sta_min_data = [abs(x - y) for x, y in zip(sta_resid_data,
#                                                   sta_min_data)]
#        sta_max_data = [abs(x - y) for x, y in zip(sta_resid_data,
#                                                   sta_max_data)]

        subfig.set_xlim(min_x, max_x)
        subfig.set_ylim(min_y, max_y)
        subfig.set_title("Period = %.3f s" % (period), size=8)
        if DIST_PERIODS.index(period) % num_columns == 0:
            subfig.set_ylabel("ln (data/model)", size=8)
        subfig.tick_params(labelsize=8)
#        subfig.plot(sta_dist_data, sta_resid_data, 'o', color='black',
#                    label='_nolegend_')
        subfig.errorbar(sta_dist_data, sta_resid_data,
                        yerr=[sta_min_data, sta_max_data],
                        fmt='o', color='black', ecolor='grey',
                        label='_nolegend_')
        subfig.plot([min_x, max_x], [0.0, 0.0],
                    color='grey', label='_nolegend_')
        if log_scale:
            subfig.set_xscale('log')
        # Only add label to last row
        if DIST_PERIODS.index(period) >= (2 * num_columns) / 2:
            subfig.set_xlabel("Distance (km)", size=8)

    fig.suptitle('%s' % (plottitle), size=12)
    print "Saving dist GoF plot to %s" % (outfile)
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

# Create combined dist gof plot
plot_combined_dist_gof(INPUT_OUTDIR, OUTPUT_DIR, OPTIONS.codebase)

# All done!
print "All Done!"

