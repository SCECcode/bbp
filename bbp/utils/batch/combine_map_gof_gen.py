#!/usr/bin/env python
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

This program created a map-based GOF, combining information from all
realizations into a single map plot where the color of each station is
the average bias from all realizations.
"""

# Import Python modules
import os
import glob
import optparse
import matplotlib
if (matplotlib.get_backend() != 'agg'):
    matplotlib.use('Agg') # Disables use of Tk/X11
import matplotlib.colors as mcolors
import matplotlib.cm as cm
from matplotlib.ticker import FormatStrFormatter
import pylab
import numpy

# Import Broadband modules
from install_cfg import InstallCfg
import PlotMap
import fault_utils
import plot_config
import bband_utils

# Constants
# Use an extra buffer to plot the region around all stations (in degrees)
BUFFER_LATITUDE = 0.25
BUFFER_LONGITUDE = 0.25
DIST_PERIODS = [0.01, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0]

# --------------------------------------------------------------------------
# Functions
# --------------------------------------------------------------------------

def set_boundaries_from_lon_lat(all_sta_lon, all_sta_lat):
    """
    This function sets the north, south, east, and west boundaries
    of the region we should plot, using the stations' locations in
    a lat/lon list
    """
    # Start without anything
    north = None
    south = None
    east = None
    west = None

    all_lon = []
    all_lat = []
    for lon_list, lat_list in zip(all_sta_lon, all_sta_lat):
        for sta_lon, sta_lat in zip(lon_list, lat_list):
            all_lon.append(sta_lon)
            all_lat.append(sta_lat)

    for lon, lat in zip(all_lon, all_lat):
        # If this is the first station, use its location
        if north is None:
            north = lat
            south = lat
            east = lon
            west = lon
            # Next station
            continue
        if lat > north:
            north = lat
        elif lat < south:
            south = lat
        if lon > east:
            east = lon
        elif lon < west:
            west = lon

    # Try to make the plot more symmetric
    lat_range = abs(north - south)
    lon_range = abs(east - west)
    if (lat_range > lon_range):
        diff = lat_range - lon_range
        diff = diff / 2.0
        east = east + diff
        west = west - diff
    elif (lon_range > lat_range):
        diff = lon_range - lat_range
        diff = diff / 2.0
        north = north + diff
        south = south - diff

    # Great, now we just add a buffer on each side
    if north < (90 - BUFFER_LATITUDE):
        north = north + BUFFER_LATITUDE
    else:
        north = 90
    if south > (-90 + BUFFER_LATITUDE):
        south = south - BUFFER_LATITUDE
    else:
        south = -90
    if east < (180 - BUFFER_LONGITUDE):
        east = east + BUFFER_LONGITUDE
    else:
        east = 180
    if west > (-180 + BUFFER_LONGITUDE):
        west = west - BUFFER_LONGITUDE
    else:
        west = -180

    return north, south, east, west

def combine_realization_data(tmpdir, period):
    """
    This function reads the resid-map files from all realizations and
    returns the combined data set
    """
    data = {}

    realizations = sorted(os.listdir(tmpdir))
    for realization in realizations:
        basedir = os.path.join(tmpdir, realization)
        resid_file = glob.glob("%s%s*-resid-map-%.3f-rotd50.txt" %
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
            if len(pieces) != 3:
                continue
            # Convert to floats
            pieces = [float(piece) for piece in pieces]
            lon = pieces[0]
            lat = pieces[1]
            val = pieces[2]
            if (lon, lat) in data:
                data[(lon, lat)].append(val)
            else:
                data[(lon, lat)] = [val]
        input_file.close()
    # Ok, processed all realizations, now combine the data
    sta_x_data = []
    sta_y_data = []
    sta_resid_data = []
    for item in data:
        sta_x_data.append(item[0])
        sta_y_data.append(item[1])
        sta_resid_data.append(numpy.mean(data[item]))
    # Return the data we found
    return sta_x_data, sta_y_data, sta_resid_data

def plot_combined_map_gof(indir, tmpdir, outdir, codebase):
    """
    This function reads data from the residuals files from multiple
    realizations and plots a map gof plot with a number of periods.
    """
    # Capture number of realizations and event label
    num_realizations = len(os.listdir(tmpdir))
    basedir = os.path.join(tmpdir, os.listdir(tmpdir)[0])
    resid_file = glob.glob("%s%s*-resid-map*rotd50.txt" %
                           (basedir, os.sep))[0]
    event_label = os.path.basename(resid_file).split("-")[0]

    # Get one trace file
    basedir = os.path.join(indir, os.listdir(indir)[0])
    trace_file = glob.glob("%s%s*.trace" % (basedir, os.sep))[0]

    # Now get the SRC (or SRF file) in order to get the hypocenter
    # location. Note that this function will look for the hypocenter
    # location from the first realization. If the simulation was
    # created using randomized hypocenter locations, the plot will
    # only display the location of the hypocenter from the first
    # realization.
    src_file = glob.glob("%s%s*.src" % (basedir, os.sep))
    if not len(src_file):
        srf_file = glob.glob("%s%s*.srf" % (basedir, os.sep))
        if not len(srf_file):
            raise bband_utils.ProcessingError("Cannot find SRC/SRF file!")
        source_file = srf_file[0]
    else:
        source_file = src_file[0]

    # Get hypo_lon, hypo_lat from src/srf file
    hypo_lon, hypo_lat = fault_utils.calculate_epicenter(source_file)

    # Collect all the data from the residuals file
    all_sta_x_data = []
    all_sta_y_data = []
    all_sta_resid_data = []
    for period in DIST_PERIODS:
        (sta_x_data,
         sta_y_data,
         sta_resid_data) = combine_realization_data(tmpdir, period)
        all_sta_x_data.append(sta_x_data)
        all_sta_y_data.append(sta_y_data)
        all_sta_resid_data.append(sta_resid_data)

    # Get plot boundaries
    (north, south,
     east, west) = set_boundaries_from_lon_lat(all_sta_x_data,
                                               all_sta_y_data)

    # Get directory names
    install = InstallCfg.getInstance()
    # Prepare to plot map GOF
    plotregion = [west, east, south, north]
    topo = os.path.join(install.A_PLOT_DATA_DIR, 'calTopo18.bf')
    coastal = os.path.join(install.A_PLOT_DATA_DIR, 'gshhs_h.txt')
    border = os.path.join(install.A_PLOT_DATA_DIR, 'wdb_borders_h.txt')

    # Now create the map GOF
    outfile = os.path.join(outdir, "gof-map-combined-%s-%s-rotd50.png" %
                           (codebase, event_label))

    create_combined_map_gof(all_sta_x_data, all_sta_y_data, all_sta_resid_data,
                            plotregion, topo, coastal, border, trace_file,
                            event_label, num_realizations, codebase, outfile,
                            hypo_lat=hypo_lat, hypo_lon=hypo_lon)

def create_combined_map_gof(all_sta_x_data, all_sta_y_data, all_sta_resid_data,
                            plotregion, topo, coastal, border, fault,
                            event_label, num_realizations, codebase, outfile,
                            hypo_lat=None, hypo_lon=None):
    """
    Creates a combined gof map plot for all the data and distances
    provided
    """
    plottitle = ("GOF Comparison for %s\n%d Realizations\n%s Method" %
                 (event_label, num_realizations, codebase.upper()))

    # Read in topo data
    topo_points = PlotMap.read_topo(topo, plotregion)

    # Read in fault data
    fault_x, fault_y = PlotMap.read_fault(fault)

    # Read coastlines
    coast_x, coast_y = PlotMap.read_coastal(coastal, plotregion)

    # Read borders
    bord_x, bord_y = PlotMap.read_coastal(border, plotregion)

    # Create figure
    num_plots = len(DIST_PERIODS)
    if len(DIST_PERIODS) % 2:
        num_plots = num_plots + 1
    num_columns = num_plots / 2
    fig, axs = pylab.plt.subplots(2, num_columns)
    fig.set_size_inches(12, 6.5)
    fig.autofmt_xdate()

    # Setup color scale
    cmap = cm.gist_gray
    norm = mcolors.Normalize(vmin=-2000.0, vmax=3000.0)

    # Convert to list
    subfigs = []
    for y_subplot in range(0, 2):
        for x_subplot in range(0, num_columns):
            subfigs.append(axs[y_subplot, x_subplot])

    # Fixed vmin and vmax for all plots
    vmin = -1.5
    vmax = 1.5

    # Good, now walk through each subfig
    for (subfig, sta_x_data, sta_y_data,
         sta_resid_data, period) in zip(subfigs, all_sta_x_data, all_sta_y_data,
                                        all_sta_resid_data, DIST_PERIODS):
        # Plot basemap
        subfig.imshow(topo_points, cmap=cmap, norm=norm,
                      extent=plotregion, interpolation='nearest')

        # Freeze the axis extents
        subfig.set_autoscale_on(False)

        # Plot coast lines
        for idx in xrange(0, len(coast_x)):
            subfig.plot(coast_x[idx], coast_y[idx], linestyle='-', color='0.75')

        # Plot borders
        for idx in xrange(0, len(bord_x)):
            subfig.plot(bord_x[idx], bord_y[idx], linestyle='-', color='0.75')

        # Plot fault trace
        subfig.plot(fault_x, fault_y, linestyle='-', color='k')

        # Plot hypocenter
        if hypo_lat is not None and hypo_lon is not None:
            hypo_lat = [hypo_lat]
            hypo_lon = [hypo_lon]
            subfig.scatter(hypo_lon, hypo_lat, marker=(5,1,0),
                           color='y', s=50)

        # Plot the stations
        im = subfig.scatter(sta_x_data, sta_y_data, s=20, c=sta_resid_data,
                            cmap=cm.jet_r, vmin=vmin, vmax=vmax, marker='o')

        # Set degree formatting of tick values
        major_formatter = FormatStrFormatter(u'%.1f\u00b0')
        subfig.xaxis.set_major_formatter(major_formatter)
        subfig.yaxis.set_major_formatter(major_formatter)

        #Disable colorbar on each plot
        #subfig.figure.colorbar(im, ax=subfig)

        # Set font size
        for tick in subfig.get_xticklabels():
            tick.set_fontsize(6)
        for tick in subfig.get_yticklabels():
            tick.set_fontsize(6)

        subfig.set_title("Period = %.3f s" % (period), size=8)

    # Slightly different values for top/bottom since the combined plot
    # has 3 title lines
    fig.subplots_adjust(left = 0.05, right = 0.91, hspace = 0.0,
                        top = 0.92, bottom = 0.02)
    colorbar_ax = fig.add_axes([0.93, 0.17, 0.02, 0.6])
    fig.colorbar(im, cax=colorbar_ax)
    fig.suptitle('%s' % (plottitle), size=12)
    print "Saving map GoF plot to %s" % (outfile)
    fig.savefig(outfile, format="png", transparent=False, dpi=plot_config.dpi)

# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

PARSER = optparse.OptionParser()
PARSER.add_option("-d", "--dir", dest="input_dir",
                  help="Input directory containing simulation results")
PARSER.add_option("-o", "--output_dir", dest="output_dir",
                  help="Directory where produced map plot will go")
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
INPUT_INDIR = os.path.join(TOP_INPUT_DIR, "Sims", "indata")

if OPTIONS.output_dir is None:
    PARSER.error("error specify output directory!")
else:
    OUTPUT_DIR = OPTIONS.output_dir
    if not os.path.isdir(OUTPUT_DIR):
        PARSER.error("Invalid output directory!")

if OPTIONS.codebase is None:
    PARSER.error("Please specify codebase!")

# Create combined map gof plot
plot_combined_map_gof(INPUT_INDIR, INPUT_OUTDIR, OUTPUT_DIR, OPTIONS.codebase)

# All done!
print "All Done!"
