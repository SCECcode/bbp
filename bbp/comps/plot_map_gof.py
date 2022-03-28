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

Creates a map gof plot for a list of periods
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys
import matplotlib as mpl
if (mpl.get_backend() != 'agg'):
    mpl.use('Agg') # Disables use of Tk/X11
import matplotlib.colors as mcolors
import matplotlib.cm as cm
from matplotlib.ticker import FormatStrFormatter
import pylab

# Import Broadband modules
from install_cfg import InstallCfg
import plot_utils
import PlotMap
import fault_utils
import plot_config

# Constants
MIN_Y_AXIS = -1.75
MAX_Y_AXIS = 1.75
COMP_EXT_RD50 = 'rotd50'
COMP_TITLE_RD50 = 'RotD50'
DIST_PERIODS = [0.01, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0]

def read_resid(resid_file, period, summary_output):
    """
    Reads the residual file resid_file and returns all data for the
    requested period
    """
    # Start empty
    sta_x_data = []
    sta_y_data = []
    sta_resid_data = []

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
        print ("Residuals file %s does not have data for period %f" %
               (resid_file, period))
        # Close input file
        input_file.close()
        # Return empty sets
        return sta_x_data, sta_y_data, sta_resid_data

    # Index #3 has lon, #4 has lat
    # Index #12 has component
    # Indexes #10 and #11 have period range for valid data

    # Read the rest of the file
    for line in input_file:
        items = line.split()
        comp = items[12]
        lon = items[3]
        lat = items[4]
        tmin = items[10]
        tmax = items[11]
        value = items[index]
        # Skip components we don't know
        if comp != COMP_EXT_RD50:
            continue
        if period >= float(tmin) and period <= float(tmax):
            # Data within range, take it
            sta_x_data.append(float(lon))
            sta_y_data.append(float(lat))
            sta_resid_data.append(float(value))

    # Done reading the file
    input_file.close()

    # Write summary output for later processing
    output_file = open(summary_output, 'w')
    for lon, lat, val in zip(sta_x_data, sta_y_data, sta_resid_data):
        output_file.write("%f %f %f\n" % (lon, lat, val))
    output_file.close()

    # Return the data we found
    return sta_x_data, sta_y_data, sta_resid_data

def plot_map_gof(r_srcfile, r_stations, resid_file, comp_label, sim_id):
    """
    Reads data from resid_file and plots a map gof plot with a number
    of periods
    """
    # Make sure we have a src or srf file
    if (r_srcfile is None or r_srcfile == "" or
        (not r_srcfile.endswith(".srf") and
         not r_srcfile.endswith(".src"))):
        # We need a SRC or SRF file to get the fault geometry
        return

    # Get directory names
    install = InstallCfg.getInstance()
    a_indir = os.path.join(install.A_IN_DATA_DIR, str(sim_id))
    a_outdir = os.path.join(install.A_OUT_DATA_DIR, str(sim_id))

    a_input_file = os.path.join(a_indir, r_srcfile)
    a_station_file = os.path.join(a_indir, r_stations)

    # Define boundaries to plot using the stations in the station file
    (north, south,
     east, west) = plot_utils.set_boundaries_from_stations(a_station_file,
                                                           a_input_file)
    trace_file = "%s.trace" % (a_input_file)
    simple_station_file = "%s.simple" % (a_station_file)

    if r_srcfile.endswith(".srf"):
        plot_utils.write_fault_trace(a_input_file, trace_file)
    else:
        plot_utils.write_simple_trace(a_input_file, trace_file)
    plot_utils.write_simple_stations(a_station_file, simple_station_file)

    # Get hypo_lon, hypo_lat from src/srf file
    hypo_lon, hypo_lat = fault_utils.calculate_epicenter(a_input_file)

    plotregion = [west, east, south, north]
    topo = os.path.join(install.A_PLOT_DATA_DIR, 'calTopo18.bf')
    coastal = os.path.join(install.A_PLOT_DATA_DIR, 'gshhs_h.txt')
    border = os.path.join(install.A_PLOT_DATA_DIR, 'wdb_borders_h.txt')

    # Collect all the data from the residuals file
    all_sta_x_data = []
    all_sta_y_data = []
    all_sta_resid_data = []
    for period in DIST_PERIODS:
        summary_output = os.path.join(a_outdir, "%s-%d-resid-map-%.3f-%s.txt" %
                                      (comp_label, sim_id,
                                       period, COMP_EXT_RD50))
        sta_x_data, sta_y_data, sta_resid_data = read_resid(resid_file,
                                                            period,
                                                            summary_output)
        all_sta_x_data.append(sta_x_data)
        all_sta_y_data.append(sta_y_data)
        all_sta_resid_data.append(sta_resid_data)

    # Now create the map GOF
    outfile = os.path.join(a_outdir, "gof-map-%s-%d-rotd50.png" %
                           (comp_label, sim_id))
    create_map_gof(all_sta_x_data, all_sta_y_data, all_sta_resid_data,
                   plotregion, topo, coastal, border, trace_file,
                   comp_label, sim_id, outfile, hypo_lat=hypo_lat,
                   hypo_lon=hypo_lon)

def create_map_gof(all_sta_x_data, all_sta_y_data, all_sta_resid_data,
                   plotregion, topo, coastal, border, fault, comp_label,
                   sim_id, outfile, hypo_lat=None, hypo_lon=None):
    """
    Creates a gof distance plots for all the data and distances
    provided
    """

    plottitle = ("GOF Comparison between %s and simulation %d" %
                 (comp_label, sim_id))

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
    num_columns = num_plots // 2
    fig, axs = pylab.plt.subplots(2, num_columns)
    fig.set_size_inches(12, 6.5)
    #fig.autofmt_xdate()

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
#     # Find vmin and vmax for all plots
#     vmin = 0.0
#     vmax = 0.0
#     for sta_resid_data in all_sta_resid_data:
#         if len(sta_resid_data):
#             vmin = min(vmin, min(sta_resid_data))
#             vmax = max(vmax, max(sta_resid_data))
#     # But make it symmetrical
#     if abs(vmax) > abs(vmin):
#         vmin = -vmax
#     else:
#         vmax = -vmin

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
        for idx in range(0, len(coast_x)):
            subfig.plot(coast_x[idx], coast_y[idx], linestyle='-', color='0.75')

        # Plot borders
        for idx in range(0, len(bord_x)):
            subfig.plot(bord_x[idx], bord_y[idx], linestyle='-', color='0.75')

        # Plot fault trace
        subfig.plot(fault_x, fault_y, linestyle='-', color='k', linewidth=1.0)

        # If we don't have at least 1 station for this period, create
        # a fake station outside of the map area so that we can still
        # create the empty plot with the colobar to the right
        if not len(sta_x_data) or not len(sta_y_data):
            sta_x_data = [1000.0]
            sta_y_data = [1000.0]
            sta_resid_data = [0.0]

        # Plot hypocenter
        if hypo_lat is not None and hypo_lon is not None:
            hypo_lat = [hypo_lat]
            hypo_lon = [hypo_lon]
            subfig.scatter(hypo_lon, hypo_lat, marker=(5, 1, 0),
                           color='y', s=50)

        # Plot the stations
        im = subfig.scatter(sta_x_data, sta_y_data, s=20, c=sta_resid_data,
                            cmap=cm.jet_r, vmin=vmin, vmax=vmax,
                            marker='o', edgecolors='k')

        # Adding colorbars to the right of each row
#        if DIST_PERIODS.index(period) % num_columns == num_columns - 1:
#            subfig.figure.colorbar(im, ax=subfig)

        # Set degree formatting of tick values
        major_formatter = FormatStrFormatter(u'%.1f\u00b0')
        subfig.xaxis.set_major_formatter(major_formatter)
        subfig.yaxis.set_major_formatter(major_formatter)

#        # Turn on ticks for both sides of axis
#        for tick in subfig.xaxis.get_major_ticks():
#            tick.label1On = True
#            tick.label2On = True
#        for tick in subfig.yaxis.get_major_ticks():
#            tick.label1On = True
#            tick.label2On = True

        # Set font size
        for tick in subfig.get_xticklabels():
            tick.set_fontsize(6)
            tick.set_ha("right")
            tick.set_rotation(30)
        for tick in subfig.get_yticklabels():
            tick.set_fontsize(6)

        subfig.set_title("Period = %.3f s" % (period), size=8)

    fig.subplots_adjust(left=0.05, right=0.91, hspace=0.0,
                        top=0.95, bottom=0.05)
    colorbar_ax = fig.add_axes([0.93, 0.20, 0.02, 0.6])
    fig.colorbar(im, cax=colorbar_ax)
    fig.suptitle('%s' % (plottitle), size=12)
    print("==> Created Map GoF plot: %s" % (outfile))
    fig.savefig(outfile, format="png", transparent=False, dpi=plot_config.dpi)
    pylab.close()

def usage():
    """
    Prints usage information
    """
    print("usage: %s <src_file> <stations> <resid_file> <label> <sim_id>" %
          (sys.argv[0]))
    return

if __name__ == '__main__':
    if len(sys.argv) < 5:
        usage()
        sys.exit(1)

    plot_map_gof(sys.argv[1], sys.argv[2], sys.argv[3],
                 sys.argv[4], int(sys.argv[5]))
    sys.exit(0)
