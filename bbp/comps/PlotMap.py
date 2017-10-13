#!/bin/env python
"""
Copyright 2010-2017 University Of Southern California

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

This module plots station map and fault trace
"""
from __future__ import division, print_function

# Import Python Modules
import sys
import struct
import numpy as np
import matplotlib as mpl
if mpl.get_backend() != 'agg':
    mpl.use('Agg') # Disables use of Tk/X11
import matplotlib.colors as mcolors
import matplotlib.cm as cm
from matplotlib.ticker import FormatStrFormatter
import pylab

# Import plot config file
import plot_config

# GMT "bf" native binary header and format
GMT_HDR_FORMAT = '<3i10d80s80s80s80s320s160s'
GMT_DATA_FORMAT = '<f'
# Note: this script uses the "bf" native C float format
# The existing GMT topo grid in NetCDF format can be
# converted to this format with the following command:
# $ grdreformat calTopo18.grd calTopo18.test=bf -V

def in_region(point, plotregion):
    """
    Returns true if point falls within the plotregion
    """

    if ((point[0] >= plotregion[0]) and
        (point[0] < plotregion[1]) and
        (point[1] >= plotregion[2]) and
        (point[1] < plotregion[3])):
        return True
    else:
        return False

def read_fault(filename):
    """
    Read in fault file
    """

    fault_x = []
    fault_y = []
    fault_file = open(filename)

    for segment in fault_file:
        x, y = segment.split()
        fault_x.append(float(x))
        fault_y.append(float(y))

    fault_file.close()

    return fault_x, fault_y

def read_stations(filename):
    """
    Read in station list
    """

    sta_x = []
    sta_y = []
    station_file = open(filename)

    for station in station_file:
        x, y, _ = station.split()
        sta_x.append(float(x))
        sta_y.append(float(y))

    station_file.close()

    return sta_x, sta_y

def read_coastal(filename, plotregion):
    """
    Read in coastal geometry as a list of lists of polygon segments
    """

    # Initialize all variables
    coast_x = []
    coast_y = []
    poly_x = []
    poly_y = []
    segnum = 0
    segments = 0

    # Read in file
    polygons = open(filename)

    # Parse polygons
    for line in polygons:
        tokens = line.split()
        if (tokens[0] == 'P') or (tokens[0] == 'L'):
            if (len(poly_x) > 0):
                coast_x.append(poly_x)
                coast_y.append(poly_y)
            poly_x = []
            poly_y = []
            segnum = 0
            segments = int(tokens[2])
        else:
            if (segnum >= segments):
                print("Invalid number of segments in " +
                      "polygon from file %s" % (file))
                return([], [])
            segnum = segnum + 1
            x = float(tokens[0])
            y = float(tokens[1])
            if (in_region([x, y], plotregion)):
                poly_x.append(x)
                poly_y.append(y)
            else:
                if (len(poly_x) > 0):
                    coast_x.append(poly_x)
                    coast_y.append(poly_y)
                poly_x = []
                poly_y = []

    # Remember to close file
    polygons.close()

    return coast_x, coast_y

def read_topo(filename, plotregion):
    """
    Reads in topo data that is saved in GMT format:
    bf   GMT native, C-binary format (float)
    Header size is 892 bytes
    """
    # Open input file
    topo_file = open(filename, 'r')

    # Parse header
    buf = topo_file.read(struct.calcsize(GMT_HDR_FORMAT))
    header = struct.unpack(GMT_HDR_FORMAT, buf)
    topo_dims = [header[0], header[1]]
    topo_region = [header[3], header[4], header[5], header[6]]

    # Read elevation values
    data = np.arange(topo_dims[0] * topo_dims[1],
                     dtype=float).reshape(topo_dims[1],
                                          topo_dims[0])

    buf = topo_file.read(topo_dims[0] * topo_dims[1] *
                         struct.calcsize(GMT_DATA_FORMAT))
    topo_file.close()

    # Data is x-fast
    for y in xrange(0, topo_dims[1]):
        for x in xrange(0, topo_dims[0]):
            offset = ((y * topo_dims[0] + x) *
                      struct.calcsize(GMT_DATA_FORMAT))
            data[y][x] = struct.unpack(GMT_DATA_FORMAT,
                                       buf[offset:offset + 4])[0]

    # Pull out sub-matrix for plotregion, and invert y-axis
    x0 = int((plotregion[0] - topo_region[0]) / header[9])
    y1 = topo_dims[1] - int((plotregion[2] -
                             topo_region[2]) / header[10])
    x1 = int((plotregion[1] - topo_region[0]) / header[9])
    y0 = topo_dims[1] - int((plotregion[3] -
                             topo_region[2]) / header[10])
    subdata = np.arange((x1 - x0) * (y1 - y0),
                        dtype=float).reshape(y1 - y0, x1 - x0)

    for y in xrange(y0, y1):
        for x in xrange(x0, x1):
            if ((y >= 0) and (y < topo_dims[1]) and
                (x >= 0) and (x < topo_dims[0])):
                subdata[y - y0][x - x0] = data[y][x]

    # Mask array to hide NaNs
    masked = np.ma.masked_invalid(subdata)

    # All done
    return masked

class PlotMap(object):
    def __init__(self):
        return

    def plot(self, plottitle, plotregion, topo, coastal, border,
             fault, sta, map_prefix, hypo_lat=None, hypo_lon=None):
        """
        Produce the plot
        """

        # Read in topo data
        topo_points = read_topo(topo, plotregion)

        # Read in fault data
        fault_x, fault_y = read_fault(fault)

        # Read in station data
        sta_x, sta_y = read_stations(sta)

        # Read coastlines
        coast_x, coast_y = read_coastal(coastal, plotregion)

        # Read borders
        bord_x, bord_y = read_coastal(border, plotregion)

        # Set plot dims
        pylab.gcf().set_size_inches(6, 6)
        pylab.gcf().clf()

        # Adjust title y-position
        t = pylab.title(plottitle, size=12)
        t.set_y(1.09)

        # Setup color scale
        cmap = cm.gist_earth
        norm = mcolors.Normalize(vmin=-1000.0, vmax=3000.0)

        # Plot basemap
        pylab.imshow(topo_points, cmap=cmap, norm=norm,
                     extent=plotregion, interpolation='nearest')

        # Freeze the axis extents
        pylab.gca().set_autoscale_on(False)

        # Plot coast lines
        for i in xrange(0, len(coast_x)):
            pylab.plot(coast_x[i], coast_y[i], linestyle='-', color='0.5')

        # Plot borders
        for i in xrange(0, len(bord_x)):
            pylab.plot(bord_x[i], bord_y[i], linestyle='-', color='0.75')

        # Plot fault trace
        pylab.plot(fault_x, fault_y, linestyle='-', color='k')

        # Plot stations
        pylab.plot(sta_x, sta_y, marker='o', color='r', linewidth=0)

        # Plot hypocenter
        if hypo_lat is not None and hypo_lon is not None:
            hypo_lat = [hypo_lat]
            hypo_lon = [hypo_lon]
            pylab.plot(hypo_lon, hypo_lat, marker='*',
                       markersize=12, color='y', linewidth=0)

        # Set degree formatting of tick values
        majorFormatter = FormatStrFormatter(u'%.1f\u00b0')
        pylab.gca().xaxis.set_major_formatter(majorFormatter)
        pylab.gca().yaxis.set_major_formatter(majorFormatter)

        # Turn on ticks for both sides of axis
        for tick in pylab.gca().xaxis.get_major_ticks():
            tick.label1On = True
            tick.label2On = True
        for tick in pylab.gca().yaxis.get_major_ticks():
            tick.label1On = True
            tick.label2On = True

        # Set font size
        for tick in pylab.gca().get_xticklabels():
            tick.set_fontsize(8)
        for tick in pylab.gca().get_yticklabels():
            tick.set_fontsize(8)

        print("==> Creating Plot: %s.png" % (map_prefix))
        pylab.savefig('%s.png' % (map_prefix), format="png",
                      transparent=False, dpi=plot_config.dpi)
        pylab.close()

def usage():
    """
    Outputs instructions on how to use the program
    """
    print("usage: %s <N bound> <S bound> <W bound> <E bound> " %
          (sys.argv[0]) +
          "<fault file> <station file> <map prefix> <topo> <coastal> <border>")

if __name__ == '__main__':
    if (len(sys.argv) != 11):
        usage()
        sys.exit(1)

    PLOTTITLE = 'Fault Trace with Stations'
    PLOTREGION = [float(sys.argv[3]), float(sys.argv[4]),
                  float(sys.argv[2]), float(sys.argv[1])]
    FAULT = sys.argv[5]
    STA = sys.argv[6]
    MAP_PREFIX = sys.argv[7]
    TOPO = sys.argv[8]
    COASTAL = sys.argv[9]
    BORDER = sys.argv[10]

    PLOTTER = PlotMap()
    PLOTTER.plot(PLOTTITLE, PLOTREGION, TOPO, COASTAL, BORDER,
                 FAULT, STA, MAP_PREFIX)
    sys.exit(0)
