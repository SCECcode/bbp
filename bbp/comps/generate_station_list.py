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

Generate Station List based on Rrup
"""
from __future__ import division, print_function

# Import Python modules
import optparse
import numpy as np

from pynga.utils import (LonLatToAngleDistance, DistanceToSimpleFaultSurface,
                         EndLocation, SimpleFaultSurface, FaultTraceGen)
#from pynga.utils import *

# Parse command-line options
parser = optparse.OptionParser()
parser.add_option("-g", "--grid-size", dest="grid_size", type="float",
                  help="specify the grid size")
parser.add_option("-n", "--number-stations", dest="number_stations",
                  type="int", help="Number of stations to generate")
parser.add_option("-e", "--equal", dest="equal_distance",
                  action="store_true",
                  help="Place n stations equally distributed around the fault")
parser.add_option("-r", "--radius", dest="station_radius", type="float",
                  help="radius from the rupture to place stations")
parser.add_option("-d", "--distance", dest="distance", type="float",
                  help="distance extend from both ends of rupture")
parser.add_option("-s", "--source", dest="source_file",
                  help="provides the source file to use")
parser.add_option("-o", "--output", dest="output_file",
                  help="output file for the station list")
parser.add_option("-v", "--vs30", dest="vs30",
                  help="vs30 for the stations")
parser.add_option("-p", "--plot_prefix", dest="plot_prefix",
                  help="specifies the plot prefix to use (default: no plots)")
parser.add_option("--hw", "--hanging-wall", dest="hanging_wall",
                  action="store_true",
                  help="place stations on the hanging wall side of the fault"
                  " (default: place on footwall)")
parser.add_option("-a", "--all_around", dest="all_around", action="store_true",
                  help="place stations 360 deg (default: only on foot wall)")
parser.add_option("-j", "--stat-prefix", dest="stat_prefix",
                  help="adds a prefix to the station names")

(options, args) = parser.parse_args()

# Make sure mandatory parameters are provided
if options.grid_size is None:
    parser.print_help()
    parser.error("Must specify grid size!")
if options.number_stations is None:
    parser.print_help()
    parser.error("Must specify number of stations!")
if options.station_radius is None:
    parser.print_help()
    parser.error("Must specify radius to place stations!")
if options.distance is None:
    parser.print_help()
    parser.error("Must specify distance parameter!")
if options.source_file is None:
    parser.print_help()
    parser.error("Must specify source file")
if options.output_file is None:
    parser.print_help()
    parser.error("Must specify output station list file")
if options.vs30 is None:
    parser.print_help()
    parser.error("Must specify Vs30 to be used in station list")

if options.plot_prefix is None:
    plot_prefix = None
else:
    plot_prefix = options.plot_prefix
if options.stat_prefix is None:
    stat_prefix = ''
else:
    stat_prefix = options.stat_prefix
if options.hanging_wall is None:
    options.hanging_wall = False
if options.all_around is None:
    options.all_around = False
if options.equal_distance is None:
    options.equal_distance = False

R0 = float(options.distance)
grid_size = float(options.grid_size)
station_radius = float(options.station_radius)
source_file = options.source_file
output_file = options.output_file
vs30 = int(float(options.vs30))

# Figure out where to place stations, start with default
stations_foot_wall = True
stations_hanging_wall = False

if options.hanging_wall:
    stations_hanging_wall = True
    stations_foot_wall = False

# all_around overrides the hanging wall option
if options.all_around:
    stations_hanging_wall = True
    stations_foot_wall = True

num_stations = int(options.number_stations)
if options.equal_distance:
    # Place stations equally around the rupture
    select_method = "e"
else:
    # Default is randomly-placed stations
    select_method = "r"

# Read source file
cfg_dict = {}
src_file = open(source_file, 'r')
for line in src_file:
    # Read file, line by line
    line = line.strip()
    if line.startswith('#'):
        # Skip comments
        continue
    ml = line.split('=')
    key = ml[0].strip().upper()
    val = ml[1].strip()
    cfg_dict[key] = float(val)
src_file.close()

# Look for the keys that we need
M = cfg_dict['MAGNITUDE']
L = cfg_dict['FAULT_LENGTH']
dl = cfg_dict['DLEN']
W = cfg_dict['FAULT_WIDTH']
dw = cfg_dict['DWID']
ztor = cfg_dict['DEPTH_TO_TOP']
strike = cfg_dict['STRIKE']
rake = cfg_dict['RAKE']
dip = cfg_dict['DIP']
lat0 = cfg_dict['LAT_TOP_CENTER']
lon0 = cfg_dict['LON_TOP_CENTER']
hypoAS = cfg_dict['HYPO_ALONG_STK']
hypoDD = cfg_dict['HYPO_DOWN_DIP']

origin = lon0, lat0    # center
Dims = L, dl, W, dw, ztor
Mech = strike, dip, rake
FaultTrace1, UpperSeisDepth, LowerSeisDepth, AveDip, dl, dw = FaultTraceGen(origin, Dims, Mech)
FaultTrace, FaultSeg, AveStrike = SimpleFaultSurface(FaultTrace1,
                                                     UpperSeisDepth,
                                                     LowerSeisDepth,
                                                     AveDip)

# generate grid of stations (for grid search)  in lon/lat
strike = strike * np.pi / 180.
loc0 = [lon0, lat0, 0.0]

vD = 0.0
dhD = grid_size   # in km (grid size)
hDx = 1.5 * station_radius
hDy = 0.5 * L + R0 + dhD * 2    # along strike

az = strike + np.pi
vector = [az, hDy, vD]
loc1 = EndLocation(loc0, vector)

az = strike + np.pi * 3. / 2.
vector = [az, hDx, vD]
loc11 = EndLocation(loc1, vector)  # new origin

Nx = int(2 * hDx / dhD + 1)
Ny = int(2 * hDy / dhD + 1)

LocY = loc11
Loc2D = []
for _ in range(Ny):
    Loc2D.append(LocY)
    az = strike + np.pi / 2.
    hD = dhD
    vector = [az, hD, vD]
    LocX = LocY
    for _ in range(Nx - 1):
        LocX = EndLocation(LocX, vector)
        Loc2D.append(LocX)
    az = strike
    hD = dhD
    vector = [az, hD, vD]
    LocY = EndLocation(LocY, vector)

Nloc = len(Loc2D)
LocS0 = np.array(Loc2D)
print('Total station locations: %d' % (Nloc))

Rrups = []
Rxs = []
for iloc in range(len(Loc2D)):
    SiteGeom = Loc2D[iloc]
    Rjb, Rrup, Rx = DistanceToSimpleFaultSurface(SiteGeom,
                                                 FaultTrace1,
                                                 UpperSeisDepth,
                                                 LowerSeisDepth,
                                                 AveDip) # in km
    Rrups.append(Rrup)
    Rxs.append(Rx)

if plot_prefix is not None:
    import matplotlib.style
    import matplotlib as mpl
    # Restore Matplotlib 1.x styles if we are using Matplotlib 2.x
    try:
        mpl.style.use('classic')
    except:
        pass
    # Turn off tick label formatting with offsets
    mpl.rcParams['axes.formatter.useoffset'] = False
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    vertsClosed = FaultTrace + [FaultTrace[0]]
    verts1 = np.array(vertsClosed)
    fig = plt.figure(1)
    fig.clf()
    ax = fig.add_subplot(111)
    img = ax.scatter(LocS0[:, 0], LocS0[:, 1],
                     c=np.array(Rrups), s=np.array(Rrups), edgecolor='none')
    ax.plot(verts1[:, 0], verts1[:, 1], 'k')
    ax.set_xlabel('lon', size=8)
    ax.set_ylabel('lat', size=8)
    ax.tick_params(labelsize=7)
    ax.set_title('All station grid and Rrup with domain\n'
                 '%s*%s km^2, and grid size: %s km' %
                 (2 * hDx, 2 * hDy, dhD), size=10)
    cbar = fig.colorbar(img)
    cbar.ax.tick_params(labelsize=8)
    fig.savefig('%s_all_station_rrup_grid_size_%s.pdf' %
                (plot_prefix, '%.2f' % grid_size), format='pdf')

    fig = plt.figure(2)
    fig.clf()
    ax = fig.add_subplot(111)
    img = ax.scatter(LocS0[:, 0], LocS0[:, 1],
                     c=np.array(Rxs), s=abs(np.array(Rxs)), edgecolor='none')
    ax.plot(verts1[:, 0], verts1[:, 1], 'k')
    ax.set_xlabel('lon', size=8)
    ax.set_ylabel('lat', size=8)
    ax.tick_params(labelsize=7)
    ax.set_title('All station grid and Rx', size=10)
    cbar = fig.colorbar(img)
    cbar.ax.tick_params(labelsize=8)
    fig.savefig('%s_all_station_rx_grid_size_%s.pdf' %
                (plot_prefix, '%.2f' % grid_size), format='pdf')

# grid search (allowing +- 0.5km range around Rrup0. Rrup0 is the
# cutoff distance to select stations from the grid)
Rrup0 = station_radius
Loc2D1 = []
Rrups1 = []
Rxs1 = []

azs1 = []   # will be used to rank the index1
for iloc in range(Nloc):
    if Rrup0-dhD / 2. <= Rrups[iloc] <= Rrup0 + dhD / 2.:
        if not stations_hanging_wall and Rxs[iloc] > 0:
            continue
        if not stations_foot_wall and Rxs[iloc] < 0:
            continue

        # compute azimuth
        tmp1, tmp2, tmp3, az1to2 = LonLatToAngleDistance(loc0, Loc2D[iloc],
                                                         CalcRadius=True,
                                                         CalcDist=True,
                                                         Fast=True,
                                                         CalcAzimuth=True,
                                                         Azimuth0to2PI=True)
        azs1.append(az1to2)

        Loc2D1.append(Loc2D[iloc])
        Rrups1.append(Rrups[iloc])
        Rxs1.append(Rxs[iloc])
    else:
        continue

Nloc1 = len(Loc2D1)
print('Stations with distance Rrup=%s:  %d' % (Rrup0, Nloc1))
if Nloc1 < num_stations:
    print('Change your grid size smaller in order to do the subsampling')
    raise ValueError

LocS1 = np.array(Loc2D1)

# Generate random integer (total number would be num_stations)
# between 0 and Nloc
LocS2 = []
index = []
Rrups2 = []
Rxs2 = []
if num_stations > Nloc1:
    num_stations = Nloc1
    LocS2 = LocS1
    Rrups2 = Rrups1
    Rxs2 = Rxs1
    num_output_stations = num_stations
else:
    # Select stations (equally or random)
    if select_method == 'e':
        # equally spaced on the radius
        # sort index based on azimuth (clockwise)
        index1sorted = np.argsort(azs1)
        # Calculate interval
        index = np.linspace(0, (Nloc1-1), num=num_stations)
        for index01 in index:
            index0 = index1sorted[index01]
            LocS2.append(LocS1[index0])
            Rrups2.append(Rrups1[index0])
            Rxs2.append(Rxs1[index0])

    if select_method == 'r':
        # randomly spaced on the radius
        i = 0
        while i < num_stations:
            index0 = np.random.randint(0, Nloc1-1)
            if not index0 in index:
                LocS2.append(LocS1[index0])
                Rrups2.append(Rrups1[index0])
                Rxs2.append(Rxs1[index0])
                index.append(index0)
                i = i + 1
            else:
                continue
    num_output_stations = len(index)

LocS2 = np.array(LocS2)

# write into file
fid = open(output_file, 'w')
fid.write('#Slon\tSlat\tRSN\tVs30(m/s)\n')
for ista in range(num_output_stations):
    fid.write('%s %s %s %d\n' %
              (LocS2[ista, 0], LocS2[ista, 1],
               "sta-%s%04d" % (stat_prefix, ista), vs30))
fid.close()

if plot_prefix is not None:
    fig = plt.figure(3)
    fig.clf()
    ax = Axes3D(fig)
    DepFactor = 1. / 6371 * 180. / np.pi
    DepFactor = 1.0
    ax.plot(verts1[:, 0], verts1[:, 1],
            -verts1[:, 2]*DepFactor, 'k--', label='Fault Surface')
    ax.plot(verts1[:, 0], verts1[:, 1],
            0.0, 'k', label='Fault Surface Projection')
    ax.plot(LocS0[:, 0], LocS0[:, 1], 0.0, 'g+', label='Station Grid')
    ax.text(verts1[3, 0], verts1[3, 1],
            -verts1[3, 2]*DepFactor, 'dip=%s' % ('%.1f' % AveDip))
    ax.set_xlabel('lon', size=8)
    ax.set_ylabel('lat', size=8)
    ax.tick_params(labelsize=7)
    ax.set_zlabel('dep (km)', size=8)
    ax.set_zlim3d([-110, 0])
    ax.plot(LocS1[:, 0], LocS1[:, 1], 0.0, 'b^',
            label='Station with Rrup=~%s km' % Rrup0)

    ax.plot(LocS2[:, 0], LocS2[:, 1], 0.0, 'r*', markersize=12,
            label='Selected Stations (Rrup=~%s km)' %
            (Rrup0))
    lg = ax.legend(loc=0)
    lg.draw_frame(False)
    ltext = plt.gca().get_legend().get_texts()
    plt.setp(ltext, fontsize=8)
    fig.savefig('%s_station_selection_grid_size_%s_Nsta_%s_RrupAt_%s.pdf' %
                (plot_prefix, '%.2f' % grid_size,
                 num_stations, Rrup0), format='pdf')
