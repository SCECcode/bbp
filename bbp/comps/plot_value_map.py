#!/usr/bin/env python
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
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys
import math
import struct
import numpy as np
import matplotlib as mpl
mpl.rcParams['font.size'] = 10.
import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import matplotlib.mlab as mlab
#from matplotlib.patches import Circle
from matplotlib import lines
from pyproj import Proj

# Import Broadband modules
import gislib as GS
import validation_cfg
from station_list import StationList
from install_cfg import InstallCfg
from Projection import Projection
import fault_utils
from extract_boundaries import ExtractBoundaries
import plot_config

# Constants
COLOR_CITY = 'r'
COLOR_HYPO = 'g'
COLOR_STA = 'c'
COLOR_TXT = 'k'

class PlotValueMap(object):

    def __init__(self, station_file, sim_id=0, hypo=None):
        self.station_file = station_file
        self.sim_id = sim_id
        self.install = InstallCfg.getInstance()
        self.coast_file = os.path.join(self.install.A_PLOT_DATA_DIR,
                                       "cali_coastline.mapgen")
        if not os.path.isfile(self.coast_file):
            self.coast_file = ""
        self.value = "GOF"
        self.stats = []
        self.dx = 500.0 #100 mts grid resolution
        self.spacing = [self.dx, self.dx]
        self.hypo = hypo
        self.dim = []
        self.rbounds = []
        self.nw = []
        self.sw = []
        self.se = []
        self.ne = []
        self.PLOT_MAP_LOC = [0.10, 0.15, 0.8, 0.8]
        self.origin = []
        self.offset = []
        self.x_invert = False
        self.y_invert = False
        self.init_dims()

    def init_dims(self):
        a_stationlist = os.path.join(self.install.A_IN_DATA_DIR,
                                     str(self.sim_id), self.station_file)
        if not os.path.isfile(a_stationlist):
            a_stationlist = os.path.join(os.getcwd(), self.station_file)
            if not os.path.isfile(a_stationlist):
                print("Error (plot_value_map): Unable to locate station file: "
                      "%s!" % self.station_file)
                sys.exit()
        self.station_file = a_stationlist
        print("Using Station File: %s" % (self.station_file))
        # a_statfile = (self.install.A_IN_DATA_DIR +
        #               "/%d/%s"%(self.sim_id,self.station_file))
        slo = StationList(self.station_file)
        site_list = slo.getStationList()
        w_lon = 0.0
        e_lon = 0.0
        n_lat = 0.0
        s_lat = 0.0
        for sites in site_list:
            slon = float(sites.lon)
            slat = float(sites.lat)
            if w_lon == 0.0:
                w_lon = slon
            elif slon < w_lon:
                w_lon = slon
            if e_lon == 0.0:
                e_lon = slon
            elif slon > e_lon:
                e_lon = slon
            if n_lat == 0.0:
                n_lat = slat
            elif slat > n_lat:
                n_lat = slat
            if s_lat == 0.0:
                s_lat = slat
            elif slat < s_lat:
                s_lat = slat
        self.rbounds = [(n_lat + 0.1), (s_lat - 0.1),
                        (e_lon + 0.1), (w_lon - 0.1)]
        print("Region Bounds: ", self.rbounds)

        self.nw = (self.rbounds[3], self.rbounds[0])
        self.sw = (self.rbounds[3], self.rbounds[1])
        self.se = (self.rbounds[2], self.rbounds[1])
        self.ne = (self.rbounds[2], self.rbounds[0])
        self.PLOT_MAP_LOC = [0.10, 0.15, 0.8, 0.8]

        self.origin = self.nw # North - West Corner
        self.x_invert = False
        self.y_invert = True

        rzone = 11
#        if self.region !=  None:
#            if self.region.getName() == "Northern California":
#                rzone = 10
#                print "Region : %s, UTM Zone: %d" % (self.region.getName(), rzone)
#        else:
        print("Region : None, UTM Zone: %d" % (rzone))

        pobj = Proj(proj='utm', zone=rzone, ellps='WGS84')
        self.offset = map(round, pobj(self.origin[0], self.origin[1]))
        #Calculate region dimension in km
        dim_y = math.ceil(GS.get_distance(self.nw, self.sw)) * (1000.0 / self.dx) #(KM*1000/dx)
        dim_x = math.ceil(GS.get_distance(self.sw, self.se)) * (1000.0 / self.dx)
        dim_z = 1.0 * (1000.0 / self.dx) #Just want to plot surface so we use 1 KM for Z
        self.dim = [int(dim_x), int(dim_y), int(dim_z)]
#               print "Self.dx, self.offset, self.dim:", self.dx, self.offset, self.dim

        self.projobj = Projection(self.dx, self.dim, self.offset, "utm", rzone)
        self.build_station_list(self.station_file)
        self.boundfile = self.build_coastline(self.coast_file, self.projobj,
                                              self.offset, self.dx, self.dim,
                                              self.x_invert, self.y_invert)

        return

    def build_coastline(self, mapfile, proj, offsets, dx, dim, x_invert, y_invert):
        if mapfile == "":
            print("Warning (plot_Value_map): "
                  "Missing coast line data! Skipping coast line plot!")
            return

        boundfile = "%s/%d/boundaries.txt" % (self.install.A_TMP_DATA_DIR, self.sim_id)
        #mapfile, outfile, proj, offsets, dx, dim
        print("Mapfile is: %s" % (mapfile))
        prog = ExtractBoundaries(mapfile, boundfile, proj, offsets,
                                 dx, dim, x_invert, y_invert)
        prog.run()
        return boundfile

    def build_station_list(self, station_file):
        work_dir = os.getcwd()
        proj = self.projobj
        sfname = os.path.splitext(os.path.basename(station_file))[0]
        fname = '%s/%s.txt' % (work_dir, sfname)
        sfile = open(fname, 'w')
        stats = []
        # a_statfile = (self.install.A_IN_DATA_DIR +
        #               "/%d/%s"%(self.sim_id,self.station_file))
        slo = StationList(self.station_file)
        site_list = slo.getStationList()

        for sites in site_list:
            slon = float(sites.lon)
            slat = float(sites.lat)
            site = sites.scode
            x, y = proj.get_xy_from_geo(slon, slat)
            # if x < 0 or y >0:
            #     print "Station oob :", slon, slat, x, y
            stat_data = (x, y)
            stats.append(stat_data)
            sfile.write("%-12s\t%f\t%f\t%f\t%f\n" % (site, slon, slat, x, y))
        self.stats = stats

        # Hypo
        if self.hypo != [] and self.hypo != None:
            self.hypo[0], self.hypo[1] = proj.get_xy_from_geo(self.hypo[0],
                                                              self.hypo[1])

        sfile.close()
        return

    def get_plot_points(self, datalist):
        x = []
        y = []
        z = []

        for sx, sy, data in datalist:
            x.append(sx)
            y.append(sy)
            z.append(data)
        # print sx, sy, data

        x = np.array(x)
        y = np.array(y)
        z = np.array(z)
        assert x.ndim == y.ndim == z.ndim == 1 and  len(x) == len(y) == len(z)

        # Define grid and interpolate the rest
        xi = np.linspace(0.0, float(self.dim[0]-1), self.dim[0])
        yi = np.linspace(0.0, float((self.dim[1]-1)*-1), self.dim[1])

        # print "Length of xi, yi:", len(xi), len(yi)
        zi = mlab.griddata(x, y, z, xi, yi)
        # print "Shape of zi:", zi.shape
        # nmask = np.ma.count_masked(zi)
        # if nmask > 0:
        #print("info: griddata: %d of %d points are masked, not interpolated" %
        #      (nmask, zi.size))
        return zi

    def get_boundary_lines(self, boundfile):
        boundaries = []
        poly = []
        readsize = struct.calcsize('iff')

        ip = open(boundfile, 'rb')
        packed = ip.read(readsize)
        while packed != '':
            data = struct.unpack('iff', packed)
            if data[0] == 0:
                if len(poly) > 1:
                    boundaries.append(poly)
                poly = []
            else:
                x = (data[1]) #*float(self.dx)/1000.0)
                y = (data[2]) #*float(self.dx)/1000.0)
                poly.append([x, y])
            packed = ip.read(readsize)

        if len(poly) > 1:
            boundaries.append(poly)
        ip.close()

        return boundaries

    def plot_grid_array(self, fig, loc, points,
                        labels, units, cmap, norm, title,
                        ticks=None, invert_y=True):

        ax = fig.add_axes(loc, frameon=True)
        ax.set_title('%s' % (title))
        ax.set_xlabel('%s (%s)' % (labels[0], units[0]))
        ax.set_ylabel('%s (%s)' % (labels[1], units[1]))

        # print "plot_grid_array: Points.shape", points.shape
        ax.imshow(points, cmap=cmap, norm=norm, interpolation='nearest', alpha=0.8)

        dims = points.shape

        # Setup custom axis
        ax.set_xlim(0, dims[1])
        ax.set_ylim(0, dims[0])
        if invert_y:
            ax.invert_yaxis()

        if ticks != None:
            plt.xticks(ticks[0][0], ticks[0][1])
            plt.yticks(ticks[1][0], ticks[1][1])

        return 0

    def plot_colorbar(self, loc, value_type, value_units, cmap, norm,
                      value_min, value_max, orient):
        cax = plt.axes(loc)
        # Compute ticks
        ticks = []
        num_ticks = 10
        diff = (value_max - value_min) / float(num_ticks)
        # print "Colorbar diff: ", diff
        for i in xrange(0, num_ticks + 1):
            # print (value_min + (i * diff))
            ticks.append(value_min + (i * diff))


        cax.set_title("%s (%s)" % (value_type, value_units))
        mpl.colorbar.ColorbarBase(cax, cmap=cmap,
                                  norm=norm,
                                  ticks=ticks,
                                  format='%2.1f',
                                  orientation=orient)
        return 0

    def plot_site(self, fig, ax, label, point1, point2, sitetype, r):
        grid_x = [point1]
        grid_y = [point2]
        if sitetype == 'sta':
            fig.plot(grid_x, grid_y, color=COLOR_STA, marker='.')
        elif sitetype == 'hypo':
        #print "Plotting hypo at %d, %d" % (grid_x[0], grid_y[0])
            fig.plot(grid_x, grid_y, color=COLOR_HYPO, marker='*',
                     markersize=10, markeredgecolor='red')

            fig.text(point1, point2+2,
                     '%s' % (label), color=COLOR_HYPO,
                     horizontalalignment='center')
                     # cir = Circle( (point1,point2), radius=r, \
                     #    edgecolor=COLOR_HYPO, facecolor=COLOR_TXT, alpha=0.5)
                     #    ax.add_patch(cir)
                     #    fig.text(point1, point2+r, \
                     #    '%s' % ('Hypo'), color=COLOR_HYPO, \
                     #    horizontalalignment='center')
        elif sitetype == 'txt':
            fig.text(point1, point2,
                     '%s' % (label), color=COLOR_TXT,
                     horizontalalignment='center')
        else:
            fig.plot(grid_x, grid_y, color=COLOR_CITY, marker='.')
            fig.text(point1, point2+2,
                     '%s' % (label), color=COLOR_CITY,
                     horizontalalignment='center')
        return

    def plot_boundaries(self, fig, boundfile):

        boundaries = []
        boundaries = self.get_boundary_lines(boundfile)

        for poly in boundaries:
            prevline = []
            for line in poly:
                # print line
                if prevline != []:
                    grid_x = [prevline[0], line[0]]
                    grid_y = [prevline[1], line[1]]
                    # print "Line:", prevline, line
                    cline = lines.Line2D(grid_x, grid_y, lw=1., color='k')
                    # fig.plot(grid_x, grid_y, '-', color='k')
                    # cline.set_clip_on(False)
                    # cline.set_axes(fig)
                    fig.add_line(cline)
                prevline = line
        return

    def do_plot(self, datalist, title, outfile):
        plot_x_size = 8
        plot_y_size = 6
        fig = plt.figure(figsize=(plot_x_size, plot_y_size))

        # Get plot points
        points = self.get_plot_points(datalist)
        if points is None:
            print("ERROR (plot_valu_map): Failed to get plot points!")
            sys.exit(-1)

        value_min = points.min()
        value_max = points.max()
        # if the values are between 0-100, used fixed color bar
        if value_min >= 0:
            value_min = 0
        if value_max <= 100:
            value_max = 100
        # print "Colorbar range: %f to %f" % (value_min, value_max)

        cmap = cm.gist_rainbow_r
        norm = mcolors.Normalize(vmin=value_min, vmax=value_max)

        # Redefine tick labels
        ticks = [[[], []], [[], []]]
        ticks[0][0] = []
        ticks[0][1] = []
        for i in xrange(0, self.dim[0] * int(self.spacing[0]) + 1, 50000):
            ticks[0][0].append(i / self.spacing[0])
            ticks[0][1].append('%d' % (i / 1000))
        ticks[1][0] = []
        ticks[1][1] = []
        for i in xrange(0, self.dim[1] * int(self.spacing[1]) + 1, 50000):
            ticks[1][0].append(i / self.spacing[1])
            ticks[1][1].append('%d' % (i / 1000))

            # y_scale = self.dim[1]/float(self.dim[0])
            # x_scale = 1.0
            # if (y_scale > 1.0):
            #     x_scale = 1/y_scale
            #     y_scale = 1.0
            # y_scale = y_scale * plot_x_size / float(plot_y_size)
            # print "X-scale: %lf" % (x_scale)
            # print "Y-scale: %lf" % (y_scale)

        origin_x = 0.05
        origin_y = 0.18 #(1.0-(0.9*y_scale))/2 + 0.08
        x_length = 0.85 # * x_scale
        y_length = 0.74 # * y_scale
        loc = [origin_x, origin_y, x_length, y_length]
        # Plot the depth map
        invert_y = True
        self.plot_grid_array(fig, loc, points,
                             ['X', 'Y'], ['km', 'km'],
                             cmap, norm, title,
                             ticks, invert_y)

        #plot Sites
        radius = 1.0 * float(self.dx) / 1000.0
        ax = fig.get_axes()
        for site in self.stats:
            self.plot_site(plt, ax[0], "", site[0], (-1*site[1]), 'sta', radius)

        #plot coast line
        self.plot_boundaries(ax[0], self.boundfile)

        #plot Hypo
        if self.hypo != [] and self.hypo is not None:
            self.plot_site(plt, ax[0], "", self.hypo[0],
                           (-1 * self.hypo[1]), 'hypo', radius * 10.0)

        # Plot colorbar
        value_units = '%'
        c_orient = "horizontal"
        if c_orient == "horizontal":
            cloc = [origin_x + (x_length * 0.1), origin_y - 0.13,
                    (x_length * 0.8), 0.02]
        # else:
        #     if invert_y:
        #         cloc = [(origin_x+x_length),
        #                  origin_y-y_length+ 0.05, 0.02, (y_length*0.8), ]
        #     else:
        #         cloc = [(origin_x+x_length+0.02),
        #                  origin_y+ 0.05, 0.02, (y_length*0.8), ]
        # print cloc

        self.plot_colorbar(cloc, self.value, value_units, cmap, norm,
                           value_min, value_max, c_orient)

        outfile.replace(' ', '_')
        outfile.replace('/', '_')
        print("Saving plot file %s" % (outfile))
        plt.savefig(outfile, dpi=plot_config.dpi)
        plt.show()

        return 0

    def run(self, data, label):
        self.label = label
        self.data = data
        # a_indir = "%s/%d" % (self.install.A_IN_DATA_DIR, self.sim_id)
        a_outdir = os.path.join(self.install.A_OUT_DATA_DIR, str(self.sim_id))
        outfile = "%s/%d_%s_map.png" % (a_outdir, self.sim_id, label)
        title = "Run %d - %s Map" % (self.sim_id, label)
        self.log = "%s/%d/%d.plot_value_map.log" % (self.install.A_OUT_LOG_DIR,
                                                    self.sim_id, self.sim_id)
        if len(data) != len(self.stats):
            print("ERROR (plot_valu_map): Data file length "
                  "(%d) is not equal to number of stations (%d)!" \
                    % (len(data), len(self.stats)))
            sys.exit(-1)
        iindex = 0
        datalist = []
        for value in self.stats:
            dval = (float(value[0]), float(value[1]), data[iindex][0])
            datalist.append(dval)
            iindex += 1
        self.do_plot(datalist, title, outfile)

if __name__ == '__main__':
    STATION_FILE = sys.argv[1]
    DATA = sys.argv[2]
    LABEL = sys.argv[3]
    VNAME = sys.argv[4]
    print("Validation Event - %s" % (VNAME))
    HYPO = []

    if VNAME is not None:
        VAL_OBJ = validation_cfg.VE_EVENTS.get_event_by_name(VNAME)
        SRF_FILE = VAL_OBJ.get_input("GP", "srf")
        HYPO = fault_utils.get_hypocenter(SRF_FILE)

    PVM = PlotValueMap(STATION_FILE, int(sys.argv[5]), HYPO)
    PVM.run(DATA, LABEL)
