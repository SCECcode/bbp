#!/usr/bin/env python3
"""
BSD 3-Clause License

Copyright (c) 2023, University of Southern California
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

# This is PYTHON port of Walter Imperatori GeoBB_srf.m matlab code
# Function to calculate stations position in CompSyn geometry convention.
# It also calculates stations and hypocenter positions in BBTool geometry.
# Then it provides some parameters useful for subsequent analysis, such
# as different type of station-to-fault distances.
# INPUT:
#
#        stat_file   - file with coordinates (lon/lat) of the stations
#        filename    - srf input filename
#        extended    - flag for extended fault computations ('y' or 'n')
#
# OUTPUT:
#        par    - array with several parameters:
#                 par.bbextent     - "far-side" for the BBTool code
#                 par.hypo         - hypocentral coordinates for BBTool
#                 par.Mw           - moment magnitude
#                 par.mecha        - general mechanism
#                 par.subfault     - subfaults coordinates
#
#
# Remarks: current version work only with positive dip
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys
import math
import numpy as np
import pyproj

# Import Broadband modules
import bband_utils
from station_list import StationList

class GeoBBSRF(object):

    def __init__(self, sim_id=0):
        """
        Initialize class variables
        """

        self.station_file = ''
        self.sim_id = sim_id
        self.mw = 0.0
        self.mecha = ''
        self.f_len = []
        self.f_width =[]
        self.f_dip = []
        self.f_strike = []
        self.f_depth = []
        self.rake = 0.0
        self.hyp = []
        self.fault_maxlon = 0.0
        self.fault_minlon = 0.0
        self.fault_maxlat = 0.0
        self.fault_minlat = 0.0
        self.corn = []
        self.projobj = None
        self.extended = 'n'
        self.n_total_cels = 0

    def geo2cart(self, lon, lat, lonmin=0, latmin=0):
        """
        # Function to convert geographic coordinates into cartesian ones
        # Transversal Mercator projection is used
        #
        # INPUT:
        # lon,lat       -> input vectors with geographic coordinates
        # lonmin,latmin -> lower-left-corner (on x-y plane)
        #                  of plot box (used for shifting)
        # OUTPUT:
        # [x,y]         -> computed cartesian coordinates (in meters)
        """

        # m.geoid= [6378135 0.08181881]
        # print "proj obj in geo2cart", self.projobj
        [a, b] = self.projobj(lon, lat)
        # print "[a,b] = self.projobj(lat,lon)", a,b
        if lonmin == 0:
            a0 = 0
            b0 = 0
        else:
            [a0, b0] = self.projobj(lonmin, latmin)

        # shift of [a0,b0]
        x = (a - a0)
        y = (b - b0)
        return [x, y]

    def translation_matrix(self, direction):
        T = np.identity(4)
        T[:3, 3] = direction[:3]
        return T

    def read_srf_line(self, srf_fh):
        """
        This function returns the next line from a SRF file,
        skipping lines with comments as needed
        """
        while True:
            line = srf_fh.readline()
            if len(line) == 0:
                # End of file
                break
            if not line.strip().startswith("#"):
                break
        return line

    def read_srf(self, srffile):
        if os.path.exists(srffile):
            srff = open(srffile, 'r')
        else:
            raise bband_utils.ParameterError("Missing SRF file (%s)!" %
                                             (srffile))

        # Check number of planes
        line = self.read_srf_line(srff).strip()
        version = float(line.split()[0])
        line = self.read_srf_line(srff).strip()
        tokens = line.split()
        # Make sure we have a valid SRF file
        if len(tokens) != 2:
            raise bband_utils.ProcessingError("Invalid SRF file (%s)!" %
                                              (srffile))
        # Number of planes
        planes = int(tokens[1])

        # Read fault data for each fault plane
        for plane in range(0, planes):

            tokens = self.read_srf_line(srff).strip().split()
            if len(tokens) != 6:
                raise bband_utils.ProcessingError("Invalid SRF file (%s)!" %
                                                  (srffile))
            self.f_len.append(float(tokens[4]))
            self.f_width.append(float(tokens[5]))
            tokens = self.read_srf_line(srff).strip().split()
            if len(tokens) != 5:
                raise bband_utils.ProcessingError("Invalid SRF file (%s)!" %
                                                  (srffile))
            self.f_strike.append(float(tokens[0]))
            self.f_dip.append(float(tokens[1]))
            self.f_depth.append(float(tokens[2]))

        # Now read data
        M0 = 0.0
        lon = []
        lat = []
        dep = []
        tinit = []
        rake = []

        self.n_total_cels = 0
        for plane in range(0, planes):
            tokens = self.read_srf_line(srff).strip().split()
            if len(tokens) != 2:
                raise bband_utils.ProcessingError("Invalid SRF file (%s)!" %
                                                  (srffile))
            n_cels = int(tokens[1])
            self.n_total_cels = self.n_total_cels + n_cels

            for i in range(0, n_cels):
                tokens = self.read_srf_line(srff).strip().split()
                if version == 1.0 and len(tokens) != 8:
                    raise bband_utils.ProcessingError("Invalid SRF version 1 "
                                                      "file (%s)!" %
                                                      (srffile))
                if version == 2.0 and len(tokens) != 10:
                    raise bband_utils.ProcessingError("Invalid SRF version 2 "
                                                      "file (%s)!" %
                                                      (srffile))
                lon.append(float(tokens[0]))
                lat.append(float(tokens[1]))
                dep.append(float(tokens[2]))
                area = float(tokens[5])
                tinit.append(float(tokens[6]))
                tokens = self.read_srf_line(srff).strip().split()
                if len(tokens) != 7:
                    raise bband_utils.ProcessingError("Invalid SRF file (%s)!" %
                                                      (srffile))
                rake.append(float(tokens[0]))
                slip1 = float(tokens[1])
                slip2 = float(tokens[3])
                slip3 = float(tokens[5])
                nt1 = float(tokens[2])
                nt2 = float(tokens[4])
                nt3 = float(tokens[6])
                if nt1 > 0:
                    for k in range(0, int(math.ceil(nt1 / 6.0))):
                        token = self.read_srf_line(srff)
                        if token == "":
                            raise bband_utils.ProcessingError("Invalid SRF "
                                                              "file (%s)!" %
                                                              (srffile))
                if nt2 > 0:
                    for k in range(0, int(math.ceil(nt2 / 6.0))):
                        token = self.read_srf_line(srff)
                        if token == "":
                            raise bband_utils.ProcessingError("Invalid SRF "
                                                              "file (%s)!" %
                                                              (srffile))
                if nt3 > 0:
                    for k in range(0, int(math.ceil(nt3 / 6.0))):
                        token = self.read_srf_line(srff)
                        if token == "":
                            raise bband_utils.ProcessingError("Invalid SRF "
                                                              "file (%s)!" %
                                                              (srffile))

                M0 = M0 + (area *
                           (math.sqrt(slip1**2 + slip2**2 + slip3**2))*3)*(10**11)

        # Close SRF file
        srff.close()

        # Figure out lowest tinit
        np_tinit = np.array(tinit)
        tinit_index = []
        [tinit_index] = np.nonzero(np_tinit == np.min(np_tinit))

        # and now find the hypocenter
        hyp = []
        hyp.append(lon[tinit_index[0]])
        hyp.append(lat[tinit_index[0]])
        hyp.append(dep[tinit_index[0]])
        self.hyp = hyp

        # Find mw
        mw = (math.log10(M0)) / 1.5 - 10.73
        self.mw = mw

        # Find rake (i.e. mechanism)
        rake_ave = np.mean(np.array(rake))
        self.rake = rake_ave
        if rake_ave > 45 and rake_ave < 135:
            mecha = 'rs'
        if rake_ave >= 135 and rake_ave < 225:
            mecha = 'ss'
        if rake_ave >= 225 and rake_ave < 315:
            mecha = 'ns'
        if rake_ave <= 45 or rake_ave >= 315:
            mecha = 'ss'
        self.mecha = mecha
        # print "mw: %f, rake_ave: %f, mecha:%s"%(mw, rake_ave, mecha)

        if self.extended == 'y':
            # Find fault corners
            i = []
            j = []
            k = []
            np_dep = np.array(dep)
            # print np.min(np_dep), len(np_dep)
            [i] = np.nonzero(np_dep == min(dep))

            if len(i) < 1:
                raise bband_utils.ProcessingError("Invalid SRF file: %s\n" %
                                                  (srffile) +
                                                  "len(i)=%d Failed to " %
                                                  len(i) + "calculate " +
                                                  "extended fault data!")
            np_lon = np.array(lon)
            [j] = np.nonzero(np_lon[i] == np.min(np_lon[i]))

            if len(j) < 1:
                raise bband_utils.ProcessingError("Invalid SRF file: %s\n" %
                                                  (srffile) +
                                                  "len(j)=%d Failed to " %
                                                  len(j) + "calculate " +
                                                  "extended fault data!")
            np_lat = np.array(lat)
            [k] = np.nonzero(np_lat[j] == np.min(np_lat[j]))

            if len(k) == 1:
                lat_index = k[0]
                self.corn = [lon[lat_index], lat[lat_index]]
                # print "Corner:", self.corn
            else:
                raise bband_utils.ProcessingError("Invalid SRF file: %s\n" %
                                                  (srffile) +
                                                  "len(k)=%d Failed to " %
                                                  len(k) + "calculate " +
                                                  "extended fault data!")

            self.fault_maxlon = np.max(np_lon)
            self.fault_minlon = np.min(np_lon)
            self.fault_maxlat = np.max(np_lat)
            self.fault_minlat = np.min(np_lat)
            #print ("fault_maxlon: %f, fault_minlon: %f, fault_maxlat: %f, " %
            #       (self.fault_maxlon,self.fault_minlon,self.fault_maxlat) +
            #       "fault_minlat: %f" % (self.fault_minlat))
            # par.mw = mw; par.mecha = mecha;
        return 0

    def write_xyz_srf(self, srf_file, xyz_srf_file, T3M):
        """
        Reads the SRF file and converts its lat/lon to XYZ format
        using self.latmin and self.lonmin as reference points. The T3M
        matrix is used to shift the coordinates to the new reference
        point calculated in the run function.
        """
        # Open files
        infile = open(srf_file, 'r')
        outfile = open(xyz_srf_file, 'w')

        # Pick up version number from SRF file
        version = float(infile.readline().strip().split()[0])
        # Now go back to the start
        infile.seek(0)

        # Copy lines
        while True:
            line = infile.readline()
            # Cannot mix for line in infile with readline...
            if line is None:
                break
            # Don't copy comments
            if line.strip().startswith("#"):
                continue
            # Until we find the plane line
            if line.find("PLANE") >= 0:
                tokens = line.strip().split()
                outfile.write("PLANE %d %d\n" %
                              (int(tokens[1]), self.n_total_cels))
                break

            outfile.write(line)

        tokens = line.strip().split()

        # Make sure we have a valid SRF file
        if len(tokens) != 2:
            raise bband_utils.ProcessingError("Invalid SRF file (%s)!" %
                                              (srf_file))
        planes = int(tokens[1])

        for plane in range(0, planes):
            line = self.read_srf_line(infile)
            tokens = line.strip().split()
            if len(tokens) != 6:
                raise bband_utils.ProcessingError("Invalid SRF file (%s)!" %
                                                  (srf_file))
            # Convert to XYZ
            lon = float(tokens[0])
            lat = float(tokens[1])
            [x_cart, y_cart] = self.geo2cart(lon, lat, self.min_lon, self.min_lat)
            tmp = T3M * np.mat([x_cart, y_cart, 0, 1]).transpose()
            tokens[0] = str(float(tmp[0]))
            tokens[1] = str(float(tmp[1]))
            outfile.write(" %s\n" % "   ".join(tokens))

            # Copy next line without any changes
            line = self.read_srf_line(infile)
            outfile.write(line)

        for plane in range(0, planes):
            line = self.read_srf_line(infile)

            # Figure out how many cells
            tokens = line.strip().split()
            if len(tokens) != 2:
                raise bband_utils.ProcessingError("Invalid SRF file (%s)!" %
                                                  (srf_file))
            outfile.write(line)
            n_cels = int(tokens[1])

            # Go through each cell
            for i in range(0, n_cels):
                tokens = self.read_srf_line(infile).strip().split()
                if version == 1.0 and len(tokens) != 8:
                    raise bband_utils.ProcessingError("Invalid SRF version 1 "
                                                      "file (%s)!" %
                                                      (srffile))
                if version == 2.0 and len(tokens) != 10:
                    raise bband_utils.ProcessingError("Invalid SRF version 2 "
                                                      "file (%s)!" %
                                                      (srffile))
                lon = float(tokens[0])
                lat = float(tokens[1])
                [x_cart, y_cart] = self.geo2cart(lon, lat,
                                                 self.min_lon, self.min_lat)
                tmp = T3M * np.mat([x_cart, y_cart, 0, 1]).transpose()
                tokens[0] = str(float(tmp[0]))
                tokens[1] = str(float(tmp[1]))
                outfile.write(" %s\n" % "   ".join(tokens))
                line = self.read_srf_line(infile)
                tokens = line.strip().split()
                if len(tokens) != 7:
                    raise bband_utils.ProcessingError("Invalid SRF file (%s)!" %
                                                      (srf_file))
                nt1 = float(tokens[2])
                nt2 = float(tokens[4])
                nt3 = float(tokens[6])
                # Write line
                outfile.write(line)
                if nt1 > 0:
                    for k in range(0, int(math.ceil(nt1 / 6.0))):
                        token = self.read_srf_line(infile)
                        if token == "":
                            raise bband_utils.ProcessingError("Invalid SRF "
                                                              "file (%s)!" %
                                                              (srf_file))
                        outfile.write(token)
                if nt2 > 0:
                    for k in range(0, int(math.ceil(nt2 / 6.0))):
                        token = self.read_srf_line(infile)
                        if token == "":
                            raise bband_utils.ProcessingError("Invalid SRF "
                                                              "file (%s)!" %
                                                              (srf_file))
                        outfile.write(token)
                if nt3 > 0:
                    for k in range(0, int(math.ceil(nt3 / 6.0))):
                        token = self.read_srf_line(infile)
                        if token == "":
                            raise bband_utils.ProcessingError("Invalid SRF "
                                                              "file (%s)!" %
                                                              (srf_file))
                        outfile.write(token)

        # All done, close input and output files
        infile.close()
        outfile.close()

    def run(self, slo, coord_out_file, fault_out_file,
            srf_file, xyz_srf_file, extended='y'):
        """
        Reads the SRF file and extracts needed parameters for BBToolbox
        """
        self.extended = extended
        # Read SRF file
        self.read_srf(srf_file)

        # Read station list
        if slo is None:
            raise bband_utils.ParameterError("Cannot open station list")

        site_list = slo.get_station_list()
        slon = []
        slat = []
        site = []
        for sites in site_list:
            slon.append(float(sites.lon))
            slat.append(float(sites.lat))
            site.append(sites.scode)
        n_stat = len(site_list)

        # find minimum/maximum values for stations latitude and longitude
        min_stat_lon = np.min(np.array(slon))
        max_stat_lon = np.max(np.array(slon))
        min_stat_lat = np.min(np.array(slat))
        max_stat_lat = np.max(np.array(slat))

        # find minimum/maximum absolute values
        self.min_lon = min([min_stat_lon, self.hyp[0]])
        self.min_lat = min([min_stat_lat, self.hyp[1]])
        self.max_lon = max([max_stat_lon, self.hyp[0]])
        self.max_lat = max([max_stat_lat, self.hyp[1]])

        if extended == 'y':
            self.min_lon = min([self.min_lon, self.fault_minlon])
            self.min_lat = min([self.min_lat, self.fault_minlat])
            self.max_lon = max([self.max_lon, self.fault_maxlon])
            self.max_lat = max([self.max_lat, self.fault_maxlat])

        # compute average values for projection transformation (center
        # of transformation)
        ave_lon = (self.max_lon + self.min_lon) / 2.0
        ave_lat = (self.max_lat + self.min_lat) / 2.0

        self.projobj = pyproj.Proj(proj='tmerc', lon_0=ave_lon, lat_0=ave_lat,
                                   k=0.001, ellps='WGS72')#, e=0.08181881)

        # transform GEOGRAPHIC coordinates into CARTESIAN ones (unit: Km)
        # (origin shifted at: min_lon/min_lat)
        x_stat = []
        y_stat = []
        for i in range(0, n_stat):
            [x, y] = self.geo2cart(slon[i], slat[i], self.min_lon, self.min_lat)
            x_stat.append(x)
            y_stat.append(y)

        [hypo_x_cart, hypo_y_cart] = self.geo2cart(self.hyp[0], self.hyp[1],
                                                   self.min_lon, self.min_lat)
        [x_min, y_min] = self.geo2cart(self.min_lon, self.min_lat,
                                       self.min_lon, self.min_lat)
        [x_max, y_max] = self.geo2cart(self.max_lon, self.max_lat,
                                       self.min_lon, self.min_lat)

        if extended == 'y':
            [corn_x, corn_y] = self.geo2cart(self.corn[0], self.corn[1],
                                             self.min_lon, self.min_lat)

        # ** objects coordinates in BROAD-BAND TOOLBOX code reference system **
        # NOTE: here X and Y are intented as standard cartesian axes
        #
        # add a "security" factor (in Km) to extremes coordinates: this
        # represents
        # ray-tracing spatial domain
        x_max = x_max + 10
        x_min = x_min - 10
        y_max = y_max + 10
        y_min = y_min - 10

        if extended == 'y':
            corn_x = corn_x + 10
            corn_y = corn_y + 10

        # % matrix to translate to new minimum values (minimum is [0 0])
        T3 = self.translation_matrix([-x_min, -y_min, 0])
        T3M = np.mat(T3)
        # % shift opposite corner in x-y
        tmp = T3M * np.mat([x_max, y_max, 0, 1]).transpose()
        bbextent = [float(tmp[0]), float(tmp[1])]

        # Write station coord_out file
        BBx_stat = []
        BBy_stat = []

        sfile = open(coord_out_file, 'w')
        for i in range(0, n_stat):
            tmp = T3M * np.mat([x_stat[i], y_stat[i], 0, 1]).transpose()
            BBx_stat.append(float(tmp[0]))
            BBy_stat.append(float(tmp[1]))
            sfile.write("%8.4f\t%8.4f\n" % (BBx_stat[i], BBy_stat[i]))
        sfile.close()

        tmp = T3M * np.mat([hypo_x_cart, hypo_y_cart, 0, 1]).transpose()
        BBx_hypo = float(tmp[0])
        BBy_hypo = float(tmp[1])

        # Write a XYZ SRF file
        self.write_xyz_srf(srf_file, xyz_srf_file, T3M)

        # Compute number of subfaults (only INTEGER numbers!!)
        if extended == 'y':
            # Open extended_fault file
            ffile = open("%s.tmp" % (fault_out_file), 'w')

            # Write number of planes and subfaults
            ffile.write("%8i     %8i\n" % (len(self.f_len), 0))

            # < put the 'where' point at the axis origin >
            T1 = self.translation_matrix([-corn_x, -corn_y, 0])
            T1M = np.mat(T1)

            # < align the fault plane to the north axis >
            c = np.cos(np.pi / 180 * self.f_strike[0])
            s = np.sin(np.pi / 180 * self.f_strike[0])
            T2 = np.array([[c, -s, 0, 0], [s, c, 0, 0],
                           [0, 0, 1, 0], [0, 0, 0, 1]])
            T2M = np.mat(T2)

            # inverse matrixes for inverse transform
            T1M_inv = np.linalg.inv(T1M)
            T2M_inv = np.linalg.inv(T2M)

            # Add subfaults from all planes
            total_subfaults = 0

            for plane in range(0, len(self.f_len)):
                num_strike = int(round(self.f_len[plane]))
                num_dip = int(round(self.f_width[plane]))

                # ...and their coordinates
                sub_coord = []
                for j in range(0, num_dip):
                    sub_coord.append([])
                    for i in range(0, num_strike):
                        x_term = (0 + 0.5 *
                                  math.cos(math.radians(self.f_dip[plane])) +
                                  j * math.cos(math.radians(self.f_dip[plane])))
                        y_term = (0 + 0.5 + i)
                        z_term = (self.f_depth[plane] +
                                  0.5 *
                                  math.sin(math.radians(self.f_dip[plane])) +
                                  j * math.sin(math.radians(self.f_dip[plane])))
                        tmp = (T1M_inv * T2M_inv *
                               np.mat([x_term, y_term, 0, 1]).transpose())
                        sub_coord[j].append([float(tmp[0]),
                                             float(tmp[1]),
                                             z_term])

                # Now let's make this into a 1km grid
                data = []
                for j in range(0, num_dip):
                    for i in range(0, num_strike):
                        point = sub_coord[j][i]
                        point[0] = round(float(point[0]))
                        point[1] = round(float(point[1]))
                        point[2] = round(float(point[2]))
                        if (point[0], point[1], point[2]) in data:
                            continue
                        data.append((point[0], point[1], point[2]))

                # Write the extended_fault file
                total_subfaults = total_subfaults + len(data)
                ffile.write('%8i     %8i\n' % (plane + 1, len(data)))
                for point in data:
                    ffile.write('%8.4f %8.4f %8.4f\n' % (point[0],
                                                         point[1],
                                                         point[2]))
            # All done, close file
            ffile.close()

            # Nowm rewrite header
            ffile_in = open("%s.tmp" % (fault_out_file), 'r')
            ffile_out = open(fault_out_file, 'w')
            tokens = ffile_in.readline().strip().split()
            ffile_out.write("%8i     %8i\n" %
                            (int(tokens[0]), total_subfaults))
            # Copy rest of file
            for line in ffile_in:
                ffile_out.write(line)
            ffile_in.close()
            ffile_out.close()

        # Write param file
        fname = '%s.param' % (coord_out_file)
        pfile = open(fname, 'w')
        pfile.write('bbextent: %f %f\n' % (bbextent[0], bbextent[1]))
        pfile.write('bbhypo: %f %f %f\n' % (BBx_hypo, BBy_hypo, self.hyp[2]))
        pfile.write('bbMw: %f\n' % (self.mw))
        pfile.write('mecha: %s\n' % (self.mecha))
        pfile.close()

        return 0

if __name__ == '__main__':
    PROG_BASE = os.path.basename(sys.argv[0])
    if len(sys.argv) != 5:
        print("Usage: %s " % (PROG_BASE) +
              "station_file coord_file fault_file srf_file")
        sys.exit(1)
    print("Testing Module: %s" % (PROG_BASE))
    STAT_LIST = StationList(sys.argv[1])
    ME = GeoBBSRF()
    ME.run(STAT_LIST, sys.argv[2], sys.argv[3], 'n', sys.argv[4], 'y')
    sys.exit(0)
