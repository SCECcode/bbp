#!/usr/bin/env python
"""
Southern California Earthquake Center Broadband Platform
Copyright 2010-2016 Southern California Earthquake Center

##########################################################
#
# Script: ExtractBoundaries
#
# Description: Parse NOAA coastline coordinate files in
# Mapgen format to produce a list of boundary lines in
# simulation coordinates. These lines get saved in a custom
# binary format of 'iff', where:
#
# Field 1: Boolean flag denoting if this is a break
# Field 2: X grid coord
# Field 3: Y grid coord
#
##########################################################
"""

# Basic modules
from __future__ import division, print_function
import sys
import struct

class ExtractBoundaries(object):

    def __init__(self, mapfile, outfile, proj, offsets,
                 dx, dim, x_invert=False, y_invert=False):

        self.valid = False
        self.mapfile = mapfile
        self.outfile = outfile
        self.proj = proj
        self.offsets = offsets
        self.dx = dx
        self.dim = dim
        print("ExtractBoundaries: off, dx, dim, xinv, yinv",
              offsets, dx, dim, x_invert, y_invert)
        self.x_invert = x_invert
        self.y_invert = y_invert

        self.valid = True

    def is_valid(self):
        return self.valid

    def cleanup(self):
        return

    def _parse_mapgen(self, proj):
        print("parsing mapgen file %s" % (self.mapfile))
        boundaries = []
        poly = []

        # Read each line
        ip = open(self.mapfile, 'r')

        for line in ip:
            if line.startswith('# -b'):
                if len(poly) > 1:
#                    print "Appending polygon"
                    boundaries.append(poly)
                poly = []
            else:
                # Convert to grid coord
                tokens = line.split()
                lon = float(tokens[0])
                lat = float(tokens[1])
                x, y = proj.get_xy_from_geo(lon, lat)
                append = False
#                print lon,lat, x, y
                # Save point if in region of interest
                if self.x_invert == True and self.y_invert == True:
                    #origin is NE corner
                    if ((x <= 0) and (abs(x) > self.dim[0]) and \
                        (y <= 0) and (abs(y) > self.dim[1])):
#                        print "Appending point NE"
                        x = -1.0 * x
                        y = -1.0 * y
                        append = True

                elif self.x_invert == True and self.y_invert == False:
                    #origin is SE corner
                    if ((x <= 0) and (abs(x) > self.dim[0]) and \
                        (y >= 0) and (y < self.dim[1])):
#                        print "Appending point SE"
                        x = -1.0 * x
                        append = True

                elif self.x_invert == False and self.y_invert == True:
                    #origin is NW corner
                    if ((x >= 0) and (x < self.dim[0]) and \
                        (y <= 0) and (abs(y) < self.dim[1])):
#                        print "Appending point NW"
#                        lo, la = proj.getGeoFromXY(x,y)
#                        print lo, la
                        y = -1.0 * y
                        append = True

                elif self.x_invert == False and self.y_invert == False:
                    #origin is SW corner
                    if ((x >= 0) and (x < self.dim[0]) and \
                        (y >= 0) and (y < self.dim[1])):
#                        print "Appending point SW"
                        append = True

                if append == True:
#                    print "Appending point"
                    poly.append([x, y])
                else:
                    # Poly is clipped
                    if len(poly) > 1:
#                        print "Appending polygon"
                        boundaries.append(poly)
                        poly = []

        if len(poly) > 1:
#            print "Appending polygon"
            boundaries.append(poly)

        ip.close()
        return boundaries

    def run(self):
        # Parse mapgen file
        boundaries = self._parse_mapgen(self.proj)
        print("Found %d polygons" % (len(boundaries)))

        # Open outfile
        outfile = open(self.outfile, 'w')

        # Write lines to file
        for poly in boundaries:
            for line in poly:
                packed = struct.pack('iff', 1, line[0], line[1])
                outfile.write(packed)

            packed = struct.pack('iff', 0, 0.0, 0.0)
            outfile.write(packed)

        # Close open files
        outfile.close()

        return 0

def usage():
    print("usage: %s <mapfile> <outfile> <proj> <offsets> <dx> <dim>" %
          (sys.argv[0]))
    return

if __name__ == '__main__':
    #mapfile, outfile, proj, offsets, dx, dim
    if len(sys.argv) != 7:
        print("ERROR! Incorrect number of arguments!")
        usage()
        sys.exit(1)

    PROG = ExtractBoundaries(sys.argv[1], sys.argv[2],
                             sys.argv[3], sys.argv[4],
                             sys.argv[5], sys.argv[6])
    sys.exit(PROG.run())
