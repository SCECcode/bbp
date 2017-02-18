#!/usr/bin/env python
"""
Southern California Earthquake Center Broadband Platform
Copyright 2010-2016 Southern California Earthquake Center

Projection package for converting GEO<->UTM coords
$Id: Projection.py 1730 2016-09-06 20:26:43Z fsilva $
"""
from __future__ import division, print_function

# Basic Python modules
import math
from pyproj import Proj

class Projection(object):
    def __init__(self, dx, dim, offsets, uproj="utm", uzone=11):

        self.valid = False
        self.proj = uproj
        self.dx = dx
        self.dim = dim

        if self.proj == 'utm':
            self.projobj = Proj(proj='utm', zone=uzone, ellps='WGS84')
            if offsets != None:
                self.xoffset = -1.0 * offsets[0]
                self.yoffset = -1.0 * offsets[1]
            else:
                self.xoffset = 0
                self.yoffset = 0

            self.valid = True
        else:
            self.projobj = None

    def is_valid(self):
        return self.valid

    def cleanup(self):
        return

    def get_raw_xy_from_geo(self, lon, lat):
        x, y = self.projobj(lon, lat)
        x = x + self.xoffset
        y = y + self.yoffset
        return (x, y)

    def get_xy_from_geo(self, lon, lat):
        if self.proj == 'utm':
            x, y = self.projobj(lon, lat)
            x = x + self.xoffset
            y = y + self.yoffset
            x = int(x / float(self.dx) + 0.5)
            y = int(y / float(self.dx) + 0.5)
            #print("x: %f, y: %f, self.dim[0]: %f, self.dim[1]: %f" %
            #      (x, y, self.dim[0], self.dim[1]))
            if (x >= self.dim[0]) or (y >= self.dim[1]):
                return (-1.0, -1.0)
            return (x, y)
        else:
            return (-1.0, -1.0)

    def get_geo_from_xy(self, x, y):
        if self.proj == 'utm':
            #print x,y
            lon, lat = self.projobj(x * self.dx - self.xoffset,
                                    y * self.dx - self.yoffset,
                                    inverse=True)
            return (lon, lat)
        else:
            return (-1.0, -1.0)

    def get_dist_km(self, lon1, lat1, lon2, lat2):
        dist = 0.0
        if self.proj == 'CMU':
            x1, y1 = self.projobj(lon1, lat1)
            x2, y2 = self.projobj(lon2, lat2)

            dist = math.sqrt((x2 - x1) * (x2 - x1) +
                             (y2 - y1) * (y2 - y1)) / 1000.0

        return dist
