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

Song RMG Rupture Generator Multi Segment
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys
import math
import time
import cmath
import random
import numpy as np
from scipy.interpolate import griddata
import warnings

# Import BBP modules
import plot_srf
import bband_utils
from song_rmg_multi_seg_cfg import SongRMGMSCfg
from install_cfg import InstallCfg

# Constants from gen_stats_inp.m
mu = 33e9 # Pa

# Constants from src_stats_mod.mat
mu_ax = [1.5329, 1.4496, 1.7063, 1.2242, 1.5395, 1.6174]
mu_az = [1.2215, 1.0458, 1.0519, 1.0540, 0.8088, 1.0153]
mu_cc = [0.5807, 0.7681, 0.6343]
mu_mu = [115.2695, 2.0327, 111.6469]
mu_sig = [58.8421, 0.7846, 78.9101]
sig_ax = [0.4864, 0.8706, 0.8875, 0.5799, 1.0582, 0.7069]
sig_az = [0.7453, 0.9564, 0.8797, 0.5549, 0.8503, 0.5774]
sig_cc = [0.1231, 0.0874, 0.1164]
sig_mu = [39.7429, 0.2514, 31.9504]
sig_sig = [22.9453, 0.1122, 27.0968]
C = [[1.0000, 0.0351, 0.6673, 0.9803, 0.2617, 0.6657, -0.1094,
      -0.0136, -0.1807, 0.0021, -0.0622, -0.2019, 0.5704, 0.3749,
      0.4942, 0.1264, 0.1764, 0.1810, 0.3585, -0.2776, 0.2169],
     [0.0351, 1.0000, 0.4690, 0.0608, 0.3916, 0.3144, 0.0745,
      0.2339, 0.0862, 0.3514, 0.3348, 0.2345, 0.0086, -0.1889,
      0.0396, -0.2999, -0.2126, -0.1075, -0.1158, 0.3768, -0.3763],
     [0.6673, 0.4690, 1.0000, 0.6760, 0.1915, 0.9213, 0.1151,
      0.3258, 0.0969, 0.4111, 0.4195, 0.0421, 0.2675, 0.0019,
      0.2679, -0.2285, 0.0404, 0.1922, 0.4571, -0.0356, 0.2327],
     [0.9803, 0.0608, 0.6760, 1.0000, 0.2641, 0.6988, -0.0773,
      -0.0106, -0.1464, 0.0092, -0.0255, -0.1432, 0.5238, 0.3484,
      0.4867, 0.1499, 0.1806, 0.2068, 0.3694, -0.1675, 0.2355],
     [0.2617, 0.3916, 0.1915, 0.2641, 1.0000, 0.2469, -0.2244,
      -0.2883, -0.3076, -0.1196, -0.1709, -0.2304, 0.2745, 0.3352,
      0.3204, 0.2161, 0.2028, 0.1481, -0.0627, 0.1207, -0.3495],
     [0.6657, 0.3144, 0.9213, 0.6988, 0.2469, 1.0000, 0.0727,
      0.1566, 0.0476, 0.2599, 0.3412, 0.0322, 0.2397, 0.0994,
      0.2816, -0.0514, 0.1232, 0.2563, 0.5068, -0.0576, 0.2777],
     [-0.1094, 0.0745, 0.1151, -0.0773, -0.2244, 0.0727, 1.0000,
      0.5633, 0.6220, 0.3922, 0.4548, 0.4867, -0.2456, -0.3080,
      -0.2000, -0.3335, -0.2204, -0.2238, 0.1223, 0.1121, 0.2360],
     [-0.0136, 0.2339, 0.3258, -0.0106, -0.2883, 0.1566, 0.5633,
      1.0000, 0.5692, 0.6432, 0.6992, 0.3994, -0.1658, -0.3600,
      -0.1948, -0.5043, -0.2639, -0.2222, 0.1836, 0.0165, 0.2772],
     [-0.1807, 0.0862, 0.0969, -0.1464, -0.3076, 0.0476, 0.6220,
      0.5692, 1.0000, 0.3915, 0.5456, 0.5740, -0.2986, -0.3650,
      -0.3755, -0.4064, -0.3091, -0.2343, -0.0003, 0.0209, 0.0343],
     [0.0021, 0.3514, 0.4111, 0.0092, -0.1196, 0.2599, 0.3922,
      0.6432, 0.3915, 1.0000, 0.7637, 0.5002, -0.0948, -0.2817,
      -0.0153, -0.4343, -0.1754, -0.0763, 0.1807, 0.0773, 0.1505],
     [-0.0622, 0.3348, 0.4195, -0.0255, -0.1709, 0.3412, 0.4548,
      0.6992, 0.5456, 0.7637, 1.0000, 0.6543, -0.2287, -0.3714,
      -0.2375, -0.4589, -0.2221, -0.1780, 0.1993, 0.0458, 0.2033],
     [-0.2019, 0.2345, 0.0421, -0.1432, -0.2304, 0.0322, 0.4867,
      0.3994, 0.5740, 0.5002, 0.6543, 1.0000, -0.2518, -0.3249,
      -0.2753, -0.3914, -0.2784, -0.2883, 0.0099, 0.1892, 0.0639],
     [0.5704, 0.0086, 0.2675, 0.5238, 0.2745, 0.2397, -0.2456,
      -0.1658, -0.2986, -0.0948, -0.2287, -0.2518, 1.0000, 0.5875,
      0.5187, 0.2251, 0.3262, 0.1950, 0.2035, -0.1137, -0.0447],
     [0.3749, -0.1889, 0.0019, 0.3484, 0.3352, 0.0994, -0.3080,
      -0.3600, -0.3650, -0.2817, -0.3714, -0.3249, 0.5875, 1.0000,
      0.4580, 0.4671, 0.4855, 0.2041, 0.0488, -0.2098, -0.0927],
     [0.4942, 0.0396, 0.2679, 0.4867, 0.3204, 0.2816, -0.2000,
      -0.1948, -0.3755, -0.0153, -0.2375, -0.2753, 0.5187, 0.4580,
      1.0000, 0.3694, 0.4056, 0.4700, 0.1531, 0.0067, 0.0377],
     [0.1264, -0.2999, -0.2285, 0.1499, 0.2161, -0.0514, -0.3335,
      -0.5043, -0.4064, -0.4343, -0.4589, -0.3914, 0.2251, 0.4671,
      0.3694, 1.0000, 0.4558, 0.3702, -0.0696, -0.0459, -0.0714],
     [0.1764, -0.2126, 0.0404, 0.1806, 0.2028, 0.1232, -0.2204,
      -0.2639, -0.3091, -0.1754, -0.2221, -0.2784, 0.3262, 0.4855,
      0.4056, 0.4558, 1.0000, 0.6706, 0.1294, -0.1370, 0.0518],
     [0.1810, -0.1075, 0.1922, 0.2068, 0.1481, 0.2563, -0.2238,
      -0.2222, -0.2343, -0.0763, -0.1780, -0.2883, 0.1950, 0.2041,
      0.4700, 0.3702, 0.6706, 1.0000, 0.1610, 0.0159, 0.0401],
     [0.3585, -0.1158, 0.4571, 0.3694, -0.0627, 0.5068, 0.1223,
      0.1836, -0.0003, 0.1807, 0.1993, 0.0099, 0.2035, 0.0488,
      0.1531, -0.0696, 0.1294, 0.1610, 1.0000, 0.1075, 0.7584],
     [-0.2776, 0.3768, -0.0356, -0.1675, 0.1207, -0.0576, 0.1121,
      0.0165, 0.0209, 0.0773, 0.0458, 0.1892, -0.1137, -0.2098,
      0.0067, -0.0459, -0.1370, 0.0159, 0.1075, 1.0000, -0.0486],
     [0.2169, -0.3763, 0.2327, 0.2355, -0.3495, 0.2777, 0.2360,
      0.2772, 0.0343, 0.1505, 0.2033, 0.0639, -0.0447, -0.0927,
      0.0377, -0.0714, 0.0518, 0.0401, 0.7584, -0.0486, 1.0000]]

# Constant from gen_rup
par_t_acc = 0.2

def get_nodes(points, dx, dy):
    """
    """
    nodes = np.array([0, 0])
    nodes[0] = math.ceil(points[0] / dx)
    nodes[1] = math.ceil(points[1] / dy)
    nodes[nodes == 0] = 1.0

    return nodes

def get_points(nodes, dx, dy):
    """
    Does the opposite of get_nodes
    """
    # Scale points
    points = np.dot(np.diag([dx, dy]), nodes)
    # Subtract dx/2
    points[0] = points[0] - (dx / 2.0)
    points[1] = points[1] - (dy / 2.0)

    return points

def sub2ind(i, j, m, n):
#    return int(np.ravel_multi_index([[i-1], [j-1]], [m, n], order='F'))
    return((i + (j - 1) * m) - 1)

def ind2sub(ind, m, n):
    return np.unravel_index(ind, (m, n), order='F')
#    ind = ind + 1
#    j = math.ceil(ind / (1.0 * m))
#    i = ind - (j - 1) * m
#    return i - 1, j - 1

def is_in_domain(i, j, m, n):
    """
    Checks if i,j inside m,n
    """
    return (i >= 0) and (i < m) and (j >= 0) and (j < n)

def is_in_domain2(i, j, m, n):
    """
    Checks if i,j inside m,n
    """
    return (i > 0) and (i < m) and (j > 0) and (j < n)

def reckon(lat0, lon0, rng, az):
    """
    Computes the latitude and longitude positions for selected ranges
    and azimuths from a starting point along a great circle path. The
    range is input as degrees of arc length on a sphere The input
    azimuth is measured clockwise from due north
    """
    # epsilon = 1.7453e-07
    lat = lat0 * math.pi / 180.0
    lon = lon0 * math.pi / 180.0
    az = az * math.pi / 180.0
    rng = rng * math.pi / 180.0

    # Compute latitude
    temp1 = math.sin(lat) * math.cos(rng)
    temp2 = math.cos(lat) * math.sin(rng) * math.cos(az)
    tmplat = math.asin(temp1 + temp2)

    # Compute longitude
    temp1 = math.sin(rng) * math.sin(az)
    temp2 = math.cos(lat) * math.cos(rng)
    temp3 = math.sin(lat) * math.sin(rng) * math.cos(az)
    tmplon = lon + math.atan2(temp1, temp2-temp3)

    # Convert back to degrees
    newlat = tmplat * 180.0 / math.pi
    newlon = (math.pi *
              ((abs(tmplon) / math.pi) -
               2 * math.ceil(((abs(tmplon) / math.pi) - 1) / 2.0)) *
              math.copysign(1, tmplon))
    newlon = newlon * 180.0 / math.pi

    return newlat, newlon

def calc_time(i, j, Fij, T, frozen, m, n, dy, dx, order):
    """
    Calculates the time-distance at (i,j) from the neighboring Frozen
    pixels. The formula is the standard quadratic equation along with
    the min-max switches.
    """
    patch = np.empty((5, 5))
    patch[:] = np.inf

    for pi in range(0, 5):
        for pj in range(0, 5):
            if ((is_in_domain(i-3+pi, j-3+pj, m, n)) and
                (frozen[i-3+pi][j-3+pj])):
                patch[pi][pj] = T[i-3+pi][j-3+pj]

    # If all the surrounding cross values are inf, set the time to
    # Inf, otherwise it would results in b^2 - 4 * a(=0) *
    # c(=1/Fij=Inf)= NaN
    if min(patch[1][2], patch[3][2], patch[2][1], patch[2][3]) == np.inf:
        time = np.inf
        e_flag = 0
        return time, e_flag

    # Get the minimum vertical neighbor or don't use vertical
    # component in gradient
    oy = 0
    ymin1 = 0
    ymin2 = 0
    if not (patch[1][2] == np.inf and patch[3][2] == np.inf):
        oy = 1
        if patch[1][2] <= patch[3][2]:
            ymin1 = patch[1][2]
            if order == 2 and (patch[0][2] <= patch[1][2]):
                ymin2 = patch[0][2]
                oy = 2
        else:
            ymin1 = patch[3][2]
            if order == 2 and (patch[4][2] <= patch[3][2]):
                ymin2 = patch[4][2]
                oy = 2

    # Get the minimum horizontal neighbor or don't use horizontal
    # component in gradient
    ox = 0
    xmin1 = 0
    xmin2 = 0
    if not (patch[2][1] == np.inf and patch[2][3] == np.inf):
        ox = 1
        if patch[2][1] <= patch[2][3]:
            xmin1 = patch[2][1]
            if order == 2 and (patch[2][0] <= patch[2][1]):
                xmin2 = patch[2][0]
                ox = 2
        else:
            xmin1 = patch[2][3]
            if order == 2 and (patch[2][4] <= patch[2][3]):
                xmin2 = patch[2][4]
                ox = 2

    # Calculate coefficients of 2nd degree polynomial corresponding to
    # the time-distance equation. Coeffs: aT^2 + bT + c = 0. The
    # (ox==1) and (ox==2) act as switches. Remove all the (ox==2) and
    # what is left are the first order scheme coefficients.
    a = ((oy == 1)+9.0/4*(oy == 2))/dy**2 + ((ox == 1)+9.0/4*(ox == 2))/dx**2
    b = (((-2.0*(oy==1)-6.0*(oy == 2))*ymin1 + 3.0/2.0*(oy == 2)*ymin2)/dy**2 +
         ((-2.0*(ox==1)-6.0*(ox==2))*xmin1 + 3.0/2.0*(ox == 2)*xmin2)/dx**2)
    c = (((oy == 1)*ymin1**2 +
          (oy == 2)*(4.0*ymin1**2 + 1.0/4.0*ymin2**2 - 2.0*ymin1*ymin2))/dy**2 +
         ((ox == 1)*xmin1**2 +
          (ox == 2)*(4.0*xmin1**2 + 1.0/4.0*xmin2**2 - 2.0*xmin1*xmin2))/dx**2 -
         (1.0/Fij)**2)
    d = b**2-4.0*a*c

    # One can analytically prove that even if causality is not
    # violated (i.e. all Frozen pixels should always have lower values
    # than all narrow band pixels), complex time-dist results may
    # still arise for 2nd order schemes. This is not the case for 1st
    # order schemes, which may only yield complex results if the
    # causality aspect of the algorithm has been violated. Higher
    # order schemes are less tolerant of extreme variations in the
    # speed-map values.

    # This implementation first attempts to use the 1st order scheme
    # if the second order scheme fails. If that also fails, meaning
    # the causality has been violated (which might happen if you're
    # very unlucky, because of previous reversions,or very extreme
    # speed-map differences around the source points) we simply add
    # 1/Fij to the smallest neighbour. In case of a reversion, an
    # persistent error flag 'eFlag' is set to 1. In case of a
    # violation, it's set to 2.

    # This is better error treatment than what is used in the 'msfm'
    # (multistencil fast-marching) code imho, where any complex
    # calculation always results in [add 1/Fij to smallest neigh.],
    # and which doesn't notify the user with any error flags.
    if d < 0 and order == 2:
        # Revert to order 1
        order = 1
        (time, tmp_flag) = calc_time(i, j, Fij, T, frozen, m, n, dy, dx, order)
        e_flag = max(1, tmp_flag)
    elif d < 0 and order == 1:
        # Add 1/Fij to smallest neighbour
        if oy == 0:
            ymin1 = np.inf
        if ox == 0:
            xmin1 = np.inf
        time = min(xmin1, ymin1) + 1.0/Fij
        e_flag = 2
    else:
        # All Good!
        # Solve quadratic equation, only use the maximum root
        time = (-b+math.sqrt(d))/(2.0*a)
        e_flag = 0

    return time, e_flag

class SongRMGMS(object):
    """
    This class implements the Song RMG Multi Segment rupture generator
    """
    # Number of stats points
    n_stats = 1000

    @staticmethod
    def fm(F, source_points, dxyz):
        """
        Implements the Fast-marching method to solve the Eikonal eqn
        in 2 dimensions
        """

        # Initialize parameters
        # ndims = 2
        m, n = F.shape
        dx = dxyz[0]
        dy = dxyz[1]

        # Check error conditions
        if not np.all(F >= 0):
            raise ValueError("F must be >= 0")

        if ((n < 1) or (m < 1)):
            raise ValueError("m and n must by > 0")

        if (dx <= 0) or (dy <= 0):
            raise ValueError("dxyz must be > 0")

        # Flip source points
        source_points = source_points[::-1]

        # Time distance matrix, initialize all points to -1
        T = np.zeros((m, n)) - 1

        # Matrix that keeps tabs on Frozen pixels
        frozen = np.zeros((m, n))

        # Construct dictionary that we use as a heap
        narrow_band = {}

        # Relative indices of the 4 neighbors to a given point
        i_offsets = [-1, 1, 0, 0]
        j_offsets = [0, 0, -1, 1]

        # Relative indices of the 8 neighbors to a given point
        i_offsets_full = [-1, -1, -1, 0, 0, 1, 1, 1]
        j_offsets_full = [-1, 0, 1, -1, 1, -1, 0, 1]

        # Initialize error flag
        e_flag = 0

        # First we calculate the time-distances to all 8 neighboring
        # pixels of the source point. This calculation is done simply
        # by taking the distance and dividing by the speed. We also
        # freeze this

        # Get node and its indice
        cp = source_points
        cp_node = get_nodes(cp, dx, dy)
        cp_ind = sub2ind(cp_node[0], cp_node[1], m, n)

        # Calculate travel time and freeze
        cp_points = get_points(cp_node, dx, dy)
        FT = T.flatten(order='F')
        FT[cp_ind] = (np.sqrt(np.dot(cp - cp_points, cp - cp_points)) /
                      F.flatten(order='F')[cp_ind])
        T = FT.reshape((m, n), order='F')
        f_frozen = frozen.flatten(order='F')
        f_frozen[cp_ind] = 1
        frozen = f_frozen.reshape((m, n), order='F')

        # For all 8 neighbors of the source point
        for neigh in range(0, 8):
            # Get index of neighbor, store as i and j
            ni = cp_node[0] + i_offsets_full[neigh]
            nj = cp_node[1] + j_offsets_full[neigh]
            n_ind = sub2ind(ni, nj, m, n)

            # Only check if is_in_domain
            if is_in_domain2(ni, nj, m, n):
                val = cp - get_points([ni, nj], dx, dy)
                time = (np.sqrt(np.dot(val, val)) /
                        F.flatten(order='F')[n_ind])
                FT = T.flatten(order='F')
                if FT[n_ind] >= 0:
                    FT[n_ind] = min(time, FT[n_ind])
                else:
                    FT[n_ind] = time
                T = FT.reshape((m, n), order='F')

                f_frozen = frozen.flatten(order='F')
                f_frozen[n_ind] = 1
                frozen = f_frozen.reshape((m, n), order='F')

        # Calculate the initial narrow band as all neighboring pixels
        # to the ones that have been frozen. Note that this time,
        # unlike the source-point loop, the henceforth in the
        # algorithm, the neighbors of a pixel only include its 4
        # non-diagonal neighbors.
        for ind in range(0, m*n):
            if frozen.flatten(order='F')[ind]:
                i, j = ind2sub(ind, m, n)
                # For all 4 neighbors of the frozen points
                for neigh in range(0, 4):
                    # Get index of neighbor, store as i and j
                    ni = i + i_offsets[neigh]
                    nj = j + j_offsets[neigh]
                    n_ind = sub2ind(ni+1, nj+1, m, n)

                    # If (i,j) valid for consideration
                    if ((is_in_domain(ni, nj, m, n)) and
                        (not frozen.flatten(order='F')[n_ind])):
                        if not n_ind in narrow_band:
                            (time,
                             flag) = calc_time(ni+1, nj+1,
                                               F.flatten(order='F')[n_ind],
                                               T, frozen, m, n, dy, dx, 2)
                            narrow_band[n_ind] = time
                            e_flag = max(flag, e_flag)

        # Main Loop
        # Now start the main loop.
        # Loop until there are no more narrow band
        # neighbours, meaning the algorithm has finished
        l_count = 0
        while len(narrow_band) > 0:
            l_count = l_count + 1

            # Get min heap element.
            # This will be the new "center pixel" (CP)
            cp = min(narrow_band, key=narrow_band.get)
            time = narrow_band[cp]
            narrow_band.pop(cp)

            i, j = ind2sub(cp, m, n)
            # Freeze and set time
            f_frozen = frozen.flatten(order='F')
            f_frozen[cp] = 1
            frozen = f_frozen.reshape((m, n), order='F')
            FT = T.flatten(order='F')
            FT[cp] = time
            T = FT.reshape((m, n), order='F')

            # For all neighbours of CP
            for neigh in range(0, 4):
                # Get index of neighbor, store as i and j
                ni = i + i_offsets[neigh]
                nj = j + j_offsets[neigh]
                if not is_in_domain(ni, nj, m, n):
                    continue
                n_ind = sub2ind(ni+1, nj+1, m, n)

                # If (i,j) valid for consideration
                if ((is_in_domain(ni, nj, m, n)) and
                    (not frozen.flatten(order='F')[n_ind])):
                    (time,
                     flag) = calc_time(ni+1, nj+1,
                                       F.flatten(order='F')[n_ind],
                                       T, frozen, m, n, dy, dx, 2)
                    narrow_band[n_ind] = time

                    e_flag = max(flag, e_flag)

        return T, e_flag

    @staticmethod
    def f_moment_n(slip_values, width, length):
        """
        Calculates the moment of an event for a given slip
        distribution (in cm)
        """
        mu = 3.3 * 1e10 # shear modulus = rigidity [N/m^2]
        slip = slip_values.mean() / 100 # Average slip over fault surface [m]
        area = width * length * 1e6 # Fault surface [m]
        moment = mu * slip * area # Moment [Nm]
        magnitude = (2.0 / 3.0) * (math.log10(moment) - 9.05)

        return moment, magnitude

    @staticmethod
    def b_taper(nz, nx, dtop):
        """
        Tapering slip values near the boundary
        """
        nz_bt = int(np.round(nz / 5))
        nx_bt = int(np.round(nx / 5))

        sfac = np.ones((nz, nx))
        x = np.linspace(0, 1, nx_bt)
        z = np.linspace(0, 1, nz_bt)

        sfac_x = np.sqrt(1 - x**2)
        sfac_z = np.sqrt(1 - z**2)

        if dtop < 0.1:
            sfac_x = np.hstack([np.tile(sfac_x[::-1], (nz, 1)),
                                np.ones((nz, nx-2*nx_bt)),
                                np.tile(sfac_x, (nz, 1))])
            sfac_z = np.vstack([np.ones((nz-nz_bt, nx)),
                                np.tile(sfac_z, (nx, 1)).T])
        else:
            sfac_x = np.hstack([np.tile(sfac_x[::-1], (nz, 1)),
                                np.ones((nz, nx-2*nx_bt)),
                                np.tile(sfac_x, (nz, 1))])
            sfac_z = np.vstack([np.tile(sfac_z[::-1], (nx, 1)).T,
                                np.ones((nz-2*nz_bt, nx)),
                                np.tile(sfac_z, (nx, 1)).T])

        sfac = sfac_x * sfac_z

        return sfac

    def __init__(self, r_velmodel, r_srcfile,
                 o_r_srffile, i_vmodel_name,
                 sim_id=0, **kwargs):
        """
        Initialize class variables
        """
        self.sim_id = sim_id
        self.r_srffile = o_r_srffile
        self.rup_stats = None
        self.config = None
        self.r_srcfile = r_srcfile
        self.rup = {}
        self.use_interpolation = True
        self.r_srcfiles = []

        # Get all src files that were passed to us
        if kwargs is not None:
            for idx in range(len(kwargs)):
                self.r_srcfiles.append(kwargs['src%d' % (idx)])
        else:
            # Not a multisegment run, just use the single src file
            self.r_srcfiles.append(r_srcfile)

    def gen_stats_inp(self):
        """
        This function generates the src_stats_inp structures
        """
        self.rup_stats = [None for _ in range(self.config.num_srcfiles)]

        for segment in range(self.config.num_srcfiles):
            # Calculate some values
            fault_area = (self.config.CFGDICT[segment]["fault_length"] *
                          self.config.CFGDICT[segment]["fault_width"] * 1e6)
            target_mo = 10**(3.0 * self.config.CFGDICT[segment]['magnitude'] /
                             2.0 + 9.05)
            seg_mo = target_mo * self.config.CFGDICT[segment]["moment_fraction"]
            mu_slip = 1.0 * seg_mo / mu / fault_area * 100.0

            L = np.linalg.cholesky(C)
            N = len(C)
            np.random.seed(int(self.config.CFGDICT[segment]['seed']))

            # Initialize arrays
            s = np.zeros((N, self.n_stats))
            sigs = np.append(sig_mu, sig_sig)
            sigs = np.append(sigs, sig_ax)
            sigs = np.append(sigs, sig_az)
            sigs = np.append(sigs, sig_cc)
            mus = np.append(mu_mu, mu_sig)
            mus = np.append(mus, mu_ax)
            mus = np.append(mus, mu_az)
            mus = np.append(mus, mu_cc)

            # Loop
            for i in range(0, self.n_stats):
                s[:, i] = np.random.randn(N)
                s[0, i] = (mu_slip - mu_mu[0]) / sig_mu[0] / L[0][0]
                s[:, i] = np.dot(L, s[:, i])
                s[:, i] = s[:, i] * sigs + mus
                s[6:18, i] = np.exp(s[6:18, i])

            # Now we transpose it
            s = s.T

            # Save copy
            self.rup_stats[segment] = s

    def gen_src(self):
        """
        This function sets up the basic input parameters to describe a
        finite source model
        """
        config = self.config
        rup = self.rup
        rup["nx"] = []
        rup["nz"] = []
        rup["lx"] = []
        rup["lz"] = []
        #rup["dis"] = []

        for segment in range(config.num_srcfiles):
            cfg_dict = config.CFGDICT[segment]
            rup["nx"].append(int(math.ceil(cfg_dict["fault_length"] /
                                           cfg_dict["dlen"])))
            rup["nz"].append(int(math.ceil(cfg_dict["fault_width"] /
                                           cfg_dict["dwid"])))

            rup["lx"].append(np.linspace(cfg_dict["dlen"] / 2,
                                         ((cfg_dict["dlen"] *
                                           rup["nx"][segment]) -
                                          cfg_dict["dlen"] / 2),
                                         num=int(math.ceil((cfg_dict["dlen"] *
                                                            rup["nx"][segment] -
                                                            cfg_dict["dlen"] / 2) /
                                                           cfg_dict["dlen"]))))
            #rup["lx"][segment] = rup["lx"][segment] - (rup["nx"][segment] *
            #                                           cfg_dict["dlen"] /
            #                                           2 + cfg_dict["hypo_along_stk"])

            rup["lz"].append(np.linspace(cfg_dict["dwid"] / 2,
                                         ((cfg_dict["dwid"] *
                                           rup["nz"][segment]) -
                                          cfg_dict["dwid"] / 2),
                                         num=int(math.ceil((cfg_dict["dwid"] *
                                                            rup["nz"][segment] -
                                                            cfg_dict["dwid"] / 2) /
                                                           cfg_dict["dwid"]))))
            #rup["lz"][segment] = rup["lz"][segment] - cfg_dict["hypo_down_dip"]

            #XX, ZZ = np.meshgrid(rup["lx"][segment], rup["lz"][segment])
            #rup["dis"].append(np.sqrt(XX*XX + ZZ*ZZ))

        # If we are not using interpolation, skip the setup below
        if not self.use_interpolation:
            return

        rup["dx1"] = []
        rup["dz1"] = []
        rup["nx1"] = []
        rup["nz1"] = []
        rup["lx1"] = []
        rup["lz1"] = []

        # Parameters for the interpolation
        for segment in range(config.num_srcfiles):
            cfg_dict = config.CFGDICT[segment]

            rup["dx1"].append(cfg_dict["dlen"] * 4.0)
            rup["dz1"].append(cfg_dict["dwid"] * 4.0)

            rup["nx1"].append(int(math.ceil(cfg_dict["fault_length"] /
                                            rup["dx1"][segment])))
            rup["nz1"].append(int(math.ceil(cfg_dict["fault_width"] /
                                            rup["dz1"][segment])))

            rup["lx1"].append(np.linspace(rup["dx1"][segment] / 2,
                                          ((rup["dx1"][segment] *
                                            rup["nx1"][segment]) -
                                           rup["dx1"][segment] / 2),
                                          num=int(math.ceil((rup["dx1"][segment] *
                                                             rup["nx1"][segment] -
                                                             rup["dx1"][segment] / 2) /
                                                            rup["dx1"][segment]))))
            #rup["lx1"][segment] = rup["lx1"][segment] - (rup["nx1"][segment] *
            #                                             rup["dx1"][segment] /
            #                                             2 + cfg_dict["hypo_along_stk"])

            rup["lz1"].append(np.linspace(rup["dz1"][segment] / 2,
                                          ((rup["dz1"][segment] *
                                            rup["nz1"][segment]) -
                                           rup["dz1"][segment] / 2),
                                          num=int(math.ceil((rup["dz1"][segment] *
                                                             rup["nz1"][segment] -
                                                             rup["dz1"][segment] / 2) /
                                                            rup["dz1"][segment]))))
            #rup["lz1"][segment] = rup["lz1"][segment] - cfg_dict["hypo_down_dip"]

            rup["lx1"][segment] = np.insert(rup["lx1"][segment], 0,
                                            rup["lx1"][segment][0] -
                                            rup["dx1"][segment])
            rup["lx1"][segment] = np.append(rup["lx1"][segment],
                                            rup["lx1"][segment][-1] +
                                            rup["dx1"][segment])

            rup["lz1"][segment] = np.insert(rup["lz1"][segment], 0,
                                            rup["lz1"][segment][0] -
                                            rup["dz1"][segment])
            rup["lz1"][segment] = np.append(rup["lz1"][segment],
                                            rup["lz1"][segment][-1] +
                                            rup["dz1"][segment])

            rup["nx1"][segment] = rup["nx1"][segment] + 2
            rup["nz1"][segment] = rup["nz1"][segment] + 2

    def gen_stats(self):
        """
        This function set up target 1-point and 2-point statistics
        """
        config = self.config
        rup = self.rup
        # Pick random number
        random.seed(config.CFGDICT[0]['seed'])
        stats_id = int(random.random() * 1000)
        seg_ihypo = config.seg_hypocenter

        # Create data structures
        rup["p1_mu"] = [None for _ in range(config.num_srcfiles)]
        rup["p1_sig"] = [None for _ in range(config.num_srcfiles)]
        rup["p1_mu1"] = [None for _ in range(config.num_srcfiles)]
        rup["p1_sig1"] = [None for _ in range(config.num_srcfiles)]
        rup["p2_tag"] = [None for _ in range(config.num_srcfiles)]
        rup["p2_ax"] = [None for _ in range(config.num_srcfiles)]
        rup["p2_az"] = [None for _ in range(config.num_srcfiles)]
        rup["p2_cc"] = [None for _ in range(config.num_srcfiles)]
        rup["p2_RDx"] = [None for _ in range(config.num_srcfiles)]
        rup["p2_RDz"] = [None for _ in range(config.num_srcfiles)]

        rup["p1_min"] = [0.0, 1.0, 1.0, 0.2]
        rup["p1_max"] = [500.0, 6.0, 500.0, 10.0]
        rup["lambda"] = 1.0
        rup["fN"] = 6

        for segment in range(config.num_srcfiles):
            # 1-point statistics
            # [slip, Vr, Vmax, risT]

            rup["p1_mu"][segment] = [self.rup_stats[seg_ihypo][stats_id, 0],
                                     self.rup_stats[seg_ihypo][stats_id, 1]*1.4,
                                     self.rup_stats[seg_ihypo][stats_id, 2]]
            rup["p1_sig"][segment] = [self.rup_stats[seg_ihypo][stats_id, 3],
                                      self.rup_stats[seg_ihypo][stats_id, 4]*0.7,
                                      self.rup_stats[seg_ihypo][stats_id, 5]]
            # To be completed later...
            rup["p1_mu1"][segment] = [None, None, None, None]
            rup["p1_sig1"][segment] = [None, None, None, None]

            # 2-point statistics
            rup["p2_tag"][segment] = "positive_eigen"

            rup["p2_ax"][segment] = [[self.rup_stats[seg_ihypo][stats_id, 6],
                                      self.rup_stats[seg_ihypo][stats_id, 7],
                                      self.rup_stats[seg_ihypo][stats_id, 8]],
                                     [np.nan,
                                      self.rup_stats[seg_ihypo][stats_id, 9],
                                      self.rup_stats[seg_ihypo][stats_id, 10]],
                                     [np.nan,
                                      np.nan,
                                      self.rup_stats[seg_ihypo][stats_id, 11]]]
            rup["p2_az"][segment] = [[self.rup_stats[seg_ihypo][stats_id, 12],
                                      self.rup_stats[seg_ihypo][stats_id, 13],
                                      self.rup_stats[seg_ihypo][stats_id, 14]],
                                     [np.nan,
                                      self.rup_stats[seg_ihypo][stats_id, 15],
                                      self.rup_stats[seg_ihypo][stats_id, 16]],
                                     [np.nan,
                                      np.nan,
                                      self.rup_stats[seg_ihypo][stats_id, 17]]]
            rup["p2_cc"][segment] = [[1.0,
                                      self.rup_stats[seg_ihypo][stats_id, 18],
                                      self.rup_stats[seg_ihypo][stats_id, 19]],
                                     [np.nan,
                                      1.0,
                                      self.rup_stats[seg_ihypo][stats_id, 20]],
                                     [np.nan,
                                      np.nan,
                                      1.0]]
            rup["p2_RDx"][segment] = [[0.0, 0.0, 0.0],
                                      [np.nan, 0.0, 0.0],
                                      [np.nan, np.nan, 0.0]]
            rup["p2_RDz"][segment] = [[0.0, 0.0, 0.0],
                                      [np.nan, 0.0, 0.0],
                                      [np.nan, np.nan, 0.0]]

    def gen_rho(self, r_x, r_z, seg):
        """
        Generate spatial coherence structure, currently exponential
        function only
        """
        rup = self.rup

        N = len(rup["p1_mu"][seg])

        rho = np.empty((N, N), dtype=object)
        for k in range(0, N):
            for i in range(k, N):
                rho[k][i] = np.exp(-np.sqrt(((r_x - rup["p2_RDx"][seg][k][i]) /
                                             rup["p2_ax"][seg][k][i]) ** 2 +
                                            ((r_z - rup["p2_RDz"][seg][k][i]) /
                                             rup["p2_az"][seg][k][i]) ** 2))
                rho[k][i] = np.dot(rup["p2_cc"][seg][k][i], rho[k][i])

        return rho

    def gen_dist(self):
        """
        Generate 2D distributions of source parameters using the
        Cholesky factorization based on 1-point and 2-point statistics
        """
        config = self.config
        rup = self.rup

        # Initialize data structures
        rup["slip_dist"] = []
        rup["vr_dist"] = []
        rup["psv_dist"] = []
        rup["slip1_dist"] = []
        rup["vr1_dist"] = []
        rup["psv1_dist"] = []

        for segment in range(config.num_srcfiles):

            print()
            print("** Segment: %d" % (segment + 1))
            print()

            t1 = time.time()

            N = int(rup["nx"][segment] * rup["nz"][segment])
            print()
            print("=> Number of subfaults (fine grid): %d" % (N))
            if self.use_interpolation:
                # Also print information about the coarse grid
                N1 = int(rup["nx1"][segment] * rup["nz1"][segment])
                print("=> Number of subfaults (coarse grid): %d" % (N1))
            print()

            if not self.use_interpolation:
                # Does not use interpolation
                XX, ZZ = np.meshgrid(rup["lx"][segment], rup["lz"][segment])
                X = XX.flatten(order="F").astype('float32')
                Z = ZZ.flatten(order="F").astype('float32')

                Z1 = np.tile(Z, (N, 1)).T
                del Z
                r_z = (Z1.T - Z1)
                del Z1

                X1 = np.tile(X, (N, 1)).T
                del X
                r_x = (X1.T - X1)
                del X1
            else:
                # Uses the interpolation code
                XX, ZZ = np.meshgrid(rup["lx1"][segment], rup["lz1"][segment])
                X = XX.flatten(order="F").astype('float32')
                Z = ZZ.flatten(order="F").astype('float32')

                Z1 = np.tile(Z, (N1, 1)).T
                del Z
                r_z = (Z1.T - Z1)
                del Z1

                X1 = np.tile(X, (N1, 1)).T
                del X
                r_x = (X1.T - X1)
                del X1

            rho = self.gen_rho(r_x, r_z, segment)

            cm = np.vstack([np.hstack([rho[0][0], rho[0][1], rho[0][2]]),
                            np.hstack([rho[0][1].T, rho[1][1], rho[1][2]]),
                            np.hstack([rho[0][2].T, rho[1][2].T, rho[2][2]])])

            # Step 1 done, look how long it took
            t2 = time.time()
            print("=> Elapsed time for constructing Cm: %10.2f" % (t2 - t1))

            t1 = time.time()

            eig_d, eig_v = np.linalg.eig(cm)
            eig_d[eig_d < 0.01] = 0.01
            cm = np.dot(np.dot(eig_v, np.diag(eig_d)), (eig_v.T))

            # Step 2 done, look how long it took
            t2 = time.time()
            print("=> Elapsed time for Eigen decomposition: %10.2f" % (t2 - t1))

            t1 = time.time()

            L = np.linalg.cholesky(cm)

            # Step 3 done, look how long it took
            t2 = time.time()
            print("=> Elapsed time for Cholesky Factorization: %10.2f" % (t2 - t1))

            t1 = time.time()

            if not self.use_interpolation:
                s0 = np.random.randn(3 * N)
                s = np.dot(L, s0)
                slip = s[:N]
                vr = s[N:2*N]
                psv = s[2*N:]

                rup["slip_dist"].append(np.reshape(slip, (rup["nz"][segment],
                                                          rup["nx"][segment]),
                                                   order='F'))
                rup["vr_dist"].append(np.reshape(vr, (rup["nz"][segment],
                                                      rup["nx"][segment]),
                                                 order='F'))
                rup["psv_dist"].append(np.reshape(psv, (rup["nz"][segment],
                                                        rup["nx"][segment]),
                                                  order='F'))
            else:
                s0 = np.random.randn(3 * N1)
                s = np.dot(L, s0)
                slip1 = s[:N1]
                vr1 = s[N1:2*N1]
                psv1 = s[2*N1:]

                rup["slip1_dist"].append(np.reshape(slip1,
                                                    (rup["nz1"][segment],
                                                     rup["nx1"][segment]),
                                                    order='F'))
                rup["vr1_dist"].append(np.reshape(vr1,
                                                  (rup["nz1"][segment],
                                                   rup["nx1"][segment]),
                                                  order='F'))
                rup["psv1_dist"].append(np.reshape(psv1,
                                                   (rup["nz1"][segment],
                                                    rup["nx1"][segment]),
                                                   order='F'))

            # Step 4 done, look how long it took
            t2 = time.time()
            print("=> Elapsed time for random sampling: %10.2f" % (t2 - t1))

    def svf_etinti(self, tau_s, tau_r, dt, nt):
        """
        Generating the shape of slip velocity function (SVF)
        based on Tinti et al. (BSSA, 2005)
        """
        t = np.linspace(0, (nt-1)*dt, int(nt))

        if tau_r < tau_s:
            raise ValueError("Tau_r should be larger than Tau_s")

        # U_tot (final slip) is one here
        K = 2.0 / (math.pi * tau_r * tau_s**2)

        # Ignore the divide by zero error
        with np.errstate(divide='ignore'):
            C1 = ((t / 2.0 + 1.0 / 4 * tau_r) *
                  np.sqrt(np.array((t * (tau_r - t)), dtype=complex)) +
                  (t * tau_r - tau_r ** 2) *
                  np.arcsin(np.sqrt(np.array((t / tau_r)), dtype=complex)) -
                  [3.0 / 4 * tau_r ** 2 * y for y in
                   [cmath.atan(x) for x in np.sqrt(np.array(((tau_r - t) / t),
                                                            dtype=complex))]])

        C2 = 3.0 / 8 * math.pi * tau_r ** 2

        C3 = ((tau_s - t - 1.0 / 2 * tau_r) *
              np.sqrt(np.array(((t - tau_s) * (tau_r - t + tau_s)),
                               dtype=complex)) +
              tau_r * (2.0 * tau_r - 2 * t + 2 * tau_s) *
              np.arcsin(np.sqrt(np.array(((t - tau_s) / tau_r),
                                         dtype=complex))) +
              [3.0 / 2 * tau_r ** 2 * y for y in
               [cmath.atan(x) for x in np.sqrt(np.array((tau_r - t + tau_s) /
                                                        (t - tau_s),
                                                        dtype=complex))]])

        C4 = ((-tau_s + 1.0 / 2.0 * t + 1.0 / 4.0 * tau_r) *
              np.sqrt(np.array(((t - 2 * tau_s) * (tau_r - t + 2 * tau_s)),
                               dtype=complex)) +
              tau_r * (-tau_r + t - 2 * tau_s) *
              np.arcsin(np.sqrt(np.array(((t - 2 * tau_s) / tau_r),
                                         dtype=complex))) -
              [3.0 / 4.0 * tau_r ** 2 * y for y in
               [cmath.atan(x) for x in np.sqrt(np.array((tau_r - t + 2*tau_s) /
                                                        (t - 2 * tau_s),
                                                        dtype=complex))]])

        C5 = math.pi / 2.0 * tau_r * (t - tau_r)

        C6 = math.pi / 2.0 * tau_r * (2.0 * tau_s - t + tau_r)

        # Filters ComplexWarning when typecasting array to float
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            if tau_r > 2 * tau_s:
                svf = (np.array(((C1 + C2) *
                                 [x1 and x2 for x1, x2 in zip(t >= 0,
                                                              t < tau_s)]),
                                dtype=float) +
                       (np.array(((C1 - C2 + C3) *
                                  [x1 and x2 for x1, x2 in zip(t >= tau_s,
                                                               t < 2*tau_s)]),
                                 dtype=float)) +
                       (np.array(((C1 + C3 + C4) *
                                  [x1 and x2 for x1, x2 in zip(t >= 2*tau_s,
                                                               t < tau_r)]),
                                 dtype=float)) +
                       (np.array(((C5 + C3 + C4) *
                                  [x1 and x2 for x1, x2 in zip(t >= tau_r,
                                                               t < tau_r +
                                                               tau_s)]),
                                 dtype=float)) +
                       (np.array(((C4 + C6) *
                                  [x1 and x2 for x1, x2 in zip(t >= (tau_r +
                                                                     tau_s),
                                                               t < (tau_r+2 *
                                                                    tau_s))]),
                                 dtype=float)))# +
#                   (np.array(((0) *
#                              [x1 for x1 in t>=(tau_r+2*tau_s)]),
#                             dtype=float)))
            else:
                svf = (np.array(((C1 + C2) *
                                 [x1 and x2 for x1, x2 in zip(t >= 0,
                                                              t < tau_s)]),
                                dtype=float) +
                       (np.array(((C1 - C2 + C3) *
                                  [x1 and x2 for x1, x2 in zip(t >= tau_s,
                                                               t < tau_r)]),
                                 dtype=float)) +
                       (np.array(((C5 + C3 - C2) *
                                  [x1 and x2 for x1, x2 in zip(t >= tau_r,
                                                               t < 2*tau_s)]),
                                 dtype=float)) +
                       (np.array(((C5 + C3 + C4) *
                                  [x1 and x2 for x1, x2 in zip(t >= 2*tau_s,
                                                               t < tau_r +
                                                               tau_s)]),
                                 dtype=float)) +
                       (np.array(((C4 + C6) *
                                  [x1 and x2 for x1, x2 in zip(t >= (tau_r +
                                                                     tau_s),
                                                               t < (tau_r+2 *
                                                                    tau_s))]),
                                 dtype=float)))# +
#                   (np.array(((0) *
#                              [x1 for x1 in t>=(tau_r+2*tau_s)]),
#                             dtype=float)))

        svf = np.dot(K, svf)

        return svf, t

    def gen_svf(self, tw, dt, nt, svf_type):
        """
        Generating triangular or rectangular shape Slip Velocity
        Function
        """
        # Implements the 'etinti' svf_type only

        # 10% of total risetime
        #t_acc = tw * 0.1
        # Fixed T_acc
        t_acc = par_t_acc

        tau_s = t_acc / 1.3
        tau_r = tw - 2.0 * tau_s
        if tau_r <= tau_s:
            tau_r = tau_s * 1.01
        svf, time_vec = self.svf_etinti(tau_s, tau_r, dt, nt)

        # All done
        return svf, time_vec

    # def calculate_interp(self, X, Z, func):
    #    """
    #    Calls function 'func' as many times as needed to calculate
    #    the interpolated values for the points in X, Z
    #    """
    #    # Create output array
    #    out_data = np.zeros(X.shape)
    #    rows = X.shape[0]
    #    columns = X.shape[1]
    #
    #    # # Loop through arrays
    #    for row in range(rows):
    #        for column in range(columns):
    #            result = func(X[row][column], Z[row][column])
    #            out_data[row][column] = result
    #
    #    # Return results
    #    return out_data

    def gen_rup(self):
        """
        This function implements the rupture model generator based on
        1-point and 2-point statistics
        """
        rup = self.rup
        config = self.config
        num_segments = config.num_srcfiles
        seg_ihypo = config.seg_hypocenter

        print()
        print("### Step I: 2-point statistics")
        print("### generating 2D distributions of source parameters")
        print("### based on covariance matrix constructed from")
        print("### auto- and cross-correlation")
        print()

        self.gen_dist()

        if self.use_interpolation:
            # Added code for the interpolation
            rup["slip_dist"] = []
            rup["vr_dist"] = []
            rup["psv_dist"] = []

            for seg in range(num_segments):
                X, Z = np.meshgrid(rup["lx"][seg],
                                   rup["lz"][seg])
                X1, Z1 = np.meshgrid(rup["lx1"][seg],
                                     rup["lz1"][seg])

                rup["slip_dist"].append(griddata((X1.flatten(),
                                                  Z1.flatten()),
                                                 rup["slip1_dist"][seg].flatten(),
                                                 (X[0][None, :],
                                                  Z.T[0][:, None])))
                rup["vr_dist"].append(griddata((X1.flatten(), Z1.flatten()),
                                                rup["vr1_dist"][seg].flatten(),
                                                (X[0][None, :],
                                                 Z.T[0][:, None])))
                rup["psv_dist"].append(griddata((X1.flatten(), Z1.flatten()),
                                                rup["psv1_dist"][seg].flatten(),
                                                (X[0][None, :],
                                                 Z.T[0][:, None])))

                #f_slip_dist = interpolate.interp2d(X1, Z1, rup["slip1_dist"])
                #rup["slip_dist"] = self.calculate_interp(X, Z, f_slip_dist)
                ## rup["slip_dist"] = f_slip_dist(X, Z)
                #f_vr_dist = interpolate.interp2d(X1, Z1, rup["vr1_dist"])
                ## rup["vr_dist"] = f_vr_dist(X, Z)
                #rup["vr_dist"] = self.calculate_interp(X, Z, f_vr_dist)
                #f_psv_dist = interpolate.interp2d(X1, Z1, rup["psv1_dist"])
                ## rup["psv_dist"] = f_psv_dist(X, Z)
                #rup["psv_dist"] = self.calculate_interp(X, Z, f_psv_dist)

        print()
        print("### Step II: 1-point statistics")
        print("### Adjust 1-point statistics, e.g., mean and sigma")
        print("### Currently assume Gaussian distribution")
        print("### But it can be transformed into non Gaussian distribution")
        print("### if necessary")
        print()

        t1 = time.time()

        rup["rist_dist"] = []
        rup["seg_mo"] = [None for _ in range(num_segments)]
        rup["seg_mw"] = [None for _ in range(num_segments)]
        rup["rupT_dist"] = [None for _ in range(num_segments)]

        for seg in range(num_segments):
            # Slip, fixing mean slip and Mo
            rup["slip_dist"][seg] = ((rup["slip_dist"][seg] -
                                      rup["slip_dist"][seg].mean()) /
                                     rup["slip_dist"][seg].std())
            rup["slip_dist"][seg] = (rup["p1_sig"][seg][0] *
                                     rup["slip_dist"][seg] +
                                     rup["p1_mu"][seg][0])
            # Boundary tapering for slip
            sfac = self.b_taper(rup["nz"][seg],
                                rup["nx"][seg],
                                config.CFGDICT[seg]["depth_to_top"])
            rup["slip_dist"][seg] = rup["slip_dist"][seg] * sfac

            rup["slip_dist"][seg] = ((rup["slip_dist"][seg] -
                                      rup["slip_dist"][seg].mean()) /
                                     rup["slip_dist"][seg].std())
            rup["slip_dist"][seg] = (rup["p1_sig"][seg][0] *
                                     rup["slip_dist"][seg] +
                                     rup["p1_mu"][seg][0])
            rup["slip_dist"][seg][rup["slip_dist"][seg] < rup["p1_min"][0]] = rup["p1_min"][0]
            rup["slip_dist"][seg][rup["slip_dist"][seg] > rup["p1_max"][0]] = rup["p1_max"][0]
            rup["p1_mu1"][seg][0] = rup["slip_dist"][seg].mean()
            rup["p1_sig1"][seg][0] = rup["slip_dist"][seg].std()

            # Rupture Velocity
            rup["vr_dist"][seg] = (rup["p1_sig"][seg][1] *
                                   rup["vr_dist"][seg] +
                                   rup["p1_mu"][seg][1])
            rup["vr_dist"][seg][rup["vr_dist"][seg] < rup["p1_min"][1]] = rup["p1_min"][1]
            rup["vr_dist"][seg][rup["vr_dist"][seg] > rup["p1_max"][1]] = rup["p1_max"][1]
            rup["p1_mu1"][seg][1] = rup["vr_dist"][seg].mean()
            rup["p1_sig1"][seg][1] = rup["vr_dist"][seg].std()

            # Peak Slip Velocity
            rup["psv_dist"][seg] = (rup["p1_sig"][seg][2] *
                                    rup["psv_dist"][seg] +
                                    rup["p1_mu"][seg][2])
            rup["psv_dist"][seg][rup["psv_dist"][seg] < rup["p1_min"][2]] = rup["p1_min"][2]
            rup["psv_dist"][seg][rup["psv_dist"][seg] > rup["p1_max"][2]] = rup["p1_max"][2]
            rup["p1_mu1"][seg][2] = rup["psv_dist"][seg].mean()
            rup["p1_sig1"][seg][2] = rup["psv_dist"][seg].std()

            # Implement etinti method
            rup["rist_dist"].append(((1.04 * rup["slip_dist"][seg] /
                                      (par_t_acc**0.54 *
                                       rup["psv_dist"][seg]))**(1.0 / 0.47)))
            rup["rist_dist"][seg][rup["rist_dist"][seg] < rup["p1_min"][3]] = rup["p1_min"][3]
            rup["rist_dist"][seg][rup["rist_dist"][seg] > rup["p1_max"][3]] = rup["p1_max"][3]
            rup["p1_mu1"][seg][3] = rup["rist_dist"][seg].mean()
            rup["p1_sig1"][seg][3] = rup["rist_dist"][seg].std()
            rup["svf_t_acc"] = par_t_acc

            (rup["seg_mo"][seg],
             rup["seg_mw"][seg]) = self.f_moment_n(rup["slip_dist"][seg],
                                                    (rup["nz"][seg] *
                                                     config.CFGDICT[seg]["dwid"]),
                                                    (rup["nx"][seg] *
                                                     config.CFGDICT[seg]["dlen"]))

        rup["mo"] = sum(rup["seg_mo"])
        rup["mw"] = (2.0 / 3.0) * (math.log10(rup["mo"]) - 9.05)

        # Done, check how long it took
        t2 = time.time()

        print("=> Elapsed time for 1-point adjustment: %10.2f" % (t2 - t1))

        t1 = time.time()

        # Step III: post-processing of generated source models
        # a)rupture time calculation from simulated local rupture velocity
        print()
        print("### Step III: Post-processing")
        print("###")

        shypo = (config.CFGDICT[seg_ihypo]["fault_length"] /
                 2 + config.CFGDICT[seg_ihypo]["hypo_along_stk"])
        dhypo = config.CFGDICT[seg_ihypo]["hypo_down_dip"]
        (rup["rupT_dist"][seg_ihypo], _) = self.fm(rup["vr_dist"][seg_ihypo],
                                                   [shypo, dhypo],
                                                   (config.CFGDICT[seg_ihypo]["dlen"],
                                                    config.CFGDICT[seg_ihypo]["dwid"]))

        last_column = list(rup["rupT_dist"][seg_ihypo][:,-1])
        iz_hypo = last_column.index(min(last_column))

        for seg in range(seg_ihypo + 1, num_segments):
            print("Running right segment: %d" % (seg))
            shypo = 0.0
            dhypo = rup["lz"][seg][iz_hypo]
            (rup["rupT_dist"][seg], _) = self.fm(rup["vr_dist"][seg],
                                                 [shypo, dhypo],
                                                 (config.CFGDICT[seg_ihypo]["dlen"],
                                                  config.CFGDICT[seg_ihypo]["dwid"]))
            last_column = rup["rupT_dist"][seg - 1][:,-1]
            rup["rupT_dist"][seg] = (rup["rupT_dist"][seg] +
                                      min(last_column))

            last_column = list(rup["rupT_dist"][seg][:,-1])
            iz_hypo = last_column.index(min(last_column))

        first_column = list(rup["rupT_dist"][seg_ihypo][:,0])
        iz_hypo = first_column.index(min(first_column))

        seg_nl = (seg_ihypo + 1) - 1

        for seg in range(0, seg_nl):
            print("Running left segment: %d" % (seg))
            seg_idx = seg_ihypo - seg - 1
            shypo = config.CFGDICT[seg_idx]["fault_length"]
            dhypo = rup["lz"][seg_idx][iz_hypo]
            (rup["rupT_dist"][seg_idx], _) = self.fm(rup["vr_dist"][seg_idx],
                                                     [shypo, dhypo],
                                                     (config.CFGDICT[seg_ihypo]["dlen"],
                                                      config.CFGDICT[seg_ihypo]["dwid"]))
            first_column = rup["rupT_dist"][seg_ihypo - seg][:,0]
            rup["rupT_dist"][seg_idx] = (rup["rupT_dist"][seg_idx] +
                                         min(first_column))

            first_column = list(rup["rupT_dist"][seg_idx][:,0])
            iz_hypo = first_column.index(min(first_column))

        # Done, check how long it took
        t2 = time.time()

        print("=> Elapsed time for computing rupture time distribution: "
              "%10.2f" % (t2 - t1))

        t1 = time.time()

        # b) generating slip velocity functions (SVF)
        # Initialize data structures
        rup["svf_nt"] = []
        #rup["psv1_dist"] = []
        rup["svf_svf"] = []
        rup["svf_time"] = []

        for seg in range(num_segments):
            rup["svf_nt"].append(np.zeros((rup["nz"][seg],
                                           rup["nx"][seg])))
            #rup["psv1_dist"].append(np.zeros((rup["nz"][seg],
            #                                  rup["nx"][seg])))
            rup["svf_svf"].append(np.empty((rup["nz"][seg],
                                            rup["nx"][seg]),
                                           dtype=object))
            rup["svf_time"].append(np.empty((rup["nz"][seg],
                                             rup["nx"][seg]),
                                            dtype=object))

        for seg in range(num_segments):
            for k in range(0, int(rup["nz"][seg])):
                for i in range(0, int(rup["nx"][seg])):
                    rup["svf_nt"][seg][k][i] = math.ceil((rup["rist_dist"][seg][k][i]+0.5) /
                                                         config.svf_dt)
                    (rup["svf_svf"][seg][k][i],
                     rup["svf_time"][seg][k][i]) = self.gen_svf(rup["rist_dist"][seg][k][i],
                                                                config.svf_dt,
                                                                rup["svf_nt"][seg][k][i],
                                                                config.svf_type)
                     #rup["psv1_dist"][k][i] = (max(rup["svf_svf"][k][i]) *
                     #                            rup["slip_dist"][k][i])

        # Done, check how long it took
        t2 = time.time()

        print("=> Elapsed time for generating slip velocity functions: "
              "%10.2f" % (t2 - t1))

    def gen_srf(self, outfl):
        """
        Writing an SRF source file (SRF: Standard Rupture Format by
        R. Graves)
        """
        rup = self.rup
        config = self.config
        deg2rad = math.pi / 180.0
        den = -1
        vs = -1
        num_segments = config.num_srcfiles

        print()
        print("### Step IV: Generating SRF files...")
        print("###")

        # Data structures
        subf = {}
        subf["lat"] = []
        subf["lon"] = []
        subf["dep"] = []
        subf["stk"] = []
        subf["dip"] = []
        subf["area"] = []
        subf["tinit"] = []
        subf["dt"] = []
        subf["rak"] = []
        subf["slip1"] = []
        subf["nt1"] = []
        subf["slip2"] = []
        subf["nt2"] = []
        subf["slip3"] = []
        subf["nt3"] = []

        for seg in range(num_segments):
            subf["lat"].append(np.zeros((rup["nz"][seg],
                                         rup["nx"][seg])))
            subf["lon"].append(np.zeros((rup["nz"][seg],
                                         rup["nx"][seg])))
            subf["dep"].append(np.zeros((rup["nz"][seg],
                                         rup["nx"][seg])))
            subf["stk"].append(np.zeros((rup["nz"][seg],
                                         rup["nx"][seg])))
            subf["dip"].append(np.zeros((rup["nz"][seg],
                                         rup["nx"][seg])))
            subf["area"].append(np.zeros((rup["nz"][seg],
                                          rup["nx"][seg])))
            subf["tinit"].append(np.zeros((rup["nz"][seg],
                                           rup["nx"][seg])))
            subf["dt"].append(np.zeros((rup["nz"][seg],
                                        rup["nx"][seg])))
            subf["rak"].append(np.zeros((rup["nz"][seg],
                                         rup["nx"][seg])))
            subf["slip1"].append(np.zeros((rup["nz"][seg],
                                           rup["nx"][seg])))
            subf["nt1"].append(np.zeros((rup["nz"][seg],
                                         rup["nx"][seg])))
            subf["slip2"].append(np.zeros((rup["nz"][seg],
                                           rup["nx"][seg])))
            subf["nt2"].append(np.zeros((rup["nz"][seg],
                                         rup["nx"][seg])))
            subf["slip3"].append(np.zeros((rup["nz"][seg],
                                           rup["nx"][seg])))
            subf["nt3"].append(np.zeros((rup["nz"][seg],
                                         rup["nx"][seg])))

        # SRF Data Block Input
        for seg in range(num_segments):
            for i in range(0, int(rup["nx"][seg])):
                for k in range(0, int(rup["nz"][seg])):

                    z_km = (config.CFGDICT[seg]["depth_to_top"] + ((k+1)-0.5) *
                            config.CFGDICT[seg]["dwid"] *
                            math.sin(deg2rad * config.CFGDICT[seg]['dip']))

                    azi = config.CFGDICT[seg]['strike'] + 90 # Fault normal direction
                    length = (((i+1) - 0.5 - rup["nx"][seg] / 2.0) *
                              config.CFGDICT[seg]["dlen"])

                    x_km = (length * math.sin(deg2rad * config.CFGDICT[seg]['strike']) +
                            ((k+1) - 0.5) * config.CFGDICT[seg]["dwid"] *
                            math.cos(deg2rad * config.CFGDICT[seg]['dip']) *
                            math.sin(deg2rad * azi))
                    y_km = (length * math.cos(deg2rad * config.CFGDICT[seg]['strike']) +
                            ((k+1) - 0.5) * config.CFGDICT[seg]["dwid"] *
                            math.cos(deg2rad * config.CFGDICT[seg]['dip']) *
                            math.cos(deg2rad * azi))

                    rng = math.sqrt(x_km**2 + y_km**2) * 180.0 / 6371.0 / math.pi
                    azi = math.atan2(y_km, x_km)
                    azi = ((azi >= 0 and azi <= math.pi/2) * (math.pi/2 - azi) +
                           (azi > math.pi/2) * (2 * math.pi + math.pi / 2 - azi) +
                           (azi < 0) * (abs(azi) + math.pi/2))
                    azi = azi * 180.0 / math.pi

                    (subf["lat"][seg][k][i],
                     subf["lon"][seg][k][i]) = reckon(config.CFGDICT[seg]["lat_top_center"],
                                                      config.CFGDICT[seg]["lon_top_center"],
                                                      rng, azi)
                    subf["dep"][seg][k][i] = z_km
                    subf["stk"][seg][k][i] = config.CFGDICT[seg]['strike']
                    subf["dip"][seg][k][i] = config.CFGDICT[seg]['dip']
                    # cm * cm
                    subf["area"][seg][k][i] = (config.CFGDICT[seg]["dlen"] * 1E5 *
                                               config.CFGDICT[seg]["dwid"] * 1E5)
                    subf["tinit"][seg][k][i] = rup["rupT_dist"][seg][k][i]
                    subf["dt"][seg][k][i] = config.svf_dt
                    subf["rak"][seg][k][i] = config.CFGDICT[seg]['rake']
                    subf["slip1"][seg][k][i] = rup["slip_dist"][seg][k][i]
                    subf["nt1"][seg][k][i] = rup["svf_nt"][seg][k][i]
                    subf["slip2"][seg][k][i] = 0.0
                    subf["nt2"][seg][k][i] = 0
                    subf["slip3"][seg][k][i] = 0.0
                    subf["nt3"][seg][k][i] = 0

        print()
        print("=> Writing SRF file : %s" % (outfl))
        print()

        outf = open(outfl, 'w')

        # Write Header Block
        outf.write('2.0\n') # Version 2.0
        outf.write('PLANE %d\n' % (num_segments)) # Currently one segment only
        for seg in range(num_segments):
            outf.write('%10.4f %10.4f %7d %7d %9.2f %9.2f\n' %
                       (config.CFGDICT[seg]['lon_top_center'],
                        config.CFGDICT[seg]['lat_top_center'],
                        rup["nx"][seg], rup["nz"][seg],
                        config.CFGDICT[seg]["fault_length"],
                        config.CFGDICT[seg]["fault_width"]))
            outf.write('%10.2f %10.2f %7.2f %7.2f %9.2f\n' %
                       (config.CFGDICT[seg]['strike'],
                        config.CFGDICT[seg]['dip'],
                        config.CFGDICT[seg]["depth_to_top"],
                        config.CFGDICT[seg]["hypo_along_stk"],
                        config.CFGDICT[seg]["hypo_down_dip"]))

        # Write Data Block
        for seg in range(num_segments):
            outf.write('POINTS %7d\n' % (rup["nx"][seg] * rup["nz"][seg]))
            for k in range(0, int(rup["nz"][seg])):
                for i in range(0, int(rup["nx"][seg])):
                    outf.write(' %10.4f %10.4f %7.2f %10.2f '
                               '%10.2f %e %7.2f %7.4f %4d %4d\n' %
                               (subf["lon"][seg][k][i], subf["lat"][seg][k][i],
                                subf["dep"][seg][k][i], subf["stk"][seg][k][i],
                                subf["dip"][seg][k][i], subf["area"][seg][k][i],
                                subf["tinit"][seg][k][i], subf["dt"][seg][k][i],
                                vs, den))
                    outf.write(' %10.2f %10.2f %6d %10.2f %6d %10.2f %6d' %
                               (subf["rak"][seg][k][i], subf["slip1"][seg][k][i],
                                subf["nt1"][seg][k][i], subf["slip2"][seg][k][i],
                                subf["nt2"][seg][k][i], subf["slip3"][seg][k][i],
                                subf["nt3"][seg][k][i]))
                    for j in range(0, int(rup["svf_nt"][seg][k][i])):
                        if j % 6 == 0:
                            outf.write('\n  ')
                        outf.write('%e' % (rup["svf_svf"][seg][k][i][j] *
                                           subf["slip1"][seg][k][i]))
                        outf.write('   ')
                    outf.write('\n')
        outf.close()

    def run(self):
        """
        This function runs the RMG method
        """
        print("SONG Multi-Segment Rupture Generator RMG".center(80, '-'))

        # Load configuration, set sim_id
        install = InstallCfg.getInstance()
        sim_id = self.sim_id

        # Build directory paths
        a_tmpdir = os.path.join(install.A_TMP_DATA_DIR, str(sim_id))
        a_indir = os.path.join(install.A_IN_DATA_DIR, str(sim_id))
        a_outdir = os.path.join(install.A_OUT_DATA_DIR, str(sim_id))
        a_logdir = os.path.join(install.A_OUT_LOG_DIR, str(sim_id))

        # Make sure the output and tmp directories exist
        bband_utils.mkdirs([a_tmpdir, a_indir, a_outdir], print_cmd=False)

        # Now, file paths
        self.log = os.path.join(a_logdir, "%d.song_rmg.log" % (sim_id))
        a_srcfiles = [os.path.join(a_indir,
                                   srcfile) for srcfile in self.r_srcfiles]

        # Path to SRF file
        a_srffile = os.path.join(a_outdir, self.r_srffile)

        # Read SRC file
        self.config = SongRMGMSCfg(a_srcfiles)

        # Run the Song RMG
        self.gen_stats_inp()
        self.gen_src()
        self.gen_stats()
        self.gen_rup()
        self.gen_srf(a_srffile)

        # Copy SRF file to other places
        progstring = "cp %s %s" % (a_srffile, os.path.join(a_tmpdir, self.r_srffile))
        bband_utils.runprog(progstring, print_cmd=False)
        progstring = "cp %s %s" % (a_srffile, os.path.join(a_indir, self.r_srffile))
        bband_utils.runprog(progstring, print_cmd=False)

        # Plot SRF
        plot_srf.run(self.r_srffile, sim_id=self.sim_id)

        print("SONG Multi-Segment Rupture Generator RMG Completed".center(80, '-'))

if __name__ == "__main__":
    print("Testing module: %s" % (os.path.basename(sys.argv[0])))
    RMG = SongRMGMS(None, sys.argv[1], sys.argv[2], None, int(sys.argv[3]))
    RMG.run()
