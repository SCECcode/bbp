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
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys
import math

import numpy as np
import matplotlib
if matplotlib.get_backend() != 'agg':
    matplotlib.use('Agg') # Disable use of Tk/X11
import pylab as py
import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import matplotlib.pyplot as plt
from pylab import arange
from scipy import integrate
from scipy import signal, stats
from scipy.signal import butter, filtfilt

# Import BBP modules
import bband_utils
import plot_config
from install_cfg import InstallCfg
from station_list import StationList

def school_round(value):
    """
    Implements the school round where round(0.5) = 1.0
    """
    if (float(value) % 1) >= 0.5:
        if value > 0:
            return float(math.ceil(value))
        elif value < 0:
            return float(math.floor(value))

    return round(value)

class AndersonGOF(object):
    """
    This class implements the Anderson GoF
    """
    g = 981.
    NMAX = 300

    B1 = [0.1, 0.2]
    B2 = [0.2, 0.5]
    B3 = [0.5, 1.0]
    B4 = [1.0, 2.0]
    B5 = [2.0, 5.0]
    B6 = [5.0, 10.0]
    B7 = [10.0, 20.0]
    B8 = [20.0, 50.0]
    B9 = [50.0, 100.0]
    B10 = [0.1, 100.0]

    B = [B1, B2, B3, B4, B5, B6, B7, B8, B9, B10]
    BNAMES = ['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B9', 'B10']
    BMAX = len(B)

    C1 = np.empty([NMAX, BMAX])
    C1[:] = np.nan
    C2 = np.empty([NMAX, BMAX])
    C2[:] = np.nan
    C3 = np.empty([NMAX, BMAX])
    C3[:] = np.nan
    C4 = np.empty([NMAX, BMAX])
    C4[:] = np.nan
    C5 = np.empty([NMAX, BMAX])
    C5[:] = np.nan
    C6 = np.empty([NMAX, BMAX])
    C6[:] = np.nan
    C7 = np.empty([NMAX, BMAX])
    C7[:] = np.nan
    C8 = np.empty([NMAX, BMAX])
    C8[:] = np.nan
    C9 = np.empty([NMAX, BMAX])
    C9[:] = np.nan
    C10 = np.empty([NMAX, BMAX])
    C10[:] = np.nan

    S1 = np.empty(NMAX)
    S1[:] = np.nan

    def __init__(self, stations, eventname, sim_id=0):
        """
        Initilize class variables
        """
        self.eventname = eventname
        self.stations = stations
        self.dt = None
        self.log = None
        self.sim_id = sim_id

    @staticmethod
    def eval_func(p1, p2):
        return 10. * np.exp(-np.power(((p1 - p2) / min(p1, p2)), 2.))

    @staticmethod
    def eval_func2(p1, p2):
        return 10. * np.exp(-np.power((np.abs(p1 - p2) / min(p1, p2)), 2.))

    @staticmethod
    def integral(x, idt):
        return integrate.simps(x, x=None, dx=idt)

    @staticmethod
    def integral2(x, idt):
        return integrate.simps(x * x, x=None, dx=idt)

    @staticmethod
    def integ(array_in, idt):
        return integrate.cumtrapz(array_in, dx=idt)

    @staticmethod
    def shift_bit_length(x):
        return 1 << (x - 1).bit_length()

    @staticmethod
    def align_seismograms(obs_times, obs_data,
                          sym_times, sym_data):
        """
        This function aligns the observed and calculated seismograms. It
        first brings both seismograms to the same dt by downsampling as
        needed, and then truncates the ground motion according to
        cross-correlation.

        If should be called as:
        (new_obs_times, new_obs_data,
         new_sym_times, new_sym_data) = align_seismograms(obs_times, obs_data,
                                                          sym_times, sym_data)
        """

        # Calculate dt for each timeseries
        tdeltat = obs_times[1] - obs_times[0]
        deltat = sym_times[1] - sym_times[0]

        # Initial timeseries from files
        tfi = obs_data
        t = obs_times
        xii = sym_data
        tt1 = sym_times

        # Now downsample as needed
        if deltat > tdeltat:
            # Filter tfi
            X5, Y5 = butter(8, (1.0 / deltat) / (1.0 / tdeltat), 'low')
            tfif = filtfilt(X5, Y5, tfi,
                            padlen=3*(max(len(X5), len(Y5))-1))
            idx = [int(school_round(x)) for x in np.arange(1,
                                                           len(tfif) + 0.0000001,
                                                           (deltat/tdeltat))]
            tf = [tfif[i-1] for i in idx]
            # Keep xii
            xi = xii
        elif deltat < tdeltat:
            # Filter xii
            X5, Y5 = butter(8, (1.0 / tdeltat) / (1.0 / deltat), 'low')
            xif = filtfilt(X5, Y5, xii,
                           padlen=3*(max(len(X5), len(Y5))-1))
            idx = [int(round(x)) for x in np.arange(1,
                                                    len(xif) + 0.0000001,
                                                    (tdeltat/deltat))]
            xi = [xif[i-1] for i in idx]
            # Keep tfi
            tf = tfi
        else:
            tf = tfi
            xi = xii

        # Truncate simulated ground motion according to cross-correlation
        tq2_integ = np.zeros((1, len(tf)))[0]
        tq2_integ[0] = 0

        for r in range(1, len(tf)):
            tq2_integ[r] = tq2_integ[r-1] + (tf[r]**2)
        tq2_integ = tq2_integ * deltat
        cav_tq2_integ_total = tq2_integ[-1]

        oq2_integ = np.zeros((1, len(xi)))[0]
        oq2_integ[0] = 0

        for r in range(1, len(xi)):
            oq2_integ[r] = oq2_integ[r-1] + (xi[r]**2)
        oq2_integ = oq2_integ * tdeltat
        cav_oq2_integ_total = oq2_integ[-1]

        for cav_tp in range(0, len(tf)):
            if tq2_integ[cav_tp] >= (0.05 * cav_tq2_integ_total):
                cav_tii = cav_tp
                break

        for cav_p in range(0, len(xi)):
            if oq2_integ[cav_p] >= (0.05 * cav_oq2_integ_total):
                cav_ii = cav_p
                break

        v_sp = cav_ii - cav_tii
        if v_sp >= 0:
            ii = v_sp
        else:
            ii = 0

        if len(xi) < ii + len(tf):
            #if len(xi) < len(tf):
            #    tf = tf[0:len(xi)]
            new_max = len(xi) - ii
            tf = tf[0:new_max]
            x = xi[ii:]
            # x = xi[0:len(tf)]
        else:
            x = xi[ii:(ii + len(tf))]

        ttt = np.linspace(0, len(tf) * deltat, len(tf))

        # Return time, recorded_data, time, simulated_data
        return ttt, tf, ttt, x

    def c13_eval(self, acc1, acc2):
        """
        ************** Arias duration and intensity ****************
        The time series must have the same duration, the same origin
         time and must have the same sample rate
        """
        IA1t = (0.5 * np.pi / self.g) * self.dt * np.cumsum(acc1*acc1)
        IA2t = (0.5 * np.pi / self.g) * self.dt * np.cumsum(acc2*acc2)
        IA1 = IA1t[(len(IA1t)-1)]
        IA2 = IA2t[(len(IA2t)-1)]
        N1 = IA1t / IA1
        N2 = IA2t / IA2
        F = np.absolute(N1 - N2)
        return max(0, 10. * (1. - np.amax(F))), self.eval_func(IA1, IA2)

    def c24_eval(self, vel1, vel2):
        """
        ********************* Energy duration *********************
        """
        IE1t = self.dt * np.cumsum(vel1*vel1)
        IE2t = self.dt * np.cumsum(vel2*vel2)
        IE1 = IE1t[(len(IE1t)-1)]
        IE2 = IE2t[(len(IE2t)-1)]
        N1 = IE1t / IE1
        N2 = IE2t / IE2
        F = np.absolute(N1 - N2)
        return max(0, 10. * (1. - np.amax(F))), self.eval_func(IE1, IE2)

    def c5_eval(self, acc1, acc2):
        """
        ********************* Peak Acceleration *********************
        """
        pga1 = np.amax(acc1)
        pga2 = np.amax(acc2)
        return self.eval_func(pga1, pga2)

    def c6_eval(self, vel1, vel2):
        """
        ********************* Peak Velocity *********************
        """
        pgv1 = np.amax(vel1)
        pgv2 = np.amax(vel2)
        return self.eval_func(pgv1, pgv2)

    def c7_eval(self, dis1, dis2):
        """
        ********************* Peak Displacement *********************
        """
        pgd1 = np.amax(dis1)
        pgd2 = np.amax(dis2)
        return self.eval_func(pgd1, pgd2)

    def c8_eval(self, rs1, rs2, period):
        """
        ********************* Response Spectra *********************
        """
        c = [self.eval_func2(rs1[i], rs2[i]) for i in range(len(period))]
        return np.nanmean(c)

    def c9_eval(self, fs1, fs2, f):
        """
        ********************* Fourier Spectra *********************
        """
        c = [self.eval_func2(fs1[i], fs2[i]) for i in range(len(f))]
        return np.nanmean(c)

    def c10_eval(self, acc1, acc2):
        """
        ********************* Cross Correlation *********************
        """
        ICC1 = np.sqrt(self.integral2(acc1, self.dt))
        ICC2 = np.sqrt(self.integral2(acc2, self.dt))
        ICC12 = self.integral(acc1 * acc2, self.dt)
        return 10. * max(ICC12 / (ICC1 * ICC2), 0.)

    def cplots(self, irec, station_name, output_file):
        """
        This function creates a per-station plot summarizing all 10
        metrics from the Anderson GoF
        """

        my_s1 = self.S1[irec]
        my_c1 = self.C1[irec, :]
        my_c2 = self.C2[irec, :]
        my_c3 = self.C3[irec, :]
        my_c4 = self.C4[irec, :]
        my_c5 = self.C5[irec, :]
        my_c6 = self.C6[irec, :]
        my_c7 = self.C7[irec, :]
        my_c8 = self.C8[irec, :]
        my_c9 = self.C9[irec, :]
        my_c10 = self.C10[irec, :]

        fig = py.figure()
        my_s1 = "{:3.1f}".format(my_s1)
        title = "%s - %s - Score S1 : %s" % (self.eventname,
                                             station_name,
                                             str(my_s1))
        fig.suptitle(title, fontsize=18)
        x_vals = np.arange(1, len(self.B) + 1)
        label_size = 9
        matplotlib.rcParams['xtick.labelsize'] = label_size
        matplotlib.rcParams['ytick.labelsize'] = label_size

        plt.subplot(5, 2, 1)
        plt.grid(True, which='both', color='0.35', ls=":")
        plt.xlim(0, 11)
        plt.xticks(arange(1, 11))
        plt.yticks(arange(0, 12, 2))
        plt.ylim(0, 10)
        plt.ylabel('C1', fontsize=12, labelpad=-1)
        plt.plot(x_vals, my_c1, 'ro', ms=8)
        plt.plot(x_vals, my_c1, 'b-')

        plt.subplot(5, 2, 2)
        plt.grid(True, which='both', color='0.35', ls=":")
        plt.xlim(0, 11)
        plt.ylim(0, 10)
        plt.xticks(arange(1, 11))
        plt.yticks(arange(0, 12, 2))
        plt.ylabel('C2', fontsize=12, labelpad=-1)
        plt.plot(x_vals, my_c2, 'ro', ms=8)
        plt.plot(x_vals, my_c2, 'b-')

        plt.subplot(5, 2, 3)
        plt.grid(True, which='both', color='0.35', ls=":")
        plt.xlim(0, 11)
        plt.ylim(0, 10)
        plt.xticks(arange(1, 11))
        plt.yticks(arange(0, 12, 2))
        plt.ylabel('C3', fontsize=12, labelpad=-1)
        plt.plot(x_vals, my_c3, 'ro', ms=8)
        plt.plot(x_vals, my_c3, 'b-')

        plt.subplot(5, 2, 4)
        plt.grid(True, which='both', color='0.35', ls=":")
        plt.xlim(0, 11)
        plt.ylim(0, 10)
        plt.xticks(arange(1, 11))
        plt.yticks(arange(0, 12, 2))
        plt.ylabel('C4', fontsize=12, labelpad=-1)
        plt.plot(x_vals, my_c4, 'ro', ms=8)
        plt.plot(x_vals, my_c4, 'b-')

        plt.subplot(5, 2, 5)
        plt.grid(True, which='both', color='0.35', ls=":")
        plt.xlim(0, 11)
        plt.ylim(0, 10)
        plt.xticks(arange(1, 11))
        plt.yticks(arange(0, 12, 2))
        plt.ylabel('C5', fontsize=12, labelpad=-1)
        plt.plot(x_vals, my_c5, 'ro', ms=8)
        plt.plot(x_vals, my_c5, 'b-')

        plt.subplot(5, 2, 6)
        plt.grid(True, which='both', color='0.35', ls=":")
        plt.xlim(0, 11)
        plt.ylim(0, 10)
        plt.xticks(arange(1, 11))
        plt.yticks(arange(0, 12, 2))
        plt.ylabel('C6', fontsize=12, labelpad=-1)
        plt.plot(x_vals, my_c6, 'ro', ms=8)
        plt.plot(x_vals, my_c6, 'b-')

        plt.subplot(5, 2, 7)
        plt.grid(True, which='both', color='0.35', ls=":")
        plt.xlim(0, 11)
        plt.ylim(0, 10)
        plt.xticks(arange(1, 11))
        plt.yticks(arange(0, 12, 2))
        plt.ylabel('C7', fontsize=12, labelpad=-1)
        plt.plot(x_vals, my_c7, 'ro', ms=8)
        plt.plot(x_vals, my_c7, 'b-')

        plt.subplot(5, 2, 8)
        plt.grid(True, which='both', color='0.35', ls=":")
        plt.xlim(0, 11)
        plt.ylim(0, 10)
        plt.xticks(arange(1, 11))
        plt.yticks(arange(0, 12, 2))
        plt.ylabel('C8', fontsize=12, labelpad=-1)
        plt.plot(x_vals, my_c8, 'ro', ms=8)
        plt.plot(x_vals, my_c8, 'b-')

        plt.subplot(5, 2, 9)
        plt.grid(True, which='both', color='0.35', ls=":")
        plt.xlim(0, 11)
        plt.ylim(0, 10)
        plt.xticks(arange(1, 11))
        plt.yticks(arange(0, 12, 2))
        plt.ylabel('C9', fontsize=12, labelpad=-1)
        plt.plot(x_vals, my_c9, 'ro', ms=8)
        plt.xlabel('Frequency Band', fontsize=12)
        plt.plot(x_vals, my_c9, 'b-')

        plt.subplot(5, 2, 10)
        plt.grid(True, which='both', color='0.35', ls=":")
        plt.xlim(0, 11)
        plt.ylim(0, 10)
        plt.xticks(arange(1, 11))
        plt.yticks(arange(0, 12, 2))
        plt.ylabel('C10', fontsize=12, labelpad=-1)
        plt.plot(x_vals, my_c10, 'ro', ms=8)
        plt.plot(x_vals, my_c10, 'b-')
        plt.xlabel('Frequency Band', fontsize=12)

    #   plt.show()
        plt.savefig(output_file, dpi=plot_config.dpi)
        plt.close(fig)

        return

# ****************************************************************************
    def fplots(self, S1, c1cf, c2cf, c3cf, c4cf, c5cf,
               c6cf, c7cf, c8cf, c9cf, c10cf, output_file):
        """
        This function creates a combined GoF plot for all stations
        """

        fig = py.figure()
        S1 = "{:3.1f}".format(S1)
        title = self.eventname + ' - Score S1 : ' + str(S1)
        fig.suptitle(title, fontsize=18)
        x_vals = np.arange(1, len(self.B) + 1)
        label_size = 9
        matplotlib.rcParams['xtick.labelsize'] = label_size
        matplotlib.rcParams['ytick.labelsize'] = label_size

        mean = [c1cf[i][0] for i in arange(self.BMAX)]
        # stdev = [c1cf[i][1] for i in arange(self.BMAX)]
        low_std = [c1cf[i][0] - c1cf[i][1] for i in arange(self.BMAX)]
        upp_std = [c1cf[i][0] + c1cf[i][1] for i in arange(self.BMAX)]
        low_conf = [c1cf[i][2] for i in arange(self.BMAX)]
        upp_conf = [c1cf[i][3] for i in arange(self.BMAX)]
        plt.subplot(5, 2, 1)
        plt.grid(True, which='both', color='0.35', ls=":")
        plt.xlim(0, 11)
        plt.xticks(arange(1, 11))
        plt.yticks(arange(0, 12, 2))
        plt.ylim(0, 10)
        plt.ylabel('C1', fontsize=12, labelpad=-1)
        plt.plot(x_vals, mean, 'ro', ms=4)
        plt.plot(x_vals, mean, 'r-', mew=4)
        plt.fill_between(x_vals, low_conf, upp_conf,
                         color='cyan', alpha='1.0')
        plt.fill_between(x_vals, low_std, upp_std,
                         color='yellow', alpha='1.0')

        mean = [c2cf[i][0] for i in arange(self.BMAX)]
        # stdev = [c2cf[i][1] for i in arange(self.BMAX)]
        low_std = [c2cf[i][0] - c2cf[i][1] for i in arange(self.BMAX)]
        upp_std = [c2cf[i][0] + c2cf[i][1] for i in arange(self.BMAX)]
        low_conf = [c2cf[i][2] for i in arange(self.BMAX)]
        upp_conf = [c2cf[i][3] for i in arange(self.BMAX)]
        plt.subplot(5, 2, 2)
        plt.grid(True, which='both', color='0.35', ls=":")
        plt.xlim(0, 11)
        plt.ylim(0, 10)
        plt.xticks(arange(1, 11))
        plt.yticks(arange(0, 12, 2))
        plt.ylabel('C2', fontsize=12, labelpad=-1)
        plt.plot(x_vals, mean, 'ro', ms=4)
        plt.plot(x_vals, mean, 'r-', mew=4)
        plt.fill_between(x_vals, low_conf, upp_conf,
                         color='cyan', alpha='1.0')
        plt.fill_between(x_vals, low_std, upp_std,
                         color='yellow', alpha='1.0')

        mean = [c3cf[i][0] for i in arange(self.BMAX)]
        # stdev = [c3cf[i][1] for i in arange(self.BMAX)]
        low_std = [c3cf[i][0] - c3cf[i][1] for i in arange(self.BMAX)]
        upp_std = [c3cf[i][0] + c3cf[i][1] for i in arange(self.BMAX)]
        low_conf = [c3cf[i][2] for i in arange(self.BMAX)]
        upp_conf = [c3cf[i][3] for i in arange(self.BMAX)]
        plt.subplot(5, 2, 3)
        plt.grid(True, which='both', color='0.35', ls=":")
        plt.xlim(0, 11)
        plt.ylim(0, 10)
        plt.xticks(arange(1, 11))
        plt.yticks(arange(0, 12, 2))
        plt.ylabel('C3', fontsize=12, labelpad=-1)
        plt.plot(x_vals, mean, 'ro', ms=4)
        plt.plot(x_vals, mean, 'r-', mew=4)
        plt.fill_between(x_vals, low_conf, upp_conf,
                         color='cyan', alpha='1.0')
        plt.fill_between(x_vals, low_std, upp_std,
                         color='yellow', alpha='1.0')

        mean = [c4cf[i][0] for i in arange(self.BMAX)]
        # stdev = [c4cf[i][1] for i in arange(self.BMAX)]
        low_std = [c4cf[i][0] - c4cf[i][1] for i in arange(self.BMAX)]
        upp_std = [c4cf[i][0] + c4cf[i][1] for i in arange(self.BMAX)]
        low_conf = [c4cf[i][2] for i in arange(self.BMAX)]
        upp_conf = [c4cf[i][3] for i in arange(self.BMAX)]
        plt.subplot(5, 2, 4)
        plt.grid(True, which='both', color='0.35', ls=":")
        plt.xlim(0, 11)
        plt.ylim(0, 10)
        plt.xticks(arange(1, 11))
        plt.yticks(arange(0, 12, 2))
        plt.ylabel('C4', fontsize=12, labelpad=-1)
        plt.plot(x_vals, mean, 'ro', ms=4)
        plt.plot(x_vals, mean, 'r-', mew=4)
        plt.fill_between(x_vals, low_conf, upp_conf,
                         color='cyan', alpha='1.0')
        plt.fill_between(x_vals, low_std, upp_std,
                         color='yellow', alpha='1.0')

        mean = [c5cf[i][0] for i in arange(self.BMAX)]
        # stdev = [c5cf[i][1] for i in arange(self.BMAX)]
        low_std = [c5cf[i][0] - c5cf[i][1] for i in arange(self.BMAX)]
        upp_std = [c5cf[i][0] + c5cf[i][1] for i in arange(self.BMAX)]
        low_conf = [c5cf[i][2] for i in arange(self.BMAX)]
        upp_conf = [c5cf[i][3] for i in arange(self.BMAX)]
        plt.subplot(5, 2, 5)
        plt.grid(True, which='both', color='0.35', ls=":")
        plt.xlim(0, 11)
        plt.ylim(0, 10)
        plt.xticks(arange(1, 11))
        plt.yticks(arange(0, 12, 2))
        plt.ylabel('C5', fontsize=12, labelpad=-1)
        plt.plot(x_vals, mean, 'ro', ms=4)
        plt.plot(x_vals, mean, 'r-', mew=4)
        plt.fill_between(x_vals, low_conf, upp_conf,
                         color='cyan', alpha='1.0')
        plt.fill_between(x_vals, low_std, upp_std,
                         color='yellow', alpha='1.0')

        mean = [c6cf[i][0] for i in arange(self.BMAX)]
        # stdev = [c6cf[i][1] for i in arange(self.BMAX)]
        low_std = [c6cf[i][0] - c6cf[i][1] for i in arange(self.BMAX)]
        upp_std = [c6cf[i][0] + c6cf[i][1] for i in arange(self.BMAX)]
        low_conf = [c6cf[i][2] for i in arange(self.BMAX)]
        upp_conf = [c6cf[i][3] for i in arange(self.BMAX)]
        plt.subplot(5, 2, 6)
        plt.grid(True, which='both', color='0.35', ls=":")
        plt.xlim(0, 11)
        plt.ylim(0, 10)
        plt.xticks(arange(1, 11))
        plt.yticks(arange(0, 12, 2))
        plt.ylabel('C6', fontsize=12, labelpad=-1)
        plt.plot(x_vals, mean, 'ro', ms=4)
        plt.plot(x_vals, mean, 'r-', mew=4)
        plt.fill_between(x_vals, low_conf, upp_conf,
                         color='cyan', alpha='1.0')
        plt.fill_between(x_vals, low_std, upp_std,
                         color='yellow', alpha='1.0')

        mean = [c7cf[i][0] for i in arange(self.BMAX)]
        # stdev = [c7cf[i][1] for i in arange(self.BMAX)]
        low_std = [c7cf[i][0] - c7cf[i][1] for i in arange(self.BMAX)]
        upp_std = [c7cf[i][0] + c7cf[i][1] for i in arange(self.BMAX)]
        low_conf = [c7cf[i][2] for i in arange(self.BMAX)]
        upp_conf = [c7cf[i][3] for i in arange(self.BMAX)]
        plt.subplot(5, 2, 7)
        plt.grid(True, which='both', color='0.35', ls=":")
        plt.xlim(0, 11)
        plt.ylim(0, 10)
        plt.xticks(arange(1, 11))
        plt.yticks(arange(0, 12, 2))
        plt.ylabel('C7', fontsize=12, labelpad=-1)
        plt.plot(x_vals, mean, 'ro', ms=4)
        plt.plot(x_vals, mean, 'r-', mew=4)
        plt.fill_between(x_vals, low_conf, upp_conf,
                         color='cyan', alpha='1.0')
        plt.fill_between(x_vals, low_std, upp_std,
                         color='yellow', alpha='1.0')

        mean = [c8cf[i][0] for i in arange(self.BMAX)]
        # stdev = [c8cf[i][1] for i in arange(self.BMAX)]
        low_std = [c8cf[i][0] - c8cf[i][1] for i in arange(self.BMAX)]
        upp_std = [c8cf[i][0] + c8cf[i][1] for i in arange(self.BMAX)]
        low_conf = [c8cf[i][2] for i in arange(self.BMAX)]
        upp_conf = [c8cf[i][3] for i in arange(self.BMAX)]
        plt.subplot(5, 2, 8)
        plt.grid(True, which='both', color='0.35', ls=":")
        plt.xlim(0, 11)
        plt.ylim(0, 10)
        plt.xticks(arange(1, 11))
        plt.yticks(arange(0, 12, 2))
        plt.ylabel('C8', fontsize=12, labelpad=-1)
        plt.plot(x_vals, mean, 'ro', ms=4)
        plt.plot(x_vals, mean, 'r-', mew=4)
        plt.fill_between(x_vals, low_conf, upp_conf,
                         color='cyan', alpha='1.0')
        plt.fill_between(x_vals, low_std, upp_std,
                         color='yellow', alpha='1.0')

        mean = [c9cf[i][0] for i in arange(self.BMAX)]
        # stdev = [c9cf[i][1] for i in arange(self.BMAX)]
        low_std = [c9cf[i][0] - c9cf[i][1] for i in arange(self.BMAX)]
        upp_std = [c9cf[i][0] + c9cf[i][1] for i in arange(self.BMAX)]
        low_conf = [c9cf[i][2] for i in arange(self.BMAX)]
        upp_conf = [c9cf[i][3] for i in arange(self.BMAX)]
        plt.subplot(5, 2, 9)
        plt.grid(True, which='both', color='0.35', ls=":")
        plt.xlim(0, 11)
        plt.ylim(0, 10)
        plt.xticks(arange(1, 11))
        plt.yticks(arange(0, 12, 2))
        plt.ylabel('C9', fontsize=12, labelpad=-1)
        plt.plot(x_vals, mean, 'ro', ms=4)
        plt.xlabel('Frequency Band', fontsize=12)
        plt.plot(x_vals, mean, 'r-', mew=4)
        plt.fill_between(x_vals, low_conf, upp_conf,
                         color='cyan', alpha='1.0')
        plt.fill_between(x_vals, low_std, upp_std,
                         color='yellow', alpha='1.0')

        mean = [c10cf[i][0] for i in arange(self.BMAX)]
        # stdev = [c10cf[i][1] for i in arange(self.BMAX)]
        low_std = [c10cf[i][0] - c10cf[i][1] for i in arange(self.BMAX)]
        upp_std = [c10cf[i][0] + c10cf[i][1] for i in arange(self.BMAX)]
        low_conf = [c10cf[i][2] for i in arange(self.BMAX)]
        upp_conf = [c10cf[i][3] for i in arange(self.BMAX)]
        plt.subplot(5, 2, 10)
        plt.grid(True, which='both', color='0.35', ls=":")
        plt.xlim(0, 11)
        plt.ylim(0, 10)
        plt.xticks(arange(1, 11))
        plt.yticks(arange(0, 12, 2))
        plt.ylabel('C10', fontsize=12, labelpad=-1)
        plt.plot(x_vals, mean, 'ro', ms=4)
        plt.plot(x_vals, mean, 'r-', mew=4)
        plt.xlabel('Frequency Band', fontsize=12)
        plt.fill_between(x_vals, low_conf, upp_conf,
                         color='cyan', alpha='1.0')
        plt.fill_between(x_vals, low_std, upp_std,
                         color='yellow', alpha='1.0')

        #plt.show()
        plt.savefig(output_file, dpi=plot_config.dpi)
        plt.close(fig)

        return

    @staticmethod
    def padts(acc1, acc2, acc3, acc4, lowcut):
        """
        This function pads the data
        """

        padlen = int(1.5 * 8 / lowcut)
        ndata1 = len(acc1)
        ndata2 = len(acc2)
        ndata3 = len(acc3)
        ndata4 = len(acc4)
        datamax = max(ndata1, ndata2, ndata3, ndata4)
    #   imax = [i for i, j in enumerate(ndata) if j == datamax]
        ndata = AndersonGOF.shift_bit_length(datamax + padlen)

        return ndata

    @staticmethod
    def smcpadf(y1, y2, y3, y4, dt, flc, flc_slope, fhc, fhc_slope,
                replace_discarded_data_with_zeros):
        """
        ********************* zero pad time series *********************
        """
        small = 1.0e-10

        npts_in1 = len(y1)
        npts_in2 = len(y2)
        npts_in3 = len(y3)
        npts_in4 = len(y4)

        if flc == 0.0 and fhc == 0.0:   # no filtering
            npad = 0
        else:                           # other filtering
            npad1 = 3. * flc_slope / (flc * dt)
            npad2 = 6. * fhc_slope / (np.abs((fhc - flc)) * dt)
            npad = max(npad1, npad2)

        nlz = int(npad/2)
        ntz = int(npad/2)
        izcl = nlz
        izct = ntz

        # Get peak motion:
        peak_mtn1 = np.amax(y1)
        peak_mtn2 = np.amax(y2)
        peak_mtn3 = np.amax(y3)
        peak_mtn4 = np.amax(y4)

        # Define fraction small of peak:
        smidge1 = small * np.abs(peak_mtn1)
        smidge2 = small * np.abs(peak_mtn2)
        smidge3 = small * np.abs(peak_mtn3)
        smidge4 = small * np.abs(peak_mtn4)

        zerocross1 = [i for i, x in enumerate(y1) if x <= smidge1]
        zerocross2 = [i for i, x in enumerate(y2) if x <= smidge2]
        zerocross3 = [i for i, x in enumerate(y3) if x <= smidge3]
        zerocross4 = [i for i, x in enumerate(y4) if x <= smidge4]
        # izcl, izct = index in input smc file of zero crossings
        izcl1 = zerocross1[0]
        izct1 = zerocross1[len(zerocross1) - 1]
        izcl2 = zerocross2[0]
        izct2 = zerocross2[len(zerocross2) - 1]
        izcl3 = zerocross3[0]
        izct3 = zerocross3[len(zerocross3) - 1]
        izcl4 = zerocross4[0]
        izct4 = zerocross4[len(zerocross4) - 1]

        # Shift the time series to make room for the leading zeros:
        npts_snip1 = izct1 - izcl1 + 1 # number of original points between zeros
        npts_snip2 = izct2 - izcl2 + 1
        npts_snip3 = izct3 - izcl3 + 1
        npts_snip4 = izct4 - izcl4 + 1

        if replace_discarded_data_with_zeros:
            nlznew1 = nlz + izcl1
            ntznew1 = ntz + npts_in1 - izct1 - 1
            nlznew2 = nlz + izcl2
            ntznew2 = ntz + npts_in2 - izct2 - 1
            nlznew3 = nlz + izcl3
            ntznew3 = ntz + npts_in3 - izct3 - 1
            nlznew4 = nlz + izcl4
            ntznew4 = ntz + npts_in4 - izct4 - 1
        else:
            nlznew1 = nlz
            ntznew1 = ntz
            nlznew2 = nlz
            ntznew2 = ntz
            nlznew3 = nlz
            ntznew3 = ntz
            nlznew4 = nlz
            ntznew4 = ntz

        # Compute new npts_out:
        npts_out1 = AndersonGOF.shift_bit_length(npts_snip1 + nlznew1 + ntznew1)
        npts_out2 = AndersonGOF.shift_bit_length(npts_snip2 + nlznew2 + ntznew2)
        npts_out3 = AndersonGOF.shift_bit_length(npts_snip3 + nlznew3 + ntznew3)
        npts_out4 = AndersonGOF.shift_bit_length(npts_snip4 + nlznew4 + ntznew4)

        npts_out = max(npts_out1, npts_out2, npts_out3, npts_out4)

        y_out1 = np.zeros(npts_out)
        y_out2 = np.zeros(npts_out)
        y_out3 = np.zeros(npts_out)
        y_out4 = np.zeros(npts_out)

        y_out1[nlznew1:nlznew1+npts_snip1 - 1] = y1[izcl1:izct1]
        y_out2[nlznew2:nlznew2+npts_snip2 - 1] = y2[izcl2:izct2]
        y_out3[nlznew3:nlznew3+npts_snip3 - 1] = y3[izcl3:izct3]
        y_out4[nlznew4:nlznew4+npts_snip4 - 1] = y4[izcl4:izct4]

        return y_out1, y_out2, y_out3, y_out4, npts_out

    @staticmethod
    def statts(x_in):
        """
        Calculates statistics
        """
        # Filters RuntimeWarning from python's stats package
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            a = np.asarray(x_in)
            s = a[~np.isnan(a)]
            if len(s) > 1:
                mean = np.mean(s)
                stdev = np.std(s, dtype=np.float32, ddof=1)
                conf = stats.t.interval(0.95, len(s), loc=mean, scale=stdev)
                m95 = conf[0]
                p95 = conf[1]
            else:
                mean = 0.
                stdev = 0.
                m95 = 0.
                p95 = 0.

        dist_stat = (mean, stdev, m95, p95)

        return dist_stat

    @staticmethod
    def butter_bandpass(lowcut, highcut, nyq, y, order):
        """
        This function implements a butter bandpass filter
        """

        low = lowcut / nyq
        high = min(highcut / nyq, 0.9999)
        b, a = signal.butter(order, [low, high], btype='bandpass')
        # y = signal.lfilter(b, a, y)
        y = signal.filtfilt(b, a, y)

        return y

    def run(self):
        """
        Runs the Anderson GoF code
        """
        print("Anderson GoF".center(80, '-'))

        # Load configuration, set sim_id
        install = InstallCfg.getInstance()
        sim_id = self.sim_id

        # Build directory paths
        a_tmpdir = os.path.join(install.A_TMP_DATA_DIR, str(sim_id))
        a_indir = os.path.join(install.A_IN_DATA_DIR, str(sim_id))
        a_outdir = os.path.join(install.A_OUT_DATA_DIR, str(sim_id))
        a_logdir = os.path.join(install.A_OUT_LOG_DIR, str(sim_id))
        a_validation_outdir = os.path.join(a_outdir,
                                           "validations",
                                           "anderson_gof")

        # Make sure the output and tmp directories exist
        bband_utils.mkdirs([a_tmpdir, a_indir, a_outdir, a_validation_outdir],
                           print_cmd=False)

        # Now the file paths
        self.log = os.path.join(a_logdir, "%d.anderson.log" % (sim_id))
        sta_file = os.path.join(a_indir, self.stations)
        sta_base = os.path.basename(os.path.splitext(self.stations)[0])
        sims_dir = a_outdir
        obs_dir = os.path.join(a_tmpdir, "obs_seis_%s" % (sta_base))

        # Start with first record
        irec = 0

        # Read station list
        slo = StationList(sta_file)
        site_list = slo.getStationList()

        # Figure out station names
        station_names = []
        for station in site_list:
            station_names.append(station.scode)

        # Loop over stations
        for site in site_list:
            station = site.scode

            print("==> Processing station: %s" % (station))

            file_sims_acc = os.path.join(sims_dir, "%d.%s.acc.bbp" %
                                         (sim_id, station))
            file_sims_rd50 = os.path.join(sims_dir, "%d.%s.rd50" %
                                          (sim_id, station))
            lowcut = site.low_freq_corner
            highcut = site.high_freq_corner
            #print(lowcut, highcut)

            (sims_acc_org_time, sims_acc_org_ns,
             sims_acc_org_ew, sims_acc_org_ver) = np.genfromtxt(file_sims_acc,
                                                                skip_header=2,
                                                                dtype='float64',
                                                                unpack='TRUE')
            (sims_perd, sims_rd50_ns,
             sims_rd50_ew, sims_rd50_ver) = np.genfromtxt(file_sims_rd50,
                                                          skip_header=2,
                                                          dtype='float64',
                                                          unpack='TRUE')

            file_obs_acc = os.path.join(obs_dir, "%s.bbp" %
                                        (station))
            file_obs_rd50 = os.path.join(obs_dir, "%s.rd50" %
                                         (station))

            (obs_acc_org_time, obs_acc_org_ns,
             obs_acc_org_ew, obs_acc_org_ver) = np.genfromtxt(file_obs_acc,
                                                              skip_header=2,
                                                              dtype='float64',
                                                              unpack='TRUE')
            (obs_perd, obs_rd50_ns,
             obs_rd50_ew, obs_rd50_ver) = np.genfromtxt(file_obs_rd50,
                                                        skip_header=2,
                                                        dtype='float64',
                                                        unpack='TRUE')
            # Intitialize the rd50 arrays
            RD50PER = len(obs_perd)
            rd1 = np.zeros(RD50PER)
            rd2 = np.zeros(RD50PER)
            rd3 = np.zeros(RD50PER)
            rd4 = np.zeros(RD50PER)

            # Resample and align the time series
            (obs_acc_time,
             obs_acc_ew,
             sims_acc_time,
             sims_acc_ew) = self.align_seismograms(obs_acc_org_time,
                                                   obs_acc_org_ew,
                                                   sims_acc_org_time,
                                                   sims_acc_org_ew)

            (obs_acc_time,
             obs_acc_ns,
             sims_acc_time,
             sims_acc_ns) = self.align_seismograms(obs_acc_org_time,
                                                   obs_acc_org_ns,
                                                   sims_acc_org_time,
                                                   sims_acc_org_ns)

            obs_org_dt = obs_acc_org_time[1] - obs_acc_org_time[0]
            sims_org_dt = sims_acc_org_time[1] - sims_acc_org_time[0]
            obs_dt = obs_acc_time[1] - obs_acc_time[0]
            sims_dt = sims_acc_time[1] - sims_acc_time[0]

            if obs_dt == sims_dt:
                self.dt = obs_dt

            fs = 1. / self.dt
            fnyq = 0.5 * fs

            # Compute the number of pads for the time series
            # to have equal number of points for the fft
            # and for criteria 1 and 2.

            (sims_acc_ns, sims_acc_ew,
             obs_acc_ns, obs_acc_ew, ndata) = self.smcpadf(sims_acc_ns,
                                                           sims_acc_ew,
                                                           obs_acc_ns,
                                                           obs_acc_ew,
                                                           self.dt, lowcut,
                                                           8, highcut, 8,
                                                           'FALSE')
            # Start the loop for the different frequency bands
            for iband in range(len(self.B)):
                f1 = self.B[iband][0]
                f2 = self.B[iband][1]
                # Do the job only if the frequency band is within
                # the filtered band and if fnyq is higher than f1
                if f1 >= lowcut and f2 <= highcut and fnyq >= f2:
                    #print("Working on Period Band :", iband + 1,
                    #      "[", 1. / f2, 1. / f1, "]")
                    T1 = 1. / f1
                    T2 = 1. / f2
                    t_tmp = sims_perd[(sims_perd <= T1) & (T2 <= sims_perd)]

                    acc_1_flt = self.butter_bandpass(f1, f2, fnyq, sims_acc_ns, 2)
                    acc_2_flt = self.butter_bandpass(f1, f2, fnyq, sims_acc_ew, 2)
                    acc_3_flt = self.butter_bandpass(f1, f2, fnyq, obs_acc_ns, 2)
                    acc_4_flt = self.butter_bandpass(f1, f2, fnyq, obs_acc_ew, 2)

                    # Work on the frequency domain

                    # Do the response spectra
                    # Save the rsp for the specific frequency band
                    rd1 = sims_rd50_ns[(sims_perd <= T1) & (T2 <= sims_perd)]
                    rd2 = sims_rd50_ew[(sims_perd <= T1) & (T2 <= sims_perd)]
                    rd3 = obs_rd50_ns[(obs_perd <= T1) & (T2 <= obs_perd)]
                    rd4 = obs_rd50_ew[(obs_perd <= T1) & (T2 <= obs_perd)]

                    self.C8[irec, iband] = np.nanmean(
                        [self.c8_eval(rd1, rd3, t_tmp),
                         self.c8_eval(rd2, rd4, t_tmp)])

                    # Now the FFT
                    # Compute the FFT frequencies
                    F = np.fft.fftfreq(ndata, self.dt)
                    # Compute the fft and the amplitudes
                    fft_1 = np.fft.fft(sims_acc_ns)
                    fft_1 = fft_1[(0. <= F) & (f1 <= F) & (F <= f2)]
                    fft_2 = np.fft.fft(sims_acc_ew)
                    fft_2 = fft_2[(0. <= F) & (f1 <= F) & (F <= f2)]
                    fft_3 = np.fft.fft(obs_acc_ns)
                    fft_3 = fft_3[(0. <= F) & (f1 <= F) & (F <= f2)]
                    fft_4 = np.fft.fft(obs_acc_ew)
                    fft_4 = fft_4[(0. <= F) & (f1 <= F) & (F <= f2)]
                    # Slice the FFT frequencies for the working frequency band
                    F = F[(f1 <= F) & (F <= f2)]

                    fs1 = abs(fft_1) / len(fft_1)
                    fs2 = abs(fft_2) / len(fft_2)
                    fs3 = abs(fft_3) / len(fft_3)
                    fs4 = abs(fft_4) / len(fft_4)

                    self.C9[irec, iband] = np.nanmean(
                        [self.c9_eval(fs1, fs3, F),
                        self.c9_eval(fs2, fs4, F)])

                    # Work on the time domain
                    """
                    # Compute the site corrected accelerograms
                    acc_3_scor = np.fft.ifft(fft_3)
                    acc_4_scor = np.fft.ifft(fft_4)

                    # These do not need filtering because I'm working
                    # in the sliced frequency domain
                    acc_3_flt = abs(acc_3_scor)
                    acc_4_flt = abs(acc_4_scor)
                    """
                    vel_1 = self.integ(acc_1_flt, self.dt)
                    vel_2 = self.integ(acc_2_flt, self.dt)
                    vel_3 = self.integ(acc_3_flt, self.dt)
                    vel_4 = self.integ(acc_4_flt, self.dt)

                    dis_1 = self.integ(vel_1, self.dt)
                    dis_2 = self.integ(vel_2, self.dt)
                    dis_3 = self.integ(vel_3, self.dt)
                    dis_4 = self.integ(vel_4, self.dt)

                    c11, c31 = self.c13_eval(acc_1_flt,
                                             acc_3_flt)
                    c12, c32 = self.c13_eval(acc_2_flt,
                                             acc_4_flt)

                    c21, c41 = self.c24_eval(vel_1,
                                             vel_3)
                    c22, c42 = self.c24_eval(vel_2,
                                             vel_4)

                    self.C1[irec, iband] = np.nanmean(
                                           np.array(c11, c12))
                    self.C2[irec, iband] = np.nanmean(
                                           np.array(c21, c22))
                    self.C3[irec, iband] = np.nanmean(
                                           np.array(c31, c32))
                    self.C4[irec, iband] = np.nanmean(
                                           np.array(c41, c42))
                    self.C5[irec, iband] = np.nanmean(
                                           [self.c5_eval(acc_1_flt,
                                                         acc_3_flt),
                                            self.c5_eval(acc_2_flt,
                                                         acc_4_flt)])
                    self.C6[irec, iband] = np.nanmean(
                                           [self.c6_eval(vel_1, vel_3),
                                            self.c6_eval(vel_2, vel_4)])
                    self.C7[irec, iband] = np.nanmean(
                                           [self.c7_eval(dis_1, dis_3),
                                            self.c7_eval(dis_2, dis_4)])
                    self.C10[irec, iband] = np.nanmean(
                                            [self.c10_eval(acc_1_flt,
                                                           acc_3_flt),
                                             self.c10_eval(acc_2_flt,
                                                           acc_4_flt)])

                    #print(self.C1[irec, iband],
                    #      self.C2[irec, iband],
                    #      self.C3[irec, iband],
                    #      self.C4[irec, iband],
                    #      self.C5[irec, iband],
                    #      self.C6[irec, iband],
                    #      self.C7[irec, iband],
                    #      self.C8[irec, iband],
                    #      self.C9[irec, iband],
                    #      self.C10[irec, iband])

            self.S1[irec] = np.nanmean(
                np.array([np.nanmean(self.C1[irec, :]),
                          np.nanmean(self.C2[irec, :]),
                          np.nanmean(self.C3[irec, :]),
                          np.nanmean(self.C4[irec, :]),
                          np.nanmean(self.C5[irec, :]),
                          np.nanmean(self.C6[irec, :]),
                          np.nanmean(self.C7[irec, :]),
                          np.nanmean(self.C8[irec, :]),
                          np.nanmean(self.C9[irec, :]),
                          np.nanmean(self.C10[irec, :])]))

            output_file = os.path.join(a_validation_outdir,
                                       "gof-%s-%d-anderson-%s.txt" %
                                       (self.eventname, self.sim_id,
                                        station))
            out_file = open(output_file, 'w')
            line = ('#%s%5s%4s%4s%4s%4s%4s%4s%4s%4s%4s\n' %
                    ('band', 'C1', 'C2', 'C3', 'C4', 'C5',
                     'C6', 'C7', 'C8', 'C9', 'C10'))
            out_file.write(line)
            for i in range(self.BMAX):
                line = ('%s %3.1f %3.1f %3.1f %3.1f %3.1f %3.1f %3.1f %3.1f %3.1f %3.1f\n' %
                        (self.BNAMES[i], self.C1[irec, i], self.C2[irec, i],
                         self.C3[irec, i], self.C4[irec, i], self.C5[irec, i],
                         self.C6[irec, i], self.C7[irec, i], self.C8[irec, i],
                         self.C9[irec, i], self.C10[irec, i]))
                out_file.write(line)
            out_file.close()

            output_file = os.path.join(a_validation_outdir,
                                       "gof-%s-%d-anderson-%s.png" %
                                       (self.eventname, self.sim_id,
                                        station))
            self.cplots(irec, station, output_file)
            print('===> Station score :', "{:3.1f}".format(self.S1[irec]))

            irec = irec + 1

        print('==> Total number of stations processed: %d' % (irec))
        self.C1 = self.C1[0:irec, :]
        self.C2 = self.C2[0:irec, :]
        self.C3 = self.C3[0:irec, :]
        self.C4 = self.C4[0:irec, :]
        self.C5 = self.C5[0:irec, :]
        self.C6 = self.C6[0:irec, :]
        self.C7 = self.C7[0:irec, :]
        self.C8 = self.C8[0:irec, :]
        self.C9 = self.C9[0:irec, :]
        self.C10 = self.C10[0:irec, :]
        self.S1 = self.S1[0:irec]

        c1conf = [self.statts(self.C1[:, i]) for i in range(self.BMAX)]
        c2conf = [self.statts(self.C2[:, i]) for i in range(self.BMAX)]
        c3conf = [self.statts(self.C3[:, i]) for i in range(self.BMAX)]
        c4conf = [self.statts(self.C4[:, i]) for i in range(self.BMAX)]
        c5conf = [self.statts(self.C5[:, i]) for i in range(self.BMAX)]
        c6conf = [self.statts(self.C6[:, i]) for i in range(self.BMAX)]
        c7conf = [self.statts(self.C7[:, i]) for i in range(self.BMAX)]
        c8conf = [self.statts(self.C8[:, i]) for i in range(self.BMAX)]
        c9conf = [self.statts(self.C9[:, i]) for i in range(self.BMAX)]
        c10conf = [self.statts(self.C10[:, i]) for i in range(self.BMAX)]
        s1_event = np.nanmean(self.S1)
        print('==> Overall event score:', "{:3.1f}".format(s1_event))

        output_file = os.path.join(a_validation_outdir,
                                   '%d.gof_anderson.%s.txt' %
                                   (self.sim_id, self.eventname))
        out_file = open(output_file, 'w')
        line = ('#%s%5s%4s%4s%4s%4s%4s%4s%4s%4s%4s\n' %
                ('band', 'C1', 'C2', 'C3', 'C4', 'C5',
                 'C6', 'C7', 'C8', 'C9', 'C10'))
        out_file.write(line)
        for i in range(self.BMAX):
            line = ('%s %3.1f %3.1f %3.1f %3.1f %3.1f %3.1f %3.1f %3.1f %3.1f %3.1f\n' %
                    (self.BNAMES[i], c1conf[i][0], c2conf[i][0], c3conf[i][0],
                     c4conf[i][0], c5conf[i][0], c6conf[i][0], c7conf[i][0],
                     c8conf[i][0], c9conf[i][0], c10conf[i][0]))
            out_file.write(line)
        out_file.close()

        output_file = os.path.join(a_validation_outdir,
                                   "gof-%s-%d-anderson-summary.png" %
                                   (self.eventname, self.sim_id))
        self.fplots(s1_event, np.asarray(c1conf), np.asarray(c2conf),
                    np.asarray(c3conf), np.asarray(c4conf), np.asarray(c5conf),
                    np.asarray(c6conf), np.asarray(c7conf), np.asarray(c8conf),
                    np.asarray(c9conf), np.asarray(c10conf), output_file)

        print("Anderson GoF Completed".center(80, '-'))

if __name__ == "__main__":
    print("Testing Module: %s" % (os.path.basename(sys.argv[0])))
    if len(sys.argv) < 5:
        print("Usage: %s " % (os.path.basename(sys.argv[0]) +
              "station_list correction_file event_name sim_id"))
        sys.exit(1)
    ANDERSON_GOF = AndersonGOF(sys.argv[1], sys.argv[2],
                               sys.argv[3], int(sys.argv[4]))
    ANDERSON_GOF.run()
