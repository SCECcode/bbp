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
from scipy import integrate
from scipy.signal import butter, filtfilt
from scipy import interpolate
import numpy as np
import matplotlib as mpl
if mpl.get_backend() != 'agg':
    mpl.use('Agg') # Disable use of Tk/X11
import pylab

# Import BBP modules
import bband_utils
import plot_config
from install_cfg import InstallCfg
from station_list import StationList

CM2G = 980.664999
RZZ_DT = 0.02

class RZZ2015(object):
    """
    This class implements the Rezaeian-Zhong-Zareian 2015 validation metrics
    """

    def __init__(self, r_stations, eventname, sim_id=0):
        """
        Initilize class variables
        """
        self.eventname = eventname
        self.stations = r_stations
        self.log = None
        self.sim_id = sim_id

    def read_bbp(self, bbp_file):
        """
        This function reads the input bbp_file and returns 4 arrays
        containing the timestamps and 3 time series. This function
        converts all BBP files from cm/2^2 to g

        """
        times = []
        comp1 = []
        comp2 = []
        comp3 = []

        ifile = open(bbp_file, 'r')
        for line in ifile:
            line = line.strip()
            # Skip comments
            if line.startswith('%') or line.startswith('#'):
                continue
            pieces = [float(x) for x in line.split()]
            times.append(pieces[0])
            comp1.append(pieces[1] / CM2G)
            comp2.append(pieces[2] / CM2G)
            comp3.append(pieces[3] / CM2G)
        # Close input file
        ifile.close()

        # All done, return arrays
        return times, comp1, comp2, comp3

    def arias(self, F, dt, percent):
        """
        For a given motion, this function will tell you at what time a
        given percentage of arias intensity is reached (if time starts
        at 0 sec)
        """
        n = len(F)

        a_i = [pow(value, 2) for value in F]
        I = integrate.cumtrapz(a_i) * dt
        # Arias Intensity
        Ia = (F[0]**2) * dt / 2.0 + I[n-2] + (F[n-1]**2) * dt / 2.0
        It = (percent / 100.0) * Ia

        if I[0] < It:
            index = len(I) - len(I[I >= It])
            if index == len(I):
                index = index - 1
        else:
            index = 0

        t = index * dt

        return t, index

    def error_calc(self, target, theory):
        """
        This function calculates the relative error difference between two
        given vectors. One is the target vector and the other one is the
        theoretical vector which needs to be matched to the target
        """
        n = len(target)
        diff = np.abs(target - theory)
        diff_intg = 0
        for r in range(1, n-1):
            diff_intg = diff_intg + (diff[r] + diff[r+1]) / 2.0

        area = 0
        for r in range(1, n-1):
            area = area + (target[r] + target[r+1]) / 2.0

        return diff_intg / area

    def new_error_calc(self, target, theory):
        """
        This function calculates the relative error difference between two
        given vectors. One is the target vector and the other one is the
        theoretical vector which needs to be matched to the target
        """
        n = len(target)
        diff = target - theory
        diff_intg = 0
        for r in range(1, n-1):
            diff_intg = diff_intg + (diff[r] + diff[r+1]) / 2.0

        area = 0
        for r in range(1, n-1):
            area = area + (target[r] + target[r+1]) / 2.0

        return diff_intg / area

    def fre_slope_mid(self, integ_zero, x, deltat, t_mid):
        """
        """
        nu = integ_zero

        # Time when 1%AI is reached
        tesf, indexf = self.arias(x, deltat, 1)
        # Time when 99%AI is reached
        tesl, indexl = self.arias(x, deltat, 99)
        nuf = nu[indexf]
        nul = nu[indexl]
        t = np.linspace(0, len(x) * deltat, len(x))

        tes2 = tesf + (tesl-tesf) / 8.0
        index2 = len(t[tes2 > t]) - 1
        nu2 = nu[index2]

        tes3 = tesf + 2 * (tesl-tesf) / 8.0
        index3 = len(t[tes3 > t]) - 1
        nu3 = nu[index3]

        tes4 = tesf + 3 * (tesl-tesf) / 8.0
        index4 = len(t[tes4 > t]) - 1
        nu4 = nu[index4]

        tes5 = tesf + 4 * (tesl-tesf) / 8.0
        index5 = len(t[tes5 > t]) - 1
        nu5 = nu[index5]

        tes6 = tesf + 5 * (tesl-tesf) / 8.0
        index6 = len(t[tes6 > t]) - 1
        nu6 = nu[index6]

        tes7 = tesf + 6 * (tesl-tesf) / 8.0
        index7 = len(t[tes7 > t]) - 1
        nu7 = nu[index7]

        tes8 = tesf + 7 * (tesl-tesf) / 8.0
        index8 = len(t[tes8 > t]) - 1
        nu8 = nu[index8]

        p = np.polyfit([tesf, tes2, tes3, tes4, tes5, tes6, tes7, tes8, tesl],
                       [nuf, nu2, nu3, nu4, nu5, nu6, nu7, nu8, nul], 2)

        # Desired data for this record
        fre_mid = 2.0 * p[0] * t_mid + p[1]
        slope_mid = 2.0 * p[0]

        return fre_mid, slope_mid

    def filter_data(self, cur_times, cur_data, deltat):
        """
        This function brings the timeseries to the desired DT
        """

        # If we already have the DT, do nothing, keep what we have!
        if deltat == RZZ_DT:
            return cur_times, cur_data, deltat

        #######################################################################
        # # Now downsample as needed                                          #
        # if deltat > obs_delta_t:                                            #
        #     # Filter obs_data                                               #
        #     X5, Y5 = butter(8, (1.0 / deltat) / (1.0 / obs_delta_t), 'low') #
        #     tfif = filtfilt(X5, Y5, obs_data,                               #
        #                     padlen=3*(max(len(X5), len(Y5))-1))             #
        #     idx = [int(round(x)) for x in np.arange(1,                      #
        #                                             len(tfif) + 0.0000001,  #
        #                                             (deltat/obs_delta_t))]  #
        #     tf = [tfif[i-1] for i in idx]                                   #
        #     # Keep sym_data intact                                          #
        #     xi = sym_data                                                   #
        # elif deltat < obs_delta_t:                                          #
        #     # Filter sym_data                                               #
        #     X5, Y5 = butter(8, (1.0 / obs_delta_t) / (1.0 / deltat), 'low') #
        #     xif = filtfilt(X5, Y5, sym_data,                                #
        #                    padlen=3*(max(len(X5), len(Y5))-1))              #
        #     idx = [int(round(x)) for x in np.arange(1,                      #
        #                                             len(xif) + 0.0000001,   #
        #                                             (obs_delta_t/deltat))]  #
        #     xi = [xif[i-1] for i in idx]                                    #
        #     # Keep obs_data intact                                          #
        #     tf = obs_data                                                   #
        # else:                                                               #
        #     tf = obs_data                                                   #
        #     xi = sym_data                                                   #
        #######################################################################

        # Check if we can decimate (preferable)
        dt_ratio = 1.0 * RZZ_DT / deltat
        if (dt_ratio > 1.0) and (dt_ratio == int(dt_ratio)):
            # Yes!
            X5, Y5 = butter(8, (1.0 / RZZ_DT) / (1.0 / deltat), 'low')
            filt_data = filtfilt(X5, Y5, cur_data,
                                 padlen=3*(max(len(X5), len(Y5))-1))
            dt_ratio = int(dt_ratio)
            new_data = filt_data[::dt_ratio]
            #new_data = cur_data[::dt_ratio]
            new_times = cur_times[::dt_ratio]
            return new_times, new_data, RZZ_DT

        # Need to interpolate
        f = interpolate.interp1d(cur_times, cur_data)

        # Figure out new times
        new_times = []
        last_time = cur_times[-1]
        cur_time = 0.0
        while True:
            new_times.append(round(cur_time, 8))
            cur_time = cur_time + RZZ_DT
            if round(cur_time, 8) > last_time:
                break

        # Now, interpolate data
        new_data = []
        for new_time in new_times:
            new_data.append(float(f(new_time)))

        return new_times, new_data, RZZ_DT

    def sim_f(self, w0, wn, zeta, n, deltat):
        """
        Generates the unmodulated simulations that are "not" post-processed yet
        w0: frequency at the begining (rad/sec = 2pi*Hz)
        wn: frequency at the end
        zeta: constant damping ratio
        n: length of acceleration vector
        deltat: time steps
        """
        t = np.linspace(0, deltat * (n-1), n)
        f = np.zeros(n)
        C = np.zeros(n)

        # Find the scale parameter, C vector (see eqn 3 in the qual. report,
        # C is the denominator). And construct the Sj vectors to evaluate f
        # vector (see equations 2 and 3 in the report)
        for j in range(0, n):
            tj = t[j]
            uj = np.random.normal(0)
            omega = w0 - (w0 - wn) * tj / t[n-1]
            omegaD = omega * math.sqrt(1 - zeta**2)
            Cj = np.zeros(n)
            Sj = np.zeros(n)
            for i in range(0, n):
                ti = t[i]
                if ti >= tj:
                    hf = (omega / math.sqrt(1 - zeta**2) *
                          math.exp(-zeta * omega * (ti-tj)) *
                          math.sin(omegaD * (ti-tj)))
                    Cj[i] = hf ** 2
                    Sj[i] = hf
                else:
                    Cj[i] = 0.0
                    Sj[i] = 0.0
            C = C + Cj
            f = f + Sj * uj
        C = C ** 0.5
        C[0] = 0.1

        # Scale f to normalize its standard deviation
        f = f / C

        return f, t

    def nmax_pmin(self, F, n):
        """
        Gives the cumulative count of negative-maximums and
        positive-minimums of F, a vector of size n, as a function of time
        """
        I = np.zeros(n-2)
        count_min = 0
        count_max = 0

        for r in range(1, n-1):
            if F[r] > 0:
                if F[r-1] > F[r]:
                    if F[r+1] > F[r]:
                        count_min = count_min + 1
            else:
                if F[r-1] < F[r]:
                    if F[r+1] < F[r]:
                        count_max = count_max + 1
            I[r-1] = count_min + count_max

        I = np.insert(I, 0, I[0], axis=0)
        I = np.insert(I, len(I), I[-1], axis=0)
        return I

    def zhatcalc(self, error, target_lp):
        """
        # Find which two z's give the smallest positive
        # and the smallest negative errors
        # Interpolate and find zhat that gives zetadiff=0
        # { z-zp=[(zp-zn)/(ep-en)](e-ep)  zhat is the value of z when e=0}
        # If z < 0.1, interpolate by assuming zetadiff=sum(target) at z=0
        # implying that when damping is 0, motion is harmonic
        """
        if error[0] > 0 and error[8] < 0:
            # interpolate
            zp = len(error[error>0]) * 0.1
            zn = zp + 0.1
            ep = error[int(np.round(zp * 10)) - 1]
            en = error[int(np.round(zn * 10)) - 1]
        elif error[8] > 0:
            zp = 0.9
            ep = np.sum(target_lp)
            zn = 1.0
            en = error[0]
        else:
            zp = 0.0
            ep = np.sum(target_lp)
            zn = 0.1
            en = error[0]

        zhat = zp - ep * ((zp - zn) / (ep - en))

        return zhat

    def damping_z(self, tf, deltat, wmid, wslope):
        """
        Calculate the bandwidth parameter zeta
        """

        # Cut the target at 5% and 9% arias intensity
        # Work with f_cut and t_cut as the target after
        # this point
        t5, index5 = self.arias(tf, deltat, 5)
        t45 = self.arias(tf, deltat, 45)[0]
        t95, index95 = self.arias(tf, deltat, 95)

        tf_cut = tf[index5:index95]
        t_cut = np.linspace(t5, t95, len(tf_cut))

        # Estimate the frequency at the begining and the end of the new
        # cut-motion having the slope and a point (t45, midfreq), get the
        # equation of the line and calculate the freq at t5 and t95

        # freq (Hz) at the begining of the motion to be simulated
        w0 = wmid + wslope * (t5 - t45)
        # freq (Hz) at the end of the motion to be simulated
        wn = wmid + wslope * (t95 - t45)

        sim_lp = {}
        for k in range(0, 9):
            sim_lp[k] = None
            for j in range(0, 20):
                # Unmodulated simulation for the cut (5-95) motion
                sim_f = self.sim_f(w0*2*np.pi, wn*2*np.pi, (k+1)*0.1,
                                  len(tf_cut), deltat)[0]
                # Cumulative nmax+pmin: (Local Peaks)
                new_row = self.nmax_pmin(sim_f, len(sim_f))
                if sim_lp[k] is None:
                    sim_lp[k] = new_row
                else:
                    sim_lp[k] = np.vstack((sim_lp[k], new_row))

        # Take average of every 20 set of simulations
        avg_lp = {}
        for k in range(0, 9):
            avg_lp[k] = sum(sim_lp[k]) / 20.0

        # Calculate errors (to see which column of avg_LP is
        # closer to the target nmax+pmin)
        zetadiff = np.zeros(9)
        tlp = self.nmax_pmin(tf_cut, len(tf_cut))
        for k in range(0, 9):
            zetadiff[k] = np.sum(tlp - avg_lp[k])

        # Find which two z's give the smallest positive
        # and the smallest negative errors
        # Interpolate and find zhat that gives zetadiff=0
        # { z-zp=[(zp-zn)/(ep-en)](e-ep)  zhat is the value of z when e=0}
        # If z < 0.1, interpolate by assuming zetadiff=sum(target) at z=0
        # implying that when damping is 0, motion is harmonic
        return self.zhatcalc(zetadiff, tlp)

    def process(self, stat, comp, output_dir,
                obs_times, obs_data, sym_times, sym_data):
        """
        This function computes the validation metrics of RZZ2015 for a
        given pair of timeseries
        """
        # Create path for our plot and output file
        out_file = os.path.join(output_dir,
                                '%d.rzz2015.%s.txt' %
                                (self.sim_id, self.eventname))
        out_plot = os.path.join(output_dir,
                                '%d.rzz2015.%s.%03d.png' %
                                (self.sim_id, stat, comp))

        # Create fig
        fig, axs = pylab.plt.subplots(2, 3)
        fig.set_size_inches(17, 8.5)
        fig.suptitle("RZZ2015 - %s - %s - %03d" %
                     (self.eventname, stat, comp), size=16)
        fig.subplots_adjust(hspace=0.4)
        fig.subplots_adjust(left=0.05)
        fig.subplots_adjust(right=0.98)

        # Plot original time series
        subfig = axs[0][0]
        subfig.set_title("Non-synched Acc time series", size=14)

        # Try to get a reasonable scale
        peak_val = max(abs(max(obs_data)),
                       abs(min(obs_data)),
                       abs(max(sym_data)),
                       abs(min(sym_data)))
        plot_offset = 1
        if peak_val <= 0.2:
            plot_offset = 0.2

        subfig.plot(obs_times, np.array(obs_data)+plot_offset,
                    color='black', label='Recorded')
        subfig.plot(sym_times, np.array(sym_data)-plot_offset,
                    color='red', label='Simulated')
        subfig.set_ylabel("Acceleration (g)", size=10)
        subfig.set_xlabel("Time (s)", size=10)
        subfig.grid(True)
        subfig.minorticks_on()

        # Calculate dt for each timeseries
        obs_delta_t = obs_times[1] - obs_times[0]
        deltat = sym_times[1] - sym_times[0]

        print("==> Processing station: "
              "%s - Comp: %s - obs dt = %4.3f - sym dt = %4.3f" %
              (stat, comp, obs_delta_t, deltat))

        # Bring obs_data to desired RZZ_DT
        print("Obs in: %d pts, %f delta t" % (len(obs_data), obs_delta_t))
        obs_times, obs_data, obs_delta_t = self.filter_data(obs_times,
                                                            obs_data,
                                                            obs_delta_t)
        print("Obs out: %d pts, %f delta t" % (len(obs_data), obs_delta_t))
        # Bring sym_data to desired RZZ_DT
        print("Sym in: %d pts, %f delta t" % (len(sym_data), deltat))
        sym_times, sym_data, deltat = self.filter_data(sym_times,
                                                       sym_data,
                                                       deltat)
        print("Sym out: %d pts, %f delta t" % (len(sym_data), deltat))

        tf = obs_data
        xi = sym_data

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
        oq2_integ = oq2_integ * obs_delta_t
        cav_oq2_integ_total = oq2_integ[-1]

        for cav_tp in range(0, len(tf)):
            if tq2_integ[cav_tp] >= (0.05 * cav_tq2_integ_total):
                cav_tii = cav_tp
                break

        for cav_p in range(0, len(xi)):
            if oq2_integ[cav_p] >= (0.05 * cav_oq2_integ_total):
                cav_ii = cav_p
                break

        if cav_ii >= cav_tii:
            ii = cav_ii - cav_tii

            if len(xi) < ii + len(tf):
                new_max = len(xi) - ii
                tf = tf[0:new_max]
                x = xi[ii:]
            else:
                x = xi[ii:(ii + len(tf))]

        else:
            ii = cav_tii - cav_ii

            if len(tf) < ii + len(xi):
                new_max = len(tf) - ii
                x = xi[0:new_max]
                tf = tf[ii:]
            else:
                tf = tf[ii:(ii + len(xi))]
                x = xi

        ttt = np.linspace(0, len(tf) * deltat, len(tf))

        # Plot synched time series
        subfig = axs[0][1]
        subfig.set_title("Synched Acc time series", size=14)

        # Try to get a reasonable scale
        peak_val = max(abs(max(obs_data)),
                       abs(min(obs_data)),
                       abs(max(sym_data)),
                       abs(min(sym_data)))
        plot_offset = 1
        if peak_val <= 0.2:
            plot_offset = 0.2

        subfig.plot(ttt, np.array(tf)+plot_offset,
                    color='black', label='Recorded')
        subfig.plot(ttt, np.array(x)-plot_offset,
                    color='red', label='Simulated')
        subfig.set_ylabel("Acceleration (g)", size=10)
        subfig.set_xlabel("Time (s)", size=10)
        subfig.grid(True)
        subfig.minorticks_on()

        # Frequency plot: Plot Fourier Amplitude Spectra

        # Recorded motion
        aa1 = tf
        point = len(aa1)
        dt = deltat
        sr = 1.0 / dt
        t = np.linspace(0, dt * (point - 1), point)

        nfft1 = point
        pxx = np.abs(np.fft.fft(aa1, nfft1))
        f1 = sr * np.array(range(0, nfft1//2)) / nfft1

        # Simulated motion
        aa2 = x
        point = len(aa2)
        dt = deltat
        sr = 1.0 / dt
        t = np.linspace(0, dt * (point - 1), point)

        nfft2 = point
        pxx2 = np.abs(np.fft.fft(aa2, nfft2))
        f2 = sr * np.array(range(0, nfft2//2)) / nfft2

        # Frequency plot
        subfig = axs[0][2]
        subfig.set_xscale("log")
        subfig.set_yscale("log")
        subfig.set_title("Fourier Amplitude Spectra", size=14)
        subfig.plot(f1, pxx[0:nfft1//2],
                    color='black', label='Recorded')
        subfig.plot(f2, pxx2[0:nfft2//2],
                    color='red', lw=0.25, label='Simulated')
        subfig.set_ylabel("FAS (g.s)", size=10)
        subfig.set_xlabel("Frequency (Hz)", size=10)
        subfig.grid(True)
        subfig.minorticks_on()

        # Calculate Validation Metric1 (evolving intensity of the
        # motion), 2 Scalar Error Terms for Metric1, and the First 3
        # Parameters

        # Taking the integral of recorded accel^2
        tq2_integ = np.zeros((1, len(tf)))[0]
        tq2_integ[0] = 0

        for r in range(1, len(tf)):
            tq2_integ[r] = tq2_integ[r-1] + (tf[r]**2)
        tq2_integ = tq2_integ * deltat

        # Taking the integral of simulated accel^2
        q2_integ = np.zeros((1, len(tf)))[0]
        for r in range(1, len(tf)):
            q2_integ[r] = q2_integ[r-1] + (x[r]**2)
        q2_integ = q2_integ * deltat
        q2_integ = q2_integ[0:len(tq2_integ)]

        # Calculate the scalar error terms
        metric_error_1 = self.error_calc(tq2_integ, q2_integ)
        metric_error_4 = self.new_error_calc(tq2_integ, q2_integ)
        metric_error_shape_1 = metric_error_4 / metric_error_1

        # Calculate the first three parameters
        result_1 = tq2_integ[-1] * np.pi / 2.0
        result_2 = self.arias(tf, deltat, 95)[0] - self.arias(tf, deltat, 5)[0]
        result_3 = self.arias(tf, deltat, 45)[0]
        result_7 = q2_integ[len(x)-1] * np.pi / 2.0
        result_8 = self.arias(x, deltat, 95)[0] - self.arias(x, deltat, 5)[0]
        result_9 = self.arias(x, deltat, 45)[0]

        # Plot validation metric 1
        subfig = axs[1][0]
        subfig.set_title("Evolution of Intensity", size=14)
        subfig.plot(ttt, tq2_integ,
                    color='black', label='Recorded')
        subfig.plot(ttt, q2_integ,
                    color='red', label='Simulated')
        subfig.set_ylabel("Cumulative Intensity (g^2.s)", size=10)
        subfig.set_xlabel("Time (s)", size=10)
        subfig.grid(True)
        texth1 = max(max(tq2_integ), max(q2_integ)) / 4.0
        texth4 = max(max(tq2_integ), max(q2_integ)) / 8.0
        textw1 = (25 / 40.0) * max(ttt)
        subfig.text(textw1, texth1,
                    r'$\epsilon = %1.2f$' % (metric_error_1),
                    fontsize=14)
        subfig.text(textw1, texth4,
                    r'$\nu = %1.2f$' % (metric_error_shape_1),
                    fontsize=14)
        subfig.minorticks_on()

        # Calculate validation metric 2 (evolving frequency of the
        # motion), 2 scalar error terms for metric2, and the 2
        # frequency parameters

        # Recorded motion
        t_integ_zero = np.zeros((1, len(tf)-1))[0]
        nstart = 0
        xx = 0
        for r in range(nstart, len(tf)-1):
            count = 0
            for i in range(nstart, r + 1):
                if tf[i] < xx:
                    if tf[i+1] >= xx:
                        count = count + 1
            t_integ_zero[r] = count

        tt2 = np.linspace(0, (len(t_integ_zero)-1) * deltat, len(t_integ_zero))

        result_4, result_5 = self.fre_slope_mid(t_integ_zero, tf, deltat,
                                                result_3)

        # Simulated motion
        integ_zero = np.zeros((1, len(x)-1))[0]
        nstart = 0
        xx = 0
        for r in range(nstart, len(x)-1):
            count = 0
            for i in range(nstart, r + 1):
                if x[i] < xx:
                    if x[i+1] >= xx:
                        count = count + 1
            integ_zero[r] = count

        t2 = np.linspace(0, (len(integ_zero)-1) * deltat, len(integ_zero))
        integ_zero = integ_zero[0:len(t_integ_zero)]

        # Calculate the scalar error terms
        metric_error_2 = self.error_calc(t_integ_zero, integ_zero)
        metric_error_5 = self.new_error_calc(t_integ_zero, integ_zero)
        metric_error_shape_2 = metric_error_5 / metric_error_2

        result_10, result_11 = self.fre_slope_mid(integ_zero, x, deltat,
                                                  result_9)

        # Plot validation metric 2
        subfig = axs[1][1]
        subfig.set_title("Evolution of Frequency", size=14)
        subfig.plot(tt2, t_integ_zero,
                    color='black', label='Recorded')
        subfig.plot(tt2, integ_zero,
                    color='red', label='Simulated')
        subfig.set_ylabel("Cumulative zero-level up-crossings", size=10)
        subfig.set_xlabel("Time (s)", size=10)
        subfig.grid(True)
        texth2 = max(max(t_integ_zero), max(integ_zero)) / 4.0
        texth5 = max(max(t_integ_zero), max(integ_zero)) / 8.0
        textw2 = (25 / 40.0) * max(tt2)
        subfig.text(textw2, texth2,
                    r'$\epsilon = %1.2f$' % (metric_error_2),
                    fontsize=14)
        subfig.text(textw2, texth5,
                    r'$\nu = %1.2f$' % (metric_error_shape_2),
                    fontsize=14)
        subfig.minorticks_on()

        # Calculate validation metric 3 (evolving bandwidth of the
        # motion), and 2 scalar error terms for metric 3

        # Recorded motion
        tcountmin = 0
        tcountmax = 0
        ttot = np.zeros((1, len(tf)-2))[0]
        for r in range(1, len(tf) - 1):
            if tf[r] > 0:
                if tf[r-1] > tf[r]:
                    if tf[r+1] > tf[r]:
                        tcountmin = tcountmin + 1
            else:
                if tf[r-1] < tf[r]:
                    if tf[r+1] < tf[r]:
                        tcountmax = tcountmax + 1
            ttot[r-1] = tcountmin + tcountmax

        tt3 = np.linspace(0, (len(ttot)-1) * deltat, len(ttot))

        # Simulated motion
        tcountmin = 0
        tcountmax = 0
        tot = np.zeros((1, len(x)-2))[0]
        for r in range(1, len(x) - 1):
            if x[r] > 0:
                if x[r-1] > x[r]:
                    if x[r+1] > x[r]:
                        tcountmin = tcountmin + 1
            else:
                if x[r-1] < x[r]:
                    if x[r+1] < x[r]:
                        tcountmax = tcountmax + 1
            tot[r-1] = tcountmin + tcountmax

        tot = tot[0:len(ttot)]

        # Calculate the scalar error terms
        metric_error_3 = self.error_calc(ttot, tot)
        metric_error_6 = self.new_error_calc(ttot, tot)
        metric_error_shape_3 = metric_error_6 / metric_error_3

        # Calculate the bandwidth parameter zeta
        result_6 = 0.0
        result_12 = 0.0
        # Disabled calculations below as processing can take a long time
        # Will be re-enabled in the future once the algorithm is revised
        # result_6 = self.damping_z(tf, deltat, result_4, result_5)
        # result_12 = self.damping_z(x, deltat, result_10, result_11)

        # Plot validation metric 3
        subfig = axs[1][2]
        subfig.set_title("Evolution of Bandwidth", size=14)
        subfig.plot(tt3, ttot,
                    color='black', label='Recorded')
        subfig.plot(tt3, tot,
                    color='red', label='Simulated')
        subfig.set_ylabel("Cumulative positive min & negative max", size=10)
        subfig.set_xlabel("Time (s)", size=10)
        subfig.grid(True)
        texth3 = max(max(ttot), max(tot)) / 4.0
        texth6 = max(max(ttot), max(tot)) / 8.0
        textw3 = (25 / 40.0) * max(tt2)
        subfig.text(textw3, texth3,
                    r'$\epsilon = %1.2f$' % (metric_error_3),
                    fontsize=14)
        subfig.text(textw3, texth6,
                    r'$\nu = %1.2f$' % (metric_error_shape_3),
                    fontsize=14)
        subfig.minorticks_on()

        # All done! Save plot!
        fig.savefig(out_plot, format='png', transparent=False,
                    dpi=plot_config.dpi)

        # Close image, to save memory
        pylab.plt.close(fig)

        # Organize results
        metric_error_avg = [metric_error_1, metric_error_2, metric_error_3]
        epsilon_a = metric_error_avg[0]
        epsilon_b = metric_error_avg[1]
        epsilon_c = metric_error_avg[2]
        nu_a = metric_error_shape_1
        nu_b = metric_error_shape_2
        nu_c = metric_error_shape_3

        # Arias Intensity
        r1_record = result_1
        r1_siml = result_7
        # Duration
        r2_record = result_2
        r2_siml = result_8
        # Arias Intensity / Duration
        r3_record = r1_record / r2_record
        r3_siml = r1_siml / r2_siml
        # wmid (frequency)
        r4_record = result_4
        r4_siml = result_10
        # w' (rate of change of frequency)
        r5_record = result_5
        r5_siml = result_11
        # Zeta parameters
        r6_record = result_6
        r6_siml = result_12

        # Write data to file
        ofile = open(out_file, 'a')
        ofile.write("%s, %03d, " % (stat, comp) +
                    "%7.4f, %7.4f, " % (epsilon_a, nu_a) +
                    "%7.4f, %7.4f, " % (epsilon_b, nu_b) +
                    "%7.4f, %7.4f, " % (epsilon_c, nu_c) +
                    "%7.4f, %7.4f, " % (r1_record, r1_siml) +
                    "%7.4f, %7.4f, " % (r2_record, r2_siml) +
                    "%7.4f, %7.4f, " % (r3_record, r3_siml) +
                    "%7.4f, %7.4f, " % (r4_record, r4_siml) +
                    "%7.4f, %7.4f, " % (r5_record, r5_siml) +
                    "%7.2f, %7.2f\n" % (r6_record, r6_siml))
        ofile.close()

    def run(self):
        """
        Runs the Anderson GoF code
        """
        print("RZZ2015".center(80, '-'))

        # Load configuration, set sim_id
        install = InstallCfg.getInstance()
        sim_id = self.sim_id

        # Build directory paths
        a_tmpdir = os.path.join(install.A_TMP_DATA_DIR, str(sim_id))
        a_indir = os.path.join(install.A_IN_DATA_DIR, str(sim_id))
        a_outdir = os.path.join(install.A_OUT_DATA_DIR, str(sim_id))
        a_logdir = os.path.join(install.A_OUT_LOG_DIR, str(sim_id))
        a_validation_outdir = os.path.join(a_outdir, "validations", "rzz2015")

        # Make sure the output and tmp directories exist
        bband_utils.mkdirs([a_tmpdir, a_indir, a_outdir, a_validation_outdir],
                           print_cmd=False)

        # Now the file paths
        self.log = os.path.join(a_logdir, "%d.rzz2015.log" % (sim_id))
        sta_file = os.path.join(a_indir, self.stations)
        sta_base = os.path.basename(os.path.splitext(self.stations)[0])
        obs_dir = os.path.join(a_tmpdir, "obs_seis_%s" % (sta_base))

        # Get station list
        slo = StationList(sta_file)
        site_list = slo.getStationList()

        # Create output file, add header
        out_file = open(os.path.join(a_validation_outdir,
                                     '%d.rzz2015.%s.txt' %
                                     (self.sim_id, self.eventname)), 'w')
        out_file.write("#station, component, epsilon_a, nu_a,"
                       " epsilon_b, nu_b, epsilon_c, nu_c,"
                       " r1_record, r1_siml, r2_record, r2_siml,"
                       " r3_record, r3_siml, r4_record, r4_siml,"
                       " r5_record, r5_siml, r6_record, r6_siml\n")
        out_file.close()

        # Go through each station
        for site in site_list:
            stat = site.scode
            r_obs_bbp = "%s.bbp" % (stat)
            a_obs_bbp = os.path.join(obs_dir, r_obs_bbp)
            r_sym_bbp = "%d.%s.acc.bbp" % (sim_id, stat)
            a_sym_bbp = os.path.join(a_outdir, r_sym_bbp)

            if not (os.path.exists(a_obs_bbp) and
                    os.path.exists(a_sym_bbp)):
                # Just skip it
                print("===> Couldn't find files "
                      "%s and %s, skipping station %s" %
                      (a_obs_bbp, a_sym_bbp, stat))
                continue

            obs_data = self.read_bbp(a_obs_bbp)
            sym_data = self.read_bbp(a_sym_bbp)

            # Process each component separately
            for comp in range(1, 3):
                self.process(stat, comp, a_validation_outdir,
                             obs_data[0], obs_data[comp],
                             sym_data[0], sym_data[comp])

        print("RZZ2015 Completed".center(80, '-'))

if __name__ == "__main__":
    print("Testing Module: %s" % (os.path.basename(sys.argv[0])))
    if len(sys.argv) < 4:
        print("Usage: %s " % (os.path.basename(sys.argv[0])) +
              "station_list event_name sim_id")
        sys.exit(1)
    RZZ2015 = RZZ2015(sys.argv[1], sys.argv[2], int(sys.argv[3]))
    RZZ2015.run()
