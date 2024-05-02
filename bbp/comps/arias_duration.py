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

Created on Aug 7, 2012
@author: maechlin
Purpose: Translation of arias duration algorithm in matlab from
         Christine Goulet for use on bbp project
"""
from __future__ import division, print_function

# G2CMSS = 980.665 # Convert g to cm/s/s
import os
import sys
import math
import scipy.integrate

# Import Broadband modules
import bband_utils
import install_cfg
import bbp_formatter
from station_list import StationList

# Converting to cm units. Use approximation to g
G_TO_CMS = 981.0 # %(cm/s)

def ad_from_acc(a_in_peer_file, a_out_ad):
    """
    This function reads a PEER formatted acceleration file and outputs
    the arias duration
    """

    acc = []
    start_samples = False

    pfx = open(a_in_peer_file, "r")
    for line in pfx:
        pieces = line.split()
        if len(pieces) and pieces[0].lower() == 'acceleration':
            start_samples = True
            continue
        if start_samples and line.find("NPTS") > 0:
            pts = int(pieces[0])
            dt = float(pieces[1])
            #print("Reading Seismogram with NPTS: %d and DT: %f" % (pts, dt))
            continue
        if start_samples:
            for val in pieces:
                acc.append(float(val))
    pfx.close()

    if not start_samples:
        print("No samples found in peer file: %s" % (a_in_peer_file))
        return -1

    #
    # Full file read. Now count samples found
    # in_vals contains samples, dt timesetp, pts number of data points
    #
    if len(acc) != pts:
        print("Number of points in header %d does not match number " % (pts) +
              "of pts read from file %d. Exiting file: %s" %
              (len(acc), a_in_peer_file))
        return -1
    #
    # Reduce the sample count to be an even number for fft
    #
    if len(acc) % 2 != 0:
        pts = pts - 1
        #print("adjusted pts: %d" % (pts))

    #
    # Create a list of dt values
    #
    dt_vals = [point * dt for point in range(pts)]

    #
    # Integrate and calculate vel and disp arrays
    #
    acc_cms = [value * G_TO_CMS for value in acc]
    velo = list(dt * scipy.integrate.cumtrapz(acc_cms))
    velo.insert(0, 0)
    disp = list(dt * scipy.integrate.cumtrapz(velo))
    disp.insert(0, 0)

    #
    # Find PGA and PGV
    #
    temp1 = max(acc)
    temp2 = min(acc)
    if temp1 > abs(temp2):
        max_acc = temp1
        max_acc_idx = acc.index(temp1)
    else:
        max_acc = abs(temp2)
        max_acc_idx = acc.index(temp2)
    #print("Peak Accel (g): %f (secs): %f" % (max_acc, (max_acc_idx * dt)))

    temp1 = max(velo)
    temp2 = min(velo)
    if temp1 > abs(temp2):
        max_vel = temp1
        max_vel_idx = velo.index(temp1)
    else:
        max_vel = abs(temp2)
        max_vel_idx = velo.index(temp2)
    #print("Peak Vel (cm/s): %f (secs): %f" % (max_vel, (max_vel_idx * dt)))

    temp1 = max(disp)
    temp2 = min(disp)
    if temp1 > abs(temp2):
        max_disp = temp1
        max_disp_idx = disp.index(temp1)
    else:
        max_disp = abs(temp2)
        max_disp_idx = disp.index(temp2)
    #print("Peak Disp (cm): %f (secs): %f" % (max_disp, (max_disp_idx * dt)))

    # Arias Intensities
    # Using the trapezoidal integration
    arias_intensity = [pow((value * G_TO_CMS), 2) for value in acc]
    tsum = list(scipy.integrate.cumtrapz(arias_intensity) * dt)
    tsum.insert(0, 0)

    arias_intensity = [num * math.pi / (2 * G_TO_CMS) for num in tsum]

    #
    # Find peak arias_intensity
    #
    temp1 = max(arias_intensity)
    temp2 = min(arias_intensity)
    if temp1 > abs(temp2):
        arias_intensity_max = temp1
        arias_intensity_index = arias_intensity.index(temp1)
    else:
        arias_intensity_max = abs(temp2)
        arias_intensity_index = arias_intensity.index(temp2)
    #print("Peak Arias Intensity (cm/sec): %f (secs): %f" %
    #      (arias_intensity_max, (arias_intensity_index * dt)))

    if arias_intensity_max:
        ia_norm = [100 * (i_acc / arias_intensity_max) for i_acc in arias_intensity]
    else:
        ia_norm = [0.0 for i_acc in arias_intensity]

    # Define the time for AI=5%, 20%, 75%, 80%, 95%
    time_ai5 = 0
    for i in range(pts):
        if ia_norm[i] >= 5:
            break
        time_ai5 = dt * i

    time_ai20 = 0
    for i in range(pts):
        if ia_norm[i] >= 20:
            break
        time_ai20 = dt * i

    time_ai75 = 0
    for i in range(pts):
        if ia_norm[i] >= 75:
            break
        time_ai75 = dt * i

    time_ai80 = 0
    for i in range(pts):
        if ia_norm[i] >= 80:
            break
        time_ai80 = dt * i

    time_ai95 = 0
    for i in range(pts):
        if ia_norm[i] >= 95:
            break
        time_ai95 = dt * i

    # Now, calculate the arias intervals 5% to 75% and 5% to 95%
    time_5_75 = time_ai75 - time_ai5
    time_5_95 = time_ai95 - time_ai5
    time_20_80 = time_ai80 - time_ai20

    #print("Arias Intervals: 5 to 75 % 6f (secs), 5 to 95 % 6f (secs) " %
    #      (time_5_75, time_5_95))

    #
    # create single component avd (accel, vel, disp) bbp file
    #
    #outfile = open("%s.avd.bbp" % (a_in_peer_file), "w")
    #outfile.write("# Acc  Vel Disp (avd) from: %s\n" % a_in_peer_file)
    #outfile.write("# Peak Accel: %f at: %f (s)\n" %
    #              (max_acc, (max_acc_idx * dt)))
    #outfile.write("# Peak Vel: %f at: %f (s)\n" %
    #              (max_vel, (max_vel_idx * dt)))
    #outfile.write("# Peak Disp: %f at: %f (s)\n" %
    #              (max_disp, (max_disp_idx * dt)))
    #outfile.write("# Seconds Accel (g) Vel (cm/s) Disp (cm)\n")
    #outfile.write("# %d %f  NPTS, DT\n" % (pts, dt))
    #for i in range(pts):
    #    outfile.write("% 8f % 8f % 8f % 8f\n" %
    #                  (dt_vals[i], acc[i], velo[i], disp[i]))
    #outfile.close()

    #
    # Create single component Arias Duration bbp file
    #
    outfile = open(a_out_ad, "w")
    outfile.write("# Arias Intensities from input accel: %s\n" % a_in_peer_file)
    outfile.write("# Peak Arias Intensity (cm/sec): %f (secs): %f\n" %
                  (arias_intensity_max, (arias_intensity_index * dt)))
    outfile.write("# Arias Intervals: T5-75 %f (s), T5-95 %f (s), T20-80 %f (s)\n" %
                  (time_5_75, time_5_95, time_20_80))
    outfile.write("# Seconds Accel (g) "
                  "Arias Intensity (cm/s) "
                  "ADNormalized (%)\n")
    outfile.write("# %d %f  NPTS, DT\n" % (pts, dt))
    for i in range(pts):
        outfile.write("% 8f % 8f % 8f % 8f\n" %
                      (dt_vals[i], acc[i], arias_intensity[i], ia_norm[i]))
    outfile.close()
    return arias_intensity_max, time_5_75, time_5_95, time_20_80;

class AriasDuration(object):
    """
    BBP module implementation of arias duration
    """

    def __init__(self, i_r_stations, sim_id=0):
        """
        Initializes class variables
        """
        self.sim_id = sim_id
        self.r_stations = i_r_stations

    def run(self):
        print("AriasDuration".center(80, '-'))
        #
        # convert input bbp acc files to peer format acc files
        #

        install = install_cfg.InstallCfg.getInstance()
        sim_id = self.sim_id
        sta_base = os.path.basename(os.path.splitext(self.r_stations)[0])
        self.log = os.path.join(install.A_OUT_LOG_DIR, str(sim_id),
                                "%d.araisduration_%s.log" % (sim_id, sta_base))
        a_statfile = os.path.join(install.A_IN_DATA_DIR,
                                  str(sim_id),
                                  self.r_stations)
        a_tmpdir = os.path.join(install.A_TMP_DATA_DIR, str(sim_id))
        a_outdir = os.path.join(install.A_OUT_DATA_DIR, str(sim_id))

        #
        # Make sure the tmp and out directories exist
        #
        bband_utils.mkdirs([a_tmpdir, a_outdir], print_cmd=False)

        slo = StationList(a_statfile)
        site_list = slo.getStationList()

        for site in site_list:
            stat = site.scode
            print("==> Processing station: %s" % (stat))
            bbpfile = os.path.join(a_outdir,
                                   "%d.%s.acc.bbp" % (sim_id, stat))

            # Now we need to convert to peer format
            out_n_acc = os.path.join(a_tmpdir,
                                     "%d.%s.peer_n.acc" % (sim_id, stat))
            out_e_acc = os.path.join(a_tmpdir,
                                     "%d.%s.peer_e.acc" % (sim_id, stat))
            out_z_acc = os.path.join(a_tmpdir,
                                     "%d.%s.peer_z.acc" % (sim_id, stat))
            bbp_formatter.bbp2peer(bbpfile, out_n_acc, out_e_acc, out_z_acc, accel=True)

            # Duration output files
            out_n_arias = os.path.join(a_outdir,
                                     "%d.%s_N.arias" % (sim_id, stat))
            out_e_arias = os.path.join(a_outdir,
                                     "%d.%s_E.arias" % (sim_id, stat))
            out_z_arias = os.path.join(a_outdir,
                                     "%d.%s_Z.arias" % (sim_id, stat))

            # compute each one
            n_max, n_time_5_75, n_time_5_95, n_time_20_80 = ad_from_acc(out_n_acc, out_n_arias)
            e_max, e_time_5_75, e_time_5_95, e_time_20_80 = ad_from_acc(out_e_acc, out_e_arias)
            z_max, z_time_5_75, z_time_5_95, z_time_20_80 = ad_from_acc(out_z_acc, out_z_arias)
            geo_max = math.sqrt(n_max*e_max)
            geo_time_5_75 = math.sqrt(n_time_5_75*e_time_5_75)
            geo_time_5_95 = math.sqrt(n_time_5_95*e_time_5_95)
            geo_time_20_80 = math.sqrt(n_time_20_80*e_time_20_80)

            # write summary file
            out_duration = os.path.join(a_outdir,
                                     "%d.%s.ard" % (sim_id, stat))
            outfile = open(out_duration, "w")
            outfile.write("# Component  Peak Arias (cm/s)  T5-75 (s)  T5-95 (s)  T20-80 (s)\n")
            outfile.write("N %f %f %f %f\n" % (n_max, n_time_5_75, n_time_5_95, n_time_20_80))
            outfile.write("E %f %f %f %f\n" % (e_max, e_time_5_75, e_time_5_95, e_time_20_80))
            outfile.write("Z %f %f %f %f\n" % (z_max, z_time_5_75, z_time_5_95, z_time_20_80))
            outfile.write("GEOM %f %f %f %f\n" % (geo_max, geo_time_5_75, geo_time_5_95, geo_time_20_80))
            outfile.close()

        # All done!
        print("AriasDuration Completed".center(80, '-'))

if __name__ == '__main__':
    print("Testing Module: %s" % (os.path.basename(sys.argv[0])))
    ME = AriasDuration(sys.argv[1], sim_id=int(sys.argv[2]))
    ME.run()
    sys.exit(0)
