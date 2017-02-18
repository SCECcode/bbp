"""
Southern California Earthquake Center Broadband Platform
Copyright 2010-2016 Southern California Earthquake Center

Created on Aug 7, 2012
@author: maechlin
Purpose: Translation of arias duration algorithm in matlab from
         Christine Goulet for use on bbp project

$Id: arias_duration.py 1801 2017-02-14 23:18:41Z fsilva $
"""
from __future__ import division, print_function

# G2CMSS = 980.665 # Convert g to cm/s/s
import math
import scipy.integrate

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

    # Define the time for AI=5%, 75%, 95%
    time_ai5 = 0
    for i in range(pts):
        if ia_norm[i] >= 5:
            break
        time_ai5 = dt * i

    time_ai75 = 0
    for i in range(pts):
        if ia_norm[i] >= 75:
            break
        time_ai75 = dt * i

    time_ai95 = 0
    for i in range(pts):
        if ia_norm[i] >= 95:
            break
        time_ai95 = dt * i

    # Now, calculate the arias intervals 5% to 75% and 5% to 95%
    time_5_75 = time_ai75 - time_ai5
    time_5_95 = time_ai95 - time_ai5

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
    outfile.write("# Arias Intervals: T5-75 %f (s), T5-95 %f (secs)\n" %
                  (time_5_75, time_5_95))
    outfile.write("# Seconds Accel (g) "
                  "Arias Intensity (cm/s) "
                  "ADNormalized (%)\n")
    outfile.write("# %d %f  NPTS, DT\n" % (pts, dt))
    for i in range(pts):
        outfile.write("% 8f % 8f % 8f % 8f\n" %
                      (dt_vals[i], acc[i], arias_intensity[i], ia_norm[i]))
    outfile.close()
    return 0
