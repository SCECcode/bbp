#!/usr/bin/env python
"""
Copyright 2010-2018 University Of Southern California

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

Utilities converting PEER seismograms into bbp format, and bbp
into PEER format.  Designed specifically for use on SCEC
validation study that needs to convert 6 header line PEER-format
acceleration time series into bbp format and bbp format to a 6
header line PEER format.

Created on Aug 17, 2012
@author: maechlin
"""
from __future__ import division, print_function

# Number of header lines in PEER files
PEER_HEADER_LINES = 6

# Import Python modules
import sys

# Import Broadband modules
import bband_utils

def peer2bbp(in_peer_n_file, in_peer_e_file, in_peer_z_file, out_bbp_file):
    """
    This function converts the 3 input peer files (N/E/Z) to a
    3-component bbp file
    """
    n_vals = []
    e_vals = []
    z_vals = []
    header_lines = []
    dt_vals = []
    start_line = 0
    start_samples = False

    # Read input files
    pfn = open(in_peer_n_file, "r")
    lines_n = pfn.readlines()
    pfn.close()
    pfe = open(in_peer_e_file, "r")
    lines_e = pfe.readlines()
    pfe.close()
    pfz = open(in_peer_z_file, "r")
    lines_z = pfz.readlines()
    pfz.close()

    #
    # Test for empty input file and return error if found
    if len(lines_n) < 1:
        raise bband_utils.ProcessingError("Input file %s is empty!" %
                                          (in_peer_n_file))
    elif len(lines_e) != len(lines_n):
        raise bband_utils.ProcessingError("N and E peer files do not have "
                                          "same number of lines!")
    elif len(lines_z) != len(lines_n):
        raise bband_utils.ProcessingError("N and Z peer files do not have "
                                          "same number of lines!")
    else:
        pass
        #print("Input PEER-format seismogram files with len: %d" %
        #      (len(lines_n)))

    while not start_samples:
        for line in lines_n:
            start_line = start_line + 1
            elems = line.split()
            if len(elems) > 0:
                if elems[0] == "ACCELERATION" or elems[0] == 'Acceleration':
                    start_samples = True
                    break
            header_lines.append(line)

    if not start_samples:
        raise bband_utils.ProcessingError("No samples found in peer file: %s" %
                                          (in_peer_n_file))
    else:
        #print "starting line: %d"%(start_line)
        cur_line = 0
        pts = 0
        dt = 0.0
        for line in lines_n:
            cur_line = cur_line + 1
            if cur_line == start_line + 1:
                pts_dt = line.split()
                pts = float(pts_dt[0])
                dt = float(pts_dt[1])
                #print("Reading Seismogram with NPTS: %d and DT: %f" %
                #      (pts, dt))
            elif cur_line > start_line + 1:
                vals = line.split()
                for x in vals:
                    n_vals.append(float(x) * bband_utils.G2CMSS)
            else:
                # This block will skip the header lines use
                # ACCELERATION TAG as starting key
                pass

    cur_line = 0
    for line in lines_e:
        cur_line = cur_line + 1
        if cur_line > start_line + 1:
            vals = line.split()
            for y in vals:
                e_vals.append(float(y) * bband_utils.G2CMSS)
            else:
                # Skip header lines
                pass

    cur_line = 0
    for line in lines_z:
        cur_line = cur_line + 1
        if cur_line > start_line + 1:
            vals = line.split()
            for z in vals:
                z_vals.append(float(z) * bband_utils.G2CMSS)
        else:
            # Skip header lines
            pass
    #
    # Populate the dt values for bbp format
    # Use number of pts in time series as counter
    #
    dt_val = float(dt)
    pts_count = len(n_vals)
    #
    # Define the first sample at time dt
    cur_dt = 0.0
    for x in range (pts_count):
        dt_vals.append(cur_dt)
        cur_dt = cur_dt + dt_val
    #
    #
    if len(dt_vals) == pts_count:
        pass
        #print("Read format consistent time series with %d samples." %
        #      (len(dt_vals)))
    else:
        print("Inconsistent time series: %d %d" % (len(dt_vals), pts_count))

    #
    # Print files in bbp
    #
    # Open file
    bbp_file = open(out_bbp_file, "w")
    # Write header
    bbp_file.write("#    time(sec)      N-S(cm/s/s)      E-W(cm/s/s)      U-D(cm/s/s)\n")
    #    for line in header_lines:
    #        bbp_file.write("# %s" % (line))
    #    bbp_file.write("# Column 1: Time (s)\n")
    #    bbp_file.write("# Column 2: North-south acceleration (cm/s/s) (+ is northward)\n")
    #    bbp_file.write("# Column 3: East-west acceleration (cm/s/s) (+ is eastward)\n")
    #    bbp_file.write("# Column 4: Up-down acceleration (cm/s/s) (+ is upward)\n#\n")
    #    bbp_file.write("# NPTS  DT  \n")
    #    bbp_file.write("# %d %s\n" % (pts_count, dt))
    # Write the data
    for x in range(pts_count):
        bbp_file.write("%7e   % 8e   % 8e   % 8e\n" % (dt_vals[x], n_vals[x],
                                                       e_vals[x], z_vals[x]))
    # Lastly, close the file
    bbp_file.close()

def bbp2peer(in_bbp_file, out_peer_n_file, out_peer_e_file, out_peer_z_file):
    """
    Convert bbp file into three peer files for use by RotD50 and
    other programs that input PEER format seismograms
    """
    num_header_lines = 0
    bbp = open(in_bbp_file, "r")
    lines = bbp.readlines()
    bbp.close()
    for line in lines:
        elems = line.split()
        if elems[0] == "#":
            num_header_lines = num_header_lines + 1
        else:
            break

    #print "Counted %d header lines" % (num_header_lines)

    cur_line = 0
    header_lines = []
    dt_vals = []
    n_vals = []
    e_vals = []
    z_vals = []
    for line in lines:
        if cur_line <= (num_header_lines - 2):
            cur_line = cur_line + 1
            header_lines.append(line)
        elif cur_line == (num_header_lines - 1):
            # Last line before start gives dt and npts
            cur_line = cur_line + 1
            # Print line
            elems = line.split()
            try:
                npts = int(elems[1])
                dt = float(elems[2])
            except (IndexError, ValueError) as err:
                # Ok, it doesn't seem we are re-converting to PEER
                # format a previously PEER-BBP converted file
                # Let's try to figure things out...
                npts = len(lines) - num_header_lines
                try:
                    time_1 = float(lines[num_header_lines].split()[0])
                    time_2 = float(lines[num_header_lines+1].split()[0])
                except ValueError:
                    print("Cannot figure out npts and dt from this bbp file!")
                    sys.exit(-1)
                dt = time_2 - time_1
            #print("Reformating BBP file with %d header lines." %
            #      (num_header_lines))
            #print("Reformating BBP file with dt: %f " % (dt))
            #print("Reformating BBP file with npts: %d" % (npts))
        else:
            elems = line.split()
            if len(elems) != 4:
                raise bband_utils.ProcessingError("Unexpected BBP time series "
                                                  "line format found."
                                                  "Error in conversion.")
            else:
                dt_vals.append(dt)
                n_vals.append(float(elems[1]) / bband_utils.G2CMSS)
                e_vals.append(float(elems[2]) / bband_utils.G2CMSS)
                z_vals.append(float(elems[3]) / bband_utils.G2CMSS)

    # Prepare to write 6 colume format
    n_file = open(out_peer_n_file, "w")
    e_file = open(out_peer_e_file, "w")
    z_file = open(out_peer_z_file, "w")

    #n_file.write("Created by: bbp2peer v12.8.0\n")
    #e_file.write("Created by: bbp2peer v12.8.0\n")
    #z_file.write("Created by: bbp2peer v12.8.0\n")

    # Adjust header lines, so we always have enough
    while len(header_lines) <= (PEER_HEADER_LINES - 2):
        header_lines.append("\n")

    for line in header_lines[0:(PEER_HEADER_LINES - 2)]:
        n_file.write(line)
        e_file.write(line)
        z_file.write(line)

    n_file.write("Acceleration in g\n")
    n_file.write("  %d   %1.6f   NPTS, DT\n" % (npts, dt))
    e_file.write("Acceleration in g\n")
    e_file.write("  %d   %1.6f   NPTS, DT\n" % (npts, dt))
    z_file.write("Acceleration in g\n")
    z_file.write("  %d   %1.6f   NPTS, DT\n" % (npts, dt))

    cur_line = 0
    for index, elem in enumerate(dt_vals):
        n_file.write("% 12.7E " % (n_vals[index]))
        e_file.write("% 12.7E " % (e_vals[index]))
        z_file.write("% 12.7E " % (z_vals[index]))
        #print "%f"%(dt_vals[index])
        #print "%e"%(n_vals[index])
        #print "%e"%(e_vals[index])
        #print "%e"%(z_vals[index])
        if (index % 5) == 4:
            n_file.write("\n")
            e_file.write("\n")
            z_file.write("\n")
        #else:
        #    n_file.write("\t")
        #    e_file.write("\t")
        #    z_file.write("\t")

    # Add newline at the end of last line to avoid issue when rotd50.f
    # reads the file (only when compiled with gfortran 4.3.3 on HPCC)
    n_file.write("\n")
    e_file.write("\n")
    z_file.write("\n")
    # Close all files
    n_file.close()
    e_file.close()
    z_file.close()

def exsim2bbp(in_exsim_n, in_exsim_e, in_exsim_z, out_bbp_file):
    """
    Converts 3 exsim acc format files into 3 component bbp file.
    Assumes all three files are for the same site.
    """

    n_vals = []
    e_vals = []
    z_vals = []
    header_lines = []
    dt_vals = []

    pfn = open(in_exsim_n, "r")
    lines_n = pfn.readlines()
    pfn.close()
    pfe = open(in_exsim_e, "r")
    lines_e = pfe.readlines()
    pfe.close()
    pfz = open(in_exsim_z, "r")
    lines_z = pfz.readlines()
    pfz.close()

    #
    # Test for empty input file and return error if found
    if len(lines_n) < 1:
        raise bband_utils.ProcessingError("Input file %s is empty!" %
                                          (in_exsim_n))
    elif len(lines_e) != len(lines_n):
        raise bband_utils.ProcessingError("N and E peer files do not have "
                                          "same number of lines!")
    elif len(lines_z) != len(lines_n):
        raise bband_utils.ProcessingError("N and Z peer files do not have "
                                          "same number of lines!")
    else:
        print("Input EXSIM-format seismogram files with len: %d" %
              (len(lines_n)))

    pts = 0
    dt = 0.0
    start_line = 0
    start_samples = False

    while not start_samples:
        for line in lines_n:
            start_line = start_line + 1
            elems = line.split()
            if len(elems) == 2 and elems[1] == "samples":
                pts = int(elems[0])
            elif len(elems) >= 2 and elems[0] == "dt:":
                dt = float(elems[1])
            elif len(elems) >= 2 and elems[0] == "time(s)":
                start_samples = True
                break
            else:
                header_lines.append(line)

    if not start_samples:
        print("No samples found in peer file: " + in_exsim_n)
        sys.exit(0)
    else:
        print("starting line: %d" % (start_line))
        print("Reading Seismogram with NPTS: %d and DT: %f" % (pts, dt))

    cur_line = 0
    for line in lines_n:
        cur_line = cur_line + 1
        if cur_line > start_line:
            vals = line.split()
            n_vals.append(float(vals[1]))

    cur_line = 0
    for line in lines_e:
        cur_line = cur_line + 1
        if cur_line > start_line:
            vals = line.split()
            e_vals.append(float(vals[1]))

    cur_line = 0
    for line in lines_z:
        cur_line = cur_line + 1
        if cur_line > start_line:
            vals = line.split()
            z_vals.append(float(vals[1]))

    #
    # Populate the dt values for bbp format
    # Use number of pts in time series as counter
    #
    dt_val = dt
    pts_count = len(n_vals)
    #
    # Define the first sample at time dt
    cur_dt = 0.0
    for x in range (pts_count):
        dt_vals.append(cur_dt)
        cur_dt = cur_dt + dt_val
    #
    #
    if len(dt_vals) == pts_count:
        print("Read format consistent time series with %d samples." %
              (len(dt_vals)))
    else:
        print("Inconsistent time series: %d %d" % (len(dt_vals), pts_count))

    #
    # Print bbp file
    #
    bbp_file = open(out_bbp_file, "w")
    # Print header
    bbp_file.write("#    time(sec)      N-S(cm/s/s)      E-W(cm/s/s)      U-D(cm/s/s)\n")
    #for line in header_lines:
    #    bbp_file.write("# %s" % (line))
    # bbp_file.write("# Column 1: Time (s)\n")
    # bbp_file.write("# Column 2: North-south acceleration (cm/s/s) (+ is northward)\n")
    # bbp_file.write("# Column 3: East-west acceleration (cm/s/s) (+ is eastward)\n")
    # bbp_file.write("# Column 4: Up-down acceleration (cm/s/s) (+ is upward)\n#\n")
    # bbp_file.write("# NPTS  DT  \n")
    # bbp_file.write("# %d %s\n" % (pts_count, dt))
    for x in range(pts_count):
        bbp_file.write("%7e   % 8e   % 8e   % 8e\n" % (dt_vals[x], n_vals[x],
                                                       e_vals[x], z_vals[x]))

    bbp_file.close()

if __name__ == '__main__':
    print("Creating BBP Formatter")
