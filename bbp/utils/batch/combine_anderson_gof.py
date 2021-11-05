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
import glob
import optparse
import numpy as np
import matplotlib
if (matplotlib.get_backend() != 'agg'):
    matplotlib.use('Agg') # Disable use of Tk/X11
import pylab as py
import matplotlib.pyplot as plt
from pylab import arange
from scipy import stats

# Import Broadband modules
import bband_utils
from station_list import StationList

BNAMES = ['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8', 'B9', 'B10']
BMAX = len(BNAMES)

# ****************************************************************************
def statts(x_in):
    """
    Calculates statistics
    """
    s = x_in[~np.isnan(x_in)]
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

# ****************************************************************************
def fplots(eventname, method, S1, c1cf, c2cf, c3cf, c4cf, c5cf,
           c6cf, c7cf, c8cf, c9cf, c10cf, output_file):
    """
    This function creates a combined GoF plot for all stations
    """

    fig = py.figure()
    S1 = "{:3.1f}".format(S1)
    title = eventname + ' - ' + method + ' - Score S1 : ' + str(S1)
    fig.suptitle(title, fontsize=18)
    x_vals = np.arange(1, BMAX + 1)
    label_size = 9
    matplotlib.rcParams['xtick.labelsize'] = label_size
    matplotlib.rcParams['ytick.labelsize'] = label_size

    mean = [c1cf[i][0] for i in arange(BMAX)]
    #stdev = [c1cf[i][1] for i in arange(BMAX)]
    low_std = [c1cf[i][0] - c1cf[i][1] for i in arange(BMAX)]
    upp_std = [c1cf[i][0] + c1cf[i][1] for i in arange(BMAX)]
    low_conf = [c1cf[i][2] for i in arange(BMAX)]
    upp_conf = [c1cf[i][3] for i in arange(BMAX)]
    plt.subplot(5, 2, 1)
    plt.grid(True, which='both', color='0.35', ls=":")
    plt.xlim(0, 11)
    plt.xticks(arange(1, 11))
    plt.yticks(arange(0, 12, 2))
    plt.ylim(0, 10)
    plt.ylabel('C1', fontsize=12, labelpad=-1)
    plt.plot(x_vals, mean, 'ko', ms=4)
    plt.plot(x_vals, mean, 'k-', mew=4)
    plt.fill_between(x_vals, low_conf, upp_conf,
                     color=(1, 0.3, 0.9), alpha='1.0')
    plt.fill_between(x_vals, low_std, upp_std,
                     color=(0.8, 0.8, 1), alpha='1.0')

    mean = [c2cf[i][0] for i in arange(BMAX)]
    #stdev = [c2cf[i][1] for i in arange(BMAX)]
    low_std = [c2cf[i][0] - c2cf[i][1] for i in arange(BMAX)]
    upp_std = [c2cf[i][0] + c2cf[i][1] for i in arange(BMAX)]
    low_conf = [c2cf[i][2] for i in arange(BMAX)]
    upp_conf = [c2cf[i][3] for i in arange(BMAX)]
    plt.subplot(5, 2, 2)
    plt.grid(True, which='both', color='0.35', ls=":")
    plt.xlim(0, 11)
    plt.ylim(0, 10)
    plt.xticks(arange(1, 11))
    plt.yticks(arange(0, 12, 2))
    plt.ylabel('C2', fontsize=12, labelpad=-1)
    plt.plot(x_vals, mean, 'ko', ms=4)
    plt.plot(x_vals, mean, 'k-', mew=4)
    plt.fill_between(x_vals, low_conf, upp_conf,
                     color=(1, 0.3, 0.9), alpha='1.0')
    plt.fill_between(x_vals, low_std, upp_std,
                     color=(0.8, 0.8, 1), alpha='1.0')

    mean = [c3cf[i][0] for i in arange(BMAX)]
    #stdev = [c3cf[i][1] for i in arange(BMAX)]
    low_std = [c3cf[i][0] - c3cf[i][1] for i in arange(BMAX)]
    upp_std = [c3cf[i][0] + c3cf[i][1] for i in arange(BMAX)]
    low_conf = [c3cf[i][2] for i in arange(BMAX)]
    upp_conf = [c3cf[i][3] for i in arange(BMAX)]
    plt.subplot(5, 2, 3)
    plt.grid(True, which='both', color='0.35', ls=":")
    plt.xlim(0, 11)
    plt.ylim(0, 10)
    plt.xticks(arange(1, 11))
    plt.yticks(arange(0, 12, 2))
    plt.ylabel('C3', fontsize=12, labelpad=-1)
    plt.plot(x_vals, mean, 'ko', ms=4)
    plt.plot(x_vals, mean, 'k-', mew=4)
    plt.fill_between(x_vals, low_conf, upp_conf,
                     color=(1, 0.3, 0.9), alpha='1.0')
    plt.fill_between(x_vals, low_std, upp_std,
                     color=(0.8, 0.8, 1), alpha='1.0')

    mean = [c4cf[i][0] for i in arange(BMAX)]
    #stdev = [c4cf[i][1] for i in arange(BMAX)]
    low_std = [c4cf[i][0] - c4cf[i][1] for i in arange(BMAX)]
    upp_std = [c4cf[i][0] + c4cf[i][1] for i in arange(BMAX)]
    low_conf = [c4cf[i][2] for i in arange(BMAX)]
    upp_conf = [c4cf[i][3] for i in arange(BMAX)]
    plt.subplot(5, 2, 4)
    plt.grid(True, which='both', color='0.35', ls=":")
    plt.xlim(0, 11)
    plt.ylim(0, 10)
    plt.xticks(arange(1, 11))
    plt.yticks(arange(0, 12, 2))
    plt.ylabel('C4', fontsize=12, labelpad=-1)
    plt.plot(x_vals, mean, 'ko', ms=4)
    plt.plot(x_vals, mean, 'k-', mew=4)
    plt.fill_between(x_vals, low_conf, upp_conf,
                     color=(1, 0.3, 0.9), alpha='1.0')
    plt.fill_between(x_vals, low_std, upp_std,
                     color=(0.8, 0.8, 1), alpha='1.0')

    mean = [c5cf[i][0] for i in arange(BMAX)]
    #stdev = [c5cf[i][1] for i in arange(BMAX)]
    low_std = [c5cf[i][0] - c5cf[i][1] for i in arange(BMAX)]
    upp_std = [c5cf[i][0] + c5cf[i][1] for i in arange(BMAX)]
    low_conf = [c5cf[i][2] for i in arange(BMAX)]
    upp_conf = [c5cf[i][3] for i in arange(BMAX)]
    plt.subplot(5, 2, 5)
    plt.grid(True, which='both', color='0.35', ls=":")
    plt.xlim(0, 11)
    plt.ylim(0, 10)
    plt.xticks(arange(1, 11))
    plt.yticks(arange(0, 12, 2))
    plt.ylabel('C5', fontsize=12, labelpad=-1)
    plt.plot(x_vals, mean, 'ko', ms=4)
    plt.plot(x_vals, mean, 'k-', mew=4)
    plt.fill_between(x_vals, low_conf, upp_conf,
                     color=(1, 0.3, 0.9), alpha='1.0')
    plt.fill_between(x_vals, low_std, upp_std,
                     color=(0.8, 0.8, 1), alpha='1.0')

    mean = [c6cf[i][0] for i in arange(BMAX)]
    #stdev = [c6cf[i][1] for i in arange(BMAX)]
    low_std = [c6cf[i][0] - c6cf[i][1] for i in arange(BMAX)]
    upp_std = [c6cf[i][0] + c6cf[i][1] for i in arange(BMAX)]
    low_conf = [c6cf[i][2] for i in arange(BMAX)]
    upp_conf = [c6cf[i][3] for i in arange(BMAX)]
    plt.subplot(5, 2, 6)
    plt.grid(True, which='both', color='0.35', ls=":")
    plt.xlim(0, 11)
    plt.ylim(0, 10)
    plt.xticks(arange(1, 11))
    plt.yticks(arange(0, 12, 2))
    plt.ylabel('C6', fontsize=12, labelpad=-1)
    plt.plot(x_vals, mean, 'ko', ms=4)
    plt.plot(x_vals, mean, 'k-', mew=4)
    plt.fill_between(x_vals, low_conf, upp_conf,
                     color=(1, 0.3, 0.9), alpha='1.0')
    plt.fill_between(x_vals, low_std, upp_std,
                     color=(0.8, 0.8, 1), alpha='1.0')

    mean = [c7cf[i][0] for i in arange(BMAX)]
    #stdev = [c7cf[i][1] for i in arange(BMAX)]
    low_std = [c7cf[i][0] - c7cf[i][1] for i in arange(BMAX)]
    upp_std = [c7cf[i][0] + c7cf[i][1] for i in arange(BMAX)]
    low_conf = [c7cf[i][2] for i in arange(BMAX)]
    upp_conf = [c7cf[i][3] for i in arange(BMAX)]
    plt.subplot(5, 2, 7)
    plt.grid(True, which='both', color='0.35', ls=":")
    plt.xlim(0, 11)
    plt.ylim(0, 10)
    plt.xticks(arange(1, 11))
    plt.yticks(arange(0, 12, 2))
    plt.ylabel('C7', fontsize=12, labelpad=-1)
    plt.plot(x_vals, mean, 'ko', ms=4)
    plt.plot(x_vals, mean, 'k-', mew=4)
    plt.fill_between(x_vals, low_conf, upp_conf,
                     color=(1, 0.3, 0.9), alpha='1.0')
    plt.fill_between(x_vals, low_std, upp_std,
                     color=(0.8, 0.8, 1), alpha='1.0')

    mean = [c8cf[i][0] for i in arange(BMAX)]
    #stdev = [c8cf[i][1] for i in arange(BMAX)]
    low_std = [c8cf[i][0] - c8cf[i][1] for i in arange(BMAX)]
    upp_std = [c8cf[i][0] + c8cf[i][1] for i in arange(BMAX)]
    low_conf = [c8cf[i][2] for i in arange(BMAX)]
    upp_conf = [c8cf[i][3] for i in arange(BMAX)]
    plt.subplot(5, 2, 8)
    plt.grid(True, which='both', color='0.35', ls=":")
    plt.xlim(0, 11)
    plt.ylim(0, 10)
    plt.xticks(arange(1, 11))
    plt.yticks(arange(0, 12, 2))
    plt.ylabel('C8', fontsize=12, labelpad=-1)
    plt.plot(x_vals, mean, 'ko', ms=4)
    plt.plot(x_vals, mean, 'k-', mew=4)
    plt.fill_between(x_vals, low_conf, upp_conf,
                     color=(1, 0.3, 0.9), alpha='1.0')
    plt.fill_between(x_vals, low_std, upp_std,
                     color=(0.8, 0.8, 1), alpha='1.0')

    mean = [c9cf[i][0] for i in arange(BMAX)]
    #stdev = [c9cf[i][1] for i in arange(BMAX)]
    low_std = [c9cf[i][0] - c9cf[i][1] for i in arange(BMAX)]
    upp_std = [c9cf[i][0] + c9cf[i][1] for i in arange(BMAX)]
    low_conf = [c9cf[i][2] for i in arange(BMAX)]
    upp_conf = [c9cf[i][3] for i in arange(BMAX)]
    plt.subplot(5, 2, 9)
    plt.grid(True, which='both', color='0.35', ls=":")
    plt.xlim(0, 11)
    plt.ylim(0, 10)
    plt.xticks(arange(1, 11))
    plt.yticks(arange(0, 12, 2))
    plt.ylabel('C9', fontsize=12, labelpad=-1)
    plt.plot(x_vals, mean, 'ko', ms=4)
    plt.xlabel('Frequency Band', fontsize=12)
    plt.plot(x_vals, mean, 'k-', mew=4)
    plt.fill_between(x_vals, low_conf, upp_conf,
                     color=(1, 0.3, 0.9), alpha='1.0')
    plt.fill_between(x_vals, low_std, upp_std,
                     color=(0.8, 0.8, 1), alpha='1.0')

    mean = [c10cf[i][0] for i in arange(BMAX)]
    #stdev = [c10cf[i][1] for i in arange(BMAX)]
    low_std = [c10cf[i][0] - c10cf[i][1] for i in arange(BMAX)]
    upp_std = [c10cf[i][0] + c10cf[i][1] for i in arange(BMAX)]
    low_conf = [c10cf[i][2] for i in arange(BMAX)]
    upp_conf = [c10cf[i][3] for i in arange(BMAX)]
    plt.subplot(5, 2, 10)
    plt.grid(True, which='both', color='0.35', ls=":")
    plt.xlim(0, 11)
    plt.ylim(0, 10)
    plt.xticks(arange(1, 11))
    plt.yticks(arange(0, 12, 2))
    plt.ylabel('C10', fontsize=12, labelpad=-1)
    plt.plot(x_vals, mean, 'ko', ms=4)
    plt.plot(x_vals, mean, 'k-', mew=4)
    plt.xlabel('Frequency Band', fontsize=12)
    plt.fill_between(x_vals, low_conf, upp_conf,
                     color=(1, 0.3, 0.9), alpha='1.0')
    plt.fill_between(x_vals, low_std, upp_std,
                     color=(0.8, 0.8, 1), alpha='1.0')

    #plt.show()
    plt.savefig(output_file, dpi=300)
    plt.close(fig)

    return

def readfmt(filein):
    """
    This function reads a data file generated by the Anderson GOF
    module and parses its data
    """
    i = 0
    acc = np.empty([10, 10])
    acc[:] = np.nan
    input_file = open(filein, 'r')

    # Read and ignore header lines
    for k in range(1):
        # Skip header
        input_file.readline()
    # Loop over lines and extract variables of interest
    for line in input_file:
        line = line.strip()
        columns = line.split()
        for j in range(len(columns)-1):
            acc[i, j] = float(columns[j+1])
        i += 1
    input_file.close()
    return acc

def load_station_data(input_outdir, data, station):
    """
    Load station data into the DATA dictionary
    """
    # Get realizations
    realizations = sorted(os.listdir(input_outdir))
    for realization in realizations:
        basedir = os.path.join(input_outdir, realization,
                               "validations", "anderson_gof")
        data_file = glob.glob("%s%sgof-*-anderson-%s.txt" % (basedir,
                                                             os.sep,
                                                             station))
        if len(data_file) != 1:
            raise bband_utils.ProcessingError("Data for station %s " %
                                              (station) +
                                              "not found for "
                                              "realization %s!" %
                                              (realization))
        data_file = data_file[0]
        stat_data = readfmt(data_file)
        # Add data to the data dict
        for c_idx in range(10):
            for b_idx in range(BMAX):
                if station not in data[c_idx][b_idx]:
                    # Add 3rd level if it doesn't exist - stations
                    data[c_idx][b_idx][station] = []
                data[c_idx][b_idx][station].append(stat_data[b_idx, c_idx])

def load_all_data(input_indir, input_outdir):
    """
    This function goes through all realizations and loads all data to
    the DATA dictionary
    """
    # First create data dictionary
    data = {}
    # First level is C1..CMAX
    for i in range(10):
        data[i] = {}
        # Second level is B1..BMAX
        for j in range(BMAX):
            data[i][j] = {}

    # Get realizations
    realizations = sorted(os.listdir(input_indir))
    one_realization = realizations[0]
    basedir = os.path.join(input_indir, one_realization)

    # Get the station list
    a_statfile = glob.glob("%s%s*.stl" % (basedir, os.sep))
    if len(a_statfile) != 1:
        raise bband_utils.ProcessingError("Cannot get station list!")
    a_statfile = a_statfile[0]
    slo = StationList(a_statfile)
    site_list = slo.getStationList()

    # Go through all stations
    for site in site_list:
        station = site.scode
        print("working on station: %s" % (station))

        # Read data for this station
        load_station_data(input_outdir, data, station)

    # Return data dictionary
    return data

def process_data(data):
    """
    This function processes the data in the data dictionary
    """
    for c_idx in range(10):
        for b_idx in range(BMAX):
            data[c_idx][b_idx]["st_means"] = []
            for station in data[c_idx][b_idx]:
                # Calculate the mean per station
                data[c_idx][b_idx]["st_means"].append(np.nanmean(data[c_idx][b_idx][station]))

    # Now, data[c_idx][b_idx]["st_means"] has an array of station means

# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

PARSER = optparse.OptionParser()
PARSER.add_option("-d", "--dir", dest="input_dir",
                  help="Input directory containing simulation results")
PARSER.add_option("-o", "--output_dir", dest="output_dir",
                  help="Output file")
PARSER.add_option("-c", "--codebase", dest="codebase",
                  help="Method used for the simulation")
(OPTIONS, ARGS) = PARSER.parse_args()


if OPTIONS.input_dir is None:
    PARSER.error("Please specify the input directory!")
TOP_INPUT_DIR = OPTIONS.input_dir
if not os.path.isdir(TOP_INPUT_DIR):
    PARSER.error("Invalid input directory!")
if not "Sims" in os.listdir(TOP_INPUT_DIR):
    PARSER.error("Please provide the top-level simulation directory!\n"
                 "This is the directory given to the cluster script")
INPUT_OUTDIR = os.path.join(TOP_INPUT_DIR, "Sims", "outdata")
INPUT_TMPDIR = os.path.join(TOP_INPUT_DIR, "Sims", "tmpdata")
INPUT_INDIR = os.path.join(TOP_INPUT_DIR, "Sims", "indata")

if OPTIONS.output_dir is None:
    PARSER.error("error specify output directory!")
else:
    OUTPUT_DIR = OPTIONS.output_dir
    if not os.path.isdir(OUTPUT_DIR):
        PARSER.error("Invalid output directory!")

if OPTIONS.codebase is None:
    PARSER.error("Please specify codebase!")

evnt = os.path.basename(os.path.normpath(OPTIONS.input_dir)).split('-')
eventcode = evnt[1]
# Create data files with both gmpe and simulation data
DATA = load_all_data(INPUT_INDIR, INPUT_OUTDIR)

# Process data
process_data(DATA)

for iB in range(BMAX):
    c1conf = []
    c2conf = []
    c3conf = []
    c4conf = []
    c5conf = []
    c6conf = []
    c7conf = []
    c8conf = []
    c9conf = []
    c10conf = []

for iB in range(BMAX):
    c1conf.append(statts(np.asarray(DATA[0][iB]["st_means"])))
    c2conf.append(statts(np.asarray(DATA[1][iB]["st_means"])))
    c3conf.append(statts(np.asarray(DATA[2][iB]["st_means"])))
    c4conf.append(statts(np.asarray(DATA[3][iB]["st_means"])))
    c5conf.append(statts(np.asarray(DATA[4][iB]["st_means"])))
    c6conf.append(statts(np.asarray(DATA[5][iB]["st_means"])))
    c7conf.append(statts(np.asarray(DATA[6][iB]["st_means"])))
    c8conf.append(statts(np.asarray(DATA[7][iB]["st_means"])))
    c9conf.append(statts(np.asarray(DATA[8][iB]["st_means"])))
    c10conf.append(statts(np.asarray(DATA[9][iB]["st_means"])))

output_file = os.path.join(OPTIONS.output_dir,
                           'gof-Anderson-%s-%s-Summary.dat' %
                           (OPTIONS.codebase, eventcode))

out_file = open(output_file, 'w')
line = ('#%s%5s%4s%4s%4s%4s%4s%4s%4s%4s%4s\n' %
        ('band', 'C1', 'C2', 'C3', 'C4', 'C5',
         'C6', 'C7', 'C8', 'C9', 'C10'))
out_file.write(line)

for i in range(BMAX):
    line = ('%s %3.1f %3.1f %3.1f %3.1f %3.1f %3.1f %3.1f %3.1f %3.1f %3.1f\n' %
            (BNAMES[i], c1conf[i][0], c2conf[i][0], c3conf[i][0],
             c4conf[i][0], c5conf[i][0], c6conf[i][0], c7conf[i][0],
             c8conf[i][0], c9conf[i][0], c10conf[i][0]))
    out_file.write(line)
out_file.close()

s1_event = np.nanmean([DATA[iC][iB]["st_means"] for iC in range(10) for iB in range(BMAX)])
print('The S1 score for this event is ', "{:3.1f}".format(s1_event))

# Plot data
output_file = os.path.join(OPTIONS.output_dir,
                           "gof-Anderson-%s-%s-Summary.png" %
                           (OPTIONS.codebase, eventcode))
fplots(eventcode, OPTIONS.codebase, s1_event,
       np.asarray(c1conf), np.asarray(c2conf),
       np.asarray(c3conf), np.asarray(c4conf),
       np.asarray(c5conf), np.asarray(c6conf),
       np.asarray(c7conf), np.asarray(c8conf),
       np.asarray(c9conf), np.asarray(c10conf),
       output_file)

print("All Done!")
# Clean-up, all done!
