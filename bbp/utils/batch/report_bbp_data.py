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
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys
import glob
import numpy
import warnings

# warnings.simplefilter('error')

# Period bins
PERIODS = [[0.01, 0.1],
           [0.11, 1.0],
           [1.1, 3.0],
           [3.5, 10.0]]

# Distance bins
DIST = [[0.0, 5.0],
        [5.0, 20.0],
        [20.0, 70.0],
        [70.0, 300.0]]

# Mechanism dictionary
MECH = [("REV", ["Niigata", "NR", "WHITTIER",
                 "SanSimeon", "Iwate",
                 "RDL1K", "Mineral", "Saguenay1k"]),
        ("ROBL", ["NORTHPS", "LOMAP", "CHINOH", "Chuetsu"]),
        ("SS", ["Tottori", "Landers", "ALUMR",
                "Parkfield", "HectorMine", "Ridgecrest19A",
                "Ridgecrest19B", "Ridgecrest19C"]),
        ("NM", [])]

# Event list
EVENTS = ["CHINOH", "ALUMR", "WHITTIER", "Parkfield",
          "NORTHPS", "Tottori", "SanSimeon", "Niigata",
          "Chuetsu", "NR", "Iwate", "LOMAP", "Landers",
          "Ridgecrest19A", "Ridgecrest19B", "Ridgecrest19C",
          "RDL1K", "Mineral", "Saguenay1k", "HectorMine"]

EVENTS_CA = ["WHITTIER", "NORTHPS", "NR", "LOMAP", "Landers",
             "CHINOH", "ALUMR", "SanSimeon", "Parkfield",
             "HectorMine", "Ridgecrest19A", "Ridgecrest19B",
             "Ridgecrest19C"]

EVENTS_CENA = ["Mineral", "RDL1K", "Saguenay1k"]

EVENTS_EQUIV = {"LandersMSDr" : "Landers",
                "LandersMS" : "Landers",
                "Ridgecrest19CMS" : "Ridgecrest19C"}

def compile_input_data(input_file):
    """
    Reads the input file and combine all the data in the format we
    need to generate the method report
    """
    # Start empty
    data = {}
    t_indexes = [[], [], [], []]
    t_periods = [[], [], [], []]
    # Create top-level dictionaries
    for event in EVENTS:
        data[event] = [[] for _ in DIST]
        for idx, _ in enumerate(DIST):
            data[event][idx] = [[] for _ in PERIODS]

    # Now, process input file
    input = open(input_file, 'r')
    # Read header
    header = input.readline()
    header = header.strip()
    items = header.split()
    for idx, item in enumerate(items):
        try:
            val = float(item)
            for bin_idx, t_range in enumerate(PERIODS):
                if val >= t_range[0] and val <= t_range[1]:
                    # Found the bin to add this index
                    t_indexes[bin_idx].append(idx)
                    t_periods[bin_idx].append(val)
                    break
        except:
            pass

    # Great, we have the indexes that we need
    for line in input:
        line = line.strip()
        items = line.split()
        event = items[0]
        # Map multi-segment events properly
        if event in EVENTS_EQUIV:
            event = EVENTS_EQUIV[event]
        dist = float(items[7])
        tmin = float(items[10])
        tmax = float(items[11])
        if event not in EVENTS:
            print("Unknown event %s, skipping..." % (event))
            continue
        # Make sure we filter the psa5e and psa5n components out
        if line.find("psa5e") > 0 or line.find("psa5n") > 0:
            continue
        dist_bin = -1
        for idx, vals in enumerate(DIST):
            if idx == 0:
                if dist >= vals[0] and dist <= vals[1]:
                    dist_bin = 0
                    break
            if dist > vals[0] and dist <= vals[1]:
                dist_bin = idx
                break
        if dist_bin < 0:
            # Skip this station
            continue
        # Now add the data to the corresponding bins
        for idx, periods in enumerate(t_indexes):
            for idx2, period in enumerate(periods):
#                print t_periods
                if t_periods[idx][idx2] >= tmin and t_periods[idx][idx2] <= tmax:
                    data[event][dist_bin][idx].append(float(items[period]))
    input.close()
    return data

def write_output_file(data):
    """
    This function writes the output data file
    """
    for idx, vals in enumerate(DIST):
        print("Rrup = %.2f-%.2f km" % (vals[0], vals[1]))
        # Initialize data for calculating mean
        all_data = [[] for _ in PERIODS]
        all_data_ca = [[] for _ in PERIODS]
        all_data_cena = [[] for _ in PERIODS]
        for event in EVENTS:
            print("%-15s" % (event), end="")
            for per_range, _ in enumerate(PERIODS):
                event_data = []
                for val in data[event][idx][per_range]:
                    event_data.append(val)
                    # Also keep track of mean across events
                    all_data[per_range].append(val)
                    if event in EVENTS_CA:
                        all_data_ca[per_range].append(val)
                    if event in EVENTS_CENA:
                        all_data_cena[per_range].append(val)
                event_data_abs = [abs(x) for x in event_data]
                if not len(event_data):
                    print("%6s %6s %6s" % ("  N/A", "  N/A", "  N/A"), end="")
                else:
                    print("%6.2f %6.2f %6d" % (numpy.mean(event_data),
                                               numpy.mean(event_data_abs),
                                               len(event_data)), end="")
            print()
        print("%-15s" % ("Average (CA)"), end="")
        for per_data in all_data_ca:
            per_data_abs = [abs(x) for x in per_data]
            if not len(per_data):
                print("%6s %6s %6s" % ("  N/A", "  N/A", "  N/A"), end="")
            else:
                print("%6.2f %6.2f %6d" % (numpy.mean(per_data),
                                           numpy.mean(per_data_abs),
                                           len(per_data)), end="")
        print()
        print("%-15s" % ("Average (CENA)"), end="")
        for per_data in all_data_cena:
            per_data_abs = [abs(x) for x in per_data]
            if not len(per_data):
                print("%6s %6s %6s" % ("  N/A", "  N/A", "  N/A"), end="")
            else:
                print("%6.2f %6.2f %6d" % (numpy.mean(per_data),
                                           numpy.mean(per_data_abs),
                                           len(per_data)), end="")
        print()
        print("%-15s" % ("Average (All)"), end="")
        for per_data in all_data:
            per_data_abs = [abs(x) for x in per_data]
            if not len(per_data):
                print("%6s %6s %6s" % ("  N/A", "  N/A", "  N/A"), end="")
            else:
                print("%6.2f %6.2f %6d" % (numpy.mean(per_data),
                                           numpy.mean(per_data_abs),
                                           len(per_data)), end="")
        print()
        print()

    print("Mechanism")
    all_data = [[] for _ in PERIODS]
    all_data_ca = [[] for _ in PERIODS]
    all_data_cena = [[] for _ in PERIODS]
    for mech in MECH:
        print("%-15s" % (mech[0]), end="")
        events = mech[1]
        for per_range, _ in enumerate(PERIODS):
            event_data = []
            for event in events:
                for dist_range, _ in enumerate(DIST):
                    for val in data[event][dist_range][per_range]:
                        event_data.append(val)
                        all_data[per_range].append(val)
                        if event in EVENTS_CA:
                            all_data_ca[per_range].append(val)
                        if event in EVENTS_CENA:
                            all_data_cena[per_range].append(val)
            event_data_abs = [abs(x) for x in event_data]
            if not len(event_data):
                print("%6s %6s %6s" % ("  N/A", "  N/A", "  N/A"), end="")
            else:
                print("%6.2f %6.2f %6d" % (numpy.mean(event_data),
                                           numpy.mean(event_data_abs),
                                           len(event_data)), end="")
        print()
    print("%-15s" % ("Average (CA)"), end="")
    for per_data in all_data_ca:
        per_data_abs = [abs(x) for x in per_data]
        print("%6.2f %6.2f %6d" % (numpy.mean(per_data),
                                   numpy.mean(per_data_abs),
                                   len(per_data)), end="")
    print()
    print("%-15s" % ("Average (CENA)"), end="")
    for per_data in all_data_cena:
        per_data_abs = [abs(x) for x in per_data]
        print("%6.2f %6.2f %6d" % (numpy.mean(per_data),
                                   numpy.mean(per_data_abs),
                                   len(per_data)), end="")
    print()
    print("%-15s" % ("Average (All)"), end="")
    for per_data in all_data:
        per_data_abs = [abs(x) for x in per_data]
        print("%6.2f %6.2f %6d" % (numpy.mean(per_data),
                                   numpy.mean(per_data_abs),
                                   len(per_data)), end="")
    print()
    print()

def main():
    """
    Get input file from the command-line
    """
    if len(sys.argv) < 2:
        print("Usage: %s input_file" % (sys.argv[0]))
        sys.exit(1)

    input_file = sys.argv[1]

    data = compile_input_data(input_file)
    write_output_file(data)

if __name__ == "__main__":
    main()
