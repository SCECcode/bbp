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

# Import Python modules
import os
import sys
import glob
import numpy

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
        ("SS", ["Tottori", "Landers", "ALUMR", "Parkfield"]),
        ("NM", [])]

# Event list
EVENTS = ["CHINOH", "ALUMR", "WHITTIER", "Parkfield",
          "NORTHPS", "Tottori", "SanSimeon", "Niigata",
          "Chuetsu", "NR", "Iwate", "LOMAP", "Landers",
          "RDL1K", "Mineral", "Saguenay1k"]

EVENTS_CA = ["WHITTIER", "NORTHPS", "NR", "LOMAP", "Landers",
             "CHINOH", "ALUMR", "SanSimeon", "Parkfield"]

EVENTS_CENA = ["Mineral", "RDL1K", "Saguenay1k"]

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
        dist = float(items[7])
        tmin = float(items[10])
        tmax = float(items[11])
        if event not in EVENTS:
            print "Unknown event %s, skipping..." % event
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
        print "Rrup = %.2f-%.2f km" % (vals[0], vals[1])
        # Initialize data for calculating mean
        all_data = [[] for _ in PERIODS]
        all_data_ca = [[] for _ in PERIODS]
        all_data_cena = [[] for _ in PERIODS]
        for event in EVENTS:
            print "%-15s" % (event),
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
                    print("%6s %6s %6s" % ("  N/A", "  N/A", "  N/A")),
                else:
                    print("%6.2f %6.2f %6d" % (numpy.mean(event_data),
                                               numpy.mean(event_data_abs),
                                               len(event_data))),
            print ""
        print "%-15s" % "Average (CA)",
        for per_data in all_data_ca:
            per_data_abs = [abs(x) for x in per_data]
            if not len(per_data):
                print("%6s %6s %6s" % ("  N/A", "  N/A", "  N/A")),
            else:
                print("%6.2f %6.2f %6d" % (numpy.mean(per_data),
                                           numpy.mean(per_data_abs),
                                           len(per_data))),
        print ""
        print "%-15s" % "Average (CENA)",
        for per_data in all_data_cena:
            per_data_abs = [abs(x) for x in per_data]
            if not len(per_data):
                print("%6s %6s %6s" % ("  N/A", "  N/A", "  N/A")),
            else:
                print("%6.2f %6.2f %6d" % (numpy.mean(per_data),
                                           numpy.mean(per_data_abs),
                                           len(per_data))),
        print ""
        print "%-15s" % "Average (All)",
        for per_data in all_data:
            per_data_abs = [abs(x) for x in per_data]
            if not len(per_data):
                print("%6s %6s %6s" % ("  N/A", "  N/A", "  N/A")),
            else:
                print("%6.2f %6.2f %6d" % (numpy.mean(per_data),
                                           numpy.mean(per_data_abs),
                                           len(per_data))),
        print ""
        print ""

    print "Mechanism"
    all_data = [[] for _ in PERIODS]
    all_data_ca = [[] for _ in PERIODS]
    all_data_cena = [[] for _ in PERIODS]
    for mech in MECH:
        print "%-15s" % (mech[0]),
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
                print("%6s %6s %6s" % ("  N/A", "  N/A", "  N/A")),
            else:
                print("%6.2f %6.2f %6d" % (numpy.mean(event_data),
                                           numpy.mean(event_data_abs),
                                           len(event_data))),
        print ""
    print "%-15s" % "Average (CA)",
    for per_data in all_data_ca:
        per_data_abs = [abs(x) for x in per_data]
        print("%6.2f %6.2f %6d" % (numpy.mean(per_data),
                                   numpy.mean(per_data_abs),
                                   len(per_data))),
    print ""
    print "%-15s" % "Average (CENA)",
    for per_data in all_data_cena:
        per_data_abs = [abs(x) for x in per_data]
        print("%6.2f %6.2f %6d" % (numpy.mean(per_data),
                                   numpy.mean(per_data_abs),
                                   len(per_data))),
    print ""
    print "%-15s" % "Average (All)",
    for per_data in all_data:
        per_data_abs = [abs(x) for x in per_data]
        print("%6.2f %6.2f %6d" % (numpy.mean(per_data),
                                   numpy.mean(per_data_abs),
                                   len(per_data))),
    print ""
    print ""

def main():
    """
    Get input file from the command-line
    """
    if len(sys.argv) < 2:
        print ("Usage: %s input_file" %
               (sys.argv[0]))
        sys.exit(1)

    input_file = sys.argv[1]

    data = compile_input_data(input_file)
    write_output_file(data)

if __name__ == "__main__":
    main()
