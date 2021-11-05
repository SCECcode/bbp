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
import scipy.stats as st
import matplotlib as mpl
if mpl.get_backend() != 'agg':
    mpl.use('Agg') # Disables use of Tk/X11
import pylab

def read_data(input_file):
    """
    This function reads the input file and loads the data into
    our data structures
    """
    rrup = None
    data = []
    ifile = open(input_file, 'r')
    for line in ifile:
        line = line.strip()
        # Skip empty lines
        if not line:
            continue
        # Skip comments
        if line.startswith("%") or line.startswith("#"):
            continue
        # Skip Average lines
        if line.startswith("Average"):
            continue
        if line.startswith("Mechanism"):
            # Done with this file!
            break
        if line.startswith("Rrup"):
            # Process Rrup line
            pieces = line.split()
            distances = pieces[2]
            pieces = [float(piece) for piece in distances.split("-")]
            rrup = np.mean(pieces)
            continue
        # Real data line, process it!
        pieces = line.split()[1:]
        pieces = [np.nan if piece == "N/A" else piece for piece in pieces]
        pieces = [float(piece) for piece in pieces]
        pieces.insert(0, rrup)
        data.append(pieces)
    ifile.close()

    # All done, return data array
    return data

def summarize_and_plot_data(data, method, output_file):
    """
    Summarized all data into the format we need for plotting
    """
    mean_data = {}
    bins = 4
    titles = ["0.01 to 0.1s",
              "0.1 to 1s",
              "1 to 3s",
              "> 3s"]
    locs = [[0,0], [0,1], [1,0], [1,1]]

    # Calculate mean_data
    start = 1
    step = 3

    # Create fig
    fig, axs = pylab.plt.subplots(2, 2)
    fig.set_size_inches(17, 8.5)
    fig.suptitle("Method: %s" % (method))
    fig.subplots_adjust(hspace=0.4)
    fig.subplots_adjust(left=0.05)
    fig.subplots_adjust(right=0.98)

    current = start
    for bin in range(0, bins):
        mean_data[bin] = {}
        mean_data[bin]['mean'] = np.array([piece[current] for piece in data])
        mean_data[bin]['n'] = np.array([piece[current+2] for piece in data])
        current = current + step

    # List of distances
    r = np.array([piece[0] for piece in data])

    # Process each bin
    for bin in range(0, bins):
        x = np.log(r[~np.isnan(mean_data[bin]['mean'])])
        y = mean_data[bin]['mean'][~np.isnan(mean_data[bin]['mean'])]
        ww = mean_data[bin]['n'][~np.isnan(mean_data[bin]['n'])]
        numdata = len(y)

        A = np.array([list(np.ones(len(x))), x])
        A = A.T
        W = np.diag(ww)
        b = np.linalg.lstsq(((A.T).dot(W)).dot(A),
                            ((A.T).dot(W)).dot(np.array(y).T))[0]
        intercept = b[0]
        slope = b[1]

        degfree = len(x) - 2
        e = y - (intercept + slope * x)
        var = np.sum(e * e) / degfree
        se_y = np.sqrt(var)
        sdev = np.sqrt(var)
        se_b = sdev / np.sqrt(np.sum((x - np.mean(x)) * (x - np.mean(x))))
        se_a = sdev * np.sqrt(1.0 / len(x) + np.mean(x) * np.mean(x) /
                              np.sum((x - np.mean(x)) * (x - np.mean(x))))

        xx = np.linspace(min(x), max(x),
                         num=(int(math.ceil((max(x) - min(x)) / 0.1))))
        yy = slope * xx + intercept

        # Calculate 95% confidence bounds
        t = st.t.ppf(1.0 - 0.05 / 2, degfree)
        b95 = se_b * t
        a95 = se_a * t
        ratio = abs(slope) / b95
        ratio_round = round(ratio * 100) / 100.0
        lower95 = yy - t * se_y * np.sqrt(1.0 /
                                          len(x) + ((xx - np.mean(x)) *
                                                    (xx - np.mean(x))) /
                                          np.sum((xx - np.mean(x)) *
                                                 (xx - np.mean(x))))
        upper95 = yy + t * se_y * np.sqrt(1.0 /
                                          len(x) + ((xx - np.mean(x)) *
                                                    (xx - np.mean(x))) /
                                          np.sum((xx - np.mean(x)) *
                                                 (xx - np.mean(x))))

        # Let's plot it
        p_x = locs[bin][0]
        p_y = locs[bin][1]
        subfig = axs[p_x][p_y]
        subfig.set_title("%s - Ratio: %.2f" % (titles[bin], ratio_round))
        subfig.plot(x, y, 'k+')
        subfig.plot(xx, yy, color='green', ls='-')
        subfig.plot(xx, lower95, 'r--', xx, upper95, 'r--')
        subfig.set_ylabel('ln(data/model)', size=10)
        subfig.set_xlabel('ln(distance(km))', size=10)
        subfig.set_xlim(0, 6)
        subfig.set_ylim(-1.5, 1.5)
        subfig.grid(True)
        subfig.minorticks_on()

    # All done, save plot!
    fig.savefig(output_file, format='png', transparent=False,
                dpi=300)

def main():
    """
    Main function
    """
    if len(sys.argv) != 2:
        print("Usage: %s input_file" % (sys.argv[0]))
        sys.exit(0)

    # Output filename
    input_file = sys.argv[1]
    output_file = "%s.png" % (os.path.splitext(input_file)[0])
    method = os.path.basename(input_file).split("-")[0].upper()

    # Read input file
    data = read_data(input_file)
    summarize_and_plot_data(data, method, output_file)

if __name__ == "__main__":
    main()
