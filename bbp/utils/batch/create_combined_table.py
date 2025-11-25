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
import sys

def write_combined_table(input_files):
    """
    This function reads the input files and writes a combined output
    """

    # Create list of methods
    methods = " ".join([token.split("-")[0] for token in input_files])

    # Open all input files
    input_fps = [open(input_file) for input_file in input_files]

    # Process files
    while True:
        lines = [input_fp.readline() for input_fp in input_fps]
        # Check if we got to the end
        if not lines[0]:
            break
        lines = [line.strip() for line in lines]
        # Skip blank lines
        if not lines[0]:
            print()
            continue
        # If this is a header line
        if lines[0].startswith("Rrup") or lines[0].startswith("Mechanism"):
            print("%s - (%s)" % (lines[0], methods))
            continue
        # Ok, this is a data line
        tokens = [line.split() for line in lines]
        # First token in the event name, print it
        if tokens[0][0] == "Average":
            head = "%s %s" % (tokens[0][0], tokens[0][1])
            # The average lines headers have 2 tokens, remove the extra one
            for idx, _ in enumerate(tokens):
                del tokens[idx][0]
        else:
            head = tokens[0][0]
        # Remove line header
        for idx, _ in enumerate(tokens):
            del tokens[idx][0]
        print("%-15s" % (head), end="")
        # Figure out how many columns we have
        columns = len(tokens[0]) // 3

        for idx in range(0, columns):
            for method in tokens:
                token0 = method[idx * 3 + 0]
                token1 = method[idx * 3 + 1]
                try:
                    token0 = float(token0)
                    token1 = float(token1)
                except:
                    token0 = None
                    token1 = None
                if token0 is None:
                    print("%6s %6s" % ("  N/A", "  N/A"), end="")
                else:
                    print("%6.2f %6.2f" % (token0, token1), end="")
            print(" | ", end="")
        print()

    # All done, close input files
    for input_file in input_fps:
        input_file.close()

def main():
    """
    Get input files from the user
    """
    if len(sys.argv) < 2:
        print("Usage: %s input_file1 [input_file2 input_file3 ...]" %
              (sys.argv[0]))
        sys.exit(1)

    input_files = sys.argv[1:]

    write_combined_table(input_files)

if __name__ == "__main__":
    main()
