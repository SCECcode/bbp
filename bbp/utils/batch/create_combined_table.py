#!/usr/bin/env python
"""
Copyright 2010-2017 University Of Southern California

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
            print
            continue
        # If this is a header line
        if lines[0].startswith("Rrup") or lines[0].startswith("Mechanism"):
            print "%s - (%s)" % (lines[0], methods)
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
        print "%-15s" % (head),
        # Figure out how many columns we have
        columns = len(tokens[0]) / 3

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
                    print("%6s %6s" % ("  N/A", "  N/A")),
                else:
                    print("%6.2f %6.2f" % (token0, token1)),
            print " | ",
        print

    # All done, close input files
    for input_file in input_fps:
        input_file.close()

def main():
    """
    Get input files from the user
    """
    if len(sys.argv) < 2:
        print ("Usage: %s input_file1 [input_file2 input_file3 ...]" %
               (sys.argv[0]))
        sys.exit(1)

    input_files = sys.argv[1:]

    write_combined_table(input_files)

if __name__ == "__main__":
    main()
