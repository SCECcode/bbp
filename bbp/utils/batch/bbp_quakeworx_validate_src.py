#!/usr/bin/env python3
"""
BSD 3-Clause License

Copyright (c) 2025, University of Southern California
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

Program to validate the BBP SRC file provided by the Quakeworx gateway
"""
from __future__ import division, print_function

# Works for both Python 2 and 3
try: input = raw_input
except NameError: pass

# Import Python modules
import os
import sys
import optparse

# Import Broadband modules
import bband_utils

SRC_MAX_FILE_SIZE_KB = 2

def main():
    """
    Parse command line options and create the needed files/directories
    """
    prog_base = os.path.basename(sys.argv[0])
    usage = "usage: %s [options]" % (prog_base)
    parser = optparse.OptionParser(usage)
    parser.add_option("-i", "--input-file", type="string", action="store",
                      dest="input_file",
                      help="Input file")
    parser.add_option("-o", "--output-file", type="string", action="store",
                      dest="output_file",
                      help="Output file")
    (options, _) = parser.parse_args()

    input_file = options.input_file
    output_file = options.output_file

    # Check file size
    try:
        file_size = os.path.getsize(input_file)
        # Convert from bytes to kb
        file_size = file_size / 1024
        if file_size > SRC_MAX_FILE_SIZE_KB:
            raise bband_utils.ParameterError("Input SRC file larger than %d Kb!" % (SRC_MAX_FILE_SIZE_KB))
    except OSError as e:
        raise bband_utils.ProcessingError("Cannot get size for file %s!" % (input_file))

    # Parse input file
    src_props = bband_utils.parse_properties(input_file)

    # Create output
    output = ""
    for key in src_props:
        output = output + "%s = %s\n" % (key.upper(), src_props[key])
    outfile = open(output_file, 'w')
    outfile.write(output)
    outfile.close()

if __name__ == "__main__":
    main()
