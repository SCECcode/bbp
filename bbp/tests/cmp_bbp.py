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

This module is used to compare certain files types
"""
from __future__ import division, print_function

# Import Python modules
import math
import sys

# Boolean flag for enforcing tolerance checks
#ENFORCE_TOLERANCE = False
ENFORCE_TOLERANCE = True

def cmp_ffsp(filename1, filename2, tolerance=0.01):
    """
    Compares two FFSP files
    """
    fp1 = open(filename1, 'r')
    fp2 = open(filename2, 'r')

    # Compare headers
    header1 = fp1.readline().strip()
    header2 = fp2.readline().strip()
    items1 = [float(item) for item in header1.strip().split()]
    items2 = [float(item) for item in header2.strip().split()]
    for i in range(len(items1)):
        if items1[i] == 0.0 and items2[i] == 0.0:
            continue
        try:
            if math.fabs((items1[i] - items2[i]) / items1[i]) <= tolerance:
                continue
        except:
            pass
        fp1.close()
        fp2.close()
        print("Header difference larger than tolerance: %f versus %f" %
              (items1[i], items2[i]))
        return 1

    with fp1, fp2:
        for line1, line2 in zip(fp1, fp2):
            items1 = [float(item) for item in line1.strip().split()]
            items2 = [float(item) for item in line2.strip().split()]
            if len(items1) != len(items2):
                fp1.close()
                fp2.close()
                print("Number of items mismatch: %d versus %d!" %
                      (len(items1), len(items2)))
                return 1
            for i in range(len(items1)):
                if items1[i] == 0.0 and items2[i] == 0.0:
                    continue
                try:
                    if math.fabs((items1[i] - items2[i]) /
                                 items1[i]) <= tolerance:
                        continue
                except:
                    pass
                fp1.close()
                fp2.close()
                print("Difference larger than tolerance: %f versus %f" %
                      (items1[i], items2[i]))
                return 1
    # Done
    fp1.close()
    fp2.close()

    # All good!
    return 0

def read_srf_line(srf_fh):
    """
    This function returns the next line from a SRF file,
    skipping lines with comments as needed
    """
    while True:
        line = srf_fh.readline()
        if len(line) == 0:
            # End of file
            break
        if not line.strip().startswith("#"):
            break
    return line

def cmp_srf(filename1, filename2, tolerance=0.0011):
    """
    Compare two SRF files.  A tolerance is accepted for floating-point
    values.  If any point parameters differ by more than epsilon
    but less than the tolerance, let it go, but don't compare the slip
    values for that point (since they may differ by quite a bit).
    Count up the number of times it happens;  if it's more than 1% of the
    number of points reject it.
    """
    return_code = 0

    fp1 = open(filename1, 'r')
    fp2 = open(filename2, 'r')

    file1_version = float(read_srf_line(fp1).strip().split()[0])
    file2_version = float(read_srf_line(fp2).strip().split()[0])

    if file1_version != file2_version:
        # SRF files have different version numbers
        print("SRF versions don't match!")
        fp1.close()
        fp2.close()
        return 1

    file1_planes = int(read_srf_line(fp1).strip().split()[1])
    file2_planes = int(read_srf_line(fp2).strip().split()[1])

    if file1_planes != file2_planes:
        # SRF files have different numbers of planes
        print("Number of planes mismatch for the two SRF files!")
        fp1.close()
        fp2.close()
        return 1

    # Compare data for each plane
    for plane in range(0, file1_planes):
        file1_plane_header = []
        file2_plane_header = []
        file1_plane_header.append(read_srf_line(fp1).strip().split())
        file1_plane_header.append(read_srf_line(fp1).strip().split())
        file2_plane_header.append(read_srf_line(fp2).strip().split())
        file2_plane_header.append(read_srf_line(fp2).strip().split())
        for plane1, plane2 in zip(file1_plane_header, file2_plane_header):
            for token1, token2 in zip(plane1, plane2):
                if token1 != token2:
                    # Token mismatch
                    print("Tokens in plane header mismatch! %s != %s" % (token1, token2))
                    fp1.close()
                    fp2.close()
                    return 1

    # Compare points
    NT_INDEX = 2
    NT_TOLERANCE = 1
    LL_TOLERANCE = 0.00011
    PP_TOLERANCE = 0.011
    RAKE_TOLERANCE = 1
    SKIP_TOLERANCE = 0.0001
    total_points = 0
    
    for plane in range(0, file1_planes):
        tokens1 = read_srf_line(fp1).strip().split()
        tokens2 = read_srf_line(fp2).strip().split()
        if tokens1[0] != tokens2[0] != "POINTS":
            # Token mismatch
            print("Cannot find matching POINTS line!")
            fp1.close()
            fp2.close()
            return 1
        if tokens1[1] != tokens2[1]:
            # Token mismatch
            print("Number of points mismatch!")
            fp1.close()
            fp2.close()
            return 1
        num_points = int(float(tokens1[1]))
        total_points = total_points + num_points
        points_skipped = 0

        for j in range(0, num_points):
            # If any of the params are different, skip the slip values
            skip_slips = False
            # each point has 2 lines + ceil(NT1/6) lines
            # will permit tolerance in lat/lon points
            # and tolerance in NT values
            pieces1 = read_srf_line(fp1).strip().split()
            pieces2 = read_srf_line(fp2).strip().split()
            print("Set 1f: ", pieces1)
            print("Set 2f: ", pieces2)

            # Check line length
            if not len(pieces1) == len(pieces2):
                print("Plane %d, point %d: files have different number of point parameters." % (plane, j))
                fp1.close()
                fp2.close()
                return 1
            # Check longitudes
            lon1 = float(pieces1[0])
            lon2 = float(pieces2[0])
            if math.fabs(lon1 - lon2) > LL_TOLERANCE:
                print("Plane %d, point %d: longitudes %f and %f differ by more than the accepted tolerance %f." %
                      (plane, j, lon1, lon2, LL_TOLERANCE))
                fp1.close()
                fp2.close()
                return 2
            # Check latitudes
            lat1 = float(pieces1[1])
            lat2 = float(pieces2[1])
            if math.fabs(lat1 - lat2) > LL_TOLERANCE:
                print("Plane %d, point %d: latitudes %f and %f differ by more than the accepted tolerance %f." %
                      (plane, j, lat1, lat2, LL_TOLERANCE))
                fp1.close()
                fp2.close()
                return 3
            # Check rest of line
            for k in range(2, len(pieces1)):
                if math.fabs(float(pieces1[k]) - float(pieces2[k])) > PP_TOLERANCE:
                    print("Plane %d, point %d: point parameters in field %d disagree." %
                          (plane, j, k + 1))
                    fp1.close()
                    fp2.close()
                    return 4
                if math.fabs(float(pieces1[k]) - float(pieces2[k])) > SKIP_TOLERANCE:
                    #don't compare the slip values, they'll be different
                    skip_slips = True

            # Second line
            pieces1 = read_srf_line(fp1).strip().split()
            pieces2 = read_srf_line(fp2).strip().split()

            print("Set 1s: ", pieces1)
            print("Set 2s: ", pieces2)
            
            # Check line length
            if not len(pieces1) == len(pieces2):
                print("Plane %d, point %d: files have different number of point parameters." % (plane, j))
                fp1.close()
                fp2.close()
                return 1
            # Check rake
            r1 = int(float(pieces1[0]))
            r2 = int(float(pieces2[0]))
            if r1 != r2:
                skip_slips = True
            if abs(r1 - r2) > RAKE_TOLERANCE:
                print("Plane %d, point %d: rakes %d and %d differ by more than the accepted tolerance of 1." %
                      (plane, j, r1, r2))
                fp1.close()
                fp2.close()
                return 4
            for k in range(1, NT_INDEX):
                if math.fabs(float(pieces1[k]) - float(pieces2[k])) > PP_TOLERANCE:
                    print("Plane %d, point %d: point parameters in field %d disagree." %
                          (plane, j, k+1))
                    fp1.close()
                    fp2.close()
                    return 4
                if math.fabs(float(pieces1[k]) - float(pieces2[k])) > SKIP_TOLERANCE:
                    # ok, but skip comparisons
                    skip_slips = True

            # Compare NTs
            nt1 = int(float(pieces1[NT_INDEX]))
            nt2 = int(float(pieces2[NT_INDEX]))

            if not nt1 == nt2:
                skip_slips = True
                if abs(nt1 - nt2) > NT_TOLERANCE:
                    print("Plane %d, point %d: NT values %d and %d differ by %f, more than the accepted tolerance %d." %
                          (plane, j, nt1, nt2, abs(nt1 - nt2), NT_TOLERANCE))
                    fp1.close()
                    fp2.close()
                    return 5

            if not int(float(pieces1[4])) == 0 or not int(float(pieces1[6])) == 0:
                print("SRF has NT2 or NT3, need to alter parser.")
                print("Set 1: ", pieces1)
                print("Set 2: ", pieces2)
                fp1.close()
                fp2.close()
                sys.exit(0)

            # How many slip lines we need to read
            num_rows = int(math.ceil(nt1 / 6.0))
            for k in range(0, num_rows):
                # Read slips
                pieces1 = read_srf_line(fp1).strip().split()
                pieces2 = read_srf_line(fp2).strip().split()

                if skip_slips == False:

                    if not len(pieces1) == len(pieces2):
                        print("Plane %d, point %d, line %d: mismatch in entries in line." %
                              (plane, j, k))
                        continue
                    if not ENFORCE_TOLERANCE:
                        continue
                    for p in range(0, len(pieces1)):
                        p1 = float(pieces1[p])
                        p2 = float(pieces2[p])
                        if p1 < 1.0 or p2 < 1.0:
                            if math.fabs(p1 - p2) > tolerance:
                                print("Plane %d, point %d, line %d: %f and %f differ by more than the accepted tolerance %f (%f)." %
                                      (plane, j, k, float(pieces1[p]), float(pieces2[p]),
                                       tolerance, math.fabs(p1 - p2)))
                                return_code = 1
                        else:
                            if math.fabs(p1-p2) / p1 > tolerance:
                                print("Plane %d, point %d, line %d: %f and %f differ by more than the accepted tolerance %f%% (%f%%)." %
                                      (plane, j, k, float(pieces1[p]), float(pieces2[p]),
                                       tolerance * 100.0,
                                       math.fabs(p1 - p2) / p1 * 100.0))
                                return_code = 1
                else:
                    points_skipped += 1

        if points_skipped > num_points / 50:
            print("Too many points with different parameters. "
                  "Of %d total points %d had different parameters." %
                  (num_points, points_skipped))
            return_code = 2

    #print("Total points compared %d!" % (total_points))

    return return_code

def cmp_resid(filename1, filename2, tolerance=0.0015):
    """
    Format is long list of heaters, then for each station
    <eq> <mag> <stat name> <lon> <lat> <stat_seq_no> <vs30> <close_dist> <Xcos> <Ycos> <T_min> <T_max> <comp> <period1> ...
    """
    fp1 = open(filename1, 'r')
    fp2 = open(filename2, 'r')
    data1 = fp1.readlines()
    data2 = fp2.readlines()
    fp1.close()
    fp2.close()

    returncode = 0

    i = 0

    line1 = data1[i]
    line2 = data2[i]
    if line1 != line2:
        print("Discrepancy in header lines.")
        returncode = -1
        return returncode
    for i in range(1, len(data1)):
        pieces1 = data1[i].split()
        pieces2 = data2[i].split()
        for j in range(0, 13):
            if pieces1[j] != pieces2[j]:
                print("Line %d: %s and %s don't agree." %
                      ((i + 1), pieces1[j], pieces2[j]))
                continue
        if not ENFORCE_TOLERANCE:
            continue
        for j in range(13, len(pieces1)):
            f1 = float(pieces1[j])
            f2 = float(pieces2[j])
            if math.fabs(f1) < 1.0:
                if math.fabs(f1 - f2) > tolerance:
                    print("Line %d: %f and %f differ by more than %f tolerance." %
                          ((i + 1), f1, f2, tolerance))
                    returncode = 1
            else:
                if math.fabs(f1 - f2) / f1 > tolerance:
                    print("Line %d:  %f and %f differ by more than %f%% tolerance." %
                          ((i + 1), f1, f2, tolerance * 100.0))
                    returncode = 1
    return returncode

def cmp_bbp(filename1, filename2, tolerance=0.0015):
    fp1 = open(filename1, 'r')
    fp2 = open(filename2, 'r')
    data1 = fp1.readlines()
    data2 = fp2.readlines()
    fp1.close()
    fp2.close()

    returncode = 0

    #need to skip comments in both files
    file1_offset = 0
    file2_offset = 0
    for line1 in data1:
        if line1.strip().startswith("#") or line1.strip().startswith("%"):
            file1_offset += 1
        else:
            break
    for line2 in data2:
        if line2.strip().startswith("#") or line2.strip().startswith("%"):
            file2_offset += 1
        else:
            break

    for i in range(0, len(data1) - file1_offset):
        line1 = data1[i + file1_offset]
        pieces1 = line1.split()
        line2 = data2[i + file2_offset]
        pieces2 = line2.split()
        if not ENFORCE_TOLERANCE:
            continue
        for j in range(0, 4):
            f1 = float(pieces1[j])
            f2 = float(pieces2[j])
            if math.fabs(f1) < 1.0 or math.fabs(f2) < 1.0:
                if math.fabs(f1 - f2) > tolerance:
                    if returncode == 0:
                        print("BBP File Comparison: %s, %s" %
                              (filename1, filename2))
                    print("Line %d: %f and %f differ by more than %f tolerance." %
                          ((i + 1), f1, f2, tolerance))
                    returncode = 1
            else:
                if math.fabs(f1 - f2) / f1 > tolerance:
                    if returncode == 0:
                        print("BBP File Comparison: %s, %s" %
                              (filename1, filename2))
                    print("Line %d: %f and %f differ by more than %f%% tolerance." %
                          ((i + 1), f1, f2, tolerance*100.0))
                    returncode = 1
        if i > 1000:
            return returncode
    return returncode

def cmp_fas(filename1, filename2, tolerance=0.0015):
    """
    Compare two fas output files
    """
    return_code = 0

    fp1 = open(filename1, 'r')
    fp2 = open(filename2, 'r')

    for line1, line2 in zip(fp1, fp2):
        line1 = line1.strip()
        line2 = line2.strip()
        if line1.startswith("#") and line2.startswith("#"):
            continue
        pieces1 = line1.split()
        pieces2 = line2.split()
        if not ENFORCE_TOLERANCE:
            continue
        for token1, token2 in zip(pieces1, pieces2):
            token1 = float(token1)
            token2 = float(token2)

            if math.fabs(token1) < 1.0 or math.fabs(token2) < 1.0:
                if math.fabs(token1 - token2) > tolerance:
                    if return_code == 0:
                        print("FAS file comparison: %s, %s" %
                              (filename1, filename2))
                    print("%f and %f differ by more than %f tolerance." %
                          (token1, token2, tolerance))
                    return_code = 1
            else:
                if math.fabs(token1 - token2) / token1 > tolerance:
                    if return_code == 0:
                        print("FAS file comparison: %s, %s" %
                              (filename1, filename2))
                    print("%f and %f differ by more than %f%% tolerance." %
                          (token1, token2, tolerance * 100.0))
                    return_code = 1

    # All done, close files
    fp1.close()
    fp2.close()

    return return_code

def cmp_bias(filename1, filename2, tolerance=0.0015):
    fp1 = open(filename1, 'r')
    fp2 = open(filename2, 'r')
    data1 = fp1.readlines()
    data2 = fp2.readlines()
    fp1.close()
    fp2.close()

    returncode = 0

    for i in range(0, len(data1)):
        pieces1 = data1[i].split()
        pieces2 = data2[i].split()
        f1 = float(pieces1[1])
        f2 = float(pieces2[1])

        if pieces1[0] != pieces2[0]:
            print("Line %d: periods %f and %f don't agree." %
                  ((i + 1), float(pieces1[0]), float(pieces2[0])))
            returncode = 1
        if not ENFORCE_TOLERANCE:
            continue
        if math.fabs(f1) < 1.0:
            if math.fabs(f1 - f2) > tolerance:
                print("Line %d: %f and %f differ by more than %f tolerance." %
                      ((i + 1), f1, f2, tolerance))
                returncode = 1
        else:
            if math.fabs(f1 - f2) / f1 > tolerance:
                print("Line %d: %f and %f differ by more than %f%% tolerance." %
                      ((i + 1), f1, f2, tolerance * 100.0))
                returncode = 1
    return returncode

def cmp_gof(filename1, filename2, col_start=0, col_end=1, tolerance=0.0015):
    fp1 = open(filename1, 'r')
    fp2 = open(filename2, 'r')
    data1 = fp1.readlines()
    data2 = fp2.readlines()
    fp1.close()
    fp2.close()

    returncode = 0

    #need to skip comments in both files
    file1_offset = 0
    file2_offset = 0
    for line1 in data1:
        if line1.strip().startswith("#") or line1.strip().startswith("%"):
            file1_offset += 1
        else:
            break
        for line2 in data2:
            if line2.strip().startswith("#") or line2.strip().startswith("%"):
                file2_offset += 1
            else:
                break

    for i in range(0, len(data1) - file1_offset):
        line1 = data1[i + file1_offset]
        pieces1 = line1.split()
        line2 = data2[i + file2_offset]
        pieces2 = line2.split()
        if not ENFORCE_TOLERANCE:
            continue
        for j in range(col_start, col_end):
            f1 = float(pieces1[j])
            f2 = float(pieces2[j])
            if math.fabs(f1) < 1.0 or math.fabs(f2) < 1.0:
                if math.fabs(f1 - f2) > tolerance:
                    print("Line %d: %f and %f differ by more than %f tolerance." %
                          ((i + 1), f1, f2, tolerance))
                    returncode = 1
            else:
                if math.fabs(f1-f2) / f1 > tolerance:
                    print("Line %d: %f and %f differ by more than %f%% tolerance." %
                          ((i + 1), f1, f2, tolerance * 100.0))
                    returncode = 1
        if i > 1000:
            return returncode
    return returncode

def cmp_files_generic(filename1, filename2, tolerance=0.0015,
                      start_col=0, sep=None):
    """
    This function compares tokens from two output files
    """
    # Start with zero return code
    return_code = 0

    with open(filename1, 'r') as file1, open(filename2, 'r') as file2:
        for line1, line2 in zip(file1, file2):
            # Skip comments
            if line1.startswith("#"):
                continue
            line1 = line1.strip()
            line2 = line2.strip()
            pieces1 = line1.split(sep)[start_col:]
            pieces2 = line2.split(sep)[start_col:]
            pieces1 = [float(piece) for piece in pieces1]
            pieces2 = [float(piece) for piece in pieces2]
            # Lines must have same number of tokens
            if len(pieces1) != len(pieces2):
                print("Line contains different number of tokens!")
                return_code = 1
                break
            for token1, token2 in zip(pieces1, pieces2):
                if math.fabs(token1 - token2) > tolerance:
                    print("Tokens are different!")
                    return_code = 2
                    break

    return return_code
