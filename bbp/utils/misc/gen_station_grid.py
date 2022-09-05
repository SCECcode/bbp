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

Program to create a list of stations for use in the Broadband Platform.
"""
from __future__ import division, print_function

# Import Python modules
import sys

def main():
    """
    Get min and max latitude and longitude values from the user, also
    read the step to be used in generating the station list. This code
    can only be used to generate a maximum of 9999 stations.
    """
    if len(sys.argv) < 6 or len(sys.argv) > 7:
        print("Usage: %s lat_min lat_max lon_min lon_max step [start_st_id]" % sys.argv[0])
        sys.exit(0)

    lat_min = float(sys.argv[1])
    lat_max = float(sys.argv[2])
    lon_min = float(sys.argv[3])
    lon_max = float(sys.argv[4])
    step = float(sys.argv[5])
    if len(sys.argv) > 6:
        station_number = int(float(sys.argv[6]))
    else:
        station_number = 1

    cur_lat = lat_min
    cur_lon = lon_min

    while cur_lat <= lat_max:
        while cur_lon <= lon_max:
            print("%2.3f    %2.3f    sta%04d   " % (cur_lon,
                                                    cur_lat,
                                                    station_number))
            station_number = station_number + 1
            cur_lon = cur_lon + step
        cur_lat = cur_lat + step
        cur_lon = lon_min

if __name__ == "__main__":
    main()
