#!/usr/bin/env python3
"""
BSD 3-Clause License

Copyright (c) 2023, University of Southern California
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

Class for managing station lists in BB Platform
"""
from __future__ import division, print_function

# Import Python modules
import sys

# Import Broadband modules
from station import Station

# Sets maximum allowed len for station name, code limits are:
# jbsim: 64 characters
# hfsims: 64 characters
# bbtoolbox: 128 characters
# b_green_99v8: 15 characters
# c_simula_v12: 15 characters
# syn1D: 256 characters
MAX_STATION_NAME_LEN = 15

class StationValidationError(Exception):
    """Custom Exception class which is raised when one or more station rows fail validation.
    It contains rich information about the error when parsing each station in the .stl file"""
    def __init__(self, errors):
        """ errors: list[(scode, field, raw_value, message)]
        scode = Station ID (station code)
        field = which field (column value) has a problem (e.g., "long", "lat", vs30", "z1pt0", etc.)
        raw_value = what input value (raw value) was found for that problematic field. 
        message = error/warning message """

        #Stores the whole list inside the exception object so it can later be inspected when the error is caught. 
        self.errors = errors 
        super().__init__(f"Station list validation failed with {len(errors)} error(s).")


class StationList(object):
    """
    Input Station List file and serve up stations infor as needed
    """
    def __init__(self, a_station_list=None, on_error="raise", LF_HF_check=True):
        """
        Pass in the absolute file name to the station list and it well be parsed and
        put into a dictionary accessible with iterators.
        
        Arguments:
            a_station_list (str): path to .stl file
            on_error (str): 'raise' (default), 'collect', or 'warn'
                - 'raise': raise StationValidationError after parsing if any errors
                - 'collect': keep errors in self.errors, do not raise
                - 'warn': print warnings, do not raise
            LF_HF_check (bool): if True, require LP_Freq <= HP_Freq
        """
        if a_station_list is None:
            raise ValueError("Error reading station list - Null Station List")

        self.a_station_filename = a_station_list
        self.on_error = on_error
        self.LF_HF_check = LF_HF_check
        self.errors = []          # list of (scode, field, raw_value, message)
        self.site_list = []       # Start with empty station list. parsed Station objects

        # Open file
        try:
            station_file = open(self.a_station_filename, "r")
        except OSError as e:
            raise OSError(f"Error opening station list file: {a_station_list}") from e

        # Helper to record an error consistently
        def add_err(scode, field, raw, msg):
            if self.on_error == "warn":
                print(f"WARNING [{scode}] {field}='{raw}': {msg}")
            self.errors.append((scode, field, raw, msg))

        # Read lines one by one
        for line in station_file:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#"):
                continue
            sta = line.split()

            # Expect at least 3 columns: Longitude, Latitude, StationId
            if len(sta) < 3:
                # No StationId available → treat as format error and skip row
                add_err("<missing-id>", "format", line, "Error: Expected at least 3 columns (Longitude, Latitude, & StationId)")
                continue

            scode = sta[2]

            # Validate StationId first (so that all later errors are tagged). 
            if not scode:
                add_err("<missing station-id>", "scode", scode, "Error: StationId is empty")
                continue
            if len(scode) > MAX_STATION_NAME_LEN:
                add_err(scode, "scode", scode, f"Error: Station name too long (>{MAX_STATION_NAME_LEN})")
                continue

            station = Station()
            station.scode = scode

            # Longitude & Latitude shoule be numeric and within bounds
            try:
                station.lon = float(sta[0])
                station.lat = float(sta[1])
            except ValueError:
                add_err(scode, "Coordinates", (sta[0], sta[1]), "Error: Longitude/Latitude must be numeric")
                continue

            # Validate values of latitude and longitudes.
            lon_ok = (-180.0 <= station.lon <= 180.0)
            lat_ok = (-90.0 <= station.lat <= 90.0)
            if not lon_ok:
                add_err(scode, "lon", sta[0], "Error: Longitude must be in [-180, 180]")
            if not lat_ok:
                add_err(scode, "lat", sta[1], "Error: Latitude must be in [-90, 90]")

            # Skip the row entirely if coords invalid
            if not (lon_ok and lat_ok):
                continue

            # 4th Column: Vs30 (invalid if <= 0 or non-numeric)
            if len(sta) >= 4:
                try:
                    vs30_val = float(sta[3])
                    if vs30_val <= 0:
                        station.vs30 = None
                        if vs30_val != -999:  # -999 is treated as a 'missing data' indicator which will not result into error or warning but sets the value to None. 
                            add_err(scode, "Vs30", sta[3], "Error: Vs30 must be > 0 or -999 for 'missing data' indicator. Setting it to None")
                    else:
                        station.vs30 = vs30_val                                 
                except ValueError:
                    add_err(scode, "Vs30", sta[3], "Error: Vs30 must be numeric > 0 or -999 for 'missing data' indicator. Setting it to None")
                    station.vs30 = None
            else:
                station.vs30 = None
                
             
            # Handle BBP v22.4.0 6-column format: Longitude Latitude Station_ID Vs30(m/s) LP_Freq(Hz)  HP_Freq(Hz) 
            if len(sta) == 6:
                # Lowpass corner frequency (LP_Freq(Hz))
                try:
                    lf = float(sta[4])
                    if lf <= 0:
                        station.low_freq_corner = 1.0e-15
                        if lf != -999:
                            add_err(scode, "low_freq_corner", sta[4], "Warning: LP_Freq must be > 0 or -999 for 'missing data' indicator, using 1e-15")
                    else:
                        station.low_freq_corner = lf
                except ValueError:
                    add_err(scode, "low_freq_corner", sta[4], "Warning: LP_Freq must be numeric > 0 or -999 for 'missing data' indicator, using 1e-15")
                    station.low_freq_corner = 1.0e-15

                # Highpass corner frequency(HP_Freq(Hz))
                try:
                    hf = float(sta[5])
                    if hf <= 0:
                        station.high_freq_corner = 1.0e+15
                        if hf != -999:
                            add_err(scode, "high_freq_corner", sta[5], "Warning: HP_Freq must be > 0 or -999 for 'missing data' indicator, using 1e+15")
                        
                    else:
                        station.high_freq_corner = hf
                except ValueError:
                    add_err(scode, "high_freq_corner", sta[5], "Warning: HP_Freq must be numeric > 0 or -999 for 'missing data' indicator, using 1e+15")
                    station.high_freq_corner = 1.0e+15

            # Extended .stl file column list: Longitude Latitude Station_ID Vs30(m/s) LP_Freq(Hz)  HP_Freq(Hz) z1pt0  z2pt5  basin_id  basin_label 
            if len(sta) > 6:
                try:
                    z1_val = float(sta[6])
                    if z1_val < 0:
                        station.z1pt0 = None
                        if z1_val != -999:
                            add_err(scode, "z1pt0", sta[6], "Error: z1pt0 must be >= 0 or -999 for 'missing data' indicator, setting it to None")
                    else:
                        station.z1pt0 = z1_val
                    
                except (TypeError, ValueError):
                    add_err(scode, "z1pt0", sta[6], "Error: z1pt0 must be numeric >= 0 or -999 for 'missing data' indicator, setting it to None")
                    station.z1pt0 = None

            if len(sta) > 7:
                try:
                    z2pt5_val = float(sta[7])
                    if z2pt5_val < 0:
                        station.z2pt5 = None
                        if z2pt5_val != -999:
                            add_err(scode, "z2pt5", sta[7], "Error: z2pt5 must be >= 0 or -999 for 'missing data' indicator, setting it to None")
                    else:
                        station.z2pt5 = z2pt5_val
                
                except (TypeError, ValueError):
                    add_err(scode, "z2pt5", sta[7], "Error: z2pt5 must be numeric >= 0 or -999 for 'missing data' indicator, setting it to None")
                    station.z2pt5 = None

            if len(sta) > 8:
                try:                  
                    basin_id_val = int(sta[8])
                    if basin_id_val not in [0,1,2,3]:
                        add_err(scode, "basin_id", sta[8], "Error: basin_id for SoCal region should be either 0, 1, 2, or 3. 0 for Mountain/Hills, 1 for Valleys, \
                        2 for Basin Edges/Transitional Zones, and 3 for Basins. Setting it to None")
                        station.basin_id = None
                    else:
                        station.basin_id = basin_id_val

                except (TypeError, ValueError):
                    add_err(scode, "basin_id", sta[8], "Error: basin_id for SoCal region should be either 0, 1, 2, or 3. 0 for Mountain/Hills, 1 for Valleys, \
                        2 for Basin Edges/Transitional Zones, and 3 for Basins. Setting it to None")
                    station.basin_id = None

            if len(sta) > 9:
                try: 
                    basin_label_value = str(sta[9])
                    if basin_label_value == "-999":
                        station.basin_label = None
                    elif basin_label_value not in ["LAB", "SFB", "SGB", "CB", "SBB", "CVB", "IVB"]:
                        add_err(scode, "basin_label", sta[9], "Error: basin_label for SoCal should be either 'LAB', 'SFB', 'SGB', 'CB', 'SBB', 'CVB', 'IVB' or '-999' for 'missing data' indicator. Setting it to None")
                        basin_label_value = None
                    else: 
                        station.basin_label = basin_label_value

                except (TypeError, ValueError):
                    add_err(scode, "basin_label", sta[9], "Error: basin_label for SoCal should be either 'LAB', 'SFB', 'SGB', 'CB', 'SBB', 'CVB', 'IVB' or '-999' for 'missing data' indicator. Setting it to None")
                    station.basin_label = None
            
            # Enforce LP <= HP
            if self.LF_HF_check and getattr(station, "low_freq_corner", None) is not None and \
               getattr(station, "high_freq_corner", None) is not None and \
               station.low_freq_corner > station.high_freq_corner:
                add_err(scode, "freq_order",
                        f"{station.low_freq_corner}>{station.high_freq_corner}",
                        "LP_Freq must be <= HP_Freq")

            self.site_list.append(station)

        # Remember to close the file
        try:
            station_file.close()
        except OSError:
            pass 

        # Error message if we weren't able to read any stations
        if len(self.site_list) == 0:
            raise ValueError(f"No stations read from station file: {a_station_list}")

        # Finalize error/warning logic/policy
        if self.errors:
            if self.on_error == "raise":
                raise StationValidationError(self.errors)
            elif self.on_error == "warn":
                # already printed in add_err
                pass

    @staticmethod
    def build(stat_list, output_file):
        """
        Writes output_file containing a list of our stations
        """
        with open( output_file, 'w') as fp:
            for stat in stat_list:
                has_vs30 = getattr(stat, "vs30", None) is not None
                has_lf   = getattr(stat, "low_freq_corner", None) is not None
                has_hf   = getattr(stat, "high_freq_corner", None) is not None
                has_z1   = getattr(stat, "z1pt0", None) is not None
                has_z2pt5  = getattr(stat, "z2pt5", None) is not None
                has_basin_id = getattr(stat, "basin_id", None) is not None
                has_basin_label = getattr(stat, "basin_label", None) is not None

                ## 10 columns
                if has_vs30 and has_lf and has_hf and has_z1 and has_z2pt5 and has_basin_id and has_basin_label:
                    fp.write("%.6f\t%.6f\t%s\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%d\t%s\n" %
                            (stat.lon, stat.lat, stat.scode, stat.vs30,
                            stat.low_freq_corner, stat.high_freq_corner, stat.z1pt0, stat.z2pt5, stat.basin_id, stat.basin_label))

                ## 9 columns 
                elif has_vs30 and has_lf and has_hf and has_z1 and has_z2pt5 and has_basin_id:
                    fp.write("%.6f\t%.6f\t%s\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%d\n" %
                            (stat.lon, stat.lat, stat.scode, stat.vs30,
                            stat.low_freq_corner, stat.high_freq_corner, stat.z1pt0, stat.z2pt5, stat.basin_id))

                ## 8 columns 
                elif has_vs30 and has_lf and has_hf and has_z1 and has_z2pt5:
                    fp.write("%.6f\t%.6f\t%s\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n" %
                            (stat.lon, stat.lat, stat.scode, stat.vs30,
                            stat.low_freq_corner, stat.high_freq_corner, stat.z1pt0, stat.z2pt5))

                ## 7 columns 
                elif has_vs30 and has_lf and has_hf and has_z1:
                    fp.write("%.6f\t%.6f\t%s\t%.6f\t%.6f\t%.6f\t%.6f\n" %
                            (stat.lon, stat.lat, stat.scode, stat.vs30,
                            stat.low_freq_corner, stat.high_freq_corner, stat.z1pt0))

                ## 6 columns (BBP v22.4.0 format)
                elif has_vs30 and has_lf and has_hf:
                    fp.write("%.6f\t%.6f\t%s\t%.6f\t%.6f\t%.6f\n" %
                            (stat.lon, stat.lat, stat.scode, stat.vs30,
                            stat.low_freq_corner, stat.high_freq_corner))

                ## 4 columns
                elif has_vs30:
                    fp.write("%.6f\t%.6f\t%s\t%.6f\n" %
                            (stat.lon, stat.lat, stat.scode, stat.vs30))
                ## 3 columns 
                else:
                    fp.write("%.6f\t%.6f\t%s\n" % (stat.lon, stat.lat, stat.scode))
                    
            fp.flush()
            fp.close()

    def get_station_list(self):
        """
        Returns the station list
        """
        return self.site_list

    def find_station(self, station_name):
        """
        Returns station object for station matching station_name, otherwise returns None. 
        """
        matching_list = [station for station in self.site_list if station.scode == station_name]

        # Not found
        if not len(matching_list):
            return None

        # Alert user that are multiple matches
        if len(matching_list) > 1:
            print("[WARNING]: multiple stations match station name %s!" % (station_name))

        # Return match
        return matching_list[0]

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Validate input variables values in a .stl file")
    parser.add_argument("station_list", help="Path to .stl file")
    parser.add_argument("--on-error", choices=["raise", "warn", "collect"], default="raise", help="Error policy: raise (default), warn, or collect.")
    parser.add_argument("--no-LF-HF-check", dest="LF_HF_check", action="store_false", default=True, help="Default: Enable LP_Freq <= HP_Freq requirement when both are present.")
    args = parser.parse_args()

    try:
        sl = StationList(args.station_list, on_error=args.on_error, LF_HF_check=args.LF_HF_check)
        for st in sl.get_station_list():
            print(st.scode)

        # If caller asked to collect/warn, optionally display a summary footer
        if sl.errors and args.on_error in ("warn", "collect"):
            print("\nValidation summary:")
            for sc, fld, val, msg in sl.errors:
                print(f"  - {sc}: {fld}={val} -> {msg}")

    except StationValidationError as ex:
        # Wiriting a clean and grouped report with exit
        print("\nError: invalid station list:")
        for sc, fld, val, msg in ex.errors:
            print(f"  - {sc}: {fld}={val} -> {msg}")
        sys.exit(1)
    except Exception as ex:
        print(f"Error: {ex}")
        sys.exit(1)