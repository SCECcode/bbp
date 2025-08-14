#! /usr/bin/env python
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

from __future__ import print_function, division

import os
import sys
import unittest
from station_list import StationList

REF_DIR = os.path.join(os.path.dirname(__file__), "ref_data", "stl_inputs")

class TestStationListInputs(unittest.TestCase):

    def load_stations(self, filename):
        full_path = os.path.join(REF_DIR, filename)
        sl = StationList(full_path, on_error="collect", strict_freq_order=False)
        # Do not raise inside tests; just collect errors
        return sl.get_station_list()

    def check_common_fields(self, stations):
        self.assertEqual(len(stations), 3)
        self.assertEqual(stations[0].scode, "STN1")
        self.assertEqual(stations[1].scode, "STN2")
        self.assertEqual(stations[2].scode, "STN3")
        self.assertEqual(stations[2].lon, -125.27)
        self.assertEqual(stations[2].lat, 33.98)

    def test_3col(self):
        stations = self.load_stations("station_3col.stl")
        self.check_common_fields(stations)
        for sta in stations:
            #print(f"DEBUG: {sta.scode}, lf={getattr(sta, 'low_freq_corner', 'MISSING')}, hf={getattr(sta, 'high_freq_corner', 'MISSING')}")
            self.assertIsNone(getattr(sta, "vs30", None))
            self.assertIsNone(getattr(sta, "z1pt0", None))
            self.assertIsNone(getattr(sta, "low_freq_corner", None))
            self.assertIsNone(getattr(sta, "high_freq_corner", None))

    def test_4col(self):
        stations = self.load_stations("station_4col.stl")
        self.check_common_fields(stations)
        self.assertEqual(stations[0].vs30, 329)
        self.assertEqual(stations[1].vs30, None)
        self.assertEqual(stations[2].vs30, 199)
        for sta in stations:
            self.assertIsNone(sta.z1pt0)
            self.assertIsNone(sta.low_freq_corner)
            self.assertIsNone(sta.high_freq_corner)

    def test_5col(self):
        stations = self.load_stations("station_5col.stl")
        self.check_common_fields(stations)
        self.assertEqual(stations[0].vs30, 329)
        self.assertEqual(stations[1].vs30, None)    
        self.assertEqual(stations[2].vs30, 199)
        self.assertIsNone(stations[0].z1pt0, None)
        self.assertIsNone(stations[1].z1pt0, None)
        #self.assertEqual(stations[1].z1pt0, 650)
        self.assertEqual(stations[2].z1pt0, 975)

    def test_6col(self):
        stations = self.load_stations("station_6col.stl")
        self.check_common_fields(stations)
        # STN1
        self.assertEqual(stations[0].vs30, 329)
        self.assertIsNone(stations[0].z1pt0)
        self.assertIsNone(stations[0].low_freq_corner)
        self.assertIsNone(stations[0].high_freq_corner)
        # STN2
        self.assertIsNone(stations[1].vs30, None)
        self.assertIsNone(stations[0].z1pt0)
        self.assertEqual(stations[1].low_freq_corner, 15)
        self.assertIsNone(stations[0].high_freq_corner)
        # STN3
        self.assertEqual(stations[2].vs30, 199)
        self.assertIsNone(stations[0].z1pt0)
        self.assertEqual(stations[2].low_freq_corner, 21)
        self.assertEqual(stations[2].high_freq_corner, 0.05)

    def test_7col(self):
        stations = self.load_stations("station_7col.stl")
        self.check_common_fields(stations)
        # STN1
        self.assertEqual(stations[0].vs30, 329)
        self.assertIsNone(stations[0].z1pt0, None)
        self.assertIsNone(stations[0].low_freq_corner)
        self.assertIsNone(stations[0].high_freq_corner)
        # STN2
        self.assertIsNone(stations[1].vs30, None)
        self.assertEqual(stations[1].z1pt0, 650)
        self.assertEqual(stations[1].low_freq_corner, 15)
        self.assertIsNone(stations[1].high_freq_corner)
        # STN3
        self.assertEqual(stations[2].vs30, 199)
        self.assertEqual(stations[2].z1pt0, 975)
        self.assertEqual(stations[2].low_freq_corner, 21)
        self.assertEqual(stations[2].high_freq_corner, 0.05)

    def test_missing_z1(self):
        stations = self.load_stations("station_missing_z1.stl")
        self.check_common_fields(stations)
        self.assertEqual(stations[0].vs30, 329)
        self.assertIsNone(stations[0].z1pt0, None)
        self.assertIsNone(stations[0].low_freq_corner)
        self.assertIsNone(stations[0].high_freq_corner)
        self.assertIsNone(stations[1].vs30, None)
        self.assertIsNone(stations[1].z1pt0, None)
        self.assertEqual(stations[1].low_freq_corner, 15)
        self.assertIsNone(stations[1].high_freq_corner)

if __name__ == "__main__":
    SUITE = unittest.TestLoader().loadTestsFromTestCase(TestStationListInputs)
    RESULT = unittest.TextTestRunner(verbosity=2).run(SUITE)
    sys.exit(not RESULT.wasSuccessful())