#!/usr/bin/env python
"""
Copyright 2010-2020 University Of Southern California

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

# Import BBP modules
import bband_utils
from install_cfg import InstallCfg
from station_list import StationList

# Import Pynga and its utilities
import pynga.utils as putils

def calculate_as16(m, r, mech, vs30, z1, cj):
    """
    This function implements the Afshari and Stewart (2016) GMPE for
    significant duration

    Input parameters:
    M: Magnitude
    R: Rupture distance (km)
    Mech: Rupture mechanism (0=unknown, 1=Normal, 2=Reverse, 3=Strike-slip)
    Vs30: Time-averaged shear wave velocity of the upper 30m of the site (m/s)
    z1: Basin depth, depth to shear wave velovity of 1000 m/s isosurface (m)
    Enter -999 if unknown
    CJ: Enter 0 for California, and 1 for Japan, Enter -999 otherwise

    Output: 5-75, 5-95, and 20-80 significant duration, 5-75, 5-95, and 20-80
    between-event standard deviation (tau), 5-75, 5-95, and 20-80 within-event
    standard deviation.
    """

    # Calculating differential basin depth
    if cj == 0:
        muz1 = math.exp(-7.15/4.0*math.log((vs30**4+570.94**4) /
                                           (1360.0**4+570.94**4)) -
                        math.log(1000))
    else:
        muz1 = math.exp(-5.23/2.0*math.log((vs30**2+412.39**2) /
                                           (1360.0**2+412.39**2)) -
                        math.log(1000))
    if z1 == -999:
        dz1 = 0.0
    else:
        dz1 = z1 - muz1

    if cj == -999:
        dz1 = 0.0

    # 5-75 significant duration
    # Parameters
    # Source
    m1 = 5.35
    m2 = 7.15
    b00 = 1.280
    b01 = 1.555
    b02 = 0.7806
    b03 = 1.279
    b10 = 5.576
    b11 = 4.992
    b12 = 7.061
    b13 = 5.578
    b2 = 0.9011
    b3 = -1.684
    mstar = 6

    # Path
    c1 = 0.1159
    rr1 = 10.0
    rr2 = 50.0
    c2 = 0.1065
    c3 = 0.0682

    # Site
    c4 = -0.2246
    vref = 368.2
    v1 = 600.0
    c5 = 0.0006
    dz1ref = 200.0

    if mech == 0:
        b1 = b10
        b0 = b00
    elif mech == 1:
        b1 = b11
        b0 = b01
    elif mech == 2:
        b1 = b12
        b0 = b02
    elif mech == 3:
        b1 = b13
        b0 = b03

    if m < m2:
        ds = math.exp(b1+b2*(m-mstar))
    else:
        ds = math.exp(b1+b2*(m2-mstar)+b3*(m-m2))

    # Stress drop term
    lnsd = ds
    lnsd = lnsd / (10**(1.5*m+16.05))
    lnsd = lnsd**(-1.0/3.0)
    source = lnsd / (3.2*4.9*1000000)

    if m < m1:
        source = b0

    # Path term
    if r < rr1:
        path = c1 * r
    elif r < rr2:
        path = c1 * rr1 + c2 * (r-rr1)
    else:
        path = c1 * rr1 + c2 * (rr2-rr1) + c3 * (r-rr2)

    # Site term
    if vs30 < v1:
        site = c4 * math.log(vs30 / vref)
    else:
        site = c4 * math.log(v1 / vref)

    if dz1 <= dz1ref:
        site = site + c5 * dz1
    else:
        site = site + c5 * dz1ref

    lnsd = math.log(source + path) + site
    sd575 = math.exp(lnsd)

    # 5-95 significant duration
    # Parameters
    # Source
    m1 = 5.2
    m2 = 7.4
    b00 = 2.182
    b01 = 2.541
    b02 = 1.612
    b03 = 2.302
    b10 = 3.628
    b11 = 3.170
    b12 = 4.536
    b13 = 3.467
    b2 = 0.9443
    b3 = -3.911
    mstar = 6
    # Path
    c1 = 0.3165
    # ee1 = 10.0
    rr1 = 10.0
    rr2 = 50.0
    c2 = 0.2539
    c3 = 0.0932
    # Site
    c4 = -0.3183
    vref = 369.9
    v1 = 600.0
    c5 = 0.0006
    dz1ref = 200.0

    if mech == 0:
        b1 = b10
        b0 = b00
    elif mech == 1:
        b1 = b11
        b0 = b01
    elif mech == 2:
        b1 = b12
        b0 = b02
    elif mech == 3:
        b1 = b13
        b0 = b03

    if m < m2:
        ds = math.exp(b1+b2*(m-mstar))
    else:
        ds = math.exp(b1+b2*(m2-mstar)+b3*(m-m2))

    # Stress drop term
    lnsd = ds
    lnsd = lnsd / (10**(1.5*m+16.05))
    lnsd = lnsd**(-1.0/3.0)
    source = lnsd / (3.2*4.9*1000000)

    if m < m1:
        source = b0

    if r < rr1:
        path = c1 * r
    elif r < rr2:
        path = c1 * rr1 + c2 * (r-rr1)
    else:
        path = c1 * rr1 + c2 * (rr2-rr1) + c3 * (r-rr2)

    # Path term
    if vs30 < v1:
        site = c4 * math.log(vs30 / vref)
    else:
        site = c4 * math.log(v1 / vref)

    if dz1 <= dz1ref:
        site = site + c5 * dz1
    else:
        site = site + c5 * dz1ref

    # Site term
    lnsd = math.log(source + path) + site
    sd595 = math.exp(lnsd)

    # 20-80 significant duration
    # Parameters
    # Source
    m1 = 5.2
    m2 = 7.4
    b00 = 0.8822
    b01 = 1.409
    b02 = 0.7729
    b03 = 0.8804
    b10 = 6.182
    b11 = 4.778
    b12 = 6.579
    b13 = 6.188
    b2 = 0.7414
    b3 = -3.164
    mstar = 6.0
    # Path
    c1 = 0.0646
    rr1 = 10.0
    rr2 = 50.0
    c2 = 0.0865
    c3 = 0.0373
    # Site
    c4 = -0.4237
    vref = 369.6
    v1 = 600.0
    c5 = 0.0005
    dz1ref = 200.0

    if mech == 0:
        b1 = b10
        b0 = b00
    elif mech == 1:
        b1 = b11
        b0 = b01
    elif mech == 2:
        b1 = b12
        b0 = b02
    elif mech == 3:
        b1 = b13
        b0 = b03

    if m < m2:
        ds = math.exp(b1+b2*(m-mstar))
    else:
        ds = math.exp(b1+b2*(m2-mstar)+b3*(m-m2))

    # Stress drop term
    lnsd = ds
    lnsd = lnsd / (10**(1.5*m+16.05))
    lnsd = lnsd**(-1.0/3.0)
    source = lnsd / (3.2*4.9*1000000)

    if m < m1:
        source = b0

    if r < rr1:
        path = c1 * r
    elif r < rr2:
        path = c1 * rr1 + c2 * (r-rr1)
    else:
        path = c1 * rr1 + c2 * (rr2-rr1) + c3 * (r-rr2)

    # Path term
    if vs30 < v1:
        site = c4 * math.log(vs30 / vref)
    else:
        site = c4 * math.log(v1 / vref)

    if dz1 <= dz1ref:
        site = site + c5 * dz1
    else:
        site = site + c5 * dz1ref

    # Site term
    lnsd = math.log(source + path) + site
    sd2080 = math.exp(lnsd)

    # within-site standard deviation
    # Parameters
    phi1575 = 0.54
    phi2575 = 0.41
    phi1595 = 0.43
    phi2595 = 0.35
    phi12080 = 0.56
    phi22080 = 0.45

    if m < 5.5:
        phi575 = phi1575
        phi595 = phi1595
        phi2080 = phi12080
    elif m < 5.75:
        phi575 = phi1575 + (phi2575-phi1575)*(m-5.5)/(5.75-5.5)
        phi595 = phi1595 + (phi2595-phi1595)*(m-5.5)/(5.75-5.5)
        phi2080 = phi12080 + (phi22080-phi12080)*(m-5.5)/(5.75-5.5)
    else:
        phi575 = phi2575
        phi595 = phi2595
        phi2080 = phi22080

    # Between-site standard deviation
    # Parameters
    tau1575 = 0.28
    tau2575 = 0.25
    tau1595 = 0.25
    tau2595 = 0.19
    tau12080 = 0.3
    tau22080 = 0.19

    if m < 6.5:
        tau575 = tau1575
        tau595 = tau1595
        tau2080 = tau12080
    elif m < 7.0:
        tau575 = tau1575 + (tau2575-tau1575)*(m-6.5)/(7.0-6.5)
        tau595 = tau1595 + (tau2595-tau1595)*(m-6.5)/(7.0-6.5)
        tau2080 = tau12080 + (tau22080-tau12080)*(m-6.5)/(7.0-6.5)
    else:
        tau575 = tau2575
        tau595 = tau2595
        tau2080 = tau22080

    return (sd575, sd595, sd2080,
            tau575, tau595, tau2080,
            phi575, phi595, phi2080)

class AS16(object):
    """
    This class implements the AS16 GMPE
    """

    def __init__(self, r_stations, r_srcfile, eventname, sim_id=0):
        """
        Initialize class variables
        """
        self.eventname = eventname
        self.stations = r_stations
        self.srcfile = r_srcfile
        self.log = None
        self.sim_id = sim_id

    def run(self):
        """
        Run the AS16 validation for all stations
        """
        print("AS2016".center(80, '-'))

        # Load configuration, set sim_id
        install = InstallCfg.getInstance()
        sim_id = self.sim_id

        # Build directory paths
        a_tmpdir = os.path.join(install.A_TMP_DATA_DIR, str(sim_id))
        a_indir = os.path.join(install.A_IN_DATA_DIR, str(sim_id))
        a_outdir = os.path.join(install.A_OUT_DATA_DIR, str(sim_id))
        a_logdir = os.path.join(install.A_OUT_LOG_DIR, str(sim_id))
        a_validation_outdir = os.path.join(a_outdir, "validations",
                                           "stewart_duration_gmpe")

        # Make sure the output and tmp directories exist
        bband_utils.mkdirs([a_tmpdir, a_indir, a_outdir, a_validation_outdir],
                           print_cmd=False)

        # Now the file paths
        self.log = os.path.join(a_logdir, "%d.as16.log" % (sim_id))
        sta_file = os.path.join(a_indir, self.stations)
        a_srcfile = os.path.join(a_indir, self.srcfile)

        # Read SRC file
        src_keys = bband_utils.parse_src_file(a_srcfile)

        # Load information from SRC file
        origin = (src_keys['lon_top_center'], src_keys['lat_top_center'])
        dims = (src_keys['fault_length'], src_keys['dlen'],
                src_keys['fault_width'], src_keys['dwid'],
                src_keys['depth_to_top'])
        mech = (src_keys['strike'], src_keys['dip'],
                src_keys['rake'])

        # Set region to be unknown -- this has no effect in the AS16
        # method as z1 is not provided and that causes dz1 to be set
        # to zero and override the cj parameter
        cj = -999

        # Figure out what mechanism to use
        # 0 = unknown
        # 1 = normal
        # 2 = reverse
        # 3 = strike-slip
        rake = src_keys['rake']
        if abs(rake) <= 30 or abs(rake) >= 150:
            mechanism = 3
        elif rake > 30 and rake < 150:
            mechanism = 2
        elif rake < -30 and rake > -150:
            mechanism = 1
        else:
            print("Warning: unknown mechanism for rake = %f" % (rake))
            mechanism = 0

        # Get station list
        slo = StationList(sta_file)
        site_list = slo.getStationList()

        # Create output file, add header
        out_file = open(os.path.join(a_validation_outdir,
                                     '%d.as16.%s.txt' %
                                     (self.sim_id, self.eventname)), 'w')
        out_file.write("#station, rrup, vs30, sd575, sd595, sd2080,"
                       " tau575, tau595, tau2080, phi575, phi595, phi2080\n")

        # Go through each station
        for site in site_list:
            stat = site.scode
            vs30 = float(site.vs30)

            # Calculate Rrup
            site_geom = [site.lon, site.lat, 0.0]
            (fault_trace1, up_seis_depth,
             low_seis_depth, ave_dip,
             dummy1, dummy2) = putils.FaultTraceGen(origin, dims, mech)
            _, rrup, _ = putils.DistanceToSimpleFaultSurface(site_geom,
                                                             fault_trace1,
                                                             up_seis_depth,
                                                             low_seis_depth,
                                                             ave_dip)

            results = calculate_as16(src_keys['magnitude'], rrup,
                                     mechanism, vs30, -999.0, cj)

            out_file.write("%s, %3.2f, %3.2f" % (stat, rrup, vs30))
            for piece in results:
                out_file.write(", %7.5f" % (piece))
            out_file.write("\n")

        # All done, close output file
        out_file.close()

        # All done!
        print("AS2016 Completed".center(80, '-'))


if __name__ == "__main__":
    print("Testing Module: %s" % (os.path.basename(sys.argv[0])))
    if len(sys.argv) == 7:
        print(calculate_as16(float(sys.argv[1]), float(sys.argv[2]),
                             float(sys.argv[3]), float(sys.argv[4]),
                             float(sys.argv[5]), float(sys.argv[6])))
    elif len(sys.argv) == 5:
        AS16 = AS16(sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]))
        AS16.run()
    else:
        print("Usage: %s " % (sys.argv[0]) +
              "mw rrup mech vs30 z1 region")
        print("=== OR ===")
        print("Usage: %s " % (sys.argv[0]) +
              "station_list source_file event_name sim_id")
        sys.exit(1)
