#!/usr/bin/env python
"""
Southern California Earthquake Center Broadband Platform
Copyright 2010-2016 Southern California Earthquake Center

$Id: lf_seismograms.py 1722 2016-08-26 21:50:46Z fsilva $
"""
from __future__ import division, print_function

# Import Python modules
import os

# Import Broadband modules
from install_cfg import InstallCfg
from station_list import StationList

class LFSeismograms(object):
    """
    This module copies pre-computed low-frequency seismograms to the
    tmpdir directory
    """

    def __init__(self, i_seis_dir, i_r_stations, sim_id=0):
        """
        Initialize class variables
        """
        self.seis_dir = i_seis_dir
        self.r_stations = i_r_stations
        self.sim_id = sim_id

    def run(self):
        """
        Goes through the station list and copy each low-frequency
        seismogram from the seis_dir to the simulation's tmpdir
        """
        install = InstallCfg.getInstance()
        sim_id = self.sim_id

        a_tmpdir = os.path.join(install.A_TMP_DATA_DIR, str(sim_id))
        a_indir = os.path.join(install.A_IN_DATA_DIR, str(sim_id))
        a_stations = os.path.join(a_indir, self.r_stations)

        print(self.seis_dir)

        slo = StationList(a_stations)
        stat_list = slo.getStationList()
        for stat in stat_list:
            # Look for bbp seismogram, copy in
            print("%s/%s-lf.bbp" % (self.seis_dir, stat.scode))
            if os.path.exists("%s/%s-lf.bbp" % (self.seis_dir, stat.scode)):
                print("Copying for site %s" % (stat.scode))
                # Need to eliminate negative times
                fp_in = open("%s/%s-lf.bbp" % (self.seis_dir, stat.scode), 'r')
                fp_out = open("%s/%d.%s-lf.bbp" %
                              (a_tmpdir, sim_id, stat.scode), 'w')
                for line in fp_in:
                    pieces = line.split()
                    try:
                        if pieces[0] == '#' or pieces[0] == '%':
                            fp_out.write(line)
                        elif float(pieces[0]) < -0.0001:
                            continue
                        elif float(pieces[0]) < 0.0001:
                            fp_out.write("0.0\t%s\t%s\t%s\n" % (pieces[1],
                                                                pieces[2],
                                                                pieces[3]))
                        else:
                            fp_out.write(line)
                    except ValueError:
                        fp_out.write(line)
                fp_in.close()
                fp_out.flush()
                fp_out.close()
            else:
                print("Could not find LF seismogram for station %s!" %
                      (stat.scode))
