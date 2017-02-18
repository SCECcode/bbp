#!/usr/bin/env python
"""
Southern California Earthquake Center Broadband Platform
Copyright 2010-2016 Southern California Earthquake Center

This module creates acceleration seismograms and copies both velocity
and acceleration files to the outdata directory
$Id: copy_seismograms.py 1768 2016-10-10 16:23:49Z fsilva $
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys
import shutil

# Import Broadband modules
import bband_utils
from install_cfg import InstallCfg
from station_list import StationList

class CopySeismograms(object):
    """
    Finalizes generation of velocity and acceleration seismograms
    """

    def __init__(self, i_r_stations, sim_id=0, hybrid=False):
        self.sim_id = sim_id
        self.r_stations = i_r_stations
        self.hybrid = hybrid
        self.log = None

    def run(self):
        """
        Go through the station list and create acceleration
        seismogram. Copy results to outdata directory
        """
        print("Copy Seismograms".center(80, '-'))

        install = InstallCfg.getInstance()
        sim_id = self.sim_id

        sta_base = os.path.basename(os.path.splitext(self.r_stations)[0])
        self.log = os.path.join(install.A_OUT_LOG_DIR,
                                str(sim_id),
                                "%d.copy_seis_%s.log" % (sim_id, sta_base))
        a_statfile = os.path.join(install.A_IN_DATA_DIR,
                                  str(sim_id),
                                  self.r_stations)

        a_tmpdir = os.path.join(install.A_TMP_DATA_DIR, str(sim_id))
        a_outdir = os.path.join(install.A_OUT_DATA_DIR, str(sim_id))

        # Make sure tmpdir, outdir exist
        dirs = [a_tmpdir, a_outdir]
        bband_utils.mkdirs(dirs, print_cmd=False)

        #
        # Read and parse the statioin list with this call
        #
        slo = StationList(a_statfile)
        site_list = slo.getStationList()

        for sits in site_list:
            site = sits.scode
            print("==> Processing station: %s" % (site))

            if self.hybrid:
                expected_file = "%d.%s.bbp" % (sim_id, site)
            else:
                expected_file = "%d.%s.vel.bbp" % (sim_id, site)

            #print("Processing velocity for station %s - %s" %
            #      (site, expected_file))
            bbpfile = os.path.join(a_tmpdir, expected_file)

            # Make sure velocity file is there, otherwise, skip this station
            if not os.path.exists(bbpfile):
                print("No velocity seismograms found for station %s" %
                      (site))
                print("Skipping this station...")
                continue

            # Copy velocity bbp file to outdir
            shutil.copy2(bbpfile,
                         os.path.join(a_outdir, "%d.%s.vel.bbp" %
                                      (sim_id, site)))

            # Create path names and check if their sizes are within bounds
            nsfile = os.path.join(a_tmpdir,
                                  "%d.%s.000" % (sim_id, site))
            ewfile = os.path.join(a_tmpdir,
                                  "%d.%s.090" % (sim_id, site))
            udfile = os.path.join(a_tmpdir,
                                  "%d.%s.ver" % (sim_id, site))

            bband_utils.check_path_lengths([nsfile, ewfile, udfile],
                                           bband_utils.GP_MAX_FILENAME)

            cmd = ("%s/wcc2bbp " % (install.A_GP_BIN_DIR) +
                   "nsfile=%s ewfile=%s udfile=%s " %
                   (nsfile, ewfile, udfile) +
                   "wcc2bbp=0 < %s >> %s 2>&1" %
                   (bbpfile, self.log))
            bband_utils.runprog(cmd, abort_on_error=True, print_cmd=False)

            for comp in ['000', '090', 'ver']:
                # Create path names and check if their sizes are
                # within bounds
                filein = os.path.join(a_tmpdir,
                                      "%d.%s.%s" %
                                      (sim_id, site, comp))
                fileout = os.path.join(a_tmpdir,
                                       "%d.%s.acc.%s" %
                                       (sim_id, site, comp))

                bband_utils.check_path_lengths([filein, fileout],
                                               bband_utils.GP_MAX_FILENAME)

                cmd = ("%s/integ_diff diff=1 filein=%s fileout=%s" %
                       (install.A_GP_BIN_DIR, filein, fileout))
                bband_utils.runprog(cmd, abort_on_error=True, print_cmd=False)

            # Create path names and check if their sizes are within bounds
            nsfile = os.path.join(a_tmpdir,
                                  "%d.%s.acc.000" % (sim_id, site))
            ewfile = os.path.join(a_tmpdir,
                                  "%d.%s.acc.090" % (sim_id, site))
            udfile = os.path.join(a_tmpdir,
                                  "%d.%s.acc.ver" % (sim_id, site))
            bbpfile = os.path.join(a_tmpdir,
                                   "%d.%s.acc.bbp" % (sim_id, site))

            bband_utils.check_path_lengths([nsfile, ewfile, udfile],
                                           bband_utils.GP_MAX_FILENAME)

            cmd = ("%s/wcc2bbp " % (install.A_GP_BIN_DIR) +
                   "nsfile=%s ewfile=%s udfile=%s " %
                   (nsfile, ewfile, udfile) +
                   "units=cm/s/s wcc2bbp=1 > %s 2>> %s" %
                   (bbpfile, self.log))
            bband_utils.runprog(cmd, abort_on_error=True, print_cmd=False)

            # Copy acceleration bbp file to outdir
            shutil.copy2(bbpfile,
                         os.path.join(a_outdir, "%d.%s.acc.bbp" %
                                      (sim_id, site)))

        print("Copy Seismograms Completed".center(80, '-'))

if __name__ == "__main__":
    print("Testing Module: %s" % (os.path.basename((sys.argv[0]))))
    ME = CopySeismograms(sys.argv[1], int(sys.argv[3]), bool(str(sys.argv[4])))
    ME.run()
