#!/usr/bin/env python
"""
Southern California Earthquake Center Broadband Platform
Copyright 2010-2016 Southern California Earthquake Center

Broadband Platform Version of Rob Graves jbsim script
Outputs velocity (cm/s)
$Id: jbsim.py 1722 2016-08-26 21:50:46Z fsilva $
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys

# Import Broadband modules
import bband_utils
import fault_utils
from jbsim_cfg import JbsimCfg
from install_cfg import InstallCfg
from station_list import StationList

class Jbsim(object):
    """
    Implement Rob Graves jbsim.csh as a python component
    """

    def __init__(self, i_r_velmodel, i_r_srcfile, i_r_srffile,
                 i_r_stations, i_vmodel_name, sim_id=0):
        self.sim_id = sim_id
        self.r_velmodel = i_r_velmodel
        self.r_srcfile = i_r_srcfile
        self.r_srffile = i_r_srffile
        self.r_stations = i_r_stations
        self.vmodel_name = i_vmodel_name

    def run(self):
        """
        Runs the GP low frequency component
        """
        print("GP Jbsim".center(80, '-'))

        install = InstallCfg.getInstance()
        config = JbsimCfg(self.vmodel_name)

        sim_id = self.sim_id
        sta_base = os.path.basename(os.path.splitext(self.r_stations)[0])
        self.log = os.path.join(install.A_OUT_LOG_DIR,
                                str(sim_id),
                                "%d.jbsim_%s.log" % (sim_id, sta_base))

        a_statfile = os.path.join(install.A_IN_DATA_DIR,
                                  str(sim_id),
                                  self.r_stations)
        a_srffile = os.path.join(install.A_IN_DATA_DIR,
                                 str(sim_id),
                                 self.r_srffile)

        # Set directories, and make sure they exist
        a_indir = os.path.join(install.A_IN_DATA_DIR, str(sim_id))
        a_veldir = os.path.join(install.A_TMP_DATA_DIR, str(sim_id))
        bband_utils.mkdirs([a_veldir, a_indir], print_cmd=False)

        #
        # Read and parse the station list with this call
        #
        slo = StationList(a_statfile)
        site_list = slo.getStationList()

        for sits in site_list:
            slon = float(sits.lon)
            slat = float(sits.lat)
            site = sits.scode

            print("==> Generating LF seismogram for station: %s" % (site))

            #
            # We have a verbose of silent invocation. This is a very
            # verbose program so our default writes to dev/null
            #
            progstring = ("%s latloncoords=1 slon=%f slat=%f " %
                          (os.path.join(install.A_GP_BIN_DIR,
                                        "jbsim"), slon, slat) +
                          "rupmodtype=SRF rupmodfile=%s " % a_srffile +
                          "moment=-1 outdir=%s stat=%s " % (a_veldir, site) +
                          "min_taper_range=0.0 max_taper_range=0.0 " +
                          "gftype=fk gflocs=%s gftimes=%s gf_swap_bytes=%d " %
                          (config.GF_LOCS, config.GF_TIMES,
                           config.GF_SWAP_BYTES) +
                          "gfpath=%s gfname=%s maxnt=%d mindt=%f " %
                          (config.A_GP_GF_DIR, config.GF_NAME,
                           config.MAX_GFNT, config.MIN_GFDT) +
                          "dtout=%f ntout=%d >> %s 2>&1" %
                          (config.DTOUT, config.NTOUT, self.log))
            bband_utils.runprog(progstring, print_cmd=False)

            #This would differentiate to accel, but we want velocity
            #for entries in config.COMPS:
            #  compo = entries
            #  progstring = ("%s/integ_diff diff=1 " % install.A_GP_BIN_DIR + \
            #       "filein=%s/%s.%s fileout=%s/%d.%s-lf.%s " %
            #        (a_veldir, site, compo, a_veldir, sim_id, site, compo))
            #  bband_utils.runprog(progstring)


            # Create path names and check if their sizes are within bounds
            nsfile = os.path.join(a_veldir,
                                  "%s.000" % (site))
            ewfile = os.path.join(a_veldir,
                                  "%s.090" % (site))
            udfile = os.path.join(a_veldir,
                                  "%s.ver" % (site))
            bbpfile = os.path.join(a_veldir,
                                   "%d.%s-lf.bbp" % (sim_id, site))
            bband_utils.check_path_lengths([nsfile, ewfile, udfile],
                                           bband_utils.GP_MAX_FILENAME)

            # Run wcc2bpp
            progstring = ("%s wcc2bbp=1 " %
                          (os.path.join(install.A_GP_BIN_DIR, "wcc2bbp")) +
                          'title="LP Sim NGAH, stat=%s" ' % site +
                          'nsfile=%s ewfile=%s udfile=%s > %s 2>> %s' %
                          (nsfile, ewfile, udfile, bbpfile, self.log))
            bband_utils.runprog(progstring, abort_on_error=True,
                                print_cmd=False)


        if self.r_srcfile == "":
            #calculate magnitude and write to file
            mag = fault_utils.get_magnitude(os.path.join(a_indir,
                                                         self.r_velmodel),
                                            os.path.join(a_indir,
                                                         self.r_srffile),
                                            sta_base)
            mag_file = open("%s" % os.path.join(a_indir,
                                                "magnitude_%s" %
                                                (sta_base)), 'w')
            mag_file.write("%.2f" % mag)
            mag_file.flush()
            mag_file.close()

        print("GP Jbsim Completed".center(80, '-'))

if __name__ == "__main__":
    print("Testing Module: %s" % os.path.basename((sys.argv[0])))
    INSTALL = InstallCfg()
    ME = Jbsim(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4],
               sys.argv[5], sim_id=int(sys.argv[6]))
    ME.run()
    sys.exit(0)
