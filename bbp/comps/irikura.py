#!/usr/bin/env python
"""
Southern California Earthquake Center Broadband Platform
Copyright 2010-2016 Southern California Earthquake Center
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys

# Import Broadband modules
import plot_srf
import bband_utils
from irikura_cfg import IrikuraCfg
from install_cfg import InstallCfg

class Irikura(object):
    """
    Implements Arben's gen_srf.csh script in Python
    """
    def __init__(self, i_r_velmodel, i_r_srcfile,
                 o_r_srffile, i_vmodel_name, sim_id=0):
        self.sim_id = sim_id
        self.r_velmodel = i_r_velmodel
        self.r_srcfile = i_r_srcfile
        self.r_srffile = o_r_srffile
        self.vmodel_name = i_vmodel_name

    def run(self):
        """
        This function prepares the parameters for Irikura's gen_srf then calls it
        """
        # Load configuration, set sim_id
        install = InstallCfg.getInstance()
        sim_id = self.sim_id

        # Build directory paths
        a_tmpdir = os.path.join(install.A_TMP_DATA_DIR, str(sim_id))
        a_indir = os.path.join(install.A_IN_DATA_DIR, str(sim_id))
        a_outdir = os.path.join(install.A_OUT_DATA_DIR, str(sim_id))
        a_logdir = os.path.join(install.A_OUT_LOG_DIR, str(sim_id))

        # Make sure the output and tmp directories exist
        bband_utils.mkdirs([a_tmpdir, a_indir, a_outdir])

        # Now, file paths
        self.log = os.path.join(a_logdir, "%d.gen_srf.log" % (sim_id))
        a_srcfile = os.path.join(a_indir, self.r_srcfile)

        # Read src file
        cfg = IrikuraCfg(a_srcfile)

        # Define location of input velocity model and output srf file
        a_srffile = os.path.join(a_indir, self.r_srffile)
        a_velmod = os.path.join(install.A_IN_DATA_DIR, str(sim_id),
                                self.r_velmodel)

        # Run in tmpdir subdir to isolate temp fortran files
        # Save cwd, change back to it at the end
        old_cwd = os.getcwd()
        os.chdir(a_tmpdir)

        #
        # Run gen_srf code
        #
        progstring = ("%s >> %s 2>&1 << END\n" %
                      (os.path.join(install.A_IRIKURA_BIN_DIR, cfg.GENSRF),
                       self.log) +
                      "%s\n" % a_srffile +
                      "%f %f %f %f %f\n" %
                      (cfg.CFGDICT["fault_length"],
                       cfg.CFGDICT["fault_width"],
                       cfg.CFGDICT['strike'],
                       cfg.CFGDICT['dip'],
                       cfg.CFGDICT['rake']) +
                      "%f %f %f\n" %
                      (cfg.CFGDICT["lon_top_center"],
                       cfg.CFGDICT["lat_top_center"],
                       cfg.CFGDICT["depth_to_top"]) +
                      "%f %f\n" % (cfg.CFGDICT["dlen"],
                                   cfg.CFGDICT["dwid"]) +
                      "%f %f %f %f\n" %
                      (cfg.CFGDICT["hypo_along_stk"],
                       cfg.CFGDICT["hypo_down_dip"],
                       cfg.DENS, cfg.VS) +
                      "%f\n" % (cfg.DT) +
                      "%d\n" % (int(cfg.CFGDICT['seed'])) +
                      "%s\n" % a_velmod +
                      "END")
        bband_utils.runprog(progstring)

        # Restore working directory
        os.chdir(old_cwd)

        #
        # mv result to outputfile
        #
        progstring = "cp %s %s" % (a_srffile,
                                   os.path.join(a_tmpdir, self.r_srffile))
        bband_utils.runprog(progstring)
        progstring = "cp %s %s" % (a_srffile,
                                   os.path.join(a_outdir, self.r_srffile))
        bband_utils.runprog(progstring)

        # Plot SRF
        plot_srf.run(self.r_srffile, sim_id=self.sim_id)

        print("Completed Irikura gen_srf .....")

if __name__ == "__main__":
    print("Testing Module: %s" % os.path.basename((sys.argv[0])))
    ME = Irikura(sys.argv[1], sys.argv[2], sys.argv[3],
                 sys.argv[4], sim_id=int(sys.argv[5]))
    ME.run()
