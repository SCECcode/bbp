#!/bin/env python
"""
Southern California Earthquake Center Broadband Platform
Copyright 2010-2016 Southern California Earthquake Center

CSM Broadband Platform module
$Id: csm.py 1730 2016-09-06 20:26:43Z fsilva $
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys
import random
#import shutil

# Import Broadband modules
import bband_utils
import validation_cfg
import velocity_models
from station_list import StationList
from install_cfg import InstallCfg
from csm_cfg import CSMCfg

class CSM(object):
    """
    This class contains the glue code needed to interface the
    Composite Source Model (CSM) and the Broadband Platform.
    """
    def __init__(self, i_r_srcfile, i_r_velmodel, i_r_stations,
                 vmodel_name, val_name=None, sim_id=0):
        """
        This function initializes basic class objects
        """
        self.sim_id = sim_id
        self.r_srcfile = i_r_srcfile
        self.r_velmodel = i_r_velmodel
        self.r_stations = i_r_stations
        self.vmodel_name = vmodel_name
        self.val_name = val_name
        self.val_obj = None
        self.stat_list = None
        self.num_stations = None
        self.install = None
        self.config = None
        self.log = None
        self.csm_dir = None

    def create_mod_files(self):
        """
        This function creates the .mod files for the CSM method, based
        on the velocity model and other station parameters
        """
        vmodel = self.config.vmodel
        cfgparams = self.config.cfgparams

        # Create a velocity model file for each station
        for station in self.stat_list.getStationList():
            vmodel_file = os.path.join(self.csm_dir,
                                       "%s.mod" % station.scode)
            outfile = open(vmodel_file, 'w')
            # First write the stored velocity model
            outfile.write("%5i%5i%5i\n" % (self.config.nlay, 0, 0))
            for val1, val2, val3, val4, val5, val6 in zip(vmodel["h"],
                                                          vmodel["vp"],
                                                          vmodel["qp"],
                                                          vmodel["vs"],
                                                          vmodel["qs"],
                                                          vmodel["rho"]):
                outfile.write("%12.5f%12.5f%12.5f%12.5f%12.5f%12.5f\n" %
                              (val1, val2, val3, val4, val5, val6))
            outfile.write("%8i%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f%8.3f\n" %
                          (cfgparams["nt"], cfgparams["twins"],
                           cfgparams["t00"], cfgparams["fmax"],
                           cfgparams["fmax1"], cfgparams["ang1"],
                           cfgparams["ang2"], cfgparams["akp"]))
            outfile.write("%8.3f%8.3f%8.3f%8.3f%8.3f\n" %
                          (cfgparams["L"], cfgparams["W"],
                           cfgparams["tdep"], cfgparams["sdept"],
                           0.0))
            outfile.write("%5i%5i%8.3f%8.3f\n" %
                          (cfgparams["nx"], cfgparams["ny"],
                           cfgparams["vssp"], cfgparams["vrup"]))
            outfile.write("%5i\n" % (cfgparams["ntc"]))
            outfile.write("%12.5f%12.5f\n" %
                          (cfgparams["tlat1"], cfgparams["tlat2"]))
            outfile.write("%12.5f%12.5f\n" %
                          (cfgparams["tlon1"], cfgparams["tlon2"]))
            outfile.write("%10.2f\n" % cfgparams["D"])
            outfile.write("%10.2f\n" % cfgparams["R"])
            outfile.close()

    def create_station_list(self):
        """
        This function creates the CSM station list
        """
        # Check if Simula can handle our station list
        if self.num_stations > self.config.SIMULA_MAX_STATIONS:
            raise bband_utils.ParameterError("Too many stations in "
                                             "the station list: %d. " %
                                             (self.num_stations) +
                                             "Maximum limit is %d." %
                                             (self.config.SIMULA_MAX_STATIONS))

        csm_stations = os.path.join(self.csm_dir,
                                    self.config.csm_stations)
        outfile = open(csm_stations, 'w')
        outfile.write("%14.4f %9.4f\n" %
                      (self.config.cfgparams["hypolat"],
                       self.config.cfgparams["hypolon"]))
        outfile.write("%5i\n" % (len(self.stat_list.site_list)))
        for station in self.stat_list.getStationList():
            outfile.write("%4s %12.5f %12.5f %20s\n" %
                          (station.scode,
                           station.lat,
                           station.lon,
                           "CSMStation.txt"))
        outfile.close()

    def create_simula_stations(self):
        """
        This function creates the simula input stations file
        """
        simula_stations = os.path.join(self.csm_dir,
                                       self.config.simula_stations)
        outfile = open(simula_stations, 'w')
        for station in self.stat_list.getStationList():
            outfile.write("%s\n" % (station.scode))
        outfile.write("ENDX00\n")
        outfile.close()

    def create_simula_in(self):
        """
        This function creates the simula.in input file
        """
        simula_in = os.path.join(self.csm_dir,
                                 self.config.simula_in)
        cfgsimula = self.config.cfgparams["simula_in"]
        outfile = open(simula_in, 'w')
        outfile.write("%10i %10.4f %10.4f\n" % (cfgsimula["mfiz"],
                                                cfgsimula["flw"],
                                                cfgsimula["fhi"]))
        outfile.write("%5i %5i %8.3f %8.3f %8.3f %8.3f %8.2f %8.2f %8i " %
                      (cfgsimula["ncoda"],
                       cfgsimula["m"],
                       cfgsimula["err"],
                       cfgsimula["gi"],
                       cfgsimula["gs"],
                       cfgsimula["fmax"],
                       cfgsimula["qf"],
                       cfgsimula["cn"],
                       cfgsimula["seed2"]))
        outfile.write("%5i %8.2f %8.2f %10.4f %10.4f\n" %
                      (cfgsimula["nscat"],
                       cfgsimula["fl"],
                       cfgsimula["flb"],
                       cfgsimula["akap1"],
                       cfgsimula["akap2"]))
        outfile.write("%10s\n" % ('y'))
        outfile.close()

    def create_nuclear_in(self):
        """
        This function creates the nuclear.in input file
        """
        nuclear_in = os.path.join(self.csm_dir,
                                  self.config.nuclear_in)
        cfgnuclear = self.config.cfgparams["nuclear_in"]
        outfile = open(nuclear_in, 'w')
        outfile.write("%5i\n" % (1)) # Only 1 realization per run!
        outfile.write("%10.4f %10.4f %10.3f\n" %
                      (cfgnuclear["perw"],
                       cfgnuclear["perf"],
                       self.config.cfgparams["vrup"]))
        outfile.close()

    def create_scat1d_in(self):
        """
        This function creates the scat1d.in input file
        """
        scat1d_in = os.path.join(self.csm_dir,
                                 self.config.scat1d_in)
        cfgscat1d = self.config.cfgparams["scat1d_in"]
        outfile = open(scat1d_in, 'w')
        outfile.write("%10i %10.4f\n" %
                      (cfgscat1d["nlay"],
                       cfgscat1d["damp"]))
        outfile.write("%10.4f %10.4f %8.3f %8.3f %8.3f %8.3f %8.2f %8i\n" %
                      (cfgscat1d["dh"],
                       cfgscat1d["ddh"],
                       cfgscat1d["va"],
                       cfgscat1d["dv"],
                       cfgscat1d["dena"],
                       cfgscat1d["dden"],
                       cfgscat1d["q"],
                       cfgscat1d["seed3"]))
        outfile.close()

    def create_simula_random(self):
        """
        This function creates the file with the 100k random numbers
        needed by simula
        """
        simula_random = os.path.join(self.csm_dir,
                                     self.config.simula_random)
        # Initialize the random number generator using a variation of
        # the SEED in the SRC file
        random.seed(self.config.SEED + 13579)
        outfile = open(simula_random, 'w')
        for _ in range(0, 20000):
            for _ in range(0, 5):
                outfile.write("   %1.7e   %1.7e   %1.7e   %1.7e   %1.7e\n" %
                              (random.random(), random.random(),
                               random.random(), random.random(),
                               random.random()))
        outfile.close()

    def create_compom_in(self):
        """
        This function creates the compom.in input file
        """
        compom_in = os.path.join(self.csm_dir,
                                 self.config.compom_in)
        cfgcompom = self.config.cfgparams["compom_in"]
        outfile = open(compom_in, 'w')
        outfile.write("%10.4f %10.4f %10.3f %15e %10.3f %10.3f %12.3e\n" %
                      (self.config.cfgparams["L"],
                       self.config.cfgparams["W"],
                       self.config.cfgparams["sdept"],
                       cfgcompom["eqmom"],
                       self.config.cfgparams["sdrp1"],
                       self.config.cfgparams["sdrp2"],
                       self.config.cfgparams["mu"]))
        outfile.write("%10.4f %10.4f %10.3f %10i %10.3f\n" %
                      (cfgcompom["rmax"],
                       cfgcompom["rmin"],
                       self.config.cfgparams["fdim"],
                       self.config.cfgparams["seed1"],
                       self.config.cfgparams["plim"]))
        outfile.close()

    def create_csevents(self):
        """
        This function creates the csevents01.dat input file
        """
        csevents_dat = os.path.join(self.csm_dir,
                                    self.config.csevents_dat)
        cfg_csevents = self.config.cfgparams["csevents_dat"]
        outfile = open(csevents_dat, 'w')
        outfile.write("   %1.7e   %1.7e   %1.7e\n" %
                      (cfg_csevents["sdrpa"],
                       cfg_csevents["mu"],
                       cfg_csevents["nsub"]))
        outfile.write("   %1.7e   %1.7e   %1.7e\n" %
                      (self.config.cfgparams["seed1"],
                       0, 0))

        for val1, val2, val3 in zip(cfg_csevents["sx"],
                                    cfg_csevents["sy"],
                                    cfg_csevents["srkm"]):
            outfile.write("   %1.7e   %1.7e   %1.7e\n" %
                          (val1, val2, val3))
        outfile.close()

    def process_seismograms(self):
        """
        This function reads the seismograms generated by CSM and
        generated the velocity and acceleration seismograms in the BBP
        format needed by the platform
        """
        for station in self.stat_list.getStationList():
            in_file = os.path.join(self.csm_dir,
                                   "%s.sim01" % (station.scode))
            out_vel_file = os.path.join(self.install.A_OUT_DATA_DIR,
                                        str(self.sim_id), "%s.%s.vel.bbp" %
                                        (str(self.sim_id), station.scode))
            out_acc_file = os.path.join(self.install.A_OUT_DATA_DIR,
                                        str(self.sim_id), "%s.%s.acc.bbp" %
                                        (str(self.sim_id), station.scode))
            in_seis = open(in_file, 'r')
            out_vel = open(out_vel_file, 'w')
            out_acc = open(out_acc_file, 'w')
            out_vel.write("# Sim stat=%s\n" % (station.scode))
            out_vel.write("#    time(sec)      N-S(cm/s)      "
                          "E-W(cm/s)      U-D(cm/s)\n")
            out_acc.write("# Sim stat=%s\n" % (station.scode))
            out_acc.write("#    time(sec)      N-S(cm/s/s)      "
                          "E-W(cm/s/s)      U-D(cm/s/s)\n")
            for line in in_seis:
                line = line.strip()
                pieces = line.split()
                # Make sure lines have 10 pieces each
                if len(pieces) != 10:
                    continue
                out_vel.write("%s    %s    %s    %s\n" %
                              (pieces[0], pieces[5], pieces[4], pieces[6]))
                out_acc.write("%s    %s    %s    %s\n" %
                              (pieces[0], pieces[8], pieces[7], pieces[9]))
            out_vel.close()
            out_acc.close()
            in_seis.close()

    def run(self):
        """
        This function prepares the parameter file for CSM, invokes
        it, and formats its output to be compatible with the Broadband
        Platform
        """
        print("UNR CSM".center(80, '-'))

        self.install = InstallCfg.getInstance()
        install = self.install
        sim_id = self.sim_id

        # Get pointer to the velocity model object
        vel_obj = velocity_models.get_velocity_model_by_name(self.vmodel_name)
        if vel_obj is None:
            raise bband_utils.ParameterError("Cannot find velocity model: %s" %
                                             (self.vmodel_name))
        vmodel_params = vel_obj.get_codebase_params('csm')

        if 'SDROP' in vmodel_params:
            sdrp = float(vmodel_params['SDROP'])
        else:
            sdrp = None

        # Find validation object if this is a validation run
        if self.val_name is not None:
            self.val_obj = validation_cfg.VE_EVENTS.get_event_by_name(self.val_name)

        sta_base = os.path.basename(os.path.splitext(self.r_stations)[0])
        self.log = os.path.join(install.A_OUT_LOG_DIR,
                                str(sim_id),
                                "%d.csm_%s.log" % (sim_id, sta_base))
        a_indir = os.path.join(install.A_IN_DATA_DIR, str(sim_id))
        a_tmpdir = os.path.join(install.A_TMP_DATA_DIR, str(sim_id))
        a_tmpdir_mod = os.path.join(install.A_TMP_DATA_DIR,
                                    str(sim_id),
                                    "csm_%s" % (sta_base))
        a_outdir = os.path.join(install.A_OUT_DATA_DIR, str(sim_id))
        #
        # Make sure the output and two tmp directories exist
        #
        bband_utils.mkdirs([a_tmpdir, a_tmpdir_mod, a_outdir],
                           print_cmd=False)

        a_velmodel = os.path.join(a_indir, self.r_velmodel)
        a_stations = os.path.join(a_indir, self.r_stations)
        self.stat_list = StationList(a_stations)
        self.num_stations = len(self.stat_list.site_list)
        self.csm_dir = a_tmpdir_mod

        # Read input files, calculate CSM parameters
        self.config = CSMCfg(os.path.join(install.A_IN_DATA_DIR,
                                          str(sim_id),
                                          self.r_srcfile))
        self.config.calculate_params(a_velmodel, sdrp)

        # Create CSM station list
        self.create_station_list()

        # Create mod files
        self.create_mod_files()

        # Create Simula's station list
        self.create_simula_stations()

        # Create nuclear.in file
        self.create_nuclear_in()

        # Create simula.in file
        self.create_simula_in()

        # Create the random number file needed by Simula
        self.create_simula_random()

        # Create scat1d.in file
        self.create_scat1d_in()

        # Create compom.in file
        self.create_compom_in()

        # Create csevents01.dat file
        self.create_csevents()

        # Run in tmpdir subdir to isolate temp fortran files
        # Save cwd, change back to it at the end
        old_cwd = os.getcwd()
        os.chdir(a_tmpdir_mod)

        # Calculate the Greens Functions
        csm_gf_bin = os.path.join(install.A_UNR_BIN_DIR,
                                  "green_99v8")
        cmd = ("%s >> %s 2>&1" % (csm_gf_bin, self.log))
        bband_utils.runprog(cmd, abort_on_error=True)

        # Now, generate the simulations
        csm_sim_bin = os.path.join(install.A_UNR_BIN_DIR,
                                   "simula")
        cmd = ("%s >> %s 2>&1" % (csm_sim_bin, self.log))
        bband_utils.runprog(cmd, abort_on_error=True)

        # Restore working directory
        os.chdir(old_cwd)

        # Need to copy and re-format output seismograms
        self.process_seismograms()

        print("UNR CSM Completed".center(80, '-'))

if __name__ == "__main__":
    print("Testing Module: %s" % (os.path.basename(sys.argv[0])))
    ME = CSM(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4],
             sim_id=int(sys.argv[5]))
    ME.run()
