#!/usr/bin/env python
"""
Southern California Earthquake Center Broadband Platform
Copyright 2010-2016 Southern California Earthquake Center

This config class will encapsulate the installation configuration
parameters and several directory paths required by the Platform. There
are no user-changeable parameters in this file.
$Id: install_cfg.py 1727 2016-09-01 20:02:22Z fsilva $
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys
import time

# Import Broadband modules
import validation_cfg
import velocity_models

instance = None

class InstallCfg(object):
    """
    Define the configuration parameters that need to be edited when the
    code is moved to a new computer or computer account
    """

    @staticmethod
    def getInstance():
        """
        This function returns an existing instance of InstallCfg. If
        one doesn't exit yet, we create one.
        """
        global instance
        if instance is None:
            instance = InstallCfg()
        return instance

    def __init__(self):
        """
        This function sets up all relative directories in the
        platform. It requires the user to set up both BBP_DIR and
        BBP_GF_DIR environment variables. BBP_VAL_DIR is optional, it
        is only required for validation runs. The BBP_DATA_DIR
        environment variable defines the location where all simulation
        results will go. The Platform will create the run, indata,
        tmpdata, outdata, logs and xml subdirectories there.
        """
        #################################################################
        #                                                               #
        #                         Attention!!!                          #
        #      Users shouldn't need to change anything in this file     #
        #                                                               #
        #################################################################

        # First we make sure we have all environment variables that we need
        if "BBP_DIR" in os.environ:
            self.A_INSTALL_ROOT = os.path.normpath(os.environ["BBP_DIR"])
        else:
            print("BBP_DIR is not set. Please configure it and try again.")
            sys.exit(-1)

        if "BBP_GF_DIR" in os.environ:
            self.A_GF_DIR = os.path.normpath(os.environ["BBP_GF_DIR"])
        else:
            print("BBP_GF_DIR is not set. Please configure it and try again.")
            sys.exit(-1)

        if "BBP_VAL_DIR" in os.environ:
            self.A_VAL_DIR = os.path.normpath(os.environ["BBP_VAL_DIR"])
        else:
            # The validation package is not required... Just keep
            # track it is missing in case the user wants to try a
            # validation run...
            self.A_VAL_DIR = None

        if "BBP_DATA_DIR" in os.environ:
            self.A_DATA_ROOT = os.path.normpath(os.environ["BBP_DATA_DIR"])
        else:
            print("BBP_DATA_DIR is not set. Please configure it and try again.")
            sys.exit(-1)

        #
        # Make sure these exist before continuing and that you have permissions
        #
        if not os.path.exists(self.A_INSTALL_ROOT):
            print("Your broadband install root directory %s doesn't exist." %
                  (self.A_INSTALL_ROOT))
            print("Please make sure to set BBP_DIR correctly to point to")
            print("your broadband install root directory and try again.")
            sys.exit(-1)
        if not os.path.exists(self.A_DATA_ROOT):
            print("Your broadband data root directory %s doesn't exist." %
                  (self.A_DATA_ROOT))
            print("Please set BBP_DATA_DIR correctly and try again.")
            sys.exit(-1)
        if not (os.access(self.A_INSTALL_ROOT, os.R_OK)):
            print("You don't have read access to %s, which is your " %
                  (self.A_INSTALL_ROOT))
            print("broadband install root directory. If this is incorrect,")
            print("please set BBP_DIR correctly to point to your")
            print("broadband install root directory and try again.")
            sys.exit(-2)
        if not (os.access(self.A_DATA_ROOT, os.W_OK)):
            print("You don't have write access to %s, which is your" %
                  (self.A_DATA_ROOT))
            print("broadband data root directory. If this is incorrect,")
            print("please set BBP_DATA_DIR correctly to point to")
            print("your broadband data root directory and try again.")
            sys.exit(-2)
        if not os.path.exists(self.A_GF_DIR):
            print("Your broadband data directory %s doesn't exist." %
                  (self.A_GF_DIR))
            print("Please make sure to set BBP_GF_DIR correctly to point to")
            print("your broadband data directory and try again.")
            sys.exit(-3)
        if not os.access(self.A_GF_DIR, os.R_OK):
            print("You don't have read access to %s," %
                  (self.A_GF_DIR))
            print("which is your broadband data directory. If this is")
            print("incorrect, please set BBP_GF_DIR correctly to point to")
            print("your broadband data directory and try again.")
            sys.exit(-4)
        if self.A_VAL_DIR is not None:
            if not os.access(self.A_VAL_DIR, os.R_OK):
                print("You don't have read access to %s," %
                      (self.A_VAL_DIR))
                print("which is your broadband validation directory. If this")
                print("is incorrect, please set BBP_VAL_DIR correctly to point")
                print("to your broadband validation directory and try again.")
                sys.exit(-4)

        #
        # Component installation info
        #
        self.A_COMP_DIR = os.path.join(self.A_INSTALL_ROOT, "comps")
        self.A_TEST_DIR = os.path.join(self.A_INSTALL_ROOT, "tests")
        self.A_SRC_DIR = os.path.join(self.A_INSTALL_ROOT, "src")

        #
        # Read Broadband version
        #
        version_file = open(os.path.join(self.A_COMP_DIR, "version.txt"), 'r')
        self.VERSION = version_file.readline().strip()
        version_file.close()

        #
        # Acceptance and Unit Test Directories
        #
        self.A_TEST_REF_DIR = os.path.join(self.A_TEST_DIR, "ref_data")

        #
        # User data info
        #
        self.A_USER_DATA_DIR = os.path.join(self.A_DATA_ROOT, "run")
        self.A_IN_DATA_DIR = os.path.join(self.A_DATA_ROOT, "indata")
        self.A_OUT_DATA_DIR = os.path.join(self.A_DATA_ROOT, "outdata")
        self.A_TMP_DATA_DIR = os.path.join(self.A_DATA_ROOT, "tmpdata")
        self.A_XML_DIR = os.path.join(self.A_DATA_ROOT, "xml")

        #
        # Log file info
        #
        self.A_OUT_LOG_DIR = os.path.join(self.A_DATA_ROOT, "logs")

        # Create these directories if they don't exist
        os.system("mkdir -p %s" % self.A_IN_DATA_DIR)
        os.system("mkdir -p %s" % self.A_TMP_DATA_DIR)
        os.system("mkdir -p %s" % self.A_OUT_DATA_DIR)
        os.system("mkdir -p %s" % self.A_OUT_LOG_DIR)
        os.system("mkdir -p %s" % self.A_XML_DIR)
        os.system("mkdir -p %s" % self.A_USER_DATA_DIR)

        #
        # GP Directories
        #
        self.A_GP_BIN_DIR = os.path.join(self.A_SRC_DIR, "gp", "bin")
        if not os.path.exists(self.A_GP_BIN_DIR):
            print("Can't find GP bin directory %s." % (self.A_GP_BIN_DIR))
            print("Did you successfully build the executables?")
            print("If not, please run make in %s." % (self.A_SRC_DIR))
            sys.exit(3)

        #
        # SDSU Directories
        #
        self.A_SDSU_BIN_DIR = os.path.join(self.A_INSTALL_ROOT,
                                           "src", "sdsu", "bin")
        self.A_SDSU_DATA_DIR = os.path.join(self.A_INSTALL_ROOT,
                                            "mod_data", "sdsu")
        if not os.path.exists(self.A_SDSU_BIN_DIR):
            print("Can't find SDSU bin directory %s." % (self.A_SDSU_BIN_DIR))
            print("Did you successfully build the executables?")
            print("If not, please run make in %s." % (self.A_SRC_DIR))
            sys.exit(3)

        #
        # UCSB Directories
        #
        self.A_UCSB_BIN_DIR = os.path.join(self.A_INSTALL_ROOT,
                                           "src", "ucsb", "bin")
        self.A_UCSB_DATA_DIR = os.path.join(self.A_INSTALL_ROOT,
                                            "mod_data", "ucsb")
        if not os.path.exists(self.A_UCSB_BIN_DIR):
            print("Can't find UCSB bin directory %s." % (self.A_UCSB_BIN_DIR))
            print("Did you successfully build the executables?")
            print("If not, please run make in %s." % (self.A_SRC_DIR))
            sys.exit(3)

        #
        # UCB Directories
        #
        self.A_UCB_BIN_DIR = os.path.join(self.A_INSTALL_ROOT,
                                          "src", "ucb", "rotd50")
        if not os.path.exists(self.A_UCSB_BIN_DIR):
            print("Can't find UCB bin directory %s." % (self.A_UCB_BIN_DIR))
            print("Did you successfully build the executables?")
            print("If not, please run make in %s." % (self.A_SRC_DIR))
            sys.exit(3)

        #
        # UWO Directories
        #
        self.A_UWO_BIN_DIR = os.path.join(self.A_INSTALL_ROOT,
                                          "src", "uwo", "bin")
        if not os.path.exists(self.A_UWO_BIN_DIR):
            print("Can't find UWO bin directory %s." % (self.A_UWO_BIN_DIR))
            print("Did you successfully build the executables?")
            print("If not, please run make in %s." % (self.A_SRC_DIR))
            sys.exit(3)

        #
        # UNR Directories
        #
        self.A_UNR_BIN_DIR = os.path.join(self.A_INSTALL_ROOT,
                                          "src", "unr", "bin")
        if not os.path.exists(self.A_UNR_BIN_DIR):
            print("Can't find UNR bin directory %s." % (self.A_UNR_BIN_DIR))
            print("Did you successfully build the executables?")
            print("If not, please run make in %s." % (self.A_SRC_DIR))
            sys.exit(3)

        #
        # Irikura Directories
        #
        self.A_IRIKURA_BIN_DIR = os.path.join(self.A_INSTALL_ROOT,
                                              "src", "irikura", "bin")
        if not os.path.exists(self.A_IRIKURA_BIN_DIR):
            print("Can't find IRIKURA bin directory %s." %
                  (self.A_IRIKURA_BIN_DIR))
            print("Did you successfully build the executables?")
            print("If not, please run make in %s." % (self.A_SRC_DIR))
            sys.exit(3)

        #
        # USGS Directories
        #
        self.A_USGS_BIN_DIR = os.path.join(self.A_INSTALL_ROOT,
                                           "src", "usgs", "bin")
        if not os.path.exists(self.A_USGS_BIN_DIR):
            print("Can't find USGS bin directory %s." %
                  (self.A_USGS_BIN_DIR))
            print("Did you successfully build the executables?")
            print("If not, please run make in %s." % (self.A_SRC_DIR))
            sys.exit(3)

        # Plot directories
        self.A_PLOT_DATA_DIR = os.path.join(self.A_INSTALL_ROOT, "plot")

        # Initialize velocity models
        velocity_models.init_velocity_models(self)

        # Validation events setup
        validation_cfg.init_validation_events(self)

        # Keep track of starting time
        self.start_time = None

    def set_start_time(self):
        """
        This function sets the simulation start time
        """
        self.start_time = time.localtime()

if __name__ == "__main__":
    print("Test Config Class: %s" % (__file__))
