#!/usr/bin/env python
"""
Southern California Earthquake Center Broadband Platform
Copyright 2010-2016 Southern California Earthquake Center

Keeps track of the velocity models available for use in the platform
$Id: velocity_models.py 1730 2016-09-06 20:26:43Z fsilva $
"""
from __future__ import division, print_function

# Import Python modules
import os
import glob

# Import Broadband modules
import bband_utils

# Global variable with a list of available velocity models
VE_MODELS = []

# Constants used in this module
VELMODEL_METAFILE_SUFFIX = "_velmodel_config.txt"
VELMODEL_NAME_PROP = "velmodel_name"
VELMODEL_CODE_PROP = "velmodel_code_"
VELMODEL_VERSION_PROP = "velmodel_version"
CODE_PROP = "codebase_"

class MissingVelocityModel(Exception):
    """
    Exception used to indicate that a requested velocity model is not
    available
    """
    pass

class VelocityModel(object):
    """
    This class stores the needed information for a given velocity model
    """
    def __init__(self, name, version, base_dir):
        """
        This function initializes the velocity model with its name
        """
        self.name = name
        self.version = version
        self.base_dir = base_dir
        self.velmodels = {}
        self.codes = []
        self.params = {}

    def add_velocity_model(self, codebase, vel_file):
        """
        Sets the velocity model file for code codebase
        """
        if codebase in self.velmodels:
            print("velocity model for codebase %s already set, ignoring %s" %
                  (codebase, vel_file))
        else:
            self.velmodels[codebase] = vel_file

    def get_velocity_model(self, codebase):
        """
        Return the velocity model file for a certain codebase
        """
        return os.path.join(self.base_dir, self.velmodels[codebase])

    def get_name(self):
        """
        Returns the velocity model name
        """
        return self.name

    def get_version(self):
        """
        Returns the velocity model version tag
        """
        return self.version

    def add_code(self, codebase):
        """
        Adds a codebase to this velocity model
        """
        if codebase in self.codes:
            print("codebase %s already defined, ignoring duplicate" %
                  (codebase))
        else:
            self.codes.append(codebase)
            # Start new empty dictionary
            self.params[codebase] = {}

    def get_codebases(self):
        """
        Returns the available codebases for this velocity model
        """
        return self.codes

    def get_codebase_params(self, codebase):
        """
        Returns the parameter dictionary associated with a certain codebase
        """
        return self.params[codebase.upper()]

    def add_codebase_param(self, codebase, param, value):
        """
        Adds a param=value to a certain codebase
        """
        self.params[codebase][param] = value

def init_velocity_models(install_obj):
    """
    This function compiles a list of available velocity models to be
    used within the Bradband platform by scanning the directories
    insite the Greens Functions directory
    """
    # First we get all subdirectories in the GF_DIR top directory
    sub_dirs = bband_utils.list_subdirs(install_obj.A_GF_DIR)

    # Now we add ones that have a proper configuration file
    for directory in sub_dirs:
        base_dir = os.path.join(install_obj.A_GF_DIR, directory)

        # Look for velocity model configuration file
        file_list = glob.glob(os.path.join(base_dir,
                                           "*%s" %
                                           (VELMODEL_METAFILE_SUFFIX)))

        if len(file_list) == 0:
            # No configuration file found, skip this directory
            continue
        if len(file_list) > 1:
            # Multiple configuration files found, warn user, and skip
            print("Multiple velocity model configuration files found!")
            print("Skipping: %s" % (base_dir))
            continue

        # Got it!
        metafile = file_list[0]

        # Parse file
        file_props = bband_utils.parse_properties(metafile)
        if VELMODEL_NAME_PROP not in file_props:
            # Doesn't have name key
            continue
        if VELMODEL_VERSION_PROP not in file_props:
            # Doesn't have version key
            continue
        ve_model = VelocityModel(file_props[VELMODEL_NAME_PROP],
                                 file_props[VELMODEL_VERSION_PROP],
                                 base_dir)
         # Dictionaries are not ordered, so we have to go over them
        # twice to first collect the available codebases and their
        # velocity model files, and later collect their
        # codebase-specific parameters
        for key in file_props:
            if key.startswith(VELMODEL_CODE_PROP):
                # This key has a velocity model file mapping
                codebase = key.split('_')[2].upper()
                ve_model.add_velocity_model(codebase, file_props[key])
                ve_model.add_code(codebase)
        for key in file_props:
            if key.startswith(CODE_PROP):
                # This key has a codebase-specific parameter
                codebase = key.split('_')[1].upper()
                # For parameter, only split the first two '_'
                param = key.split('_', 2)[2].upper()
                if codebase not in ve_model.get_codebases():
                    print("Codebase %s not defined, ignoring parameter %s" %
                          (codebase, param))
                    continue
                # Add parameter to codebase
                ve_model.add_codebase_param(codebase, param, file_props[key])
        VE_MODELS.append(ve_model)

def get_velocity_model_by_name(name):
    """
    This function returns a velocity_model object given its name, it
    returns None, if no velocity model object exists with that name
    """
    for vmodel_obj in VE_MODELS:
        if vmodel_obj.get_name() == name:
            # Found it!
            return vmodel_obj
    return None

def get_all_names():
    """
    This function returns an array containing all the names of the
    available velocity models
    """
    models = []
    for vmodel_obj in VE_MODELS:
        models.append(vmodel_obj.get_name())
    return models
