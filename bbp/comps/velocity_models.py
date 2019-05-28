#!/usr/bin/env python
"""
Copyright 2010-2018 University Of Southern California

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

Keeps track of the velocity models available for use in the platform
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
VELMODEL_FORMAT_PROP = "velmodel_format"
VELMODEL_ACTIVE_REGION_PROP = "velmodel_active_region"
CODE_PROP = "codebase_"

# Velocity model version needed by this version of the BBP
VELMODEL_FORMAT = "2"

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
    def __init__(self, name, version,
                 active_region, base_dir):
        """
        This function initializes the velocity model with its name
        """
        self.name = name
        self.version = version
        self.active_region = active_region
        self.base_dir = base_dir
        self.velmodels = {}
        self.codes = []
        self.params = {}

    def add_velocity_model(self, codebase, vel_file):
        """
        Sets the velocity model file for code codebase
        """
        codebase = codebase.upper()

        if codebase in self.velmodels:
            print("velocity model for codebase %s already set, ignoring %s" %
                  (codebase, vel_file))
        else:
            self.velmodels[codebase] = vel_file

    def get_velocity_model(self, codebase):
        """
        Return the velocity model file for a certain codebase
        """
        return os.path.join(self.base_dir, self.velmodels[codebase.upper()])

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

    def is_active_region(self):
        """
        Returns True if this is an active tectonic region
        """
        return self.active_region

    def add_code(self, codebase):
        """
        Adds a codebase to this velocity model
        """
        codebase = codebase.upper()

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
        self.params[codebase.upper()][param] = value

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
            print("Cannot find velocity model name attribute!")
            print("Skipping: %s" % (base_dir))
            continue
        if VELMODEL_FORMAT_PROP not in file_props:
            # Doesn't have a format key
            print("Cannot find velocity model format attribute!")
            print("Skipping: %s" % (base_dir))
            continue
        if VELMODEL_VERSION_PROP not in file_props:
            # Doesn't have version key
            print("Cannot find velocity model version attribute!")
            print("Skipping: %s" % (base_dir))
            continue
        if VELMODEL_ACTIVE_REGION_PROP not in file_props:
            # Missing active region flag
            print("Cannot find active region attribute!")
            print("Skipping: %s" % (base_dir))
            continue
        # Check if velocity model format matching this version of the BBP
        velmodel_format = file_props[VELMODEL_FORMAT_PROP]
        if velmodel_format.split(".")[0] != VELMODEL_FORMAT:
            # Major velocity model format version has to match
            print("Cannot use velocity model format %s!" % (velmodel_format))
            print("Skipping: %s" % (base_dir))
            continue
        if file_props[VELMODEL_ACTIVE_REGION_PROP].lower() in ["1", "true"]:
            active_region = True
        elif file_props[VELMODEL_ACTIVE_REGION_PROP].lower() in ["0", "false"]:
            active_region = False
        else:
            # Cannot parse active region flag
            print("Cannot parse velocity model active region flag!")
            print("Skipping: %s" % (base_dir))
            continue
        ve_model = VelocityModel(file_props[VELMODEL_NAME_PROP],
                                 file_props[VELMODEL_VERSION_PROP],
                                 active_region,
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
