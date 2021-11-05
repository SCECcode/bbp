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

This module defines validation events and their associated input files
"""
from __future__ import division, print_function

# Import Python modules
import os
import glob

# Import Broadband modules
import bband_utils

# Global variable where validation events go (initialized in the
# function init_validation_events below)
VE_EVENTS = None

# Constants used in this module
DEFAULT_EVENT_TYPE = "validation"
DEFAULT_EVENT_GOF_PLOT_LIMIT = "0.01"
VALIDATION_METAFILE_SUFFIX = "_validation_config.txt"
PACKAGE_VERSION_PROP = "package_version"
EVENT_TYPE_PROP = "event_type"
EVENT_NAME_PROP = "event_name"
EVENT_PRINTNAME_PROP = "event_printname"
EVENT_CUTOFF_PROP = "event_cutoff"
EVENT_OBS_PATH_PROP = "event_obs_path"
EVENT_OBS_FORMAT_PROP = "event_obs_format"
EVENT_OBS_CORR_PROP = "event_obs_corrections"
EVENT_VELOCITY_MODEL_PROP = "event_velocity_model"
EVENT_GMPE_SET_PROP = "event_gmpe_set"
EVENT_GOF_PLOT_LIMIT_PROP = "event_gof_plot_limit"
PATH_CODE_PROP = "p_codebase_"
CODE_PROP = "codebase_"

class ValidationEvents(object):
    """
    This class keep track of all validation events
    """
    def __init__(self):
        """
        This function initializes the validation events class
        """
        self.events = []

    def add_event(self, event):
        """
        This function adds a validation event to our list
        """
        self.events.append(event)
        setattr(self, event.get_validation_name(), event)

    def get_event(self, i):
        """
        This function returns the validation event object whose index
        corresponds to i
        """
        return self.events[i]

    def get_event_by_name(self, event_name):
        """
        Given a validation event name, this function returns the
        corresponding validation event object
        """
        for event in self.events:
            if event.get_validation_name() == event_name:
                return event
        return None

    def get_event_by_print_name(self, event_print_name):
        """
        Given a validation event print name, this function returns the
        corresponding validation event object
        """
        for event in self.events:
            if (event.get_print_name()).upper() == event_print_name.upper():
                return event
        return None

    def get_num_events(self):
        """
        This function returns the number of validation events available
        """
        return len(self.events)

    def get_all_names(self):
        """
        This function returns an array containing all the names of the
        available validation events
        """
        event_names = []
        for event_obj in self.events:
            event_names.append(event_obj.get_print_name())
        return event_names

class ValidationEvent(object):
    """
    This class defines a validation event and contains all the input
    files and other data associated with it.
    """
    def __init__(self, name, print_name, base_dir, event_type, version):
        """
        This function initialized the basic parameters of the
        Validation Event class
        """
        self.name = name
        self.print_name = print_name
        self.base_dir = base_dir
        self.version = version
        self.event_type = event_type
        self.args = {}
        self.obs_path = None
        self.obs_format = None
        self.obs_corrections = None
        self.mag = None
        self.cutoff = None
        self.gof_plot_limit = None
        self.vmodel_name = None
        self.gmpe_set = None

    def get_validation_name(self):
        """
        This function returns the name of a validation event
        """
        return self.name

    def get_print_name(self):
        """
        This function returns the print name of this validation event
        """
        return self.print_name

    def get_event_type(self):
        """
        This function returns the event type for this validation event
        """
        return self.event_type

    def get_version(self):
        """
        This function returns the version of this validation package
        """
        return self.version

    def set_input(self, codebase, arg, value):
        """
        This function set a codebase parameter for this validation event
        """
        self.args["%s_%s" % (codebase.upper(), arg.upper())] = value

    def set_input_path(self, codebase, arg, filename):
        """
        This function set a path codebase parameter for this validation event
        """
        if filename == "":
            self.args["%s_%s" % (codebase.upper(), arg.upper())] = ""
            return

        # Make list
        pieces = filename.split(",")
        pieces = [piece.strip() for piece in pieces]
        pieces = [os.path.join(self.base_dir, piece) for piece in pieces]

        # Check if single file or list of files
        if len(pieces) == 1:
            pieces = pieces[0]

        self.args["%s_%s" % (codebase.upper(), arg.upper())] = pieces

    def get_input(self, codebase, arg):
        """
        Returns an input parameter for a certain codebase
        """
        # A codebase-specific parameter always overrides a generic parameter
        param = None
        if ("%s_%s" % (codebase.upper(), arg.upper())) in self.args:
            param = self.args["%s_%s" % (codebase.upper(), arg.upper())]
        elif ("ALL_%s" % (arg.upper())) in self.args:
            param = self.args["ALL_%s" % (arg.upper())]
        return param

    def set_velocity_model(self, vmodel_name):
        """
        Sets the velocity model to use with this validation event
        """
        self.vmodel_name = vmodel_name

    def get_velocity_model(self):
        """
        Returns the velocity model used with this validation event
        """
        return self.vmodel_name

    def set_obs_format(self, obs_format):
        """
        This function sets the format of the observation files
        """
        self.obs_format = obs_format

    def get_obs_format(self):
        """
        Returns the format of the observation files
        """
        return self.obs_format

    def set_obs_path(self, path):
        """
        Sets the path where the observed seismograms are stored for this event
        """
        self.obs_path = path

    def get_obs_path(self):
        """
        Returns the path to the observed seismograms
        """
        return os.path.join(self.base_dir, self.obs_path)

    def set_obs_corrections(self, obs_corrections):
        """
        This function sets the path for the observation corrections
        """
        self.obs_corrections = obs_corrections.strip()

    def get_obs_corrections(self):
        """
        Returns the path to the observation corrections file
        """
        if not self.obs_corrections:
            return ''
        return os.path.join(self.base_dir, self.obs_corrections)

    def set_cutoff(self, cutoff):
        """
        Sets the cutoff for this event
        """
        self.cutoff = cutoff

    def get_cutoff(self):
        """
        Returns the cutoff for this event
        """
        return self.cutoff

    def set_gof_plot_limit(self, gof_plot_limit):
        """
        Sets the GoF plot limit for this event
        """
        self.gof_plot_limit = gof_plot_limit

    def get_gof_plot_limit(self):
        """
        Returns the GoF plot limit for this event
        """
        return self.gof_plot_limit

    def set_gmpe_set(self, gmpe_set):
        """
        Sets the gmpe set for this event
        """
        self.gmpe_set = gmpe_set

    def get_gmpe_set(self):
        """
        Returns the gmpe_set for this event
        """
        return self.gmpe_set

def init_validation_events(install_obj):
    """
    This function creates the Validation Events class and adds a
    validation event for each event supported in the Broadband
    platform
    """
    global VE_EVENTS
    VE_EVENTS = ValidationEvents()

    # First we get all subdirectories in the VAL_DIR top directory
    sub_dirs = bband_utils.list_subdirs(install_obj.A_VAL_DIR)

    # Now we add ones that have a proper configuration file
    for directory in sub_dirs:
        base_dir = os.path.join(install_obj.A_VAL_DIR, directory)

        # Look for the validation configuration file
        file_list = glob.glob(os.path.join(base_dir,
                                           "*%s" %
                                           (VALIDATION_METAFILE_SUFFIX)))
        if len(file_list) == 0:
            # No configuration file found, skip this directory
            continue
        if len(file_list) > 1:
            # Multiple configuration files found, warn user, and skip
            print("Multiple validation configuration files found!")
            print("Skipping: %s" % (base_dir))
            continue

        # Got it!
        metafile = file_list[0]

        # Parse file
        file_props = bband_utils.parse_properties(metafile)

        # Check for the required event properties
        req_props = [EVENT_NAME_PROP, EVENT_PRINTNAME_PROP,
                     EVENT_CUTOFF_PROP, EVENT_OBS_PATH_PROP,
                     EVENT_OBS_FORMAT_PROP, EVENT_OBS_CORR_PROP,
                     EVENT_VELOCITY_MODEL_PROP, PACKAGE_VERSION_PROP,
                     EVENT_GMPE_SET_PROP]
        all_properties = True
        for prop in req_props:
            if prop not in file_props:
                # Validation event doesn't have a required property
                all_properties = False
                break
        # Missing required properties, skip this event
        if not all_properties:
            continue

        # Read the event_type property, if specified
        event_type = DEFAULT_EVENT_TYPE
        if EVENT_TYPE_PROP in file_props:
            event_type = file_props[EVENT_TYPE_PROP]

        # Read GoF limit, if specified
        event_gof_plot_limit = DEFAULT_EVENT_GOF_PLOT_LIMIT
        if EVENT_GOF_PLOT_LIMIT_PROP in file_props:
            event_gof_plot_limit = file_props[EVENT_GOF_PLOT_LIMIT_PROP]

        # All properties are there, create event
        val_event = ValidationEvent(file_props[EVENT_NAME_PROP],
                                    file_props[EVENT_PRINTNAME_PROP],
                                    base_dir, event_type,
                                    file_props[PACKAGE_VERSION_PROP])
        # Set basic event parameters
        val_event.set_cutoff(int(file_props[EVENT_CUTOFF_PROP]))
        val_event.set_gof_plot_limit(float(event_gof_plot_limit))
        val_event.set_obs_path(file_props[EVENT_OBS_PATH_PROP])
        val_event.set_obs_format(file_props[EVENT_OBS_FORMAT_PROP])
        val_event.set_obs_corrections(file_props[EVENT_OBS_CORR_PROP])
        val_event.set_velocity_model(file_props[EVENT_VELOCITY_MODEL_PROP])
        val_event.set_gmpe_set(file_props[EVENT_GMPE_SET_PROP])

        # Go through the list of parsed properties and populate
        for key in file_props:
            if key.startswith(PATH_CODE_PROP):
                # This key has a codebase parameters that is a path to a file
                codebase = key.split('_')[2].upper()
                param = key.split('_', 3)[3].upper()
                value = file_props[key]
                val_event.set_input_path(codebase, param, value)
            if key.startswith(CODE_PROP):
                # This key has a regular codebase parameter
                codebase = key.split('_')[1].upper()
                param = key.split('_', 2)[2].upper()
                value = file_props[key]
                val_event.set_input(codebase, param, value)
        VE_EVENTS.add_event(val_event)
