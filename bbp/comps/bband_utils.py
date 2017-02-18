#!/usr/bin/env python
"""
Southern California Earthquake Center Broadband Platform
Copyright 2010-2016 Southern California Earthquake Center

Utility classes for the SCEC Broadband Platform
$Id: bband_utils.py 1730 2016-09-06 20:26:43Z fsilva $
"""
from __future__ import division, print_function

# Import Python modules
import os
import re
import sys
import traceback
import subprocess

# Compile regular expressions
re_parse_property = re.compile(r'([^:= \t]+)\s*[:=]?\s*(.*)')

# Constants used by several Python scripts

# This is used to convert from accel in g to accel in cm/s/s
G2CMSS = 980.665 # Convert g to cm/s/s

# Set to the maximum allows filename in the GP codebase
GP_MAX_FILENAME = 256
# Set to the maximum allowed filename in the SDSU codebase
SDSU_MAX_FILENAME = 256

class BroadbandExternalError(Exception):
    """
    Exception when an external program invoked by the Broadband
    platform fails
    """
    pass

class ParameterError(Exception):
    """
    Exception when a parameter provided to a module is deemed invalid
    """
    pass

class ProcessingError(Exception):
    """
    Exception raised when a Broadband module finds an unrecoverable
    error during processing and cannot go any further
    """
    pass

def runprog(cmd, print_cmd=True, abort_on_error=False):
    """
    Run a program on the command line and capture the output and print
    the output to stdout
    """
    # Check if we have a binary to run
    if not os.access(cmd.split()[0], os.X_OK) and cmd.startswith("/"):
        raise BroadbandExternalError("%s does not seem an executable path!" %
                                     (cmd.split()[0]))

    try:
        if print_cmd:
            print("Running: %s" % (cmd))
        proc = subprocess.Popen(cmd, shell=True)
        proc.wait()
    except KeyboardInterrupt:
        print("Interrupted!")
        sys.exit(1)
    except:
        print("Unexpected error returned from Subprocess call: ",
              sys.exc_info()[0])

    if abort_on_error:
        # If we got a non-zero exit code, abort
        if proc.returncode != 0:
            # Check if interrupted
            if proc.returncode is None:
                raise BroadbandExternalError("%s\n" %
                                             (traceback.format_exc()) +
                                             "%s failed!" %
                                             (cmd))
            raise BroadbandExternalError("%s\n" %
                                         (traceback.format_exc()) +
                                         "%s returned %d" %
                                         (cmd, proc.returncode))

    return proc.returncode

def get_command_output(cmd, output_on_stderr=False, abort_on_error=False):
    """
    Get the output of the command from the shell. Adapter from CSEP's
    commandOutput function in Environment.py
    """
    # Execute command using the UNIX shell
    child = subprocess.Popen(cmd,
                             shell=True,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
    child_data, child_error = child.communicate()

    if child_error and output_on_stderr is False:
        if abort_on_error:
            error_msg = ("Child process '%s' failed with error code %s" %
                         (cmd, child_error))
            raise BroadbandExternalError("%s\n" %
                                         (traceback.format_exc()) +
                                         "%s" % (error_msg))
        else:
            return ""

    # Check for non-empty result string from the command
    if ((child_data is None or len(child_data) == 0) and
        output_on_stderr is False):
        if abort_on_error:
            error_msg = "Child process '%s' returned no data!" % (cmd)
            raise BroadbandExternalError("%s\n" %
                                         (traceback.format_exc()) +
                                         "%s" % (error_msg))
        else:
            return ""

    # If command output is on stderr
    if output_on_stderr is True:
        child_data = child_error

    return child_data

def mkdirs(list_of_dirs, print_cmd=True):
    """
    Creates all directories specified in the list_of_dirs
    """
    for my_dir in list_of_dirs:
        cmd = "mkdir -p %s" % (my_dir)
        runprog(cmd, print_cmd=print_cmd, abort_on_error=True)

def relpath(path, start=os.curdir):
    """
    Return a relative version of a path
    (from Python 2.6 os.path.relpath() implementation)
    """
    sep = os.sep
    if not path:
        raise ValueError("no path specified")

    start_list = os.path.abspath(start).split(sep)
    path_list = os.path.abspath(path).split(sep)

    # Work out how much of the filepath is shared by start and path.
    i = len(os.path.commonprefix([start_list, path_list]))

    rel_list = [os.pardir] * (len(start_list) - i) + path_list[i:]
    if not rel_list:
        return '.'
    return os.path.join(*rel_list)

def check_path_lengths(variables, max_length):
    """
    This function checks each variable in the variables list and makes
    sure their path lenghts are less than max_length. It raises a
    ValueError exception otherwise.
    """
    for var in variables:
        if len(var) > max_length:
            raise ValueError("Path len for %s " % (var) +
                             " is %d characters long, maximum is %d" %
                             (len(var), max_length))

def list_subdirs(d):
    """
    This function returns all subdirectories inside the directory d
    """
    # Return empty array if d is None
    if d is None:
        return []
    # Use list comprehension
    return [sub for sub in os.listdir(d) if os.path.isdir(os.path.join(d, sub))]

def parse_properties(filename):
    """
    This function reads all properties from filename and returns a
    dictionary containing all key=value pairs found in the file
    """
    my_file = open(filename, 'r')
    props = {}

    for line in my_file:
        # Strip tabs, spaces and newlines from both ends
        line = line.strip(' \t\n')
        # Skip comments
        if line.startswith('#'):
            continue
        # Remove inline comments
        line = line.split('#')[0]
        # Skip empty lines
        if len(line) == 0:
            continue
        result = re_parse_property.search(line)
        if result:
            # Property parsing successful
            key = result.group(1)
            val = result.group(2)
            # Make key lowercase
            key = key.lower()
            props[key] = val

    # Don't forget to close the file
    my_file.close()

    # All done!
    return props

def parse_src_file(a_srcfile):
    """
    Function parses the SRC file and checks for needed keys. It
    returns a dictionary containing the keys found in the src file.
    """
    src_keys = parse_properties(a_srcfile)
    required_keys = ["magnitude", "fault_length", "fault_width", "dlen",
                     "dwid", "depth_to_top", "strike", "rake", "dip",
                     "lat_top_center", "lon_top_center"]
    for key in required_keys:
        if key not in src_keys:
            raise ParameterError("key %s missing in src file" % (key))
    # Convert keys to floats
    for key in src_keys:
        src_keys[key] = float(src_keys[key])

    return src_keys

def count_header_lines(a_bbpfile):
    """
    Function counts and returns the number of header lines in a BBP file
    """
    header_lines = 0

    my_file = open(a_bbpfile, 'r')
    for line in my_file:
        line = line.strip()
        # Check for empty lines, we count them too
        if not line:
            header_lines = header_lines + 1
            continue
        # Check for comments
        if line.startswith('%') or line.startswith('#'):
            header_lines = header_lines + 1
            continue
        # Reached non header line
        break
    my_file.close()

    return header_lines

if __name__ == "__main__":
    print("Testing: %s" % (sys.argv))
    CMD = "/bin/date"
    RESULT = runprog(CMD)
    if RESULT != 0:
        print("Error running cmd: %s" % (CMD))
    else:
        print("Success!")
