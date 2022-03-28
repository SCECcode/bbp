#!/usr/bin/env python
"""
BSD 3-Clause License

Copyright (c) 2021, University of Southern California
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Module BBP Status (adapted from CSEP Status)
"""
from __future__ import division, print_function

# Import Python modules first
import os
import sys
import time

# Import Broadband modules now
from install_cfg import InstallCfg
import bband_utils

class BBPStatus(object):
    """
    This class is designed to acquire the current BBP system and
    software status
    """
    # Static data members

    # Name of the file with system status information
    system_type = "system_status"

    # Name of the file with software status information
    software_type = "software_status"

    # List of commands that used to capture external software version
    # and flag if command output is on stderr (True). If command
    # output is redirected to the stderr, the flag should be set to
    # 'True' we do not trigger it as a failure
    __all_packages = [["GCC Version", "cat $BBP_DIR/src/gcc.version", False],
                      ["GFORTRAN Version",
                       "cat $BBP_DIR/src/gfortran.version", False],
                      ["NumPy Version",
                       "python -c 'from __future__ import print_function; import numpy; print(numpy.__version__);'",
                       False],
                      ["SciPy Version",
                       "python -c 'from __future__ import print_function; import scipy; print(scipy.__version__);'",
                       False],
                      ["Matplotlib Version",
                       "python -c 'from __future__ import print_function; import matplotlib; print(matplotlib.__version__);'",
                       False],
                      ["PyProj Version",
                       "python -c 'from __future__ import print_function; import pyproj; print(pyproj.__version__);'",
                       False]]

    # Leave some envrionment variables out, to avoid capturing too much
    # personal information (and others that capture too much useless
    # stuff!)
    __exclude_env = ["LS_COLORS",
                     "SSH_CONNECTION",
                     "SSH_TTY",
                     "SSH_CLIENT",
                     "SSH_ASKPASS",
                     "MAIL"]

    #--------------------------------------------------------------------
    #
    # Initialization
    #
    # Input:
    #        options - Command-line options including defaults ones used
    #                  by caller program. Default is None.
    #
    def __init__(self, sim_id=0, options=None):
        """
        Initialization for BBPStatus class
        """
        self.sim_id = sim_id
        self.__options = options
        self.install = InstallCfg.getInstance()
        self.outdir = os.path.join(self.install.A_OUT_DATA_DIR, str(sim_id))

    #--------------------------------------------------------------------
    #
    # Get names of the files for system status
    #
    # Input: None
    #
    # Output: A tuple of data and corresponding metadata filenames
    #
    def system_filename(self):
        """
        Get the name of the system status file
        """
        filename = os.path.join(self.outdir,
                                "%s-%d.txt" % (self.system_type, self.sim_id))
        return filename

    #--------------------------------------------------------------------
    #
    # Get names of the files for system status
    #
    # Input: None
    #
    # Output: A list of data filename and corresponding metadata filename
    #
    def software_filename(self):
        """
        Get the name of the software status file
        """
        filename = os.path.join(self.outdir,
                                "%s-%d.txt" % (self.software_type, self.sim_id))
        return filename

    #--------------------------------------------------------------------
    #
    # Capture status of the system
    #
    # Input:  datafile - Name of the file to capture status information to
    #                    Default is None
    #
    # Output: None
    #
    def system(self, datafile=None):
        """
        Capture system status to the file
        """
        if datafile is None:
            datafile = self.system_filename()

        # Create data file
        try:
            fhandle = open(datafile, 'a')
        except IOError:
            print("Unable to open output data file: %s!" % (datafile))
            sys.exit(1)

        fhandle.write("\n%s\n" % ('='*80))
        fhandle.write("%s" %
                      (time.strftime("%a %d %b %Y %X %Z", time.localtime())))
        fhandle.write("\n%s\n\n" % ('='*80))

        # Store host info
        fhandle.write("%s: %s\n\n" % ("uname", os.uname()))

        # Store user info
        command = "id"
        fhandle.write("%s: %s\n" %
                      (command,
                       bband_utils.get_command_output(command,
                                                      abort_on_error=True)))

        # Store executable and command-line options
        fhandle.write("%s: %s\n\n" % ("argv", sys.argv))

        # Store executable and command-line options including the default ones
        if self.__options is not None:
            fhandle.write("%s: %s\n\n" % ("command-line options", self.__options))

        # Store environment variables
        fhandle.write("Environment Variables\n")
        for item in os.environ:
            if item in BBPStatus.__exclude_env:
                continue
            fhandle.write("%s = %s\n" % (item, os.environ[item]))

        # Store shell resource information
        fhandle.write("\nResource Information (ulimit -a)\n")
        command = "ulimit -a"
        fhandle.write("%s\n" %
                      (bband_utils.get_command_output(command,
                                                      abort_on_error=True)))

        # Close the file
        fhandle.close()

    #--------------------------------------------------------------------
    #
    # Capture status of the software used by the system
    #
    # Input:  program_name - Name of the calling program
    #         program_version - Version of the calling program
    #         datafile - Names of the file and metadata file to capture status
    #                    information to. Default is None
    #
    # Output: None
    #
    def software(self, datafile=None):
        """
        Capture software status to the file
        """
        # Make sure we have the data filename
        if datafile is None:
            datafile = self.software_filename()

        # Create data file
        try:
            fhandle = open(datafile, 'a')
        except IOError:
            print("Unable to open output data file: %s!" % (datafile))
            sys.exit(1)

        fhandle.write("\n%s\n" % ('='*80))
        fhandle.write("%s" %
                      (time.strftime("%a %d %b %Y %X %Z", time.localtime())))
        fhandle.write("\n%s\n\n" % ('='*80))

        # Store version of calling program
        fhandle.write("Broadband version: %s\n" % (self.install.VERSION))
        # Store python version
        fhandle.write("%s: %s\n\n" % ("Python version", sys.version))

        for package in BBPStatus.__all_packages:
            name = package[0]
            command = package[1]
            output_on_stderr = package[2]
            if isinstance(output_on_stderr, bool) is True:
                # Version command is provided
                msg = bband_utils.get_command_output(command,
                                                     output_on_stderr,
                                                     abort_on_error=True)
                fhandle.write("%s: %s\n" % (name, msg))
            else:
                # Version string is provided in 'output_on_stderr'
                # local variable
                fhandle.write("%s: %s\n" % (name, output_on_stderr))

        # Close the file
        fhandle.close()

    #--------------------------------------------------------------------
    #
    # Get user name of the running process
    #
    # Input: None
    #
    # Output: A username
    #
    def user_name():
        """
        Get the user name of the running process
        """
        name = bband_utils.get_command_output("whoami",
                                              abort_on_error=True)
        # Strip newline if any
        return name.replace("\n", "")

    user_name = staticmethod(user_name)

# Invoke the module
if __name__ == '__main__':
    STATUS = BBPStatus()
    # System status
    FILENAME = STATUS.system_filename()
    STATUS.system(FILENAME)

    # Software status
    FILENAME = STATUS.software_filename()
    STATUS.software(FILENAME)
# end of main
