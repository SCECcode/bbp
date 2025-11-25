#!/usr/bin/env python3
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

SCEC Broadband Platform
"""
from __future__ import division, print_function

# Works for both Python 2 and 3
try: input = raw_input
except NameError: pass

# Import Python modules first
import os
import sys
import atexit
import shutil
import optparse
import tempfile
import traceback

# Function that will remove the matplotlib tempdir
def cleanup(dir_name):
    """
    This function removes the temporary matplotlib directory
    """
    shutil.rmtree(os.environ["MPLCONFIGDIR"])

# Create a temporary directory for matplotlib
os.environ["MPLCONFIGDIR"] = tempfile.mkdtemp()
# Now register function to remove it at the end
atexit.register(cleanup, os.environ["MPLCONFIGDIR"])

# Import Broadband modules now
from install_cfg import InstallCfg
import bband_utils
import seqnum
import bbp_status
import xml_handler
import build_workflow
from optfile import OptFile

# Here some global variables
install = InstallCfg.getInstance()

class Logger(object):
    """
    Class used to re-direct stdout and stderr to a logfile
    """
    def __init__(self, filename):
        """
        Saves log filename and opens the file for output
        """
        self.filename = filename
        self.fp = open(self.filename, 'w')

    def write(self, string):
        """
        Write output to the logfile
        """
        self.fp.write(string)

    def flush(self):
        """
        Flush output
        """
        self.fp.flush()

    def close(self):
        """
        Closes and flushes the logfile
        """
        self.fp.flush()
        self.fp.close()

# -----------------------------------------------------------------------------
# Main starts here
# -----------------------------------------------------------------------------

parser = optparse.OptionParser()
parser.add_option("-x", "--xml-file", dest="xmlFile",
                  help="Run using XML description of workflow",
                  metavar="XML_FILE")
parser.add_option("-s", "--sim-id", dest="inputSimID", type="int",
                  help="Force a sim id", metavar="SIM_ID")
parser.add_option("-o", "--option-file", dest="optFile",
                  help=("File containing responses to interactive "
                        "platform prompts"))
parser.add_option("-v", "--version", action="store_true", dest="version",
                  help="Broadband platform version")
parser.add_option("-g", "--generate-only", action="store_true", dest="generate",
                  help=("Generates the XML description but does "
                        "not run the platform"))
parser.add_option("-l", "--log", dest="logFile",
                  help="Directs output to a file, use to run BBP in background",
                  metavar="LOG_FILE")
parser.add_option("-m", "--no-xml", action="store_true", dest="noxml",
                  help="Do not generate xml")
parser.add_option("-r", "--resume", dest="resume_module",
                  help="Resume workflow from a certain module")
parser.add_option("-e", "--end", dest="end_module",
                  help="End workflow after a certain module")
parser.add_option("--expert", action="store_true", dest="expert",
                  help="Turn on expert mode")

(options, args) = parser.parse_args()

if options.version == True:
    print("Version: %s" % (install.VERSION))
    sys.exit(0)

# Resolve conflicting command-line options
if options.xmlFile is not None and options.optFile is not None:
    parser.error("Options -x and -o are mutually exclusive.")

if (options.generate == True) and (options.noxml == True):
    parser.error("Options -g and -m are mutually exclusive.")

if options.resume_module is not None:
    if options.xmlFile is None or options.inputSimID is None:
        parser.error("Resume option requires both -x and -s options.")

# Check if user specified sim_id
if options.inputSimID:
    try:
        sim_id = int(options.inputSimID)
        if sim_id < 0:
            print("SimID must be a nonnegative integer.")
            sys.exit(2)
    except ValueError:
        print("SimID must be a nonnegative integer.")
        sys.exit(3)
else:
    # Otherwise, get sim_id from timestamp
    sim_id = int(seqnum.get_seq_num())

a_indir = os.path.join(install.A_IN_DATA_DIR, str(sim_id))
a_tmpdir = os.path.join(install.A_TMP_DATA_DIR, str(sim_id))
a_outdir = os.path.join(install.A_OUT_DATA_DIR, str(sim_id))
a_logdir = os.path.join(install.A_OUT_LOG_DIR, str(sim_id))

# Create indir, tmpdir, outdir, and logdir, if they don't exist
if not os.path.exists(a_indir):
    cmd = "mkdir %s" % (a_indir)
    rc = bband_utils.runprog(cmd, False)
    if rc != 0:
        print("Failed to create in directory %s, aborting." % (a_indir))
        sys.exit(1)
if not os.path.exists(a_tmpdir):
    cmd = "mkdir %s" % (a_tmpdir)
    rc = bband_utils.runprog(cmd, False)
    if rc != 0:
        print("Failed to create tmp directory %s, aborting." % (a_tmpdir))
        sys.exit(1)
if not os.path.exists(a_outdir):
    cmd = "mkdir %s" % (a_outdir)
    rc = bband_utils.runprog(cmd, False)
    if rc != 0:
        print("Failed to create out directory %s, aborting." % (a_outdir))
        sys.exit(1)
if not os.path.exists(a_logdir):
    cmd = "mkdir %s" % (a_logdir)
    rc = bband_utils.runprog(cmd, False)
    if rc != 0:
        print("Failed to create log directory %s, aborting." % (a_logdir))
        sys.exit(1)

if options.xmlFile is not None:
    # If we have a XML File, use it
    stations = ""
    val_obj = None
    input_file = options.xmlFile
    workflow_obj = xml_handler.parse_xml(input_file)
    if workflow_obj.val_obj is not None:
        val_obj = workflow_obj.val_obj
    stations = workflow_obj.station_file
    workflow = workflow_obj.workflow
    workflow_obj = xml_handler.Workflow(workflow, val_obj,
                                        stations, install.VERSION)
else:
    # Not starting from an XML file
    if options.optFile is not None:
        # If we have an option file, we get our answers from there
        choice_obj = OptFile(options.optFile)
        # First option selects validation (y) or simulation (n)
        val_choice = choice_obj.get_next_option()
        if val_choice.lower() != 'y' and val_choice.lower() != 'n':
            print("Invalid input.")
            sys.exit(1)
    else:
        # Otherwise, we go interactive...
        choice_obj = None
        print("Welcome to the SCEC Broadband Platform version %s." %
              (install.VERSION))
        print("=" * 80)
        print()
        print("Please select the Broadband Platform mode of operation:")
        print("   * Validation - Simulates a historical event")
        print("   * Scenario   - Runs a user-defined hypothetical event")
        print()
        while True:
            try:
                val_choice = input("Do you want to perform a "
                                   "validation simulation (y/n)? ")
                if (val_choice.lower() == 'y' or val_choice.lower() == 'n' or
                    val_choice.lower() == 'yes' or val_choice.lower() == 'no'):
                    break
                else:
                    print("Invalid input.")
            except KeyboardInterrupt:
                print("\nAborting...")
                sys.exit(1)

    expert_mode = False
    # Check for expert mode
    if options.expert == True:
        expert_mode = True

    # Initialize workflow builder class
    bband_workflow = build_workflow.WorkflowBuilder(sim_id,
                                                    expert_mode,
                                                    choice_obj)

    # Create validation or scenario simulation workflow
    try:
        if val_choice.lower() == 'y' or val_choice.lower() == 'yes':
            bband_workflow.do_validation()
        else:
            bband_workflow.do_scenario()
    except KeyboardInterrupt:
        print("\nAborting...")
        sys.exit(1)

    # Set workflow variable
    workflow_obj = xml_handler.Workflow(bband_workflow.workflow,
                                        bband_workflow.val_obj,
                                        bband_workflow.stations,
                                        install.VERSION)
    workflow = bband_workflow.workflow

if options.noxml is None:
    xml_handler.write_workflow(workflow_obj, sim_id)

if options.generate == True:
    print("XML file %s/%d.xml generated, not running workflow." %
          (install.A_XML_DIR, sim_id))
    sys.exit(0)

# If requested by user, direct stdout and stderr to logfile
if options.logFile is not None:
    log = Logger(options.logFile)
    sys.stdout = log
    sys.stderr = log

# Save system and software information
status = bbp_status.BBPStatus(sim_id=sim_id)
status.system()
status.software()

# Set the start time
install.set_start_time()

if options.noxml is None:
    # Before we start, let's copy the workflow xml description to the
    # output directory
    run_xml = os.path.join(install.A_XML_DIR, "%s.xml" % (sim_id))
    out_xml = os.path.join(a_outdir, "%s.xml" % (sim_id))
    shutil.copy2(run_xml, out_xml)

# Run xml Broadband workflow
last_module = None
for item in workflow:
    # Check if we need to stop at a certain module
    if options.end_module is not None:
        if last_module == options.end_module.lower():
            # We just completed the module we wanted, stop here!
            break
        else:
            # Update last module with the current one
            last_module = item.getName().lower()
    # Skip modules if resume_module is set
    if options.resume_module is not None:
        if item.getName().lower() != options.resume_module.lower():
            print("==> Skipping: %s" % (item.getName()))
            # Skip modules until we find the one we want
            continue
        # Ok, we got a match, run everything from now on...
        options.resume_module = None
    try:
        # print("Running %s" % (item.getName()))
        item.stage(a_indir)
        obj = item.instantiate(sim_id)
        obj.run()
    except Exception as err:
        fp = open("%s/fatal_error.log" % a_logdir, 'w')
        fp.write("%s\n" % err)
        fp.write("%s\n" % traceback.format_exc())
        fp.flush()
        fp.close()
        traceback.print_exc()
        sys.exit(-3)

# All done!
print()
print("SCEC Broadband Platform run completed.")
print("You can find results in %s" % (a_outdir))
