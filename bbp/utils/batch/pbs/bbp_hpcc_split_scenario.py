#!/usr/bin/env python
"""
Copyright 2010-2017 University Of Southern California

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

Program to set up a full BBP scenario simulation run on HPCC
"""

# Import Python modules
import os
import sys
import math
import shutil
import optparse
import tempfile

# Import Broadband modules
from install_cfg import InstallCfg
import velocity_models
import bband_utils

# Constants
BATCH_SIM_FILE = "batch_run_bbp_sims.log"
CORES_PER_NODE = 8
CORES_PER_NODE_NEW = 16
MAX_SIMULATIONS = 9999
CODEBASES = ["gp", "ucsb", "sdsu", "exsim", "csm", "irikura"]

def write_sta_file(sta_list_fn, stations):
    """
    Writes "stations" to the sta_list_fn station list file
    """
    sta_out_file = open(sta_list_fn, 'w')
    for station in stations:
        sta_out_file.write("%s\n" % (station))
    sta_out_file.close()

def generate_stl_files(station_list, numsta, stldir):
    """
    Generates the station files using at most numsta stations per
    simulation. Returns the number of simulations needed to run all
    stations.
    """
    numsim = 0
    stations = []
    sta_list_base = os.path.splitext(os.path.basename(station_list))[0]
    sta_list_base = os.path.join(stldir, sta_list_base)
    sta_file = open(station_list)
    for line in sta_file:
        line = line.strip()
        # Skip blank lines
        if not line:
            continue
        # Skip comments
        if line.startswith("#") or line.startswith("%"):
            continue
        # Add this station to the list
        stations.append(line)
        if len(stations) == numsta:
            # Write this station list
            sta_list_fn = "%s-%04d.stl" % (sta_list_base, numsim)
            write_sta_file(sta_list_fn, stations)
            # Start again, empty list and increament numsim
            numsim = numsim + 1
            stations = []
    sta_file.close()
    # Make sure we don't have an extra file to write
    if len(stations) > 0:
        # Write this station list
        sta_list_fn = "%s-%04d.stl" % (sta_list_base, numsim)
        write_sta_file(sta_list_fn, stations)
        numsim = numsim + 1

    # Return number of simulations
    return numsim, sta_list_base

def generate_src_files(numsim, source_file, srcdir, prefix, new_seed=None):
    """
    Generates num_sim source files in the srcdir using different
    random seeds
    """
    src_props = bband_utils.parse_properties(source_file)
    # Delete "seed" from the property set
    if "seed" in src_props:
        # Keep track of SRC file seed value
        seed = int(src_props["seed"])
        src_props.pop("seed")
        # But use new_seed if provided by the user
        if new_seed is not None:
            seed = new_seed
    else:
        if new_seed is None:
            raise bband_utils.ParameterError("Please specify a seed for"
                                             " this simulation!")
        seed = new_seed

    # Create common list of keys for all files
    output = ""
    for key in src_props:
        output = output + "%s = %s\n" % (key.upper(), src_props[key])
    for sim in range(0, numsim):
        srcfile = os.path.join(srcdir, "%s-%04d.src" % (prefix, sim))
        outfile = open(srcfile, 'w')
        outfile.write(output)
        outfile.write("SEED = %d\n" % (seed))
        outfile.close()

def generate_xml(install, numsim, srcdir,
                 xmldir, logdir, vmodel,
                 codebase, prefix, stlbase,
                 nosite):
    """
    Generates xml files in the xmldir for numsim simulations whose
    source files are in the srcdir and station lists are in stldir
    using the velocity model and codebase specified
    """
    tmpdir = tempfile.mkdtemp(prefix="bbp-")
    bbproot = "%s/run_bbp.py" % (install.A_COMP_DIR)
    bfn = os.path.join(xmldir, BATCH_SIM_FILE)
    batchfile = open(bfn, 'w')
    for sim in range(0, numsim):
        srcfile = os.path.join(srcdir, "%s-%04d.src" % (prefix, sim))
        stlfile = "%s-%04d.stl" % (stlbase, sim)
        ofn = os.path.join(tmpdir, "bbp.optfile")
        optfile = open(ofn, 'w')
        optfile.write('n\n') # Scenario
        optfile.write('%s\n' % (vmodel)) # Velocity model
        optfile.write('%s\n' % (codebase)) # Codebase to use
        if codebase != "exsim" and codebase != "csm":
            optfile.write('y\n') # Run rupture generator
        optfile.write('2\n') # Enter path to source file
        optfile.write('%s\n' % (srcfile)) # Source file
        optfile.write('2\n') # Enter path to station list
        optfile.write('%s\n' % (stlfile)) # Station list
        if codebase == "exsim":
            # Don't specify custom ExSIM template file
            optfile.write('n\n')
        # Select site response option
        if (codebase != "exsim" and
            codebase != "csm" and
            codebase != "irikura"):
            if nosite:
                # Skip site response
                optfile.write('n\n')
            else:
                # Run site response
                optfile.write('y\n')
        optfile.write('y\n') # Plot velocity seismograms
        optfile.write('y\n') # Plot acceleration seismograms
        optfile.flush()
        optfile.close()
        # Run BBP and generate the xml file
        bband_utils.runprog("export BBP_DATA_DIR=%s; %s -o %s -s %d -g" %
                            (tmpdir, bbproot, ofn, sim+1))
        # Copy the xml file
        srcxml = os.path.join(tmpdir, "xml", "%d.xml" % (sim+1))
        dstxml = os.path.join(xmldir, "%s-%04d.xml" % (prefix, sim))
        shutil.copy2(srcxml, dstxml)
        # Add entry to the batch file
        bbp_sim = 10000000 + sim
        logbase = os.path.join(logdir, str(bbp_sim))
        logfile = os.path.join(logbase, "%d_%s.log" %
                               (bbp_sim, prefix))
        # Make sure logdir exists
        os.makedirs(logbase)
        batchfile.write("%s -x %s -s %d -l %s\n" %
                        (bbproot, dstxml, bbp_sim, logfile))
    # Close batch file
    batchfile.flush()
    batchfile.close()
    # Clean-up
    shutil.rmtree(tmpdir)

def write_pbs(install, numsim, simdir, xmldir, email, prefix, newnodes):
    """
    Write the pbs script
    """
    # Calculate how many nodes we need
    if newnodes:
        nodes = int(math.ceil(1.0 * numsim / CORES_PER_NODE_NEW))
    else:
        nodes = int(math.ceil(1.0 * numsim / CORES_PER_NODE))
    # Some path names
    outfile = os.path.join(simdir, "%s.out" % (prefix))
    errfile = os.path.join(simdir, "%s.err" % (prefix))
    bfn = os.path.join(xmldir, BATCH_SIM_FILE)
    # Let's create the pbs file
    pbsfn = os.path.join(simdir, "%s.pbs" % (prefix))
    pbsfile = open(pbsfn, 'w')
    pbsfile.write("#!/bin/bash\n")
    pbsfile.write("\n")
    if not newnodes:
        pbsfile.write("#PBS -q nbns\n")
        pbsfile.write("#PBS -l arch=x86_64,pmem=1000mb,pvmem=2000mb,walltime=300:00:00,nodes=%d:ppn=%d\n" %
                      (nodes, CORES_PER_NODE))
    else:
        pbsfile.write("#PBS -l arch=x86_64,pmem=1000mb,pvmem=2000mb,walltime=24:00:00,nodes=%d:ppn=%d:gpus=2\n" %
                      (nodes, CORES_PER_NODE_NEW))
    pbsfile.write("#PBS -V\n")
    pbsfile.write("#PBS -m abe -M %s\n" % (email))
    pbsfile.write("#PBS -e %s\n" % (errfile))
    pbsfile.write("#PBS -o %s\n" % (outfile))
    pbsfile.write("\n")
    pbsfile.write("BBP_DIR=%s\n" % (install.A_INSTALL_ROOT))
    pbsfile.write("PYTHONPATH=%s\n" % (install.A_COMP_DIR))
    pbsfile.write("BBP_DATA_DIR=$TMPDIR/bbpruns\n")
    pbsfile.write("BBP_BASE_DIR=$TMPDIR\n")
    pbsfile.write("HOME=%s\n" % (simdir))
    pbsfile.write("\n")
    pbsfile.write("mkdir -p $BBP_DATA_DIR\n")
    pbsfile.write("mkdir -p $HOME/Sims/indata\n")
    pbsfile.write("mkdir -p $HOME/Sims/logs\n")
    pbsfile.write("mkdir -p $HOME/Sims/outdata\n")
    pbsfile.write("mkdir -p $HOME/Sims/tmpdata\n")
    pbsfile.write("\n")
    pbsfile.write('echo "Jobs start"\n')
    pbsfile.write("date\n")
    pbsfile.write('echo "BBP_DATA_DIR = $BBP_DATA_DIR"\n')
    pbsfile.write("\n")
    pbsfile.write("cd $HOME\n")
    pbsfile.write("\n")
    pbsfile.write("python $BBP_DIR/utils/batch/run_parallel.py $BBP_DIR/utils/batch/setup_bbp_env.sh %s $PBS_NODEFILE 1\n" %
                  (bfn))
    pbsfile.write("\n")
    pbsfile.write('echo "Processing end"\n')
    pbsfile.write("date\n")
    pbsfile.write("\n")
    for dir_to_copy in ['outdata', 'indata', 'logs']:
        pbsfile.write('python $BBP_DIR/utils/batch/command_parallel.py $BBP_DIR/utils/batch/setup_bbp_env.sh "cp -frp $BBP_DATA_DIR/%s/* $HOME/Sims/%s/." $PBS_NODEFILE\n' %
                      (dir_to_copy, dir_to_copy))
    pbsfile.write("\n")
    pbsfile.write('echo "Jobs end"\n')
    pbsfile.write("date\n")
    pbsfile.flush()
    pbsfile.close()

    # All done!
    print
    print "Validation run is set up on: %s" % (simdir)
    print
    print "To start the validation run, just type: "
    print "$ qsub %s" % (pbsfn)
    print
    if newnodes:
        print "Please note that the maximum walltime has been set to 24 hours!"
        print "Jobs running longer than that will be terminated at 24 hours!"
        print

def main():
    """
    Parse command line options and create the needed files/directories
    """
    # Detect BBP installation
    bbp_install = InstallCfg.getInstance()

    prog_base = os.path.basename(sys.argv[0])
    usage = "usage: %s [options]" % (prog_base)
    parser = optparse.OptionParser(usage)
    parser.add_option("-c", "--codebase", type="string", action="store",
                      dest="codebase",
                      help="Codebase for the simulation: %s" %
                      (CODEBASES))
    parser.add_option("-v", "--velocity-model", type="string", action="store",
                      dest="vmodel",
                      help="Velocity model (region) for this simulation")
    parser.add_option("--src", "--source", type="string", action="store",
                      dest="source",
                      help="Source description file for the simulation")
    parser.add_option("--stl", "--station-list", type="string", action="store",
                      dest="station_list",
                      help="Station list file for the simulation")
    parser.add_option("-d", "--dir", type="string", action="store",
                      dest="simdir",
                      help="Simulation directory")
    parser.add_option("-n", "--num-stations", type="int", action="store",
                      dest="numsta", help="Number of stations per run")
    parser.add_option("--seed", type="int", action="store",
                      dest="new_seed", help="Overrides seed in SRC file")
    parser.add_option("--email", type="string", action="store",
                      dest="email", help="Email for job notifications")
    parser.add_option("--new-nodes", action="store_true", dest="newnodes",
                      help="Schedule the job in the new HPCC nodes")
    parser.add_option("--no-site-response", action="store_true", dest="nosite",
                      help="Disable the site response module for GP/SDSU/UCSB")
    (options, _) = parser.parse_args()

    # Check if using new HPCC nodes
    if options.newnodes:
        newnodes = True
    else:
        newnodes = False

    # Check if not using site response
    if options.nosite:
        nosite = True
    else:
        nosite = False

    # Validate codebase to use
    codebase = options.codebase
    if codebase is None:
        print "Please specify a codebase!"
        sys.exit(1)
    codebase = codebase.lower()
    if codebase not in CODEBASES:
        print "Codebase needs to be one of: %s" % (CODEBASES)

    # Check for velocity model
    vmodel_names = velocity_models.get_all_names()
    vmodel = options.vmodel
    if vmodel is None:
        print "Please provide a velocity model (region) for this simulation!"
        print "Available options are: %s" % (vmodel_names)
        sys.exit(1)
    vmodels = [v_model.lower() for v_model in vmodel_names]
    if vmodel.lower() not in vmodels:
        print ("Velocity model %s does not appear to be available on BBP" %
               (vmodel))
        print ("Available options are: %s" % (vmodel_names))
        print "Please provide another velocity model or check your BBP installation."
        sys.exit(1)
    # Now get the name with the correct case
    vmodel = vmodel_names[vmodels.index(vmodel.lower())]

    # Get the source file
    source_file = options.source
    if source_file is None:
        print "Please provide a source description (src file)!"
        sys.exit(1)
    # Make it a full path
    source_file = os.path.realpath(source_file)
    # Make sure source file is in the rcf-104 filesystem
    if not "rcf-104" in source_file:
        print "Source file should be in the rcf-104 filesystem!"
        sys.exit(1)
    # Make sure source file exists and is readable
    if not os.path.isfile(source_file) or not os.access(source_file, os.R_OK):
        print "Source file does not seem to be accessible!"
        sys.exit(1)

    # Get the station list
    station_list = options.station_list
    if station_list is None:
        print "Please provide a station list (stl file)!"
        sys.exit(1)
    # Make it a full path
    station_list = os.path.realpath(station_list)
    # Make sure station list is in the rcf-104 filesystem
    if not "rcf-104" in station_list:
        print "Station list should be in the rcf-104 filesystem!"
        sys.exit(1)
    # Make sure station list exists and is readable
    if not os.path.isfile(station_list) or not os.access(station_list, os.R_OK):
        print "Station list foes not seem to be accessible!"
        sys.exit(1)

    # Check for the simulation directory
    simdir = options.simdir
    if simdir is None:
        print "Please provide a simulation directory!"
        sys.exit(1)
    simdir = os.path.abspath(simdir)
    if os.path.exists(simdir):
        print "Simulation directory exists: %s" % (simdir)
        opt = raw_input("Do you want to delete its contents (y/n)? ")
        if opt.lower() != "y":
            print "Please provide another simulation directory!"
            sys.exit(1)
        opt = raw_input("ARE YOU SURE (y/n)? ")
        if opt.lower() != "y":
            print "Please provide another simulation directory!"
            sys.exit(1)
        # Delete existing directory (we already asked the user twice!!!)
        shutil.rmtree(simdir)

    # Pick up number of simulations to run
    numsta = options.numsta
    if numsta < 1:
        print ("Number of stations should be greater than 0")
        sys.exit(1)

    # Check for user-provided seed for this simulation
    new_seed = options.new_seed

    # Check for e-mail address
    email = options.email
    if email is None:
        print "Please provide an e-mail address for job notifications"
        sys.exit(1)

    # Make sure user has configured the setup_bbp_env.sh script
    setup_bbp_env = os.path.join(bbp_install.A_INSTALL_ROOT,
                                 "utils/batch/setup_bbp_env.sh")
    if not os.path.exists(setup_bbp_env):
        print ("Cannot find setup_bbp_env.sh script!")
        print ("Expected at: %s" % (setup_bbp_env))
        sys.exit(1)
    # Create simulation directories
    prefix = "%s-%s" % (os.path.splitext(os.path.basename(source_file))[0],
                        codebase.lower())
    # Make sure we remove spaces from prefix
    prefix = prefix.replace(" ", '')
    os.makedirs(simdir)
    indir = os.path.join(simdir, "Sims", "indata")
    outdir = os.path.join(simdir, "Sims", "outdata")
    tmpdir = os.path.join(simdir, "Sims", "tmpdata")
    logsdir = os.path.join(simdir, "Sims", "logs")
    xmldir = os.path.join(simdir, "Xml")
    srcdir = os.path.join(simdir, "Src")
    stldir = os.path.join(simdir, "Stl")
    for mdir in [indir, outdir, tmpdir, logsdir, xmldir, srcdir, stldir]:
        os.makedirs(mdir)
    # Generate station lists
    numsim, stlbase = generate_stl_files(station_list, numsta, stldir)
    if numsim > MAX_SIMULATIONS:
        print "Too many simulations requested!"
        print "Maximum number allowed is %d!" % (MAX_SIMULATIONS)
        print "Try requesting more stations per simulation..."
        sys.exit(1)
    # Generate source files
    generate_src_files(numsim, source_file, srcdir, prefix, new_seed)
    # Generate xml files
    generate_xml(bbp_install, numsim, srcdir,
                 xmldir, logsdir, vmodel,
                 codebase, prefix, stlbase,
                 nosite)
    # Write pbs file
    write_pbs(bbp_install, numsim, simdir, xmldir, email, prefix, newnodes)

if __name__ == "__main__":
    main()
