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

Program to set up a full BBP scenario simulation run on the Epicenter
cluster
"""
# Import Python modules
import os
import sys
import math
import random
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
MAX_SIMULATIONS = 200
CODEBASES = ["gp", "ucsb", "sdsu", "exsim", "csm", "irikura"]

def generate_src_files(numsim, source_file, srcdir, prefix, hypo_rand):
    """
    Generates num_sim source files in the srcdir using different
    random seeds
    """
    src_props = bband_utils.parse_properties(source_file)
    # Delete "seed" from the property set
    if "seed" in src_props:
        src_props.pop("seed")
    # Get FAULT_LENGTH and FAULT_WIDTH from the SRC file
    try:
        flen = float(src_props["fault_length"])
        fwid = float(src_props["fault_width"])
    except KeyError:
        raise bband_utils.ParameterError("Cannot read fault_length/fault_width"
                                         " parameters from SRC file!")
    if hypo_rand:
        # Delete HYPO_ALONG_STK and HYPO_DOWN_DIP
        if "hypo_along_stk" in src_props:
            src_props.pop("hypo_along_stk")
        if "hypo_down_dip" in src_props:
            src_props.pop("hypo_down_dip")
    # Create common list of keys for all files
    output = ""
    for key in src_props:
        output = output + "%s = %s\n" % (key.upper(), src_props[key])
    for sim in range(0, numsim):
        random.seed(sim + 1)
        seed = int(math.exp(7 * math.log(10.0)) * random.random())
        hypo_along_stk = flen * (0.2 + 0.6 * random.random() - 0.5)
        hypo_down_dip = fwid * (0.2 + 0.6 * random.random())
        srcfile = os.path.join(srcdir, "%s-%04d.src" % (prefix, sim))
        outfile = open(srcfile, 'w')
        outfile.write(output)
        if hypo_rand:
            outfile.write("HYPO_ALONG_STK = %.2f\n" % (hypo_along_stk))
            outfile.write("HYPO_DOWN_DIP = %.2f\n" % (hypo_down_dip))
        outfile.write("SEED = %d\n" % (seed))
        outfile.close()

def generate_xml(install, numsim, srcdir, xmldir,
                 logdir, vmodel, codebase, prefix,
                 station_list):
    """
    Generates xml files in the xmldir for numsim simulations whose
    source files are in the srcdir using the validation event and
    codebase specified
    """
    tmpdir = tempfile.mkdtemp(prefix="bbp-")
    bbproot = "%s/run_bbp.py" % (install.A_COMP_DIR)
    bfn = os.path.join(xmldir, BATCH_SIM_FILE)
    batchfile = open(bfn, 'w')
    for sim in range(0, numsim):
        srcfile = os.path.join(srcdir, "%s-%04d.src" % (prefix, sim))
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
        optfile.write('%s\n' % (station_list)) # Station list
        if codebase == "exsim":
            # Don't specify custom ExSIM template file
            optfile.write('n\n')
        if (codebase != "exsim" and
            codebase != "csm" and
            codebase != "irikura"):
            # Skip site response
            optfile.write('n\n')
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

def write_pbs(install, numsim, simdir, xmldir, email, prefix):
    """
    Write the pbs script
    """
    # Calculate how many nodes we need
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
    pbsfile.write("#PBS -l walltime=300:00:00,nodes=%d:ppn=%d\n" %
                  (nodes, CORES_PER_NODE))
    pbsfile.write("#PBS -V\n")
    pbsfile.write("#PBS -m abe -M %s\n" % (email))
    pbsfile.write("#PBS -e %s\n" % (errfile))
    pbsfile.write("#PBS -o %s\n" % (outfile))
    pbsfile.write("\n")
    pbsfile.write("BBP_DIR=%s\n" % (install.A_INSTALL_ROOT))
    pbsfile.write("PYTHONPATH=%s\n" % (install.A_COMP_DIR))
    pbsfile.write("BBP_DATA_DIR=/scratch/$PBS_JOBID/bbpruns\n")
    pbsfile.write("BBP_BASE_DIR=/scratch/$PBS_JOBID\n")
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
    pbsfile.write("python $BBP_DIR/utils/batch/run_parallel.py $BBP_DIR/utils/batch/setup_bbp_epicenter_env.sh %s $PBS_NODEFILE 1\n" %
                  (bfn))
    pbsfile.write("\n")
    pbsfile.write('echo "Processing end"\n')
    pbsfile.write("date\n")
    pbsfile.write("\n")
    for dir_to_copy in ['outdata', 'indata', 'logs', 'tmpdata']:
        pbsfile.write('python $BBP_DIR/utils/batch/command_parallel.py $BBP_DIR/utils/batch/setup_bbp_epicenter_env.sh "cp -frp $BBP_DATA_DIR/%s/* $HOME/Sims/%s/." $PBS_NODEFILE\n' %
                      (dir_to_copy, dir_to_copy))
    pbsfile.write('python $BBP_DIR/utils/batch/command_parallel.py $BBP_DIR/utils/batch/setup_bbp_epicenter_env.sh "rm -rf $BBP_BASE_DIR" $PBS_NODEFILE\n')
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
    parser.add_option("--hypo-rand", action="store_true", dest="hyporand",
                      help="Enables hypocenter randomization")
    parser.add_option("--no-hypo-rand", action="store_false", dest="hyporand",
                      help="Disables hypocenter randomization")
    parser.add_option("-n", "--num-simulations", type="int", action="store",
                      dest="numsim", help="Number of simulations to run")
    parser.add_option("--email", type="string", action="store",
                      dest="email", help="Email for job notifications")
    (options, _) = parser.parse_args()

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

    # Check for hypocenter randomization
    if options.hyporand is None:
        print "Please specify --hypo-rand or --no-hypo-rand!"
        sys.exit(1)

    if options.hyporand:
        hypo_rand = True
    else:
        hypo_rand = False

    # Get the source file
    source_file = options.source
    if source_file is None:
        print "Please provide a source description (src file)!"
        sys.exit(1)
    # Make it a full path
    source_file = os.path.realpath(source_file)
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
    numsim = options.numsim
    if numsim < 1 or numsim > MAX_SIMULATIONS:
        print ("Number of simulations should be between 1 and %d" %
               (MAX_SIMULATIONS))
        sys.exit(1)

    # Check for e-mail address
    email = options.email
    if email is None:
        print "Please provide an e-mail address for job notifications"
        sys.exit(1)

    # Make sure user has configured the setup_bbp_epicenter_env.sh script
    setup_bbp_env = os.path.join(bbp_install.A_INSTALL_ROOT,
                                 "utils/batch/setup_bbp_epicenter_env.sh")
    if not os.path.exists(setup_bbp_env):
        print ("Cannot find setup_bbp_epicenter_env.sh script!")
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
    for mdir in [indir, outdir, tmpdir, logsdir, xmldir, srcdir]:
        os.makedirs(mdir)
    # Generate source files
    generate_src_files(numsim, source_file, srcdir, prefix, hypo_rand)
    # Generate xml files
    generate_xml(bbp_install, numsim, srcdir, xmldir,
                 logsdir, vmodel, codebase, prefix,
                 station_list)
    # Write pbs file
    write_pbs(bbp_install, numsim, simdir, xmldir, email, prefix)

if __name__ == "__main__":
    main()
