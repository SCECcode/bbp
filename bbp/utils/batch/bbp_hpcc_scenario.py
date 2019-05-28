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

Program to set up a full BBP scenario simulation run on HPCC
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
CORES_PER_NODE_NEW = 16
MAX_SIMULATIONS = 200
CODEBASES = ["gp", "ucsb", "sdsu", "exsim", "csm", "song", "irikura"]
CODEBASES_SITE = ["gp", "sdsu", "song", "irikura", "exsim", "ucsb"]

def generate_src_files(numsim, source_file, srcdir,
                       prefix, hypo_rand, hypo_area,
                       variation, multiseg, first_seg_dir):
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
    # Figure out are where hypocenter should go, by default we want to
    # cut 20% in each side and place the hypocenter in the center area
    if hypo_area["hdd_min"] is None:
        hypo_area["hdd_min"] = 0.2 * fwid
    if hypo_area["hdd_max"] is None:
        hypo_area["hdd_max"] = 0.8 * fwid
    if hypo_area["has_min"] is None:
        hypo_area["has_min"] = -(0.5 - 0.2) * flen
    if hypo_area["has_max"] is None:
        hypo_area["has_max"] = (0.5 - 0.2) * flen
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
    common_seeds = []
    # Check if we are doing a multi-segment run
    if multiseg and first_seg_dir is not None:
        # Read common seeds from seed file
        seed_file = open(os.path.join(first_seg_dir, "Src", "seeds.txt"), 'r')
        first_seg_sims = int(seed_file.readline().strip())
        if first_seg_sims != numsim:
            print("ERROR: Number of simulations must match across segments!")
            sys.exit(1)
        for line in seed_file:
            common_seeds.append(int(line.strip()))
        seed_file.close()
    # Generate the numsim SRC files
    all_seeds = []
    for sim in range(0, numsim):
        random.seed((sim + 1) + (variation - 1) * 500)
        seed = int(math.exp(7 * math.log(10.0)) * random.random())
        all_seeds.append(seed)
        hypo_along_stk = ((hypo_area["has_max"] - hypo_area["has_min"]) *
                          random.random()) + hypo_area["has_min"]
        hypo_down_dip = ((hypo_area["hdd_max"] - hypo_area["hdd_min"]) *
                         random.random()) + hypo_area["hdd_min"]
        srcfile = os.path.join(srcdir, "%s-%04d.src" % (prefix, sim))
        outfile = open(srcfile, 'w')
        outfile.write(output)
        if hypo_rand:
            outfile.write("HYPO_ALONG_STK = %.2f\n" % (hypo_along_stk))
            outfile.write("HYPO_DOWN_DIP = %.2f\n" % (hypo_down_dip))
        outfile.write("SEED = %d\n" % (seed))
        if multiseg and first_seg_dir is not None:
            outfile.write("COMMON_SEED = %d\n" % (common_seeds[sim]))
        outfile.close()
    # Check if we need to write file with all seeds
    if multiseg and first_seg_dir is None:
        # This is the first segment, write seeds file
        seed_file = open(os.path.join(srcdir, "seeds.txt"), 'w')
        seed_file.write("%d\n" % (numsim))
        for seed in all_seeds:
            seed_file.write("%d\n" % (seed))
        seed_file.close()

def generate_xml(install, numsim, srcdir, xmldir,
                 logdir, vmodel, codebase, prefix,
                 station_list, only_rup, srf_prefix,
                 site_response):
    """
    Generates xml files in the xmldir for numsim simulations whose
    source files are in the srcdir using the validation event and
    codebase specified
    """
    tmpdir = tempfile.mkdtemp(prefix="bbp-")
    bbproot = os.path.join(install.A_COMP_DIR, "run_bbp.py")
    bfn = os.path.join(xmldir, BATCH_SIM_FILE)
    batchfile = open(bfn, 'w')
    for sim in range(0, numsim):
        if srf_prefix is None:
            srcfile = os.path.join(srcdir, "%s-%04d.src" % (prefix, sim))
        else:
            srffile = "%s-%04d.srf" % (srf_prefix, sim)
            # Make sure srf file exists and is readable
            if (not os.path.isfile(srffile) or
                not os.access(srffile, os.R_OK)):
                print "SRF file %s does not seem to be accessible!" % (srffile)
                sys.exit(1)
        ofn = os.path.join(tmpdir, "bbp.optfile")
        optfile = open(ofn, 'w')
        optfile.write('n\n') # Scenario
        optfile.write('%s\n' % (vmodel)) # Velocity model
        optfile.write('%s\n' % (codebase)) # Codebase to use
        optfile.write('2\n') # Enter path to src file
        optfile.write('%s\n' % (srcfile)) # Source file
        if codebase != "exsim" and codebase != "csm":
            if srf_prefix is None:
                optfile.write('y\n') # Run rupture generator
            else:
                optfile.write('n\n') # Skip rupture generator
                optfile.write('2\n') # Enter path to srf file
                optfile.write('%s\n' % (srffile)) # SRF file
        optfile.write('2\n') # Enter path to station list
        optfile.write('%s\n' % (station_list)) # Station list
        if codebase == "exsim":
            # Don't specify custom ExSIM template file
            optfile.write('n\n')
        if codebase != "csm":
            if site_response:
                # Use site response
                optfile.write('y\n')
            else:
                # Skip site response
                optfile.write('n\n')
        optfile.write('y\n') # Plot velocity seismograms
        optfile.write('y\n') # Plot acceleration seismograms
        optfile.flush()
        optfile.close()
        # Run BBP and generate the xml file
        bband_utils.runprog("export BBP_DATA_DIR=%s; "
                            "%s --expert -o %s -s %d -g" %
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
        # Write this run in the batch file
        if only_rup:
            if codebase == "gp" or codebase == "sdsu":
                end_module = " -e Genslip "
            elif codebase == "ucsb":
                end_module = " -e UCrmg "
        else:
            end_module = ""
        batchfile.write("%s -x %s -s %d %s-l %s\n" %
                        (bbproot, dstxml, bbp_sim, end_module, logfile))
    # Close batch file
    batchfile.flush()
    batchfile.close()
    # Clean-up
    shutil.rmtree(tmpdir)

def write_slurm(install, numsim, simdir, xmldir, email,
                prefix, newnodes, walltime, savetemp,
                codebase):
    """
    Write the slurm script
    """
    # Calculate how many nodes we need
    if newnodes:
        nodes = int(math.ceil(1.0 * numsim / CORES_PER_NODE_NEW))
        cores_per_node = CORES_PER_NODE_NEW
    else:
        nodes = int(math.ceil(1.0 * numsim / CORES_PER_NODE))
        cores_per_node = CORES_PER_NODE
    # Some path names
    outfile = os.path.join(simdir, "%s.out" % (prefix))
    errfile = os.path.join(simdir, "%s.err" % (prefix))
    bfn = os.path.join(xmldir, BATCH_SIM_FILE)
    # Let's create the pbs file
    slurmfn = os.path.join(simdir, "%s.slurm" % (prefix))
    slurmfile = open(slurmfn, 'w')
    slurmfile.write("#!/bin/bash\n")
    slurmfile.write("\n")
    if not newnodes:
        slurmfile.write("#SBATCH --partition=scec\n")
    else:
        #slurmfile.write("#SBATCH --partition=main\n")
        slurmfile.write("#SBATCH --gres=gpu:2\n")
    #slurmfile.write("#SBATCH --mem-per-cpu=2GB\n")
    slurmfile.write("#SBATCH --mem=0\n")
    slurmfile.write("#SBATCH --nodes=%d\n" % (nodes))
    slurmfile.write("#SBATCH --time=%d:00:00\n" % (walltime))
    slurmfile.write("#SBATCH --export=all\n")
    slurmfile.write('#SBATCH --job-name="BBP-%s"\n' % (codebase))
    slurmfile.write("#SBATCH --mail-user=%s\n" % (email))
    slurmfile.write("#SBATCH --mail-type=BEGIN,END,ALL\n")
    slurmfile.write("#SBATCH --error=%s\n" % (errfile))
    slurmfile.write("#SBATCH --output=%s\n" % (outfile))
    slurmfile.write("\n")
    slurmfile.write("BBP_DIR=%s\n" % (install.A_INSTALL_ROOT))
    slurmfile.write("PYTHONPATH=%s\n" % (install.A_COMP_DIR))
    slurmfile.write("BBP_DATA_DIR=$TMPDIR/bbpruns\n")
    slurmfile.write("BBP_BASE_DIR=$TMPDIR\n")
    slurmfile.write("HOME=%s\n" % (simdir))
    slurmfile.write("SLURM_NODES=`scontrol show hostname $SLURM_JOB_NODELIST | paste -d, -s`\n")
    slurmfile.write("\n")
    slurmfile.write("mkdir -p $BBP_DATA_DIR\n")
    slurmfile.write("mkdir -p $HOME/Sims/indata\n")
    slurmfile.write("mkdir -p $HOME/Sims/logs\n")
    slurmfile.write("mkdir -p $HOME/Sims/outdata\n")
    slurmfile.write("mkdir -p $HOME/Sims/tmpdata\n")
    slurmfile.write("\n")
    slurmfile.write('echo "Jobs start"\n')
    slurmfile.write("date\n")
    slurmfile.write('echo "BBP_DATA_DIR = $BBP_DATA_DIR"\n')
    slurmfile.write("\n")
    slurmfile.write("cd $HOME\n")
    slurmfile.write("\n")
    slurmfile.write("python $BBP_DIR/utils/batch/run_parallel.py $BBP_DIR/utils/batch/setup_bbp_env.sh %s $SLURM_NODES %d\n" %
                    (bfn, cores_per_node))
    slurmfile.write("\n")
    slurmfile.write('echo "Processing end"\n')
    slurmfile.write("date\n")
    slurmfile.write("\n")
    if savetemp:
        for dir_to_copy in ['outdata', 'indata', 'logs', 'tmpdata']:
            slurmfile.write('python $BBP_DIR/utils/batch/command_parallel.py $BBP_DIR/utils/batch/setup_bbp_env.sh "cp -frp $BBP_DATA_DIR/%s/* $HOME/Sims/%s/." $SLURM_NODES\n' %
                            (dir_to_copy, dir_to_copy))
    else:
        for dir_to_copy in ['outdata', 'indata', 'logs']:
            slurmfile.write('python $BBP_DIR/utils/batch/command_parallel.py $BBP_DIR/utils/batch/setup_bbp_env.sh "cp -frp $BBP_DATA_DIR/%s/* $HOME/Sims/%s/." $SLURM_NODES\n' %
                            (dir_to_copy, dir_to_copy))
    slurmfile.write("\n")
    slurmfile.write('echo "Jobs end"\n')
    slurmfile.write("date\n")
    slurmfile.flush()
    slurmfile.close()

    # All done!
    print
    print "Validation run is set up on: %s" % (simdir)
    print
    print "To start the validation run, just type: "
    print "$ sbatch %s" % (slurmfn)
    print
    if newnodes:
        print "Please note that the maximum walltime has been set to %d hours!" % (walltime)
        print "Jobs running longer than that will be terminated at %d hours!" % (walltime)
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
    parser.add_option("--srf", "--srf-prefix", type="string", action="store",
                      dest="srf_prefix",
                      help="Prefix of SRF files to use, "
                      "only for GP, SDSU and UCSB methods. "
                      "Simulations begin after the rupture generator.")
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
    parser.add_option("-w", "--walltime", type="int", action="store",
                      dest="walltime", help="Number of hours for walltime")
    parser.add_option("--new-nodes", action="store_true", dest="newnodes",
                      help="Schedule the job in the new HPCC nodes")
    parser.add_option("--save-tmpdata", action="store_true", dest="savetemp",
                      help="Save the contents of the tmpdata directory")
    parser.add_option("--hdd-min", type="float",
                      action="store", dest="hdd_min",
                      help="Min value for hypo down dip in randomization")
    parser.add_option("--hdd-max", type="float",
                      action="store", dest="hdd_max",
                      help="Max value for hypo down dip in randomization")
    parser.add_option("--has-min", type="float",
                      action="store", dest="has_min",
                      help="Min value for hypo along strike in randomization")
    parser.add_option("--has-max", type="float",
                      action="store", dest="has_max",
                      help="Max value for hypo along strike in randomization")
    parser.add_option("--only-rup", action="store_true", dest="only_rup",
                      help="Only runs the rupture generator")
    parser.add_option("--var", "--variation", type="int", action="store",
                      dest="variation", help="seed variation (default 1)")
    parser.add_option("--multiseg", action="store_true", dest="multiseg",
                      help="Indicates simulation part of multiseg run")
    parser.add_option("--first-seg-dir", type="string", action="store",
                      dest="first_seg_dir",
                      help="required for multi-segment segments 2..n")
    parser.add_option("-s", "--site", action="store_true",
                      dest="site_response", help="Use site response module")

    (options, _) = parser.parse_args()

    # Check if using new HPCC nodes
    if options.newnodes:
        newnodes = True
    else:
        newnodes = False

    # Check if multi-segment simulation
    if options.multiseg:
        multiseg = True
    else:
        multiseg = False

    # Check for first segment directory
    if options.first_seg_dir is not None:
        first_seg_dir = os.path.abspath(options.first_seg_dir)
        if not os.path.exists(first_seg_dir):
            print("First segment directory for exists: %s" %
                  (first_seg_dir))
            sys.exit(1)
    else:
        first_seg_dir = None

    # Check if user specified custom walltime
    if options.walltime:
        if options.walltime < 1:
            print("Walltime must be at least 1 hour!")
            sys.exit(1)
        walltime = options.walltime
    else:
        if newnodes:
            walltime = 24
        else:
            walltime = 300

    # Check for variation sequence
    if options.variation:
        variation = options.variation
    else:
        variation = 1

    # Check if user wants to save the contents of tmpdata
    if options.savetemp:
        savetemp = True
    else:
        savetemp = False

    # Check if user wants to only run the rupture generator
    if options.only_rup:
        only_rup = True
    else:
        only_rup = False

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

    # Check if users wants to run site response module
    if options.site_response:
        site_response = True
        if codebase not in CODEBASES_SITE:
            print "Cannot use site response with method: %s" % (codebase)
            sys.exit(1)
    else:
        site_response = False

    # Check for hypocenter randomization
    if options.hyporand is None:
        print "Please specify --hypo-rand or --no-hypo-rand!"
        sys.exit(1)

    if options.hyporand:
        hypo_rand = True
    else:
        hypo_rand = False

    # Define area where hypocenter will be randomized
    hypo_area = {}
    hypo_area["hdd_min"] = options.hdd_min
    hypo_area["hdd_max"] = options.hdd_max
    hypo_area["has_min"] = options.has_min
    hypo_area["has_max"] = options.has_max

    # Get the source file (SRC or SRFs)
    source_file = options.source
    srf_prefix = options.srf_prefix

    if source_file is None and srf_prefix is None:
        print ("Please provide either source description (src file) "
               "or a srf prefix!")
        sys.exit(1)
    # If user specified both a source file and a srf prefix, we abort!
    if source_file is not None and srf_prefix is not None:
        print "Cannot specify both srf_prefic and source_file!"
        sys.exit(1)
    # If user specified a source file
    if source_file is not None:
        # Make it a full path
        source_file = os.path.realpath(source_file)
        # Make sure source file is in the rcf-104 / scec-00 filesystem
        if not "rcf-104" in source_file and not "scec-00" in source_file:
            print "Source file should be in the rcf-104 / scec-00 filesystems!"
            sys.exit(1)
        # Make sure source file exists and is readable
            if (not os.path.isfile(source_file) or
                not os.access(source_file, os.R_OK)):
                print "Source file does not seem to be accessible!"
                sys.exit(1)
        # Create a prefix
        prefix = ("%s-%s" %
                  (os.path.splitext(os.path.basename(source_file))[0],
                   codebase.lower()))
    # If user specified a SRF prefix
    if srf_prefix is not None:
        # Make it a full path
        srf_prefix = os.path.realpath(srf_prefix)
        # Make sure source file is in the rcf-104 or scec-00 filesystems
        if not "rcf-104" in srf_prefix and not "scec-00" in srf_prefix:
            print "SRF files should be in the rcf-104 / scec-00 filesystems!"
            sys.exit(1)
        # Create a prefix
        prefix = os.path.splitext(os.path.basename(srf_prefix))[0]
        #prefix = prefix.rsplit("-", 1)[0]
    # Make sure we remove spaces from prefix
    prefix = prefix.replace(" ", '')

    # Get the station list
    station_list = options.station_list
    if station_list is None:
        print "Please provide a station list (stl file)!"
        sys.exit(1)
    # Make it a full path
    station_list = os.path.realpath(station_list)
    # Make sure station list is in the rcf-104 or scec-00 filesystems
    if not "rcf-104" in station_list and not "scec-00" in station_list:
        print "Station list should be in the rcf-104 / scec-00 filesystems!"
        sys.exit(1)
    # Make sure station list exists and is readable
    if (not os.path.isfile(station_list) or
        not os.access(station_list, os.R_OK)):
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

    # Make sure user has configured the setup_bbp_env.sh script
    setup_bbp_env = os.path.join(bbp_install.A_INSTALL_ROOT,
                                 "utils/batch/setup_bbp_env.sh")
    if not os.path.exists(setup_bbp_env):
        print ("Cannot find setup_bbp_env.sh script!")
        print ("Expected at: %s" % (setup_bbp_env))
        sys.exit(1)
    # Create simulation directories
    os.makedirs(simdir)
    indir = os.path.join(simdir, "Sims", "indata")
    outdir = os.path.join(simdir, "Sims", "outdata")
    tmpdir = os.path.join(simdir, "Sims", "tmpdata")
    logsdir = os.path.join(simdir, "Sims", "logs")
    xmldir = os.path.join(simdir, "Xml")
    srcdir = os.path.join(simdir, "Src")
    for mdir in [indir, outdir, tmpdir, logsdir, xmldir, srcdir]:
        os.makedirs(mdir)
    if srf_prefix is None:
        # Generate source files
        generate_src_files(numsim, source_file, srcdir,
                           prefix, hypo_rand, hypo_area,
                           variation, multiseg, first_seg_dir)
    # Generate xml files
    generate_xml(bbp_install, numsim, srcdir, xmldir,
                 logsdir, vmodel, codebase, prefix,
                 station_list, only_rup, srf_prefix,
                 site_response)
    # Write pbs file
    write_slurm(bbp_install, numsim, simdir, xmldir,
                email, prefix, newnodes, walltime,
                savetemp, codebase)

    # Write .info file
    info_file = open(os.path.join(simdir, "%s.info" % (prefix)), 'w')
    info_file.write("# %s\n" % (" ".join(sys.argv)))
    info_file.close()

if __name__ == "__main__":
    main()
