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

Program to set up a full validation run on HPCC
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
import validation_cfg
import bband_utils
import gmpe_config

# Constants
BATCH_SIM_FILE = "batch_run_bbp_sims.log"
CORES_PER_NODE = 16
CORES_PER_NODE_NEW = 16
MAX_SIMULATIONS = 200
CODEBASES = ["gp", "ucsb", "sdsu", "exsim", "csm", "song", "irikura"]
CODEBASES_SITE = ["gp", "sdsu", "song", "irikura", "exsim", "ucsb"]
CODEBASES_SRF = ["gp", "sdsu", "song", "ucsb"]

def generate_src_files(numsim, source_file, srcdir,
                       prefix, hypo_rand,
                       variation, multiseg,
                       segment, first_seg_dir):
    """
    Generates num_sim source files in the srcdir using different
    random seeds
    """
    src_props = bband_utils.parse_properties(source_file)
    # Delete "seed" and "common_seed" from the property set
    if "seed" in src_props:
        src_props.pop("seed")
    if "common_seed" in src_props:
        src_props.pop("common_seed")
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
        hypo_along_stk = flen * (0.2 + 0.6 * random.random() - 0.5)
        hypo_down_dip = fwid * (0.2 + 0.6 * random.random())
        if multiseg:
            srcfile = os.path.join(srcdir,
                                   "%s-%04d_seg%02d.src" %
                                   (prefix, sim, segment))
        else:
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
                 logdir, event, codebase, prefix,
                 skip_rupgen, srf_prefix, only_rup,
                 gmpe_group_name, allmetrics,
                 site_response, multiseg, segment):
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
        if multiseg:
            srcfile = os.path.join(srcdir,
                                   "%s-%04d_seg%02d.src" %
                                   (prefix, sim, segment))
        else:
            srcfile = os.path.join(srcdir, "%s-%04d.src" % (prefix, sim))
        if srf_prefix is not None:
            srffile = "%s-%04d.srf" % (srf_prefix, sim)
            # Make sure srf file exists and is readable
            if (not os.path.isfile(srffile) or
                not os.access(srffile, os.R_OK)):
                print("SRF file %s does not seem to be accessible!" % (srffile))
                sys.exit(1)
        ofn = os.path.join(tmpdir, "bbp.optfile")
        optfile = open(ofn, 'w')
        optfile.write('y\n') # Validation
        optfile.write('%s\n' % (event)) # Validation event
        optfile.write('%s\n' % (codebase)) # Codebase to use
        # Source file, all codes need it!
        optfile.write('y\n') # Provide custom source file
        optfile.write('2\n') # Enter path to source file
        optfile.write('%s\n' % (srcfile)) # Source file
        if codebase != "exsim" and codebase != "csm":
            # Run rupture generator?
            if skip_rupgen:
                optfile.write('n\n') # Skip rupture generator
            else:
                optfile.write('y\n') # Run rupture generator
        optfile.write('1\n') # All validation stations
        if skip_rupgen:
            # We should provide SRF file
            optfile.write('2\n') # Enter path to src/srf file
            optfile.write('%s\n' % (srffile)) # SRF file
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
        if gmpe_group_name is None:
            optfile.write('n\n') # Do not generate GMPE comparison plot
        else:
            optfile.write('y\n') # Generate GMPE comparison plot
            optfile.write('%s\n' % (gmpe_group_name))
        optfile.write('y\n') # Yes for GOF
        optfile.write('1\n') # Run GP_GOF
        if not allmetrics:
            optfile.write('n\n') # No additional metrics
        else:
            optfile.write('y\n') # We want more metrics
            optfile.write('y\n') # RZZ2015
            optfile.write('y\n') # FAS
            optfile.write('y\n') # AS2016
            optfile.write('y\n') # RotD100
            optfile.write('y\n') # AndersonGOF
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
            elif codebase == "song":
                end_module = " -e RMG "
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

def write_pbs(install, numsim, simdir, xmldir, email,
              prefix, newnodes, walltime, savetemp):
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
        pbsfile.write("#PBS -q scec\n")
        pbsfile.write("#PBS -l arch=x86_64,pmem=1000mb,pvmem=2000mb,walltime=%d:00:00,nodes=%d:ppn=%d\n" %
                      (walltime, nodes, CORES_PER_NODE))
    else:
        pbsfile.write("#PBS -l arch=x86_64,pmem=1000mb,pvmem=2000mb,walltime=%d:00:00,nodes=%d:ppn=%d:gpus=2\n" %
                      (walltime, nodes, CORES_PER_NODE_NEW))
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
    if savetemp:
        for dir_to_copy in ['outdata', 'indata', 'logs', 'tmpdata']:
            pbsfile.write('python $BBP_DIR/utils/batch/command_parallel.py $BBP_DIR/utils/batch/setup_bbp_env.sh "cp -frp $BBP_DATA_DIR/%s/* $HOME/Sims/%s/." $PBS_NODEFILE\n' %
                          (dir_to_copy, dir_to_copy))
    else:
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

def main():
    """
    Parse command line options and create the needed files/directories
    """
    # Detect BBP installation
    bbp_install = InstallCfg.getInstance()

    # Get GMPE group names
    gmpe_groups_available = gmpe_config.GMPES.keys()
    gmpe_groups_available_lc = [gmpe.lower() for gmpe in gmpe_groups_available]

    prog_base = os.path.basename(sys.argv[0])
    usage = "usage: %s [options]" % (prog_base)
    parser = optparse.OptionParser(usage)
    parser.add_option("-c", "--codebase", type="string", action="store",
                      dest="codebase",
                      help="Codebase for the simulation: %s" %
                      (CODEBASES))
    parser.add_option("-e", "--event", type="string", action="store",
                      dest="event",
                      help="Validation event (should be configured in BBP)")
    parser.add_option("-d", "--dir", type="string", action="store",
                      dest="simdir",
                      help="Simulation directory")
    parser.add_option("--skip-rupgen", action="store_true", dest="skiprupgen",
                      help="Skip the rupture generator, run only 1 simulation,"
                      " unless the --srf option below is provided.")
    parser.add_option("--srf", "--srf-prefix", type="string", action="store",
                      dest="srf_prefix",
                      help="Prefix of SRF files to use, "
                      "only for GP, SDSU and UCSB methods. "
                      "Simulations begin after the rupture generator.")
    parser.add_option("--hypo-rand", action="store_true", dest="hyporand",
                      help="Enables hypocenter randomization")
    parser.add_option("--no-hypo-rand", action="store_false", dest="hyporand",
                      help="Disables hypocenter randomization")
    parser.add_option("-n", "--num-simulations", type="int", action="store",
                      dest="numsim", help="Number of simulations to run")
    parser.add_option("-w", "--walltime", type="int", action="store",
                      dest="walltime", help="Number of hours for walltime")
    parser.add_option("--email", type="string", action="store",
                      dest="email", help="Email for job notifications")
    parser.add_option("--new-nodes", action="store_true", dest="newnodes",
                      help="Schedule the job in the new HPCC nodes")
    parser.add_option("--save-tmpdata", action="store_true", dest="savetemp",
                      help="Save the contents of the tmpdata directory")
    parser.add_option("--only-rup", action="store_true", dest="only_rup",
                      help="Only runs the rupture generator")
    parser.add_option("-g", "--gmpe-group", type="string", action="store",
                      dest="gmpe_group_name",
                      help="GMPE group: %s" % (gmpe_groups_available_lc))
    parser.add_option("-a", "--all-metrics", action="store_true",
                      dest="allmetrics", help="Calculate all metrics")
    parser.add_option("--var", "--variation", type="int", action="store",
                      dest="variation", help="seed variation (default 1)")
    parser.add_option("--seg", "--segment", type="int",
                      action="store", dest="segment",
                      help="Indicates simulation part of multiseg run")
    parser.add_option("--first-seg-dir", type="string", action="store",
                      dest="first_seg_dir",
                      help="required for multi-segment segments 2..n")
    parser.add_option("-s", "--site", action="store_true",
                      dest="site_response", help="Use site response module")

    (options, args) = parser.parse_args()

    # Check if using new HPCC nodes
    if options.newnodes:
        newnodes = True
    else:
        newnodes = False

    # Check if multi-segment simulation
    if options.segment:
        multiseg = True
        segment = options.segment
    else:
        multiseg = False
        segment = None

    # Check for first segment directory
    if options.first_seg_dir is not None:
        first_seg_dir = os.path.abspath(options.first_seg_dir)
        if not os.path.exists(first_seg_dir):
            print("First segment directory does not exist: %s" %
                  (first_seg_dir))
            sys.exit(1)
    else:
        first_seg_dir = None
        if multiseg and segment > 1:
            print("Must specify first segment directory!")
            sys.exit(1)

    # Check for variation sequence
    if options.variation:
        variation = options.variation
    else:
        if multiseg:
            # If a multisegment run, variation defaults to the segment number
            variation = segment
        else:
            # Otherwise, we use 1
            variation = 1

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

    # Check if we need to calculate extra metrics
    if options.allmetrics:
        allmetrics = True
    else:
        allmetrics = False

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
        print("Please specify a codebase!")
        sys.exit(1)
    codebase = codebase.lower()
    if codebase not in CODEBASES:
        print("Codebase needs to be one of: %s" % (CODEBASES))
        sys.exit(1)

    # Check if users wants to run site response module
    if options.site_response:
        site_response = True
        if codebase not in CODEBASES_SITE:
            print("Cannot use site response with method: %s" % (codebase))
            sys.exit(1)
    else:
        site_response = False

    # Check for event
    event = options.event
    if event is None:
        print("Please provide a validation event!")
        sys.exit(1)
    event_names = validation_cfg.VE_EVENTS.get_all_names()
    events = [v_event.lower() for v_event in event_names]
    if event.lower() not in events:
        print("Event %s does not appear to be properly configured on BBP" %
              (event))
        print("Available options are: %s" % (event_names))
        print("Please provide another event or check your BBP installation.")
        sys.exit(1)
    val_obj = validation_cfg.VE_EVENTS.get_event_by_print_name(event)

    # Check if we want to run the rupture generator
    skip_rupgen = options.skiprupgen
    srf_prefix = options.srf_prefix
    if skip_rupgen is not None:
        if codebase not in CODEBASES_SRF:
            print("Cannot use SRF files with method: %s" % (codebase))
            sys.exit(1)

    # Check for hypocenter randomization
    if options.hyporand is None:
        print("Please specify --hypo-rand or --no-hypo-rand!")
        sys.exit(1)

    if options.hyporand:
        hypo_rand = True
    else:
        hypo_rand = False

    try:
        source_file = val_obj.get_input(codebase, "source")
    except KeyError:
        print("Unable to get source file for event %s, codebase %s!" %
              (event, codebase))
        sys.exit(1)
    if not source_file:
        print("Source file for event %s, codebase %s not specified!" %
              (event, codebase))
        sys.exit(1)
    # Check for multisegment events
    if isinstance(source_file, str):
        source_file = source_file.strip()
        if multiseg:
            print("This event doesn't look like a multisegment event!")
            sys.exit(1)
    else:
        # Multisegment event
        if not multiseg:
            print("This is a multisegment event! Please specify segment!")
            sys.exit(1)
        source_file = source_file[segment - 1]

    if skip_rupgen:
        if not srf_prefix:
            try:
                srf_file = val_obj.get_input(codebase, "srf").strip()
            except KeyError:
                print ("Event %s does not have a srf file for codebase %s!" %
                       (event, codebase))
                sys.exit(1)
            except AttributeError:
                print ("Event %s does not have a srf file for codebase %s!" %
                       (event, codebase))
                sys.exit(1)
            if not srf_file:
                print ("Event %s does not have a srf file for codebase %s!" %
                       (event, codebase))
                sys.exit(1)
            # Force number of simulations to 1
            options.numsim = 1
        else:
            # Use the SRF prefix to find the SRF files we need
            # Make it a full path
            srf_prefix = os.path.realpath(srf_prefix)
            # Make sure source file is in the rcf-104 or scec-00 filesystems
            if not "rcf-104" in srf_prefix and not "scec-00" in srf_prefix:
                print "SRF files should be in the rcf-104 / scec-00 filesystems!"
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

    # Figure out which gmpe group to use
    gmpe_group_name = options.gmpe_group_name
    if gmpe_group_name is not None:
        if not gmpe_group_name.lower() in gmpe_groups_available_lc:
            print "Invalid gmpe group name!"
            print "Options are: %s" % (gmpe_groups_available_lc)
            sys.exit(1)
        gmpe_group_index = gmpe_groups_available_lc.index(gmpe_group_name.lower())
        gmpe_group_name = gmpe_groups_available[gmpe_group_index]

    # Make sure user has configured the setup_bbp_env.sh script
    setup_bbp_env = os.path.join(bbp_install.A_INSTALL_ROOT,
                                 "utils/batch/setup_bbp_env.sh")
    if not os.path.exists(setup_bbp_env):
        print ("Cannot find setup_bbp_env.sh script!")
        print ("Expected at: %s" % (setup_bbp_env))
        sys.exit(1)
    # Create simulation directories
    prefix = "%s-%s" % (event.lower(), codebase.lower())
    # Make sure we remove spaces from prefix (e.g. for the "Loma Prieta" event)
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
    # Generate source files if needed
    if source_file is not None:
        generate_src_files(numsim, source_file, srcdir,
                           prefix, hypo_rand,
                           variation, multiseg,
                           segment, first_seg_dir)
    # Generate xml files
    generate_xml(bbp_install, numsim, srcdir, xmldir,
                 logsdir, event, codebase, prefix,
                 skip_rupgen, srf_prefix, only_rup, gmpe_group_name,
                 allmetrics, site_response, multiseg, segment)
    # Write pbs file
    write_pbs(bbp_install, numsim, simdir, xmldir,
              email, prefix, newnodes, walltime, savetemp)

    # Write .info file
    info_file = open(os.path.join(simdir, "%s.info" % (prefix)), 'w')
    info_file.write("# %s\n" % (" ".join(sys.argv)))
    info_file.close()

if __name__ == "__main__":
    main()
