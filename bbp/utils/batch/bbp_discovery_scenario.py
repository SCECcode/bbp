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

Program to set up a full BBP scenario simulation run on HPCC
"""
from __future__ import division, print_function

# Works for both Python 2 and 3
try: input = raw_input
except NameError: pass

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
CORES_PER_NODE = 16
CORES_PER_NODE_SCEC = 16
MAX_SIMULATIONS = 500
CODEBASES = ["gp", "ucsb", "sdsu", "exsim", "csm", "song", "irikura1", "irikura2"]
CODEBASES_SITE = ["gp", "sdsu", "song", "irikura1", "irikura2", "exsim", "ucsb"]
CODEBASES_SRF = ["gp", "sdsu", "song", "ucsb"]
SITE_MODULES = ["gp", "seismosoil"]

def generate_src_files(numsim, source_file, srcdir,
                       prefix, hypo_rand, hypo_area,
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

def generate_multi_segment_src_files(numsim, source_files,
                                     srcdir, prefix):
    """                                                                                                               
    Generates num_sim source files in the srcdir using different                                                      
    random seeds, make sure we keep the common_seed constant for                                                      
    all indivual segments in the multi-segment run                                                                    
    """
    src_props = []
    all_seeds = []
    output = []

    for source_file in source_files:
        src_props.append(bband_utils.parse_properties(source_file))

    # Delete "seed" and "common_seed" from the property set                                                           
    for props in src_props:
        if "seed" in props:
            props.pop("seed")
        if "common_seed" in props:
            props.pop("common_seed")

    # Start preparing output                                                                                          
    for idx, props in enumerate(src_props):
        # Create common list of keys for all realizations                                                             
        output.append("")
        for key in props:
            output[idx] = output[idx] + "%s = %s\n" % (key.upper(), props[key])

    # Generate seeds for all source files and realizations                                                            
    for variation in range(1, (len(source_files) + 1)):
        seeds = []
        for sim in range(0, numsim):
            random.seed((sim + 1) + (variation - 1) * 500)
            seeds.append(int(math.exp(7 * math.log(10.0)) * random.random()))
        all_seeds.append(seeds)

    # Write seeds file                                                                                                
    seed_file = open(os.path.join(srcdir, "seeds.txt"), 'w')
    seed_file.write("%d\n" % (numsim))
    for seed in all_seeds[0]:
        seed_file.write("%d\n" % (seed))
    seed_file.close()

    # Generate the numsim SRC files                                                                                   
    for sim in range(0, numsim):
        for segment in range(0, len(source_files)):
            srcfile = os.path.join(srcdir,
                                   "%s-%04d_seg%02d.src" %
                                    (prefix, sim, (segment + 1)))
            outfile = open(srcfile, 'w')
            # Write common part to all realizations                                                                   
            outfile.write(output[segment])
            # Write seed                                                                                              
            outfile.write("SEED = %d\n" % (all_seeds[segment][sim]))
            # Write common seed except for segment 1                                                                  
            if segment > 0:
                outfile.write("COMMON_SEED = %d\n" % (all_seeds[0][sim]))
            outfile.close()

def generate_xml(install, numsim, srcdir, xmldir,
                 logdir, vmodel, codebase, prefix,
                 station_list, only_rup, srf_prefix,
                 site_response, multiseg, segment,
                 source_file):
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
        # Figure out if this is a multi-segment run with                                                              
        # multiple SRC files                                                                                          
        if source_file is not None:
            if not isinstance(source_file, str):
                # Yes!                                                                                                
                srcfile = []
                for seg in range(0, len(source_file)):
                    srcfile.append(os.path.join(srcdir,
                                                "%s-%04d_seg%02d.src" %
                                                (prefix, sim, (seg + 1))))
                # Now create a string representation of the list                                                      
                srcfile = str(srcfile)
            else:
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
                optfile.write('%s\n' % (site_response))
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
            elif codebase == "song":
                if source_file is not None:
                    if not isinstance(source_file, str):
                        end_module = " -e SongRMGMS "
                    else:
                        end_module = " -e SongRMGSS "
                else:
                    end_module = " -e SongRMGSS "
            elif codebase == "irikura1" or codebase == "irikura2":
                end_module = " -e IrikuraGenSrf "
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
                prefix, scecnodes, walltime, savetemp,
                codebase):
    """
    Write the slurm script
    """
    # Calculate how many nodes we need
    if scecnodes:
        nodes = int(math.ceil(1.0 * numsim / CORES_PER_NODE_SCEC))
        cores_per_node = CORES_PER_NODE_SCEC
    else:
        nodes = int(math.ceil(1.0 * numsim / CORES_PER_NODE))
        cores_per_node = CORES_PER_NODE
    # Some path names
    outfile = os.path.join(simdir, "%s.out" % (prefix))
    errfile = os.path.join(simdir, "%s.err" % (prefix))
    bfn = os.path.join(xmldir, BATCH_SIM_FILE)
    # Let's create the slurm file
    slurmfn = os.path.join(simdir, "%s.slurm" % (prefix))
    slurmfile = open(slurmfn, 'w')
    slurmfile.write("#!/bin/bash\n")
    slurmfile.write("\n")
    if scecnodes:
        slurmfile.write("#SBATCH --partition=scec\n")
    else:
        #slurmfile.write("#SBATCH --partition=main\n")
        #slurmfile.write("#SBATCH --gres=gpu:2\n")
        pass
    slurmfile.write("#SBATCH --mem=0\n")
    slurmfile.write("#SBATCH --nodes=%d\n" % (nodes))
    slurmfile.write("#SBATCH --ntasks=%d\n" % (nodes))
    slurmfile.write("#SBATCH --cpus-per-task=%d\n" % (cores_per_node))
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
    slurmfile.write("BBP_DATA_DIR=%s/$SLURM_JOB_ID\n" %
                    (os.environ["BBP_SIMS_DIR"]))
    slurmfile.write("BBP_BASE_DIR=%s/$SLURM_JOB_ID\n" %
                    (os.environ["BBP_SIMS_DIR"]))
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
#            slurmfile.write('python $BBP_DIR/utils/batch/command_parallel.py $BBP_DIR/utils/batch/setup_bbp_env.sh "cp -frp $BBP_DATA_DIR/%s/* $HOME/Sims/%s/." $SLURM_NODES\n' %
#                            (dir_to_copy, dir_to_copy))
            slurmfile.write('cp -frp $BBP_DATA_DIR/%s/* $HOME/Sims/%s/.\n' %
                            (dir_to_copy, dir_to_copy))
    else:
        for dir_to_copy in ['outdata', 'indata', 'logs']:
#            slurmfile.write('python $BBP_DIR/utils/batch/command_parallel.py $BBP_DIR/utils/batch/setup_bbp_env.sh "cp -frp $BBP_DATA_DIR/%s/* $HOME/Sims/%s/." $SLURM_NODES\n' %
#                            (dir_to_copy, dir_to_copy))
            slurmfile.write('cp -frp $BBP_DATA_DIR/%s/* $HOME/Sims/%s/.\n' %
                            (dir_to_copy, dir_to_copy))
    slurmfile.write("\n")
    slurmfile.write('echo "Jobs end"\n')
    slurmfile.write("date\n")
    slurmfile.flush()
    slurmfile.close()

    # All done!
    print()
    print("Validation run is set up on: %s" % (simdir))
    print()
    print("To start the validation run, just type: ")
    print("$ sbatch %s" % (slurmfn))
    print()

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
    parser.add_option("--skip-rupgen", action="store_true", dest="skiprupgen",
                      help="Skip the rupture generator, run only 1 simulation,"
                      " unless the --srf option below is provided.")
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
    parser.add_option("--scec-nodes", action="store_true", dest="scecnodes",
                      help="Schedule the job in the SCEC nodes")
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
    parser.add_option("--seg", "--segment", type="int",
                      action="store", dest="segment",
                      help="Indicates simulation part of multiseg run")
    parser.add_option("--first-seg-dir", type="string", action="store",
                      dest="first_seg_dir",
                      help="required for multi-segment segments 2..n")
    parser.add_option("-s", "--site", type="string", action="store",
                      dest="site_response",
                      help="Use a site response module: %s" % (SITE_MODULES))

    (options, _) = parser.parse_args()

    # Check if using SCEC nodes
    if options.scecnodes:
        scecnodes = True
    else:
        scecnodes = False

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
            print("First segment directory for exists: %s" %
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
        if scecnodes:
            walltime = 300
        else:
            walltime = 24

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

    # Check for velocity model
    vmodel_names = velocity_models.get_all_names()
    vmodel = options.vmodel
    if vmodel is None:
        print("Please provide a velocity model (region) for this simulation!")
        print("Available options are: %s" % (vmodel_names))
        sys.exit(1)
    vmodels = [v_model.lower() for v_model in vmodel_names]
    if vmodel.lower() not in vmodels:
        print("Velocity model %s does not appear to be available on BBP" %
              (vmodel))
        print("Available options are: %s" % (vmodel_names))
        print("Please provide another velocity model or check your BBP installation.")
        sys.exit(1)
    # Now get the name with the correct case
    vmodel = vmodel_names[vmodels.index(vmodel.lower())]

    # Check if users wants to run site response module
    site_response = options.site_response
    if site_response is not None:
        if codebase not in CODEBASES_SITE:
            print("Cannot use site response with method: %s" % (codebase))
            sys.exit(1)
        site_response = site_response.lower()
        if site_response not in SITE_MODULES:
            print("Site response needs to be one of: %s" % (SITE_MODULES))
            sys.exit(1)
    else:
        site_response = False

    # Check for hypocenter randomization
    if options.hyporand is None:
        print("Please specify --hypo-rand or --no-hypo-rand!")
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

    # Get the source file
    source_file = options.source
    if source_file is None:
        print("Please provide a source description (src file)!")
        sys.exit(1)

    # Make list
    pieces = source_file.split(",")
    pieces = [piece.strip() for piece in pieces]

    # Check if single file or list of files
    if len(pieces) == 1:
        source_file = pieces[0]
    else:
        source_file = pieces

    skip_rupgen = options.skiprupgen
    srf_prefix = options.srf_prefix
    if skip_rupgen is not None:
        if codebase not in CODEBASES_SRF:
            print("Cannot use SRF files with method: %s" % (codebase))
            sys.exit(1)

    # Check for multisegment events                                                                                   
    if isinstance(source_file, str):
        source_file = source_file.strip()
        if multiseg:
            print("This event doesn't look like a multisegment event!")
            sys.exit(1)
        else:
            prefix = ("%s-%s" %
                      (os.path.splitext(os.path.basename(source_file))[0],
                       codebase.lower()))
    else:
        # Multisegment event                                                                                          
        if multiseg:
            source_file = source_file[segment - 1]
            # Create a prefix
            prefix = ("%s-%s" %
                      (os.path.splitext(os.path.basename(source_file))[0],
                       codebase.lower()))
        else:
            if codebase == "ucsb" or codebase == "exsim":
                print("This is a multisegment event! Please specify segment!")
                sys.exit(1)
            else:
                # Other methods can work with multisegment files directly
                # Create a prefix
                prefix = ("%s-%s" %
                          (os.path.splitext(os.path.basename(source_file[0]))[0],
                           codebase.lower()))

    # If user specified a SRF prefix
    if srf_prefix is not None:
        # Make it a full path
        srf_prefix = os.path.realpath(srf_prefix)
        # Create a prefix
        prefix = os.path.splitext(os.path.basename(srf_prefix))[0]

    # Make sure we remove spaces from prefix
    prefix = prefix.replace(" ", '')

    # Get the station list
    station_list = options.station_list
    if station_list is None:
        print("Please provide a station list (stl file)!")
        sys.exit(1)
    # Make it a full path
    station_list = os.path.realpath(station_list)
    # Make sure station list exists and is readable
    if (not os.path.isfile(station_list) or
        not os.access(station_list, os.R_OK)):
        print("Station list foes not seem to be accessible!")
        sys.exit(1)

    # Check for the simulation directory
    simdir = options.simdir
    if simdir is None:
        print("Please provide a simulation directory!")
        sys.exit(1)
    simdir = os.path.abspath(simdir)
    if os.path.exists(simdir):
        print("Simulation directory exists: %s" % (simdir))
        opt = input("Do you want to delete its contents (y/n)? ")
        if opt.lower() != "y":
            print("Please provide another simulation directory!")
            sys.exit(1)
        opt = input("ARE YOU SURE (y/n)? ")
        if opt.lower() != "y":
            print("Please provide another simulation directory!")
            sys.exit(1)
        # Delete existing directory (we already asked the user twice!!!)
        shutil.rmtree(simdir)

    # Pick up number of simulations to run
    numsim = options.numsim
    if numsim is None:
        print("Please provide the number of simulations to run!")
        sys.exit(1)
    if numsim < 1 or numsim > MAX_SIMULATIONS:
        print("Number of simulations should be between 1 and %d" %
              (MAX_SIMULATIONS))
        sys.exit(1)

    # Check for e-mail address
    email = options.email
    if email is None:
        print("Please provide an e-mail address for job notifications")
        sys.exit(1)

    # Make sure user has configured the setup_bbp_env.sh script
    setup_bbp_env = os.path.join(bbp_install.A_INSTALL_ROOT,
                                 "utils/batch/setup_bbp_env.sh")
    if not os.path.exists(setup_bbp_env):
        print("Cannot find setup_bbp_env.sh script!")
        print("Expected at: %s" % (setup_bbp_env))
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
    if isinstance(source_file, str):
        # Generate source files
        generate_src_files(numsim, source_file, srcdir,
                           prefix, hypo_rand, hypo_area,
                           variation, multiseg,
                           segment, first_seg_dir)
    else:
        generate_multi_segment_src_files(numsim, source_file,
                                         srcdir, prefix)
    # Generate xml files
    generate_xml(bbp_install, numsim, srcdir, xmldir,
                 logsdir, vmodel, codebase, prefix,
                 station_list, only_rup, srf_prefix,
                 site_response, multiseg, segment,
                 source_file)
    # Write slurm file
    write_slurm(bbp_install, numsim, simdir, xmldir,
                email, prefix, scecnodes, walltime,
                savetemp, codebase)

    # Write .info file
    info_file = open(os.path.join(simdir, "%s.info" % (prefix)), 'w')
    info_file.write("# %s\n" % (" ".join(sys.argv)))
    info_file.close()

if __name__ == "__main__":
    main()
