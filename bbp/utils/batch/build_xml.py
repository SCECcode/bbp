#!/usr/bin/env python

#build_xml.py
#v1.1, 20110817
#Generate XML files formatted for BBP
#Input: Run description files (txt)
#Output: BBP XML file

#Sample input run description file:
#m6.00_d20_r90_z0_0100.txt
#RUN_TAG = 0010100
#VALIDATION_RUN = n
#SOURCE_DESCRIPTION_FILE = /home/scec-00/rgraves/NgaW2/FwHw/FaultInfo/Inputs/m6.00_d20_r90_z0.src
#STATION_LIST_FILE = /home/scec-00/rgraves/NgaW2/FwHw/StatInfo/rv01-m6.00_stats.stl
#RUPTURE_GENERATOR = URS
#LOW_FREQUENCY_MODULE = URS
#HIGH_FREQUENCY_MODULE = URS
#SITE_RESPONSE_MODULE = URS
#PLOT_VEL = y
#PLOT_ACC = y
#RUN_GOF = n

import os
import shutil
import sys
import time
import optparse

from install_cfg import InstallCfg
import bband_utils
import velocity_models

def main():
    parser = optparse.OptionParser()
    parser.add_option("-x", "--xml-dir", dest="xmlDir", help="Output folder for XML files", metavar="XML_DIR")
    parser.add_option("-i", "--in-dir", dest="inDir", help="Input folder with run description files", metavar="IN_DIR")

    (options, args) = parser.parse_args()

    if not options.inDir:
        parser.error("Folder with run description files is required!")

    if os.path.exists(options.inDir):
        files = os.listdir(options.inDir)
    else:
        print "Invalid input dir: %s" % options.inDir
        sys.exit(1)

    if not files:
        print "No run description files were found in input dir: %s" % inDir
        sys.exit(1)

    if not options.xmlDir:
        cur_dir= os.getcwd()
        options.xmlDir = os.path.join(cur_dir, "run_xml")
        print ("Note: Output XML directory was not specified!\n"
               "XMLs will be written to %s" % options.xmlDir)

    if not os.path.exists(options.xmlDir):
        os.mkdir(options.xmlDir)

    required_keys = ["RUN_TAG", "VALIDATION_RUN", "VMODEL_NAME",
                     "SOURCE_DESCRIPTION_FILE", "STATION_LIST_FILE",
                     "RUPTURE_GENERATOR", "LOW_FREQUENCY_MODULE",
                     "HIGH_FREQUENCY_MODULE", "SITE_RESPONSE_MODULE",
                     "PLOT_VEL", "PLOT_ACC", "RUN_GOF", "GEN_HTML"]

    cfg = InstallCfg()

    valid_event = ['y', 'n']
    available_vmodels = velocity_models.get_all_names()
    rupture_gen = ['y', 'n']
    rup = ['GP', 'UCSB']
    lf = ['GP', 'UCSB']
    hf = ['GP', 'UCSB', 'SDSU']
    site = ['GP', 'UCSB', 'SDSU', "None"]
    plotVel = ['y', 'n']
    plotAcc = ['y', 'n']
    doGof = ['y', 'n']
    gof = ['GP', 'SDSU']
    genHtml = ['y', 'n']

    # Move out all the files we don't need from start dir
    for file in os.listdir(cfg.A_USER_DATA_DIR):
        if file.endswith(".srf") or file.endswith(".stl") or file.endswith(".src"):
            if not os.path.exists("%s/tmp" % cfg.A_USER_DATA_DIR):
                os.mkdir("%s/tmp" % cfg.A_USER_DATA_DIR)
            shutil.move("%s/%s" % (cfg.A_USER_DATA_DIR, file), "%s/tmp/%s" % (cfg.A_USER_DATA_DIR, file))

    filecount = 0

    for file in files:
        if file.endswith(".txt"):
            file_base = file[0:file.find(".txt")]
            #pieces = file_base.split("_")
            input_config = {}
            output_files = {}
            in_file = "%s/%s" % (options.inDir, file)
            fn = open(in_file, 'r')
            lines = fn.readlines()

            if not lines:
                print "Warning: Skipping empty input file: %s"% file
                continue

            for line in lines:
                key,value = line.strip().split("=")
                input_config[key.strip()] = value.strip()

            if not input_config:
                print "Warning: Skipping malformed input file: %s" % file
                continue

            input_keys = input_config.keys()

            if not set(input_keys) == set(required_keys):
                print "Warning: Skipping malformed input file: %s" % file
                print " Missing keys:", list(set(input_keys).symmetric_difference(set(required_keys)))
                continue

            valid = True
            try:
                simID = input_config["RUN_TAG"]
                if simID[0] == "0":
                    print "Warning: RUN_TAG cannot start with '0'! RUN_TAG will be left-padded with '1'!"
                    simID ="1%s" % simID
                simID = int(simID)
            except ValueError:
                print "Invalid value for RUN_TAG: %s, Integer value expected!" % input_config["RUN_TAG"]
                valid = False

            validation = input_config["VALIDATION_RUN"]
            if not validation in valid_event:
                print ("Invalid option for VALIDATION_RUN: %s, Expected: %s" %
                       (validation, str(valid_event)))
                valid = False

            vel_model = input_config["VMODEL_NAME"]
            if not vel_model in available_vmodels:
                print ("Velocity model %s not available in the platform!" %
                       (vel_model))
                print available_vmodels
                valid = False

            if validation =='n':
                r_gen = input_config["RUPTURE_GENERATOR"]
                if not r_gen in rup:
                    print "Invalid option for RUPTURE_GENERATOR: %s, Expected: %s" % (r_gen, str(rup))
                    valid = False
                    rupture_generator = 'n'
                else:
                    r_gen = rup.index(r_gen) +1
                    rupture_generator = 'y'

            lfm = input_config["LOW_FREQUENCY_MODULE"]
            if not lfm in lf:
                print "Invalid option for LOW_FREQUENCY_MODULE: %s, Expected: %s" % (lfm, str(lf))
                valid = False
            else:
                lfm = lf.index(lfm) +1

            hfm = input_config["HIGH_FREQUENCY_MODULE"]
            if not hfm in hf:
                print "Invalid option for HIGH_FREQUENCY_MODULE: %s, Expected: %s" % (hfm, str(hf))
                valid = False
            else:
                hfm = hf.index(hfm) +1

            srm = input_config["SITE_RESPONSE_MODULE"]
            if not srm in site:
                print "Invalid option for SITE_RESPONSE_MODULE: %s, Expected: %s"%(srm, str(site))
                valid = False
            else:
                srm = site.index(srm) +1

            plt_vel = input_config["PLOT_VEL"]
            if not validation in plotVel:
                print "Invalid option for PLOT_VEL: %s, Expected: %s"% (plt_vel, str(plotVel))
                valid = False

            plt_acc= input_config["PLOT_ACC"]
            if not validation in plotAcc:
                print ("Invalid option for PLOT_ACC: %s, Expected: %s" %
                       (plt_acc, str(plotAcc)))
                valid = False

            gof_opt = input_config["RUN_GOF"]
            if not validation in doGof:
                print ("Invalid option for RUN_GOF: %s, Expected: %s" %
                       (gof_opt, str(doGof)))
                valid = False

            gen_html = input_config["GEN_HTML"]
            if not gen_html in genHtml:
                print ("Invalid option for GEN_HTML: %s, Expected: %s" %
                       (gen_html, str(genHtml)))
                valid = False

            src_file = os.path.abspath(input_config["SOURCE_DESCRIPTION_FILE"])
            if not os.path.exists(src_file):
                print "Unable to locate specified source file: %s" % src_file
                valid = False

            stat_file = os.path.abspath(input_config["STATION_LIST_FILE"])
            if not os.path.exists(stat_file):
                print "Unable to locate stations file: %s" % stat_file
                valid = False

            if not valid:
                print "Skipping input file %s due to previous errors!" % file
                continue

            # Generate an options file
            optfilename = "%s/%s_%s.txt" % (options.xmlDir, str(simID), file_base)

            print "Generating options file %s" % optfilename
            fp = open(optfilename, 'w')
            fp.write("%s\n" % validation)
            fp.write("%s\n" % vel_model)
            fp.write("%s\n" % rupture_generator)
            if rupture_generator == 'y':
                fp.write("%d\n" % r_gen)
                fp.write("2\n") # Enter path to source file
                fp.write("%s\n" % (src_file))
            else:
                # Need to write a path to a srf file
                pass
            fp.write("%d\n" % lfm)
            fp.write("2\n") # Enter path to station list
            fp.write("%s\n" % (stat_file))
            fp.write("%d\n" % hfm)
            fp.write("%d\n" % srm)
            fp.write("%s\n" % plt_vel)
            fp.write("%s\n" % plt_acc)
            fp.write("%s\n" % gof_opt)
            if gof_opt == 'y':
                fp.write("1\n") #defualt to GP GOF module for now!
            fp.write("%s\n" % gen_html)
            fp.flush()
            fp.close()

            #copy src and stl files to start dir
            #shutil.copy2(src_file, "%s" % cfg.A_USER_DATA_DIR)
            #shutil.copy2(stat_file, "%s" % cfg.A_USER_DATA_DIR)

            #move bpp dirs with simID
            indatadir = "%s/%d" % (cfg.A_IN_DATA_DIR, simID)
            outdatadir = "%s/%d" % (cfg.A_OUT_DATA_DIR, simID)
            tmpdatadir = "%s/%d" % (cfg.A_TMP_DATA_DIR, simID)
            logdir = "%s/%d" % (cfg.A_OUT_LOG_DIR, simID)

            if os.path.exists(indatadir):
                shutil.move(indatadir, "%s_%s" % (indatadir, "bkp"))
            if os.path.exists(tmpdatadir):
                shutil.move(tmpdatadir, "%s_%s" % (tmpdatadir, "bkp"))
            if os.path.exists(outdatadir):
                shutil.move(outdatadir, "%s_%s" % (outdatadir, "bkp"))
            if os.path.exists(logdir):
                shutil.move(logdir, "%s_%s" % (logdir, "bkp"))

            #Generate XML
            #bband_utils.runprog("%s/run_bbp.py -o %s -s %d" % (cfg.A_COMP_DIR, filename, simID))
            bband_utils.runprog("%s/run_bbp.py -o %s -s %d -g" % (cfg.A_COMP_DIR, optfilename, simID))
            outdir = "%s/%d" % (cfg.A_OUT_DATA_DIR, simID)

            #Remove option file
            os.remove(optfilename)

            #Remove src and stl files from start dir
            #os.remove("%s/%s" % (cfg.A_USER_DATA_DIR, os.path.basename(src_file)))
            #os.remove("%s/%s" % (cfg.A_USER_DATA_DIR, os.path.basename(stat_file)))

            #move back bpp dirs with simID
            if os.path.exists("%s_%s" % (indatadir, "bkp")):
                shutil.move("%s_%s" % (indatadir, "bkp"), indatadir)
            else:
                shutil.rmtree(indatadir)
            if os.path.exists("%s_%s" % (tmpdatadir, "bkp")):
                shutil.move("%s_%s" % (tmpdatadir, "bkp"), tmpdatadir)
            else:
                shutil.rmtree(tmpdatadir)
            if os.path.exists("%s_%s" % (outdatadir, "bkp")):
                shutil.move("%s_%s" % (outdatadir, "bkp"), outdatadir)
            else:
                shutil.rmtree(outdatadir)
            if os.path.exists("%s_%s" % (logdir, "bkp")):
                shutil.move("%s_%s" % (logdir, "bkp"), logdir)
            else:
                shutil.rmtree(logdir)

            #Copy xml to xmlDir
            xmlfilename = "%s/%d_%s.xml" % (options.xmlDir, simID, file_base)
            if not os.path.exists("%s/%d.xml" % (cfg.A_XML_DIR, simID)):
                print "Failed to generate XML file for %s" % file
                continue
            shutil.copy2("%s/%d.xml" % (cfg.A_XML_DIR, simID), xmlfilename)

            #Fix path in XML file
#                       xfn = open(xmlfilename, 'r')
#                       xlines = xfn.readlines()
#                       xfn.close()
#
#                       xfn=open(xmlfilename, 'w')
#                       src_file_base = os.path.basename(src_file)
#                       stat_file_base = os.path.basename(stat_file)
#                       old_src_line = "$BBP_INSTALL/start/%s" % src_file_base
#                       old_stat_line = "$BBP_INSTALL/start/%s" % stat_file_base
#                       for xline in xlines:
#                               if xline.strip().startswith(old_src_line):
#                                       xfn.write(xline.replace(old_src_line,src_file))
#                               elif xline.strip().startswith(old_stat_line):
#                                       xfn.write(xline.replace(old_stat_line,stat_file))
#                               else:
#                                       xfn.write(xline)
#
#                       xfn.flush()
#                       xfn.close()

            filecount +=1

    print "Processed %d files, skipped %d files!" % (filecount, len(files) - filecount)

if __name__ == "__main__":
    main()
