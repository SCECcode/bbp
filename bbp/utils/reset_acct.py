#!/usr/bin/env python
"""
This program is used to reset an user broadband installation by
removing data from all previous simulations.
$Id: reset_acct.py 929 2012-08-09 23:37:42Z maechlin $
"""
# Import Python modules
import sys

# Import Broadband modules
import bband_utils
from install_cfg import InstallCfg

def do_reset():
    """
    This function resets an user account by deleting all simulation
    data from previous runs. For each of outdata, tmpdata, indata, and
    logs, it confirms if the user wants to delete the data.
    """

    install = InstallCfg()

    # Print warning message
    print "WARNING".center(80, '=')
    print ("  You are about to permanently delete Broadband simulation data")
    print ("  from all previous runs!")
    print ("  Data in indata, outdata, tmpdata, xml, and logs directories ")
    print ("  will be permanently deleted!")
    print '='*80
    print
    print (" The BBP data directory is set to: %s" % (install.A_DATA_ROOT))
    print

    while True:
        proceed = raw_input("Do you want to proceed? (yes/no) ")
        if str.lower(proceed) == "yes":
            print "Deleting simulation data directories..."
            break
        elif str.lower(proceed) == "no":
            sys.exit(0)

    #proceed = raw_input("Delete 'outdata'? (yes/no) ")
    if str.lower(proceed) == "yes":
        prog_string = "rm -r %s" % (install.A_OUT_DATA_DIR)
        bband_utils.runprog(prog_string)
        print "Deleted 'outdata'!"
    else:
        print "Skipped 'outdata'"

    #proceed = raw_input("Delete 'tmpdata'? (yes/no) ")
    if str.lower(proceed) == "yes":
        prog_string = "rm -r %s" % (install.A_TMP_DATA_DIR)
        bband_utils.runprog(prog_string)
        print "Deleted 'tmpdata'!"
    else:
        print "Skipped 'tmpdata'"

    #proceed = raw_input("Delete 'indata'? (yes/no) ")
    if str.lower(proceed) == "yes":
        prog_string = "rm -r %s" % (install.A_IN_DATA_DIR)
        bband_utils.runprog(prog_string)
        print "Deleted 'indata'!"
    else:
        print "Skipped 'indata'"

    #proceed = raw_input("Delete 'logs'? (yes/no) ")
    if str.lower(proceed) == "yes":
        prog_string = "rm -r %s" % (install.A_OUT_LOG_DIR)
        bband_utils.runprog(prog_string)
        print "Deleted 'logs'!"
    else:
        print "Skipped 'logs'"

    #proceed = raw_input("Delete 'xml'? (yes/no) ")
    if str.lower(proceed) == "yes":
        prog_string = "rm -r %s" % (install.A_XML_DIR)
        bband_utils.runprog(prog_string)
        print "Deleted 'xml'!"
    else:
        print "Skipped 'xml'"

    #
    # delete *.pyc files in comps and in test
    if str.lower(proceed) == "yes":
        prog_string = "rm %s/comps/*.pyc" % (install.A_INSTALL_ROOT)
        bband_utils.runprog(prog_string)
        print "Deleted *.pyc files in comps"
        prog_string = "rm %s/tests/*.pyc" % (install.A_INSTALL_ROOT)
        bband_utils.runprog(prog_string)
        print "Deleted *.pyc files in test"

    #
    # Make clean to remove executables
    #
    if str.lower(proceed) == "yes":
        prog_string = "cd %s/src;make clean" % (install.A_INSTALL_ROOT)
        bband_utils.runprog(prog_string)
        print "src/make clean removed executables in src directory"


if __name__ == "__main__":
    do_reset()
    # All done!
    sys.exit(0)
