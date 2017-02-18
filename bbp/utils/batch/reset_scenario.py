#!/usr/bin/env python

import sys
import os


ScenarioID = "Scenario106"
print "WARNING! This program will delete sim data from previous runs!"
print "Data in "../%s/Logs1,../%s/RunDesc1,../%s/Sims,../%s/Xml1"%(ScenarioID,ScenarioID,ScenarioID,ScenarioID)

proceed = raw_input("Do you want to proceed? (yes/no)")

if str.lower(proceed) != "yes":
  sys.exit(0)
else:
  print "Deleting simulation data directories..."
  progstring = "rm -rf ../%s/Logs1"%(ScenarioID)
  os.system(progstring)
  progstring = "rm -rf ../%s/RunDesc1"%(ScenarioID)
  os.system(progstring)
  progstring = "rm -rf ../%s/Sims"%(ScenarioID)
  os.system(progstring)
  progstring = "rm -rf ../%s/Xml1"%(ScenarioID)
  os.system(progstring)


print "For complete reset, please reset_acct.py including xml directory"
sys.exit(0)
