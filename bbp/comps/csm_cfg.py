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

This module defines the configuration parameters for the CSM module
"""
from __future__ import division, print_function

# Import Python modules
import os
import sys
import math
import numpy
import random

# Import Broadband modules
import bband_utils

class CSMCfg(object):
    """
    Define the configuration parameters for the CSM program
    """
    SIMULA_MAX_STATIONS = 300
    simula_stations = "c2_input_stn"
    simula_random = "c1_rn.dat"
    csm_stations = "stat.dat"
    nuclear_in = "nuclear.in"
    simula_in = "simula.in"
    scat1d_in = "scat1d.in"
    compom_in = "compom.in"
    csevents_dat = "csevents01.dat"
    cfgdict = {}
    cfgparams = {}
    nlay = 0
    vmodel = {}

    def getval(self, attr):
        try:
            val = self.cfgdict[attr]
        except KeyError:
            print("Invalid Source File - Missing attribute: %s" % (attr))
            print("Exiting")
            sys.exit(1)
        return val

    def parse_src(self, a_srcfile):
        """
        This function calls bband_utils' parse property file function
        to get a dictionary of key, value pairs and then looks for the
        parameters needed by CSM
        """
        self.cfgdict = bband_utils.parse_properties(a_srcfile)

        val = self.getval("magnitude")
        self.MAGNITUDE = float(val)

        val = self.getval("fault_length")
        self.LENGTH = float(val)

        val = self.getval("fault_width")
        self.WIDTH = float(val)

        val = self.getval("depth_to_top")
        self.DEPTH_TO_TOP = float(val)

        val = self.getval("strike")
        self.STRIKE = float(val)

        val = self.getval("rake")
        self.RAKE = float(val)

        val = self.getval("dip")
        self.DIP = float(val)

        val = self.getval("lat_top_center")
        self.LAT_TOP_CENTER = float(val)

        val = self.getval("lon_top_center")
        self.LON_TOP_CENTER = float(val)

        val = self.getval("hypo_along_stk")
        self.HYPO_ALONG_STK = float(val)

        val = self.getval("hypo_down_dip")
        self.HYPO_DOWN_DIP = float(val)

        val = self.getval("seed")
        self.SEED = int(val)

    def parse_velmodel(self, a_velmodel):
        """
        This function parses the velocity model file and stores its
        content in the self.vmodel variable
        """
        # Initialize velocity model structure
        self.vmodel = {'h': [],
                       'vp': [],
                       'qp': [],
                       'vs': [],
                       'qs': [],
                       'rho': []}

        vel_file = open(a_velmodel, 'r')
        for line in vel_file:
            line = line.strip()
            pieces = line.split()
            if len(pieces) == 3:
                self.nlay = float(pieces[0])
                continue
            # Skip lines without the 6 values
            if len(pieces) != 6:
                continue
            pieces = [float(piece) for piece in pieces]
            self.vmodel['h'].append(pieces[0])
            self.vmodel['vp'].append(pieces[1])
            self.vmodel['qp'].append(pieces[2])
            self.vmodel['vs'].append(pieces[3])
            self.vmodel['qs'].append(pieces[4])
            self.vmodel['rho'].append(pieces[5])
        vel_file.close()

    def zetaetaxi2xyz(self, zeta, eta, xi, strike, dip, zref):
        """
        A low level function used to transform coordinates in the
        Composite Source Model setup.

        % function zetaetaxi2xyz
        % Usage: [x,y,z,etah]=zetaetaxi2xyz(zeta,eta,xi,strike,dip,zref);
        %
        % Input:
        %   zeta, eta, xi = location in the coordinates fixed to the fault.
        %   strike, dip = strike and dip of the fault
        %   zref = reference depth
        %
        % Output
        %   x, y, z = Cartesian coordinates, with x = east, y = north, z = up
        %   etah = Cartesian horizontal distance normal to the fault.
        %
        """
        cosdip = math.cos(math.radians(dip))
        sindip = math.sin(math.radians(dip))
        cosstk = math.cos(math.radians(strike))
        sinstk = math.sin(math.radians(strike))

        etah = [el1 * cosdip - el2 * sindip for el1, el2 in zip(eta, xi)]
        z = [zref + el1 * sindip + el2 * cosdip for el1, el2 in zip(eta, xi)]
        x = [el1 * sinstk + el2 * cosstk for el1, el2 in zip(zeta, etah)]
        y = [el1 * cosstk - el2 * sinstk for el1, el2 in zip(zeta, etah)]

#         etah = (eta * math.cos(math.radians(dip)) -
#                 xi * math.sin(math.radians(dip)))
#         z = (zref + eta * math.sin(math.radians(dip)) +
#              xi * math.cos(math.radians(dip)))

#         x = (zeta * math.sin(math.radians(strike)) +
#              etah * math.cos(math.radians(strike)))
#         y = (zeta * math.cos(math.radians(strike)) -
#              etah * math.sin(math.radians(strike)))

        return x, y, z, etah

    def xyz2latlong(self, x, y, z, reflat, reflon):
        """
        Convert cartesian coordinates to map (from John Anderson's
        xyz2latlong.m program

        % function xyz2latlong
        % usage: [lat,long,depth]=xyz2latlong(x,y,z,reflat,reflong);
        %
        % Input:
        %  x = east (km)
        %  y = north (km)
        %  z = down (km)
        %  reflat, reflong = origin of the Cartesian system.
        %   So note that the xyz coordinate system is left handed.
        %
        % Output:
        %  lat = latitude (north positive, south negative)
        %  long = longitude (east positive, west negative)
        %  depth (in km)
        %
        % Limitations:
        %  mapping is only valid to a distance of less than a few hundred km.
        %  mapping will perform very badly near poles.
        %
        """
        s_y = 6371 * math.pi / 180
        s_x = s_y * math.cos(math.radians(reflat))

        lats = [reflat + elem / s_y for elem in y]
        lons = [reflon + elem / s_x for elem in x]
        return lats, lons, z

    def calculate_params(self, a_velmodel, user_sdrp=None):
        """
        This function calculates a number of parameters for the
        rupture read in the src file. The results are stored in the
        cfgparams dictionary.
        """
        cfgparams = self.cfgparams

        # Random seeds
        random.seed(self.SEED)
        cfgparams["seed1"] = int(random.random() * 10000)
        cfgparams["seed2"] = int(random.random() * 10000)
        cfgparams["seed3"] = int(random.random() * 10000)
        # cfgparams["seed1"] = 1170
        # cfgparams["seed2"] = 2000
        # cfgparams["seed3"] = 3000

        # Constants
        cfgparams["t00"] = 0.0
        cfgparams["ang1"] = 90.0
        cfgparams["ang2"] = 0.0
        cfgparams["akp"] = 0.001
        cfgparams["ntc"] = 2

        # From the Excel file
        # cfgparams["sdrp1"] = 100
        # cfgparams["sdrp2"] = 100
        cfgparams["mu"] = 3.4E+11
        cfgparams["fdim"] = 2
        cfgparams["plim"] = 1
        cfgparams["npole"] = 2
        cfgparams["flowf"] = 0.1
        cfgparams["fhighf"] = 25

        if user_sdrp is not None:
            cfgparams["sdrp1"] = user_sdrp
            cfgparams["sdrp2"] = user_sdrp
        else:
            # sdrp is calculated from Rake
            if self.RAKE >= -180 and self.RAKE <= 30:
                tmp_rake = 125
            elif self.RAKE > 30 and self.RAKE <= 60:
                tmp_rake = 125 + (200.0-125.0) / (60.0-30.0) * (self.RAKE - 30)
            elif self.RAKE > 60 and self.RAKE <= 120:
                tmp_rake = 200
            elif self.RAKE > 120 and self.RAKE <= 150:
                tmp_rake = 200 + (125.0-200.0)/(150.0-120.0) * (self.RAKE - 120)
            elif self.RAKE > 150 and self.RAKE <= 360:
                tmp_rake = 125
            else:
                raise bband_utils.ParameterError("Rake is out of range!")

            cfgparams["sdrp1"] = tmp_rake
            cfgparams["sdrp2"] = tmp_rake

        # User parameters used in the mod files
        cfgparams["dx"] = 10
        cfgparams["dy"] = 10
        cfgparams["vrup"] = 2.0
        cfgparams["dt"] = 0.02
        cfgparams["twin"] = 81.92

        # User parameters used in simula.in
        cfgparams["simula_in"] = {}
        cfgsimula = cfgparams["simula_in"]
        cfgsimula["mfiz"] = cfgparams["npole"]
        cfgsimula["flw"] = cfgparams["flowf"]
        cfgsimula["fhi"] = cfgparams["fhighf"]
        cfgsimula["ncoda"] = 400
        cfgsimula["m"] = 512
        cfgsimula["err"] = 0.001
        cfgsimula["gi"] = 0.01
        cfgsimula["gs"] = 0.005
        cfgsimula["fmax"] = cfgsimula["fhi"]
        cfgsimula["qf"] = 150.0
        cfgsimula["cn"] = 0.5
        cfgsimula["nscat"] = 200
        cfgsimula["fl"] = 1.5
        cfgsimula["flb"] = 0.5
        cfgsimula["akap1"] = 0.02
        cfgsimula["akap2"] = 0.0
        cfgsimula["seed2"] = cfgparams["seed2"]

        # User parameters used in scat1d.in
        cfgparams["scat1d_in"] = {}
        cfgscat1d = cfgparams["scat1d_in"]
        cfgscat1d["nlay"] = 300
        cfgscat1d["damp"] = 1.0
        cfgscat1d["dh"] = 0.003
        cfgscat1d["ddh"] = 0.0026
        cfgscat1d["va"] = 2.5
        cfgscat1d["dv"] = 0.2
        cfgscat1d["dena"] = 2.8
        cfgscat1d["dden"] = 0.2
        cfgscat1d["q"] = 100.0
        cfgscat1d["mscat"] = 16
        cfgscat1d["seed3"] = cfgparams["seed3"]

        # Now calculate others
        cfgparams["L"] = self.LENGTH
        cfgparams["W"] = self.WIDTH
        cfgparams["D"] = self.DIP
        cfgparams["R"] = self.RAKE

        # Set nt, dt, twins, fmax, and fmax1, fnyq
        cfgparams["nt"] = int(math.ceil(cfgparams["twin"] / cfgparams["dt"]))
        # Make sure the max number of points we have is 8192. If that
        # is not the case, adjust nt and dt accordingly
        if cfgparams["nt"] > 8192:
            cfgparams["nt"] = 8192
            cfgparams["dt"] = cfgparams["twin"] / cfgparams["nt"]
        cfgparams["twins"] = cfgparams["twin"]
        cfgparams["fnyq"] = 0.5 / cfgparams["dt"]
        cfgparams["fmax"] = cfgparams["fnyq"]
        cfgparams["fmax1"] = cfgparams["fnyq"]

        # Now, calculate nx and ny
        zref = self.DEPTH_TO_TOP
        depth_tor = self.DEPTH_TO_TOP
        wem = (depth_tor - zref) / math.sin(math.radians(self.DIP))
        wep = self.WIDTH + wem
        # Reference point is middle of the fault, so we have 1/2
        # length in each direction
        lep = self.LENGTH / 2.0
        lem = -1 * (self.LENGTH / 2.0)

        # Using the hypocenter along stk as reference point, which is
        # incorrect as per JA: 24-Apr-2013
        #lep = (self.LENGTH / 2) - self.HYPO_ALONG_STK
        #lem = lep - self.LENGTH
        kx1 = int(round(lem / cfgparams["dx"]))
        kx2 = int(round(lep / cfgparams["dx"]))
        ky1 = int(round(wem / cfgparams["dy"]))
        ky2 = int(round(wep / cfgparams["dy"]))
        cfgparams["nx"] = len(numpy.arange(kx1 * cfgparams["dx"],
                                           (kx2 * cfgparams["dx"] + 0.000001),
                                           cfgparams["dx"]))
        cfgparams["ny"] = len(numpy.arange(ky1 * cfgparams["dy"],
                                           (ky2 * cfgparams["dy"] + 0.000001),
                                           cfgparams["dy"]))

        # Calculate latftop, longftop, depthftop
        reflat = self.LAT_TOP_CENTER
        reflon = self.LON_TOP_CENTER
        zetaftopdc = [lem, lep]
        etaftopdc = [wem, wem]
        xiftopdc = [0, 0]
        xftopdc, yftopdc, zftopdc, _ = self.zetaetaxi2xyz(zetaftopdc,
                                                          etaftopdc,
                                                          xiftopdc,
                                                          self.STRIKE,
                                                          self.DIP,
                                                          zref)
        latftop, lonftop, depthftop = self.xyz2latlong(xftopdc,
                                                       yftopdc,
                                                       zftopdc,
                                                       reflat,
                                                       reflon)
        cfgparams["tlat1"] = latftop[0]
        cfgparams["tlat2"] = latftop[1]
        cfgparams["tlon1"] = lonftop[0]
        cfgparams["tlon2"] = lonftop[1]
        cfgparams["tdep"] = depthftop[0]
        cfgparams["sdept"] = (zref +
                              self.HYPO_DOWN_DIP *
                              math.sin(math.radians(self.DIP)))

        # Calculate hypolat, hypolon, hypodepth
        hypoxi = 0
        hypox, hypoy, hypoz, _ = self.zetaetaxi2xyz([self.HYPO_ALONG_STK],
                                                    [self.HYPO_DOWN_DIP],
                                                    [hypoxi],
                                                    self.STRIKE, self.DIP,
                                                    zref)
        hypolat, hypolon, hypodepth = self.xyz2latlong(hypox, hypoy, hypoz,
                                                       reflat, reflon)
        cfgparams["hypolat"] = hypolat[0]
        cfgparams["hypolon"] = hypolon[0]
        cfgparams["hypodepth"] = hypodepth[0]

        # Calculate perw and perf (for nuclear.in)
        cfgparams["nuclear_in"] = {}
        cfgnuclear = cfgparams["nuclear_in"]
        cfgnuclear["perw"] = (self.HYPO_ALONG_STK / self.LENGTH) + 0.5
        cfgnuclear["perf"] = self.HYPO_DOWN_DIP / self.WIDTH

        # Calculate realw, rmax, rmin, rmaxcm, rmincm (for compom.in)
        cfgparams["compom_in"] = {}
        cfgcompom = cfgparams["compom_in"]

        # Calculate moment
        # cfgcompom["eqmom"] = 10**(1.5 * self.MAGNITUDE + 16.1)

        # Calculate moment using new equation
        # Units N.m
        cfgcompom["eqmom"] = 10**(1.5 * self.MAGNITUDE + 16.05 - 7)

        # convert from N.m to dyne.cm
        cfgcompom["eqmom"] = cfgcompom["eqmom"] * 10**7

        # Calculate rmax and rmin
        realw = min(self.LENGTH, self.WIDTH)
        cfgcompom["rmax"] = realw / 2.0
        cfgcompom["rmin"] = cfgcompom["rmax"] / 20.0

        # Read velocity model file
        self.parse_velmodel(a_velmodel)
        if self.nlay == 0:
            raise bband_utils.ParameterError("Unable to parse velocity model!")

        # Calculate vssp = mean(fcz(ddd, depths, vssm))

        # Calculate ddd
        zetafodc = [lep, lep, lem, lem, lep]
        etafodc = [wem, wep, wep, wem, wem]
        xifodc = [0, 0, 0, 0, 0]
        xfodc, yfodc, zfodc, _ = self.zetaetaxi2xyz(zetafodc, etafodc, xifodc,
                                                    self.STRIKE, self.DIP,
                                                    zref)
        #latfo, longfo, depthfo = self.xyz2latlong(xfodc, yfodc, zfodc,
        #                                          reflat, reflon)
        _, _, depthfo = self.xyz2latlong(xfodc, yfodc, zfodc,
                                         reflat, reflon)
        d1 = min(depthfo)
        d2 = max(depthfo)
        dd = (d2 - d1) / 20
        ddd = numpy.arange(d1, d2 + 0.0001, dd)

        # Calculate depths
        zabove = [0 for _ in self.vmodel['h']]
        zbelow = [0 for _ in self.vmodel['h']]
        for idx, val in enumerate(self.vmodel['h']):
            if idx == 0:
                zabove[0] = val
                continue
            zabove[idx] = zabove[idx - 1] + val
        zbelow[0] = 0
        zbelow[1:] = zabove[:-1]
        zabove = [val * 0.995 for val in zabove]
        zbelow = [val * 1.005 for val in zbelow]
        depths = []
        for zab, zbe in zip(zabove, zbelow):
            depths.append(zbe)
            depths.append(zab)

        # Calculate vss
        vss_low = [vel * 0.995 for vel in self.vmodel['vs']]
        vss_hi = [vel * 1.005 for vel in self.vmodel['vs']]
        vss = []
        for vlow, vhi in zip(vss_low, vss_hi):
            vss.append(vlow)
            vss.append(vhi)

        # Now, finally get vssp
        cfgparams["vssp"] = numpy.mean(list(numpy.interp(ddd,
                                                         depths,
                                                         vss)))

        # Now work on the csevents file
        fdim = cfgparams["fdim"]
        plim = cfgparams["plim"]
        sdrp = (cfgparams["sdrp1"] +
                (cfgparams["sdrp2"] - cfgparams["sdrp1"]) * random.random())
        rmaxcm = cfgcompom["rmax"] * (10**5)
        rmincm = cfgcompom["rmin"] * (10**5)
        if fdim <= 3.1 and fdim >= 2.9:
            p = ((7.0 / 16) *
                 cfgcompom["eqmom"] /
                 (sdrp * (10 ** 6)) *
                 (cfgcompom["rmin"] / cfgcompom["rmax"]))
        else:
            cfdim = 3.0 - fdim
            p = ((7.0 / 16) *
                 cfgcompom["eqmom"] /
                 (sdrp * (10 ** 6)) *
                 (cfdim / (rmaxcm**cfdim - rmincm**cfdim)))
        nsub = int(math.floor((p / fdim) *
                              (rmincm**(-fdim) -
                               rmaxcm**(-fdim))))
        if nsub > 40000:
            raise bband_utils.ParameterError("Too many subevents: %d" %
                                             nsub)
        # Create subevents
        nc = [nsub * random.random() for _ in range(0, nsub)]
        srcm = [(fdim * val / p + rmaxcm**(-fdim))**(-1.0/fdim) for val in nc]
        srkm = [val / 10**5 for val in srcm]
        sx = [plim * val +
              (cfgparams["L"] - 2 * plim * val) *
              random.random() for val in srkm]
        sy = [val + (cfgparams["W"] - (1 + plim) * val) *
              random.random() for val in srkm]
        submom = [(16.0 / 7) * sdrp * 10**6 * val ** 3 for val in srcm]
        # submw = [(2.0 / 3) * (math.log10(val) - 16.1) for val in submom]
        csmom = sum(submom)
        ratiom = csmom / cfgcompom["eqmom"]
        sdrpa = sdrp / ratiom
        submom = [(16.0 / 7) * sdrpa * 10**6 * val ** 3 for val in srcm]
        submw = [(2.0 / 3) * (math.log10(val) - 16.1) for val in submom]
        cfgparams["csevents_dat"] = {}
        cfg_csevents = cfgparams["csevents_dat"]
        cfg_csevents["sdrpa"] = sdrpa
        cfg_csevents["mu"] = cfgparams["mu"]
        cfg_csevents["nsub"] = nsub
        cfg_csevents["sx"] = sx
        cfg_csevents["sy"] = sy
        cfg_csevents["srkm"] = srkm
        cfg_csevents["submw"] = submw
        cfg_csevents["zetafo"] = zetafodc
        cfg_csevents["etafo"] = etafodc

    def __init__(self, a_srcname=None):
        """
        Set up parameters for the CSM method
        """
        if a_srcname:
            self.parse_src(a_srcname)

if __name__ == "__main__":
    CSM_CFG = CSMCfg(sys.argv[1])
    print("Created Test Config Class: %s" % (os.path.basename(sys.argv[0])))
