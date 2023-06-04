#!/usr/bin/env python
# This is the main module
# when you import pynga
# what it does is to do the following statements

# Note: NGA08 provides GMRotI50, while NGA14 provides RotD50, so
# before do the comparison, do the conversion

from __future__ import division, print_function

import os
import numpy as np

# GMPE Package content
import pynga.CB08 as CB08
import pynga.BA08 as BA08
import pynga.CY08 as CY08
import pynga.AS08 as AS08
import pynga.SC08 as SC08

import pynga.BSSA14 as BSSA14
import pynga.CB14 as CB14
import pynga.CY14 as CY14
import pynga.ASK14 as ASK14

# CENA GROUP1 GMPEs
import pynga.PZT11 as PZT11
import pynga.A0811E as A0811E
import pynga.S03SCVS as S03SCVS

# PyNGA
from pynga.utils import mapfunc, logline

# NGA08 Period list (available for each NGA models)
# -1.0: PGA; -2.0: PGV
TsDict = {'BA': [0.01, 0.02, 0.03, 0.05, 0.075, 0.10, 0.15, 0.20, 0.25,
                 0.30, 0.40, 0.50, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0,
                 7.5, 10.0, -1, -2],
          'CB': [0.01, 0.02, 0.03, 0.05, 0.075, 0.10, 0.15, 0.20, 0.25,
                 0.30, 0.40, 0.50, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0,
                 7.5, 10.0, -1, -2],
          'CY': [0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.10, 0.15, 0.20,
                 0.25, 0.30, 0.40, 0.50, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0,
                 5.0, 7.5, 10.0, -1, -2],
          'AS': [0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.10, 0.15, 0.20,
                 0.25, 0.30, 0.40, 0.50, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0,
                 5.0, 7.5, 10.0, -1, -2],}

# ============================================
# Integrated function for CENA GROUP 1 models
# ============================================
def CENA1(model_name, Mw, Rjb, Rrup, period):
    """
    Combined function to calculate CENA1 models
    """
    if model_name == "PZT11":
        model = PZT11.PZT11()
        dist = Rrup
    elif model_name == "A0811E":
        model = A0811E.A0811E()
        dist = Rjb
    elif model_name == "S03SCVS":
        model = S03SCVS.S03SCVS()
        dist = Rjb
    else:
        print("Invalid CENA1 model_name")
        raise ValueError

    # Check if period requested is supported
    if period not in model.periods:
        print("Period requested is invalid")
        raise ValueError

    # Calculate median
    value = model(Mw, dist, period)

    return value

# ============================================
# Integrated function for NGA 2008 models
# ============================================
def NGA08(model_name, Mw, Rjb, Vs30, period, epislon=0, NGAs=None,
          rake=None, Mech=0, Ftype=None, Fnm=None, Frv=None,
          dip=None, W=None, Ztor=None, Zhypo=None, Fas=0,
          Rrup=None, Rx=None, Fhw=None, azimuth=None,
          VsFlag=0, Z25=None, Z15=None, Z10=None,
          AS09=None, AB11=None, ArbCB=0):
    """
    Combined function to compute median and standard deviation

    Arguments (has to be specified)
    ----------
    model_name : choose NGA model you want to use (AS,BA,CB,CY)
    Mw : moment magnitude
    Rjb: Joyner-Boore distance in km
         defined as the shortest distance from a site to the
         surface projection of the rupture surface
    Vs30: The average shear-wave velocity between 0 and 30-meters
          depth (site condition) in m/s
    period: period at which you want to use NGA
            This function allow to use periods that are not in
            the available periods (refer to TsDict)

    Keywords
    --------
    [*] shows the default value

    # ================
    # General Keywords
    # ================
    epislon : deviation from the median value [0]
    NGAs : dictionary to select terms in NGA models and use updated coefficents
              default:
                 {'CB':{'NewCoefs':None,'terms':(1,1,1,1,1,1)},
                  'BA':{'NewCoefs':None,'terms':(1,1,1)},
                  'CY':{'NewCoefs':None,'terms':(1,1,1,1,1,1)},
                  'AS':{'NewCoefs':None,'terms':(1,1,1,1,1,1,1)}}

    # ===============
    # Source Keywords
    # ===============
    rake: rake angle (in degree) [None]
          used to determine the fault type
    Mech: Used in BA model [3]
          (0:Strike-slip, 1:Normal, 2:Reverse, 3:Unknown
    Ftype: fault type string [None]
          'SS': Strike-slip, 'NM': Normal, 'RV': Reverse,
           'U': Unknown (unknown is only used in BA model)
    Fnm : 0: not a normal fault; 1: Normal [None]
          default: None
    Frv : 0: not a reverse fault; 1: reverse [None]
          default: None
    dip : dip angle of the fault plane [None]
          default: None
    W : Rupture width (down-dip) [None]
    Ztor : depth to the top of rupture [None]
    Zhypo: depth to the hypocenter location [None]
    Fas : Aftershock flag [None]
          0: Mainshock; 1: Aftershock

    # ================
    # Path Keywords
    # ================
    Rrup: Rupture distance in km [None]
          defined as the distance from a site the to the fault plane
          For simple fault geometry, function calc_Rrup in utils.py can
          be used to compute Rrup, otherwise
          use DistanceToEvenlyGriddedSurface function in utils.py
          to compute given fault geometry and site location
    Rx :  horizontal distance between a site and fault trace, in km [None]
          defined by extending the fault trace (or the top edge of
          the rupture) to infinity in both directions.
          For simple fault geometry, function calc_Rx in utils.py
          can be used to compute Rrup, otherwise,
          use DistanceToEvenlyGriddedSurface function in utils.py to
          compute given fault geometry and site location
    Fhw : hanging wall flag [None]
          0: in footwall; 1: in hanging wall
    azimuth: source-to-site azimuth [None]
           defined as the angle between the positive fault strike
           direction and the line connecting
           a site to the closet point on the surface projection
           of the top edge of rupture (clockwise)
           (used in simple fault geometry)

    # =================
    # Site Keywords
    # =================
    VsFlag : Vs30 inferred or measured flag [0]
            0: inferred Vs30; 1: measured Vs30
    Z25: basin depth to S wave velocity equal to 2.5 km/s [None], in km
         Z25 could be estimated by using calc_Z25 function
         in utils.py given Vs30
    Z15: basin depth to S wave velocity equal to 1.5 km/s [None], in km
         used to estimate Z2.5 when Z2.5 = None
    Z10: basin depth to S wave velocity equal to 1.0 km/s [None], in meter
         Z10 could be estimated by using calc_Z1 function in
         utils.py given Vs30

    # =================
    # Updated models
    # =================
    AS09 : Abrahamson and Silva 2009 updated model (taper5
           hanging wall effect) [None]
    AB11 : Atkinson and Boore 2011 updated model with correction
           term (after more small magnitude events recordings)

    # =================
    # Other Keywords
    # =================
    ArbCB: Campbell and Bozorgnia 2008 model standard deviation [0]
           0: output total standard deviation is for GMRotIpp
              intensity measures (rotation-independent)
           1: output total standard deviation is for arbitrary
              horizontal component

    """

    if NGAs == None:
        NGAs = {'CB':{'NewCoefs':None, 'terms':(1, 1, 1, 1, 1, 1)},
                'BA':{'NewCoefs':None, 'terms':(1, 1, 1)},
                'CY':{'NewCoefs':None, 'terms':(1, 1, 1, 1, 1, 1)},
                'AS':{'NewCoefs':None, 'terms':(1, 1, 1, 1, 1, 1, 1)}}

    dict1 = NGAs
    itmp = 0

    # check the input period
    if period > 10.0 or 0 < period < 0.01:
        print("Positive period value should be within [0.01,10] "
              "for SA at corresponding periods")
        raise ValueError
    if period < 0 and period not in [-1, -2]:
        print('negative period should be -1,-2 for PGA and PGV')
        raise ValueError

    if model_name == 'BA':
        ngaM = BA08.BA08_nga()
        kwds = {'Mech':Mech, 'Ftype':Ftype, 'AB11':AB11,
                'CoefTerms':dict1[model_name]} # OpenSHA doesn't have this
    if model_name == 'CB':
        ngaM = CB08.CB08_nga()
        kwds = {'Ftype':Ftype, 'Rrup':Rrup, 'Ztor':Ztor, 'dip':dip,
                'Z25':Z25, 'W':W, 'Zhypo':Zhypo, 'azimuth':azimuth,
                'Fhw':Fhw, 'Z10':Z10, 'Z15':Z15, 'Arb':ArbCB,
                'CoefTerms':dict1[model_name]}
    if model_name == 'CY':
        ngaM = CY08.CY08_nga()
        kwds = {'Ftype':Ftype, 'Rrup':Rrup, 'Rx':Rx, 'Ztor':Ztor,
                'dip':dip, 'W':W, 'Zhypo':Zhypo, 'azimuth':azimuth,
                'Fhw':Fhw, 'Z10':Z10, 'AS':Fas, 'VsFlag':VsFlag,
                'CoefTerms':dict1[model_name]}
    if model_name == 'AS':
        ngaM = AS08.AS08_nga()
        kwds = {'Ftype':Ftype, 'Rrup':Rrup, 'Rx':Rx, 'Ztor':Ztor,
                'dip':dip, 'W':W, 'Zhypo':Zhypo, 'azimuth':azimuth,
                'Fhw':Fhw, 'Z10':Z10, 'Fas':Fas, 'VsFlag':VsFlag,
                'CoefTerms':dict1[model_name]}

    # Common interpolation and calculation for all models
    periods = np.array(ngaM.periods)
    for ip in range(len(periods)):
        if abs(period-periods[ip]) < 0.0001:
            # period is within the periods list
            itmp = 1
            break

    if itmp == 1:
        # compute median, std directly for the existing period in the
        # period list of the NGA model
        values = mapfunc(ngaM, Mw, Rjb, Vs30, period, rake, **kwds)
        values = np.array(values, dtype=object)

    if itmp == 0:
        # do the interpolation for periods that is not in the period
        # list of the NGA model
        ind_low = (periods < period).nonzero()[0]
        ind_high = (periods > period).nonzero()[0]
        period_low = max(periods[ind_low])
        period_high = min(periods[ind_high])
        values_low = np.array(mapfunc(ngaM, Mw, Rjb, Vs30,
                                      period_low, rake, **kwds))
        values_high = np.array(mapfunc(ngaM, Mw, Rjb, Vs30,
                                       period_high, rake, **kwds))
        N1, N2 = np.array(values_low).shape
        values = np.zeros((N1, N2))
        for icmp in range(N2):
            if icmp != 0:
                # stardand values are in ln (g)
                values[:, icmp] = logline(np.log(period_low),
                                          np.log(period_high),
                                          values_low[:, icmp],
                                          values_high[:, icmp],
                                          np.log(period))
            else:
                # median value is in g
                values[:, icmp] = logline(np.log(period_low),
                                          np.log(period_high),
                                          np.log(values_low[:, icmp]),
                                          np.log(values_high[:, icmp]),
                                          np.log(period))
                # change the median into g unit (logline gives the
                # result in ln(g))
                values[:, icmp] = np.exp(values[:, icmp])

    # outputs
    NGAsigmaT = np.array(values[:, 1]).astype(float)
    NGAtau = np.array(values[:, 2]).astype(float)
    NGAsigma = np.array(values[:, 3]).astype(float)
    if epislon:
        NGAmedian = np.exp(np.log(values[:, 0]) + epislon * NGAsigmaT)
    else:
        NGAmedian = values[:, 0]

    # returned quantities are all in g, not in log(g), event for the
    # standard deviations
    # all in g, include the standard deviation
    return NGAmedian, np.exp(NGAsigmaT), np.exp(NGAtau), np.exp(NGAsigma)

def BA08Test(T):
    # to reproduce BA model (shown in Earthquake Spectra 2008)
    import matplotlib.pyplot as plt
    NGAs = {'CB':{'NewCoefs':None, 'terms':(1, 1, 1, 1, 1, 1)},
            'BA':{'NewCoefs':None, 'terms':(1, 1, 1)},
            'CY':{'NewCoefs':None, 'terms':(1, 1, 1, 1, 1, 1)},
            'AS':{'NewCoefs':None, 'terms':(1, 1, 1, 1, 1, 1, 1)}}

    # validation with BA
    nga = 'BA'
    Mws = [5, 6, 7, 8]
    Mws = [4,]
    Vs30 = 760
    FT = 'U'
    Rjb = np.arange(0.1, 100, 0.5)

    fig = plt.figure(1)
    ax = fig.add_subplot(111)
    lines = []
    for Mw in Mws:
        median, std, tau, sigma = NGA08(nga, Mw, Rjb, Vs30, T,
                                        Mech=1, NGAs=NGAs)
        line = ax.loglog(Rjb, median * 100 * 9.8)
        lines.append(line)
    ax.legend(lines, ('M=5', 'M=6', 'M=7', 'M=8'), loc=0)
    ax.set_title(r"T=%s, $V_{S30}$ = 760 m/s, mech='SS'" % ('%.2f' % T))
    ax.set_xlabel(r'$R_{JB}$ (km)')
    ax.set_ylabel(r'5%-damped PSA (cm/s)')
    plt.show()

def NGA08test(nga):
    # simple test comparing with file:
    # ./Validation/NGAmodelsTestFiles/nga_Sa_v19a.xls
    M = 6.93
    Ztor = 3
    Ftype = 'RV'
    W = 3.85
    dip = 70
    Rrup = Rjb = Rx = 30
    Fhw = 0
    Vs30 = 760
    Z10 = 0.024 * 1000   # in meter
    Z25 = 2.974    # in km
    VsFlag = 0

    periods = TsDict[nga]
    NT = len(periods)
    Medians = []; SigmaTs = []
    for ip in range(NT):
        Ti = periods[ip]
        median, std, tau, sigma = NGA08(nga, M, Rjb, Vs30, Ti, Ftype=Ftype,
                                        W=W, Ztor=Ztor, dip=dip, Rrup=Rrup,
                                        Rx=Rx, Fhw=Fhw, Z10=Z10, Z25=Z25,
                                        VsFlag=VsFlag)
        Medians.append(median)
        SigmaTs.append(np.log(std))
    output = np.c_[np.array(periods), np.array(Medians), np.array(SigmaTs)]
    pth = './tmp'
    if not os.path.exists(pth):
        os.mkdir(pth)
    np.savetxt(pth + '/NGA08_SimpleTest%s.txt' % nga, output)

    print(output)

# NGA 14 period list
# -1: PGA; -2: PGV
TsDict14 = {'BSSA': [0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.10, 0.15,
                     0.20, 0.25, 0.30, 0.40, 0.50, 0.75, 1.0, 1.5, 2.0,
                     3.0, 4.0, 5.0, 7.5, 10.0, -1, -2],
            'CB': [0.01, 0.02, 0.03, 0.05, 0.075, 0.10, 0.15, 0.20, 0.25,
                   0.30, 0.40, 0.50, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0,
                   7.5, 10.0, -1, -2],
            'CY': [0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.10, 0.12, 0.15,
                   0.17, 0.20, 0.25, 0.30, 0.40, 0.50, 0.75, 1.0, 1.5,
                   2.0, 3.0, 4.0, 5.0, 7.5, 10.0, -1],
            'ASK': [0.01, 0.02, 0.03, 0.05, 0.075, 0.10, 0.15, 0.20, 0.25,
                    0.30, 0.40, 0.50, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0,
                    5.0, 6.0, 7.5, 10.0, -1, -2],}

# ============================================
# Integrated function for NGA 2014 models
# =============================================
def NGA14(model_name, Mw, Rjb, Vs30, period, epislon=0, NGAs=None,
          rake=None, Mech=3, Ftype=None, Fnm=None, Frv=None,
          dip=None, W=None, Ztor=None, Zhypo=None, Fas=0,
          Rrup=None, Rx=None, Fhw=None, azimuth=None,
          VsFlag=0, Z25=None, Z15=None, Z10=None,
          ArbCB=0, SJ=0,
          country='California', region='CA',
          Dregion='GlobalCATW',
          CRjb=15, Ry0=None,
          D_DPP=0):

    if NGAs == None:
        NGAs = {'CB':{'NewCoefs':None, 'terms':(1, 1, 1, 1, 1, 1, 1, 1, 1)},
                'BSSA':{'NewCoefs':None, 'terms':(1, 1, 1)},
                'CY':{'NewCoefs':None, 'terms':(1, 1, 1, 1, 1, 1, 1)},
                'ASK':{'NewCoefs':None, 'terms':(1, 1, 1, 1, 1, 1, 1)}}

    dict1 = NGAs
    itmp = 0

    # check the input period
    # Note: this function is better used at a
    # given period with a set of other parameters (not with a set of
    # periods)
    if period > 10.0 or 0 < period < 0.01:
        print("Positive period value should be within [0.01,10] "
              "for SA at corresponding periods")
        raise ValueError
    if period < 0 and period not in [-1, -2]:
        print('negative period should be -1,-2 for PGA and PGV')
        raise ValueError

    if model_name == 'BSSA':
        ngaM = BSSA14.BSSA14_nga()
        kwds = {'Mech':Mech, 'Ftype':Ftype, 'Z10':Z10,
                'Dregion':Dregion, 'country':country,
                'CoefTerms':dict1[model_name]}

    if model_name == 'CB':
        ngaM = CB14.CB14_nga()
        kwds = {'Ftype':Ftype, 'Rrup':Rrup, 'Ztor':Ztor, 'dip':dip,
                'Z25':Z25, 'W':W, 'Zhypo':Zhypo, 'azimuth':azimuth,
                'Fhw':Fhw, 'Z10':Z10, 'Z15':Z15, 'Arb':ArbCB,
                'SJ':SJ, 'region':region, 'CoefTerms':dict1[model_name]}

    if model_name == 'CY':
        ngaM = CY14.CY14_nga()
        kwds = {'Ftype':Ftype, 'Rrup':Rrup, 'Rx':Rx, 'Ztor':Ztor,
                'dip':dip, 'W':W, 'Zhypo':Zhypo, 'azimuth':azimuth,
                'Fhw':Fhw, 'Z10':Z10, 'AS':Fas, 'VsFlag':VsFlag,
                'country':country, 'D_DPP':D_DPP,
                'CoefTerms':dict1[model_name]}
        # the new CY model treat PGA = SA(0.01)
        if period == -1:
            period = 0.01

    if model_name == 'ASK':
        ngaM = ASK14.ASK14_nga()
        kwds = {'Ftype':Ftype, 'Rrup':Rrup, 'Rx':Rx, 'Ztor':Ztor,
                'dip':dip, 'W':W, 'Zhypo':Zhypo, 'azimuth':azimuth,
                'Fhw':Fhw, 'Z10':Z10, 'Fas':Fas, 'CRjb':CRjb, 'Ry0':Ry0,
                'region':region, 'country':country, 'VsFlag':VsFlag,
                'CoefTerms':dict1[model_name]}

    # common interpolate for all models
    periods = np.array(ngaM.periods)
    for ip in range(len(periods)):
        if abs(period-periods[ip]) < 0.0001:
            # period is within the periods list
            itmp = 1
            break

    if itmp == 1:
        # compute median, std directly for the existing period in the
        # period list of the NGA model
        values = mapfunc(ngaM, Mw, Rjb, Vs30, period, rake, **kwds)
        values = np.array(values)

    if itmp == 0:
        #print 'do the interpolation for periods that is not in the
        #period list of the NGA model'
        ind_low = (periods <= period * 1.0).nonzero()[0]
        ind_high = (periods >= period * 1.0).nonzero()[0]

        period_low = max(periods[ind_low])
        period_high = min(periods[ind_high])
        values_low = np.array(mapfunc(ngaM, Mw, Rjb, Vs30,
                                      period_low, rake, **kwds))
        values_high = np.array(mapfunc(ngaM, Mw, Rjb, Vs30,
                                       period_high, rake, **kwds))
        N1, N2 = np.array(values_low).shape
        values = np.zeros((N1, N2))
        for icmp in range(N2):
            if icmp != 0:
                # stardand values are in ln (g)
                values[:, icmp] = logline(np.log(period_low),
                                          np.log(period_high),
                                          values_low[:, icmp],
                                          values_high[:, icmp],
                                          np.log(period))
            else:
                # median value is in g
                values[:, icmp] = logline(np.log(period_low),
                                          np.log(period_high),
                                          np.log(values_low[:, icmp]),
                                          np.log(values_high[:, icmp]),
                                          np.log(period))
                # change the median into g unit (logline gives the
                # result in ln(g))
                values[:, icmp] = np.exp(values[:, icmp])

    # outputs
    NGAsigmaT = values[:, 1]
    NGAtau = values[:, 2]
    NGAsigma = values[:, 3]

    if epislon:
        NGAmedian = np.exp(np.log(values[:, 0]) + epislon * NGAsigmaT)
    else:
        NGAmedian = values[:, 0]

    # returned quantities are all in g, not in log(g), event for the
    # standard deviations
    # all in g, include the standard deviation
    return NGAmedian, np.exp(NGAsigmaT), np.exp(NGAtau), np.exp(NGAsigma)

def BSSA14_validation(infile, outfile, iset):
    import matplotlib.pyplot as plt
    # read in files (mainly parameters for run using pynga)
    hdrs = open(infile, 'r').readlines()[3].strip().split()
    inputs = {}
    data = np.loadtxt(infile, skiprows=4)
    for ih in range(len(hdrs)):
        hdr = hdrs[ih]
        inputs[hdr] = data[:, ih]
    regionDict = {'0':'GlobalCATW', '1':'GlobalCATW',
                  '2':'ChinaTurkey', '3':'ItalyJapan'}

    # calculate and save to file (or plot directly) (following the same format)
    #fid = open(outfile,'w')
    BSSAnga = BSSA14.BSSA14_nga()
    BAnga = BA08.BA08_nga()
    Nl = len(inputs['T'])
    Y = []
    sig_lnY = []
    tau = []
    sigma = []
    Y1 = []
    sig_lnY1 = []
    tau1 = []
    sigma1 = []
    Nls1 = []
    Nls2 = []
    Nls3 = []
    rake = None
    for il in range(Nl):
        for key in ['T', 'M', 'Rjb', 'V30', 'mech', 'iregion', 'z1']:
            cmd = "%s = inputs['%s'][%d]" % (key, key, il)
            exec(cmd)
        Dregion = regionDict[str(int(iregion))]
        if z1 == -1.0: Z10 = None
        if z1 != -1.0: Z10 = z1
        kwds14 = {'Mech':int(mech), 'Dregion':Dregion, 'Z10':Z10}
        kwds08 = {'Mech':int(mech)}
        if T == -1.0:
            T = -2
        if T == 0.0:
            T = -1
        if T not in TsDict14['BSSA']:
            pass
        else:
            Nls1.append(il)
            Y0, sT, tau0, sigma0 = BSSAnga(M, Rjb, V30, T, rake, **kwds14)
            Y.append(Y0)
            sig_lnY.append(sT)
            tau.append(tau0)
            sigma.append(sigma0)
        if T not in TsDict['BA']:
            pass
        else:
            Nls3.append(il)
            Y0, sT, tau0, sigma0 = BAnga(M, Rjb, V30, T, rake, **kwds08)
            Y1.append(Y0)
            sig_lnY1.append(sT)
            tau1.append(tau0)
            sigma1.append(sigma0)
        Nls2.append(il)
    pyNGAs = [Y, sig_lnY, tau, sigma]
    pyNGAs08 = [Y1, sig_lnY1, tau1, sigma1]
    pyNGAs = np.array(pyNGAs)
    pyNGAs08 = np.array(pyNGAs08)
    ftNGAs = np.array([inputs['Y(g)'], inputs['sigma'],
                       inputs['tau'], inputs['phi']])
    # plot
    fig = plt.figure(1)
    texts = ['IM', r'$\sigma_T$', r'$\tau$', r'$\sigma$']
    for iax in range(len(pyNGAs)):
        ax = fig.add_subplot(2, 2, iax + 1)
        ax.plot(Nls1, pyNGAs[iax], 'bx', label='pyNGA14')
        ax.plot(Nls3, pyNGAs08[iax], 'r+', label='pyNGA08')
        ax.plot(Nls2, ftNGAs[iax], 'k.', label='orgNGA14')
        ax.set_xlabel('points')
        ax.set_ylabel('values')
        ax.legend(loc=0)
        ax.text(0.9, 0.9, texts[iax], transform=ax.transAxes)
    pltpth = './NGA_west2/validation/BSSA14/outputs'
    pltnam = pltpth + '/validation_BSSA14_set%s.png' % iset
    fig.savefig(pltnam)
    plt.show()

def NGA_Test():
    # common set to test and compare
    M = 6.93
    Ztor = 3
    Ftype = 'RV'
    Mech = 3
    W = 3.85
    rake = 90
    dip = 70
    Rrup = Rjb = Rx = 30
    Fhw = 0
    Vs30 = 760
    Vs30 = 128.
    Z10 = 0.024 * 1000   # in meter
    Z25 = 2.974    # in km
    VsFlag = 0

    # for NGA 08
    for nga in ['BA', 'CB', 'CY', 'AS']:
        periods = TsDict[nga]
        NT = len(periods)
        Medians = []
        SigmaTs = []
        for ip in range(NT):
            Ti = periods[ip]
            median, std, tau, sigma = NGA08(nga, M, Rjb, Vs30, Ti,
                                            Ftype=Ftype, W=W, Ztor=Ztor,
                                            dip=dip, Rrup=Rrup, Rx=Rx,
                                            Fhw=Fhw, Z10=Z10, Z25=Z25,
                                            VsFlag=VsFlag)
            Medians.append(median)
            SigmaTs.append(np.log(std))
        output = np.c_[np.array(periods), np.array(Medians), np.array(SigmaTs)]
        pth = './tmp'
        if not os.path.exists(pth):
            os.mkdir(pth)
        np.savetxt(pth + '/NGA08_SimpleTest%s.txt' % nga, output)

    # NGA 14
    for nga in ['BSSA', 'CB', 'CY', 'ASK']:
        periods = TsDict14[nga]
        NT = len(periods)
        Medians = []
        SigmaTs = []
        for ip in range(NT):
            Ti = periods[ip]
            median, std, tau, sigma = NGA14(nga, M, Rjb, Vs30, Ti,
                                            Ftype=Ftype, Mech=Mech,
                                            rake=rake, W=W, Ztor=Ztor,
                                            dip=dip, Rrup=Rrup, Rx=Rx,
                                            Fhw=Fhw, Z10=Z10, Z25=Z25,
                                            VsFlag=VsFlag)
            Medians.append(median)
            SigmaTs.append(np.log(std))
        output = np.c_[np.array(periods), np.array(Medians), np.array(SigmaTs)]
        pth = './tmp'
        if not os.path.exists(pth):
            os.mkdir(pth)
        np.savetxt(pth + '/NGA14_SimpleTest%s.txt' % nga, output)
        #print output

def PlotTest():
    import matplotlib.pyplot as plt
    # Debug the period for CY
    pth = './tmp'
    nga1 = ['BA', 'CB', 'CY', 'AS']
    nga2 = ['BSSA', 'CB', 'CY', 'ASK']
    fig = plt.figure(1)
    for i in range(4):
        ax = fig.add_subplot(2, 2, i + 1)
        inputs = np.loadtxt(pth + '/NGA08_SimpleTest%s.txt' % nga1[i])
        inputs1 = np.loadtxt(pth + '/NGA14_SimpleTEst%s.txt' % nga2[i])
        Ts = inputs[:-1, 0]
        values = inputs[:-1, 1]
        ax.semilogx(Ts, values, 'b+', label='%s08' % nga1[i])
        Ts = inputs1[:-1, 0]
        values = inputs1[:-1, 1]
        ax.semilogx(Ts, values, 'rx', label='%s14' % nga2[i])
        ax.legend(loc=0)
        ax.set_xlabel('period')
        ax.set_ylabel('SA (g)')
    fig.savefig(pth+'/ComparisonsNGA08_NGA14.png')

def BSSA14_test(Ti):
    # simple test comparing with file:
    # ./Validation/NGAmodelsTestFiles/nga_Sa_v19a.xls
    Rjb = Rrup = 20.
    Vs30 = 760.
    Mw = 6
    rake = 0.
    Ftype = 'SS'
    Mech = 1
    CoefTerms = {'terms':(1, 1, 1), 'NewCoefs':None}
    kwds = {'Mech':Mech, 'Ftype':Ftype, 'Z10':None,
            'Dregion':'GlobalCATW', 'country':'California',
            'CoefTerms':CoefTerms}
    BSSAnga = BSSA14.BSSA14_nga()    # BA08nga instance
    # debug mode (show each term)

    IM, sigmaT, tau, sigma = BSSAnga(Mw, Rjb, Vs30, Ti, rake, **kwds)
    print(Ti, 'BSSA14:', IM, sigmaT, tau, sigma)
    VsFlag = 0
    Z10 = Z25 = None
    dip = 90
    Ztor = 3
    W = 10
    Fhw = 0
    median, std, tau, sigma = NGA14('BSSA', Mw, Rjb, Vs30, Ti,
                                    Ftype=Ftype, Mech=Mech,
                                    rake=rake, W=W, Ztor=Ztor,
                                    dip=dip, Rrup=Rrup, Fhw=Fhw,
                                    Z10=Z10, Z25=Z25, VsFlag=VsFlag)
    print(Ti, 'NGA14_BSSA:', median, std, tau, sigma)

# ====================
# self_application
# ====================
if __name__ == '__main__':

    import sys
    opt = sys.argv[1]
    if opt == 'BA08':
        BA08Test(0.3)

    if opt == 'NGA08':
        nga = sys.argv[2]   # choose one NGA model in NGA08
        NGA08test(nga)

    if opt == 'BSSA14':
        # validation of code with BSSA outputs
        opt1 = sys.argv[3]   # 1, 2, 3 to choose reference files
        wrkpth = r'H:\local\pylib\pynga\NGA_west2\validation\BSSA14'
        inpth = wrkpth + r'\inputs'
        outpth = wrkpth + r'\outputs'
        if opt1 == '1':
            # set 1:
            file0 = r'\bssa14_vs_period_r_20_v30_760_mech_1.out'
        if opt1 == '2':
            # set 2:
            file0 = r'\bssa14_vs_period_r_20_v30_200_mech_1.out'   # (period, magnitude, distance)
        if opt1 == '3':
            # set 3:
            file0 = r'\bssa14_vs_rjb_m_4_5_6_7_8_8.5.vs30_760_mech_1.out'   # (period, magnitude)
        infile = inpth + file0
        outfile = outpth + file0
        BSSA14_validation(infile, outfile, int(opt1))

    if opt == 'NGAComparison':
        NGA_Test()

    if opt == 'PlotTest':
        PlotTest()

    if opt == 'BSSA':
        BSSA14_test(0.5)
        BSSA14_test(0.75)
