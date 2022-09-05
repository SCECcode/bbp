#!/usr/bin/env python
"""
Utilities used in NGA classes
"""
from __future__ import division, print_function

import os
import sys
import time
import numpy as np

# ===================
# General Functions
# ===================
def cpt_sqrt(a, b):
    a = np.array(a)
    b = np.array(b)
    return np.sqrt(a**2 + b**2)

def RMScalc(V1, V2, Ratio=True):
    N = len(V1)
    if len(V2) != N:
        print('length of array2 should be the same as array1')
        raise ValueError
    else:
        if 0 in V2:
            print('elements in array2 should not be zeros')
            raise ValueError
        if Ratio:
            Error = (V1 - V2) / V2
        else:
            Error = (V1 - V2)
        RMS = np.sqrt(sum(Error**2) / N)
    return RMS

def logline(x1, x2, y1, y2, x):
    # linear interpolation
    k = (y2 - y1) / (x2 - x1)
    C = y1 - k * x1
    y = k * x + C
    return y

def GetKey(key):
    return '%.3f' % (key)

def HourMinSecToSec(BlockName=None):
    hour, minute, sec = time.localtime()[3:6]
    sec1 = hour * 60 * 60 + minute * 60 + sec
    if  BlockName != None:
        print('%s' % BlockName)
    return sec1

def SecToHourMinSec(sec1, BlockName=None):
    hour = sec1 // 3600
    minute = (sec1 - hour * 3600) // 60
    sec = sec1 - hour * 3600 - minute * 60
    if BlockName == None:
        BlockName = 'the above block'
    print('Time cost of %s is %s hr %s min %s sec' %
          (BlockName, hour, minute, sec))
    return hour, minute, sec

# save and load (useful metadata operations
class namespace(object):
    """
    Namespace with object attributes initialized from a dict
    Turn Dictionary keys into object attributes
    d['KeyName'] -> d.KeyName
    d is a dictioanry
    """
    def __init__(self, d):
        self.__dict__.update(d)

# load dictionary saved in the python files
def load(f, d=None):
    """
    load variables from Python source files
    Input:
        f: file object (fid)
           or filename with full path (string)
	d: dictionary ( None default )
    Output:
        namespace(d) : updated dictionary
    """
    if type(f) is not file:
        f = open(os.path.expanduser(f))   # get the file object
    if d is None:
        d = {}
    exec(f) in d
    return namespace(d)

# save dictionary into python file (used in metadata manipulation)
def save(fd, d, expand=None, keep=None, header=''):
    """
    Write variables from a dict into a Python source file.
    """
    if type(d) is not dict:
        d = d.__dict__

    if expand is None:
        expand = []
    out = header
    for k in sorted(d):
        if k not in expand and (keep is None or k in keep):
            out += '%s = %r\n' % (k, d[k])

    for k in expand:
        print('test expand')
        if k in d:
            if type(d[k]) is tuple:
                out += k + ' = (\n'
                for item in d[k]:
                    out += '    %r,\n' % (item,)
                out += ')\n'
            elif type(d[k]) is list:
                out += k + ' = [\n'
                for item in d[k]:
                    out += '    %r,\n' % (item,)
                out += ']\n'
            elif type(d[k]) is dict:
                out += k + ' = {\n'
                for item in sorted(d[k]):
                    out += '    %r: %r,\n' % (item, d[k][item])
                out += '}\n'
            else:
                sys.exit('Cannot expand %s type %s' % (k, type(d[k])))
    if fd is not None:
        if type(fd) is not file:
            fd = open(os.path.expanduser(fd), 'w')
        fd.write(out)
    return out

# this function will be used a lot (general map function as in Python)
def mapfunc(func, *args, **kwds):
    """
    Modified function map tool
    Account for the single argument and multiple arguments
    Account for the keywords input
    """
    # arguments
    na = len(args)    # number of arguments
    args0 = {}
    for ina in range(na):
        key1 = '%s' % ina
        args0[key1] = {}
        try:
            tmp = len(args[ina])
            if tmp == 1:
                # [1,], 'a', ['AB',],{}
                key2 = '%s' % (tmp - 1)
                args0[key1][key2] = args[ina][tmp - 1]
            else:
                if isinstance(args[ina], str):
                    # 'AB' case
                    key2 = '%s' % 0
                    args0[key1][key2] = args[ina]
                else:
                    # [1,2,...],['a','b',...],['AB','CD',...] case
                    for il in range(tmp):
                        key2 = '%s' % il
                        args0[key1][key2] = args[ina][il]
        except:
            # single number as input
            key2 = '%s' % 0
            args0[key1][key2] = args[ina]

    del args

    # keywords
    keys = list(kwds.keys())
    nk = len(keys)

    if nk != 0:
        # Has keywords input
        kwds0 = {}
        for ink, key1 in enumerate(keys):
            kwds0[key1] = {}
            try:
                tmp = len(kwds[key1])   # elements each keyword has
                if tmp == 1:
                    # [1,], 'a', ['AB',]
                    key2 = '%s' % (tmp - 1)
                    kwds0[key1][key2] = kwds[key1][tmp - 1]
                else:
                    if isinstance(kwds[key1], str):
                        # 'AB' case
                        key2 = '%s' % 0
                        kwds0[key1][key2] = kwds[key1]
                    else:
                        # [1,2,...],['a','b',...],['AB','CD',...] case
                        for il in range(tmp):
                            key2 = '%s' % il
                            kwds0[key1][key2] = kwds[key1][il]
            except:
                # single number as input
                key2 = '%s' % 0
                kwds0[key1][key2] = kwds[key1]
        del kwds

        # get the maximum list length
        nl = 0; nla = 0; nlk = 0
        for ina in range(na):
            key1 = '%s' % ina
            nla0 = len(args0[key1].keys())
            if nla0 >= nla:
                nla = nla0
        for ink in range(nk):
            key1 = keys[ink]
            nlk0 = len(kwds0[key1].keys())
            if nlk0 >= nlk:
                nlk = nlk0
        nl = max(nlk, nla)

        # check input args and kwds
        for ina in range(na):
            key1 = '%s' % ina
            nl0 = len(args0[key1].keys())
            if nl0 != 1 and nl0 < nl:
                print('input argument length error!')
                raise ValueError

        for ink in range(nk):
            key1 = keys[ink]
            nl0k = len(kwds0[key1].keys())  # number of elements for each arguments (for map)
            if nl0k != 1 and nl0k < nl:
                print('input kwds element length error!')
                raise ValueError

        # map function
        value = []
        for il in range(nl):
            arg0 = []; kwd0 = {}
            for ina in range(na):
                key1 = '%s' % ina
                nl0 = len(args0[key1].keys())
                if nl0 == 1:
                    element = args0[key1]['0']
                else:
                    key2 = '%s' % il
                    element = args0[key1][key2]
                arg0.append(element)
            for ink in range(nk):
                nlk = len(kwds0[keys[ink]])  # number of elements for each arguments (for map)
                key1 = keys[ink]
                if nlk == 1:
                    kwd0[key1] = kwds0[key1]['0']
                else:
                    key2 = '%s' % il
                    kwd0[key1] = kwds0[key1][key2]
            value.append(func(*arg0, **kwd0))

    else:
        # No keywords input (use the default of the original function)
        nl = 0
        for ina in range(na):
            key1 = '%s' % ina
            nl0 = len(args0[key1].keys())
            if nl0 >= nl:
                nl = nl0

        # check input args
        for ina in range(na):
            key1 = '%s' % ina
            nl0 = len(args0[key1].keys())
            if nl0 != 1 and nl0 < nl:
                print('input argument length error!')
                raise ValueError

        # map function
        value = []
        for il in range(nl):
            arg0 = []; kwd0 = {}
            for ina in range(na):
                key1 = '%s' % ina
                nl0 = len(args0[key1].keys())
                if nl0 == 1:
                    element = args0[key1]['0']
                else:
                    key2 = '%s' % il
                    element = args0[key1][key2]
                arg0.append(element)

            value.append(func(*arg0))

    # return is a list even taking single number for each input (attention)
    return value

# geometrical projection (general)
def projection(x, y, **kwds):
    """
    Projection of lon/lat to UTM or reverse direction
    input:
    x,y ( lon/lat or x/y )
    kwds: zone, origin, rot, inverse

    output:
    x,y ( x/y or lon/lat )
    """
    import pyproj

    zone = kwds['zone']
    origin = kwds['origin']
    rot = kwds['rot']
    inverse = kwds['inverse']

    if origin == None:
        return

    # geometrical origin
    x0 = origin[0]; y0 = origin[1]

    rot = rot * np.pi / 180.
    c, s = np.cos(rot), np.sin(rot)

    x = np.array(x, 'f')
    y = np.array(y, 'f')

    # you can use other projections (modify here)
    proj = pyproj.Proj(proj='utm', zone=zone, ellps='WGS84')

    if inverse:
        x0, y0 = proj(x0, y0, inverse=False)
        x, y = c * x - s * y, s * x + c * y
        x, y = x + x0, y + y0
        x, y = proj(x, y, inverse=True)
    else:
        x0, y0 = proj(x0, y0, inverse=False)
        x, y = proj(x, y, inverse=False)
        x, y = x - x0, y - y0
        x, y = x * c + y * s, -s * x + c * y

    return x, y

# ===========================================================================
# NGA database related
# ===========================================================================
def RakeBin(rakes):
    # rakes in degree, list

    # rake: [-180,180]
    # 0: strike-slip, [-180,-150], [-30,30], [150,180]
    # 1: normal, [-120,-60]
    # 2: reverse, [60,120]
    # 3: reverse-oblique, [30,60], [120, 150]
    # 4: Normal-oblique, [-150,-120], [-60, -30]
    # These rules come from NGA flatfile
    group = {}
    groupnames = {'U':['k', 'Unknown'], 'SS':['r', 'Strike-Slip'],
                  'NM':['g', 'Normal'], 'RV':['b', 'Reverse'],
                  'NO':['#808080', 'Normal-Oblique'],
                  'RO':['m', 'Reverse-Oblique']}
    for ig, groupname in enumerate(groupnames.keys()):
        group[groupname] = []

    for ir in range(len(rakes)):
        rake = rakes[ir]

        if rake > 180. or rake < -180. or rake == None:
            group['U'].append(rake)

        if -180 <= rake <= -150 or -30 <= rake <= 30 or 150 <= rake <= 180:
            group['SS'].append(rake)

        if -120 <= rake <= -60:
            group['NM'].append(rake)

        if 60 <= rake <= 120:
            group['RV'].append(rake)

        if 30 < rake < 60 or 120 < rake < 150:
            group['RO'].append(rake)

        if -150 < rake < -120 or -60 < rake < -30:
            group['NO'].append(rake)

    return group, groupnames

def Vs30Bin(Vs30s):
    # Vs30s in m/s, list

    # A: => 1500
    # B: [760, 1500)
    # C: [360, 760)
    # D: [180, 360)
    # E: < 180
    # This rules come from NGA flatfile
    group = {}
    groupnames = {'A':['k', 'Hard Rock'], 'B':['r', 'Rock'],
                  'C':['g', 'Dense Soil and Soft Rock'],
                  'D':['b', 'Stiff Soil'], 'E':['m', 'Soft Soil']}
    for ikey, key in enumerate(groupnames.keys()):
        group[key] = []

    for iv in range(len(Vs30s)):
        Vs30 = Vs30s[iv]
        if Vs30 >= 1500.:
            group['A'].append(Vs30)
        if 760. <= Vs30 < 1500.:
            group['B'].append(Vs30)
        if 360. <= Vs30 < 760.:
            group['C'].append(Vs30)
        if 180. <= Vs30 < 360.:
            group['D'].append(Vs30)
        if Vs30 < 180.:
            group['E'].append(Vs30)
    return group, groupnames

# ===========================================================================
# Functions to compute exploratory variables used in GMPEs
# ===========================================================================

# ==============
# Fault type
# ==============
def rake2ftype_BA(rake):
    if rake == None:
        ftype = 'U'
    if rake != None:
        if -30. <= rake <= 30. or 150. <= rake <= 180. or -180. <= rake <= -150.:
            ftype = 'SS' # strike-slip
        elif 30. < rake < 150.:
            ftype = 'RV' # reverse
        elif -150. <= rake <= -30.:
            ftype = 'NM' # normal
        else:
            print('Wrong rake angle!')
            raise ValueError
    return ftype

def rake2ftype_CB(rake):
    Frv = 0; Fnm = 0
    if 30 < rake < 150:
        Frv = 1
    if -150 < rake < -30:
        Fnm = 1
    return Frv, Fnm

def rake2ftype_CY(rake):
    Frv, Fnm = 0, 0
    if 30 <= rake <= 150:
        Frv = 1
    elif -120 <= rake <= -60:
        Fnm = 1
    return Frv, Fnm

def rake2ftype_AS(rake):
    Frv, Fnm = 0, 0
    if 30 <= rake <= 150:
        Frv = 1
    elif -120 <= rake <= -60:
        Fnm = 1
    return Frv, Fnm

# ===================
# Fault geometry
# ===================
def calc_dip(rake):
    """
    Empirical determination of dip angle from the faulting style
    Input:
        rake in degree (-180<=rake<=180)
    Output:
        dip angle in degree
    """
    if abs(rake) > 180:
        print('rake angle should be within -180 and 180')
        raise ValueError

    if abs(rake) <= 30 or abs(rake) >= 150:
        dip = 90
    elif -150 < rake < -30:
        dip = 50
    elif 30 < rake < 150:
        dip = 40
    return dip

def calc_Zhypo(M, rake):
    """
    Compute Zhypo from empirical relations
    When Ztor is unknown from input models
    """
    if M < 0:
        print('Magnitude should be larger than 0')
        raise ValueError
    if abs(rake) > 180:
        print('rake angle should be within -180 and 180')
        raise ValueError

    if abs(rake) < 30 or abs(rake) > 150:
        # strike-slip
        Zhypo = 5.63 + 0.68 * M
    else:
        Zhypo = 11.24 - 0.2 * M
    return Zhypo

def calc_W(M, rake):
    """
    Compute fault width when not specified by input
    """
    if M < 0:
        print('Magnitude should be larger than 0')
        raise ValueError

    # In R
    if abs(rake) > 180:
        print('rake angle should be within -180 and 180')
        raise ValueError
    if abs(rake) < 30 or abs(rake) > 150:
        W = 10 ** (-0.76 + 0.27 * M)
    elif -150 <= rake <= -30:
        W = 10 ** (-1.14 + 0.35 * M)
    elif 30 <= rake <= 150:
        W = 10 ** (-1.61 + 0.41 * M)

    # In Matlab
    #W = 10**(-1.01+0.32*M)

    return W

def calc_Ztor(W, dip, Zhypo):
    """
    Compute Ztor if not specified by input
    dip should be in degree
    """
    if dip <= 0 or dip > 90:
        print('dip angle should be with in (0,90]')
        raise ValueError
    if W <= 0:
        print('Fault width should be larger than 0')
        raise ValueError
    if Zhypo < 0:
        print('Zhypo should be larger than 0')
        raise ValueError
    Ztor = max(Zhypo - 0.6 * W * np.sin(dip * np.pi / 180), 0)
    return Ztor

# ======================================
# Fault-Site distances (Rjb, Rx, Rrup)
# ======================================
def calc_Rx(Rjb, Ztor, W, dip, azimuth, Rrup=None):
    """
    Compute distance parameter Rx from other inputs
    """
    if Rjb < 0:
        print('Joyer-Boore distance Rjb should be larger than 0')
        raise ValueError
    if Ztor < 0:
        print('Ztor should be larger than 0')
        raise ValueError
    if W <= 0:
        print('Fault width should be larger than 0')
        raise ValueError
    if dip <= 0 or dip > 90:
        print('dip angle should be (0,90]')
        raise ValueError
    if abs(azimuth) > 180.0:
        print('azimuth should be width in -180.0 and 180.0')
        raise ValueError

    d = dip * np.pi / 180.0 # degree to radius
    a = azimuth * np.pi / 180.0
    if dip != 90:
        if azimuth > 0:
            if azimuth == 90:
                if Rjb == 0:
                    if Rrup != None:
                        if Rrup < Ztor / np.cos(d):
                            Rx = np.sqrt(Rrup**2 - Ztor**2)
                        else:
                            Rx = Rrup / np.sin(d) - Ztor / np.tan(d)
                    else:
                        Rx = W * np.cos(d) / 2.
                        # empirical relation  (Rrup is easy to compute)
                        # assume that the site is located at the
                        # center of the surface projection of the
                        # rupture plane
                else:
                    Rx = Rjb + W * np.cos(d)
            else:
                if Rjb * abs(np.tan(a)) <= W * np.cos(d):
                    Rx = Rjb * abs(np.tan(a))
                else:
                    Rx = (Rjb * np.tan(a) *
                          np.cos(a - np.arcsin(W * np.cos(d) *
                                               np.cos(a) / Rjb)))
        else:
            Rx = Rjb * np.sin(a)
    else:
        Rx = Rjb * np.sin(a)

    return Rx

def calc_Rrup(Rx, Ztor, W, dip, azimuth, Rjb=None):
    """
    Compute the closest distance from site the the surface the fault
    """
    if Ztor < 0:
        print('Ztor should be larger than 0')
        raise ValueError
    if W <= 0:
        print('Fault width should be larger than 0')
        raise ValueError
    if dip <= 0 or dip > 90:
        print('dip angle should be (0,90]')
        raise ValueError
    if abs(azimuth) > 180:
        print('azimuth should be width in -180 and 180')
        raise ValueError

    d = dip * np.pi / 180 # degree to radius
    a = azimuth * np.pi / 180

    if dip == 90 and Rjb != None:
        Rrup = np.sqrt(Rjb**2 + Ztor**2)
        return Rrup

    if dip != 90:
        if Rx < Ztor * np.tan(d):
            Rrup1 = np.sqrt(Rx**2 + Ztor**2)
        elif Rx >= Ztor * np.tan(d) and Rx <= Ztor * np.tan(d) + W / np.cos(d):
            Rrup1 = Rx * np.sin(d) + Ztor * np.cos(d)
        elif Rx > Ztor * np.tan(d) + W / np.cos(d):
            Rrup1 = np.sqrt((Rx - W * np.cos(d))**2 +
                            (Ztor + W * np.sin(d))**2)
    elif dip == 90:
        Rrup1 = np.sqrt(Rx**2 + Ztor**2)

    if azimuth == 90 or azimuth == -90:
        Ry = 0
    elif azimuth == 0 or azimuth == 180 or azimuth == -180:
        if Rjb == None:
            print("Rjb cannot be None in the case azimuth == 0 or "
                  "azimuth == 180 or azimuth == -180")
            raise ValueError
        else:
            Ry = Rjb
    else:
        Ry = abs(Rx / np.tan(a))
    Rrup = np.sqrt(Rrup1**2 + Ry**2)
    return Rrup

# One option of doing Fling and BBP distance calculation
def calc_distances(SiteGeo, Dims, Mech, ProjDict, Rrup=False, Rx=False):
    """
    Compute Rjb, Rrup, Rx implicitly given fault geometry and site
    location (in lon/lat) The grid generation is needed to get the
    explicit fault geometry for further calculation For Fling study
    and broadband platform (*.src file)
    """

    # UTM projection property
    lon0, lat0 = ProjDict['origin']   # projection origin
    ProjDict['inverse'] = False  # from ll to xy
    kwds = ProjDict

    # sites and compute the azimuth
    rlon, rlat = SiteGeo
    rx, ry = projection(rlon, rlat, **kwds)

    # Fault dimension and focal mechanism
    Fl, dfl, Fw, dfw, ztor = Dims # fault length along strike, fault width down dip
    strike, dip, rake = Mech # fault mechanism

    # create fault mesh:
    ProjDict['inverse'] = True  # from xy to ll
    kwds = ProjDict
    # along strike direction: y-axis
    # along dip direction: x-axis
    fx = np.arange(0, Fw + dfw, dfw) # along dip (x direction)
    fy = np.arange(0, Fl + dfl, dfl) - Fl / 2 # along strike (y direction)
    fx, fy = fx * 1000, fy * 1000
    fxx, fyy = np.meshgrid(fx, fy)
    fzz = fxx * np.sin(dip * np.pi / 180.) # in meter
    sdep2d = fzz / 1000 # in km
    slon2d, slat2d = projection(fxx, fyy, **kwds)

    # surface projection (change in dip direction)
    fxS = fx * np.cos(dip * np.pi / 180.)
    fxxS, fyy = np.meshgrid(fxS, fy)
    slon2dS, slat2dS = projection(fxxS, fyy, **kwds)

    Nlat, Nlon = slon2d.shape
    Nloc = Nlat * Nlon
    fxxS1d = fxxS.reshape(Nloc) / 1000.
    fxx1d = fxx.reshape(Nloc) /1000.
    fyy1d = fyy.reshape(Nloc) /1000.
    fzz1d = fzz.reshape(Nloc) /1000.

    # compute Rjb using fault and site locations and get the azimuth
    # of all sites for later use to compute Rx
    Nsta = len(rlon)
    Rjb = []; Rrup = [];  azimuth = []
    print('Computing Rjb, and Rrup, and azimuth...')
    for ista in range(Nsta):
        rx0 = rx[ista] / 1000.
        ry0 = ry[ista] / 1000.

        # Rjb
        if ((fxxS1d.min() <= rx0 <= fxxS1d.max()) and
            (fyy1d.min() <= ry0 <= fyy1d.max())):
            Rjb.append(0)
        else:
            distS = []
            for iloc in range(Nloc):
                dx = fxxS1d[iloc] - rx0
                dy = fyy1d[iloc] - ry0
                distS.append(np.sqrt(dx**2 + dy**2))
            Rjb.append(min(distS))

        # Rrup
        dist = []
        for iloc in range(Nloc):
            dx = fxx1d[iloc] - rx0
            dy = fyy1d[iloc] - ry0
            dist.append(np.sqrt(dx**2 + dy**2 + fzz1d[iloc]**2))

        Rrup.append(min(dist))

        # compute azimuth (different from common used)
        # refer to James's definition:
        # The angle between positive fault strike direction and line
        # connecting a site to the closest point on the surface
        # projection of the TOP EDGE of rupture, which clockwise
        # assumed positive !
        # fy: the surface projection of the top edge of rupture
        # different cases
        fymin = fy.min() / 1000.; fymax = fy.max() / 1000.
        if fymin <= ry0 <= fymax:
            azimuth0 = -np.pi / 2. * (rx0 < 0.0) + np.pi / 2 * (rx0 >= 0.0)

        if ry0 > fymax:
            dx = rx0 - 0.0
            dy = ry0 - fymax
            if rx0 > 0.0:
                azimuth0 = np.arctan(dx / dy)
            elif rx0 < 0.0:
                azimuth0 = -np.arctan(-dx / dy)
            elif rx0 == 0.0:
                azimuth0 = 0.0
        if ry0 < fymin:
            dx = rx0 - 0.0
            dy = fymin - ry0
            if rx0 > 0.0:
                azimuth0 = np.pi - np.arctan(dx / dy)
            elif rx0 < 0.0:
                azimuth0 = np.arctan(-dx / dy) - np.pi
            elif rx0 == 0.0:
                azimuth0 = np.pi

        azimuth.append(azimuth0*180./np.pi)

    # Compute Rx from Rjb and Rrup
    Rx = mapfunc(calc_Rx, Rjb, ztor, Fw, dip, azimuth, Rrup=Rrup)

    OutputDict = {}
    OutputDict['Rjb'] = Rjb
    OutputDict['Rx'] = Rx
    OutputDict['Rrup'] = Rrup
    OutputDict['azimuth'] = azimuth

    return OutputDict

# ============================================================================
# Utilities for general distance calculations (in spherical
# coordinates and earth-flatten)
# ============================================================================
R = 6371. # Earth radius (km)
DegToKm = np.pi/180. * R # 1 degree to km
tol = 1e-10

def LonLatToAngleDistance(loc1, loc2, CalcRadius=True,
                          CalcDist=True, Fast=True,
                          CalcAzimuth=True, Azimuth0to2PI=False):
    """
    Convert lon/lat to radius (earth center) and/or azimuth and great
    circle Distances between points on the spherical surface

    Inputs:
	Loc1: lon1, lat1 in degree, dep1 in km
              (Starting point in azimuth calculation)
	Loc2: lon2, lat2 in degree, dep2 in km
    Outputs:
	Radius: angle between two points (central angle)
	Azimuth1to2: azimuth (relative to north pole) from point 1 to point2,
                     default: [-pi,pi]
	horzDistance: distance between two points (great circle along
                      spherical surface, km)
	vertDistance: distance between two points vertically

    Angles in radius

    # converted from OpenSHA jave to python (credit to OpenSHA
    developers) org.opensha.commons.geo.Location* and in the website:
    http://en.wikipedia.org/wiki/Haversine_formula
    http://www.movable-type.co.uk/scripts/latlong.html
    """
    loc1 = np.array(loc1)
    loc2 = np.array(loc2)

    lon1, lat1 = loc1[:2] * np.pi / 180.
    lon2, lat2 = loc2[:2] * np.pi / 180.

    # initialization
    Radius = 0.
    horzDistance = 0.
    verzDistance = 0.
    Azimuth1to2 = 0.

    if CalcRadius:
        sinDlatBy2 = np.sin((lat2 - lat1) / 2.0)
        sinDlonBy2 = np.sin((lon2 - lon1) / 2.0)
        c = (sinDlatBy2**2) + np.cos(lat1) * np.cos(lat2) * sinDlonBy2**2
        # in rad (to keep angle in between pi and -pi)
        Radius = 2.0 * np.arctan2(np.sqrt(c), np.sqrt(1 - c))
        # central angle from point1 to point2

    if CalcDist:
        if Fast:
            dlat = lat1 - lat2
            dlon = (lon1 - lon2) * np.cos((lat1 + lat2) / 2.0)
            horzDistance = R * np.sqrt(dlat**2 + dlon**2)
        else:
            sinDlatBy2 = np.sin((lat2 - lat1) / 2.0)
            sinDlonBy2 = np.sin((lon2 - lon1) / 2.0)
            c = (sinDlatBy2**2) + np.cos(lat1)*np.cos(lat2) * sinDlonBy2**2
            Radius = 2.0 * np.arctan2(np.sqrt(c), np.sqrt(1 - c))
            horzDistance = R * Radius

        dep1, dep2 = loc1[2], loc2[2]
        verzDistance = dep1 - dep2

    if CalcAzimuth:
        # calculate azimuth (p1 to p2 vector relative to p1 to north pole)
        dlon = lon2 - lon1
        cosLat2 = np.cos(lat2)
        y1 = np.sin(dlon) * cosLat2
        y2 = (np.cos(lat1) * np.sin(lat2) -
              np.sin(lat1) * cosLat2 * np.cos(dlon))
        Azimuth1to2 = np.arctan2(y1, y2)
        if Azimuth0to2PI:
            Azimuth1to2 = (Azimuth1to2+2*np.pi)%(2*np.pi)

    return Radius, horzDistance, verzDistance, Azimuth1to2

def EndLocation(loc1, vector):
    """
    Given Vector and its starting point, find the end point of the
    vector, where vector has information (azimuth,horizontal distance
    and vertical distance)
    # converted from OpenSHA jave to python (credit to OpenSHA developers)
    org.opensha.commons.geo.Location*
    Line extension (limited, not to infinite)
    """
    loc1 = np.array(loc1)
    loc1[:2] = loc1[:2] * np.pi / 180.

    lon1, lat1, dep1 = loc1
    az, DH, DV = vector    # az in radius

    # refer to http://williams.best.vwh.net/avform.htm#LL
    sinLat1 = np.sin(lat1)
    cosLat1 = np.cos(lat1)
    ad = DH / R
    sinD = np.sin(ad)
    cosD = np.cos(ad)

    # compute the location information for the end point loc2
    lat2 = np.arcsin(sinLat1 * cosD + cosLat1 * sinD * np.cos(az))
    dlon = np.arctan2(np.sin(az) * sinD * cosLat1,
                      cosD - sinLat1 * np.sin(lat2))
    lon2 = lon1 + dlon
    dep2 = dep1 + DV

    lat2 = lat2 * 180. / np.pi
    lon2 = lon2 * 180. / np.pi

    return [lon2, lat2, dep2]

def CheckPointInPolygon(point, verts):
    """
    check whether a point is inside a polygon (general coordiate and
    convex or non-convex)

    refer to:
        http://www.ecse.rpi.edu/Homepages/wrf/Research/Short_Notes/pnpoly.html
    inputs:
        point: test point (x,y,[z]), z could be None
        verts: points that define the polygon shape
               (the last one and the first one are not the same!)
    return: True or False
    """
    Ndim = len(point)
    verts = np.array(verts)
    dim = verts.shape[1]
    if Ndim != dim:
        print('point and shape should be defined with two coordinates')
        raise ValueError

    # test point
    testx = point[0]
    testy = point[1]

    vertx = verts[:, 0]
    verty = verts[:, 1]
    nvert = len(vertx)
    check = False
    j = nvert - 1
    for i in range(nvert):
        c1 = verty[i] > testy
        c2 = verty[j] > testy
        factor = ((vertx[j] - vertx[i]) *
                  (testy - verty[i]) /
                  (verty[j] - verty[i]) + vertx[i])
        if c1 != c2 and testx < factor:
            check = not check
        j = i    # edge is defined from j to i

    return check

# point to 2D line
def ptToLine2D(x1, y1, x2, y2, px, py):
    """
    Compute the point (px,py) to line (x1,y1)=>(x2,y2) allowing
    infinitely-extending of the line
    Not used
    """
    # get the projected point
    p1 = np.array([x1, y1])
    p2 = np.array([x2, y2])
    p = np.array([px, py])
    v = p2 - p1
    n = v / np.sqrt(sum(v * v)) # unit vector to show the line direction
    w10 = p1 - p
    w10n = sum(w10 * n) * n
    Ploc = w10 - w10n + p

    #dist1 = np.sqrt(sum((p-Ploc)*(p-Ploc))) # should be the same as the below
    #print dist1

    # compute distance
    # adjust vectors relative to point (x1,y1)
    x2 -= x1 * 1.0
    y2 -= y1 * 1.0
    px -= x1 * 1.0
    py -= y1 * 1.0

    # 1. projection using dot production of adjusted vector (px,py) and (x2,y2)
    dotprod = (px * x2 + py * y2)
    # length of the vector (x1,y1)=>projected point of (px,py) on the line
    projLenSq = dotprod * dotprod / (x2 * x2 + y2 * y2)

    # 2. subtraction to get the closet distance (length of the vector)
    lenSq = px * px + py * py - projLenSq
    if lenSq < 0:
        # (px,py) is in the line specified by (x1,y1) and (x2,y2)
        lenSq = 0

    dist = np.sqrt(lenSq)
    #return ((-x2*py-(-px)*y2) * 1.0 /np.sqrt(x2**2+y2**2))
    return dist, Ploc

# great circle distance (point to 2D line)
def distToLine2D(loc1, loc2, loc3, Fast=False):
    """
    Compute the shortest distance between a point (loc3) and a line
    (great circle) that extends infinitely in both directiions. Depth
    is ignored.  refer to: http://williams.best.vwh.net/avform.htm#XTE
    # this distance could be postive or negative
    """
    loc1 = np.array(loc1)
    loc2 = np.array(loc2)
    loc3 = np.array(loc3)
    lon1, lat1 = loc1[:2] * np.pi / 180.
    lon2, lat2 = loc2[:2] * np.pi / 180.
    lon3, lat3 = loc3[:2] * np.pi / 180.

    if Fast:
        # with earth-flatten approximation (faster)
        # used in shorter distance (<=200km)
        lonScale = np.cos(0.5 * lat3 +
                          0.25 * lat1 +
                          0.25 * lat2)   # earth-flatten approximation factor
        x2 = (lon2 - lon1) * lonScale
        y2 = lat2 - lat1
        x3 = (lon3 - lon1) * lonScale
        y3 = lat3 - lat1
        # x1=y1 = 0
        Term1 = (x2 * (-y3) - (-x3) * y2) / np.sqrt(x2**2 + y2**2)
        # originally, Term1 =
        # abs( (x3-x1)*(y2-y1) - (y3-y1)*(x2-x1) ) / np.sqrt((x2-x1)**2+(y2-y1)**2.)
        # for x1=y1=0, Term1 = abs( x3*y2 - y3*x2 ) / [] = abs( -y3*x2 - (-x3)*y2 ) / []
        # but here, Term1 has sign which indicates which side of
        # point3 is located relative to the line vector (1to2)
        # +: right; -:left
        return Term1 * R

    else:
        # orignial method to compute the distance from a point to a line (spherical trigonometry)
        # sin(A)/sin(a) = sin(B)/sin(b) = sin(C)/sin(c)
        # A, B, and C: angle between surface circle
        # a, b, and c: center angle of surface circle (great)
        a13, hD, vD, az13 = LonLatToAngleDistance(loc1, loc3,
                                                  CalcRadius=True,
                                                  CalcDist=False,
                                                  CalcAzimuth=True)
        a12, hD, vD, az12 = LonLatToAngleDistance(loc1, loc2,
                                                  CalcRadius=True,
                                                  CalcDist=False,
                                                  CalcAzimuth=True)
        Daz13az12 = az13 - az12
        xtd = np.arcsin(np.sin(a13) * np.sin(Daz13az12))
        if abs(xtd) < tol:
            return 0.0   # point3 is on the line
        else:
            return xtd * R
        # xtd could >0 or <0 to identify the
        # location of the point3 relative to the line you could use
        # this to compute Rx without extend your fault trace to
        # infinite, but it takes time

# deal with line segments
def minDistToLine2D(loc, segs, Fast=False, Debug=False):
    """
    Compute minimum distance between loc and a line made of segments
    Segments in line are contrained by two points loc1, loc2
    """
    Npoints = len(segs)
    minDist0 = 1000.
    for iseg in range(1, Npoints):
        p1 = segs[iseg - 1]
        p2 = segs[iseg]
        dist = distToLine2D(p1, p2, loc, Fast=Fast)
        dist0 = abs(dist)
        if dist0 <= minDist0:
            minDist0 = dist0
            minDist = dist
    return minDist

# point to line segments
# Used by Rjb and Rx
def ptToLineSeg2D(x1, y1, x2, y2, px, py):
    """
    Compute the point (px,py) to line (x1,y1)=>(x2,y2) without infinitely-extending of the line
    Distance measured is the distance between the specified point and the closest point between
    the specified end points (x1,y1) and (x2,y2)
    """

    # adjust vectors relative to point (x1,y1)
    x2 -= x1 * 1.0
    y2 -= y1 * 1.0
    px -= x1 * 1.0
    py -= y1 * 1.0

    # 1. projection using dot production of adjusted vector (px,py) and (x2,y2)
    dotprod = (px * x2 + py * y2)
    if dotprod <= 0.0:
        # (px,py) is on the side of (x1,y1) (projection of the point (px,py) is not on the segment)
        #print 'beyond point1'
        projLenSq = 0.0
    else:
        # check the other side relationship
        px = x2 - px
        py = y2 - py
        dotprod = px * x2 + py * y2
        if dotprod <= 0.0:
            #   print 'beyond point2'
            # (px,py) is on the side of (x2,y2) (projection of the
            # point (px,py) is not on the segment)
            projLenSq = 0.0
        else:
            # point (px,py) 's projection is in between (x1,y1) and (x2,y2)
            # same as the ptLineDist function
            projLenSq = dotprod * dotprod / (x2 * x2 + y2 * y2)

    # 2. subtraction to get the closet distance (length of the vector)
    # if projLenSq = 0.0, then the distance would be either original
    # (px,py) to (x1,y1) or (x2,y2)
    lenSq = px * px + py * py - projLenSq
    if lenSq < 0:
        # (px,py) is in the line specified by (x1,y1) and (x2,y2)
        lenSq = 0.0
    return np.sqrt(lenSq)

# great circle
def distToLineSeg2D(loc1, loc2, loc3, Fast=False):
    """
    Compute distance between point3 and line defined by point1 and point2
    loc1, loc2, loc3 are list with 3 elements
    There are three cases: the projection point of loc3 on:
    the loc1 side, on the loc2 side, in between loc1 and loc2
    2D
    """
    loc1 = np.array(loc1)
    loc2 = np.array(loc2)
    loc3 = np.array(loc3)
    lon1, lat1 = loc1[:2] * np.pi / 180.
    lon2, lat2 = loc2[:2] * np.pi / 180.
    lon3, lat3 = loc3[:2] * np.pi / 180.

    if Fast:
        lonScale = np.cos(0.5 * lat3 + 0.25 * lat1 + 0.25 * lat2)
        x2 = (lon2 - lon1) * lonScale
        y2 = lat2 - lat1
        x3 = (lon3 - lon1) * lonScale
        y3 = lat3 - lat1
        Term1 = ptToLineSeg2D(0, 0, x2, y2, x3, y3) # always positive
        return Term1 * R
    else:
        # use cos(c) = cos(a)cos(b) + sin(a)sin(b)cos(C) to get
        a13, hD13, vD13, az13 = LonLatToAngleDistance(loc1, loc3, Fast=Fast)
        a12, hD12, vD12, az12 = LonLatToAngleDistance(loc1, loc2, Fast=Fast)
        Daz13az12 = az13 - az12

        # cross-track distance (in radius)
        xtd = np.arcsin(np.sin(a13) * np.sin(Daz13az12))

        # along-track distance (in km)
        atd = np.arccos(np.cos(a13) / np.cos(xtd)) * R
        a23, hD23, vD23, az23 = LonLatToAngleDistance(loc2, loc3,
                                                      CalcRadius=False,
                                                      CalcDist=True,
                                                      Fast=Fast,
                                                      CalcAzimuth=False)

        # check if beyond p3 (should be p2?) (different from the original Rx definition?)
        if atd > hD12:
            #print 'Beyond p2'
            return hD23

        # check if beyond p1
        if np.cos(Daz13az12) < 0:
            #print 'Beyond p1'
            return hD13

        # projection of the point is within the two points
        if abs(xtd) < tol:
            return 0.0   # point3 is on the line
        else:
            return abs(xtd) * R

# segments
def minDistToLineSeg2D(loc, segs, Fast=False, Debug=False):
    """
    Compute minimum distance between loc and a line made of segments
    Segments in line are contrained by two points loc1, loc2
    """
    Npoints = len(segs)
    minDist = 1000.
    if Debug:
        print('minDistToLineSeg2D Debug')
    for iseg in range(1, Npoints):
        p1 = segs[iseg - 1]
        p2 = segs[iseg]
        dist = abs(distToLineSeg2D(p1, p2, loc, Fast=Fast))
        if Debug:
            print('LineSeg %s' % iseg, dist)

        if dist <= minDist:
            minDist = dist
    return minDist

# Used by Rrup based on just corner points
def ptToLineSeg3D(point, point1, point2, Rscale=1.0):
    """
    get your coordinate as float number
    default is in Cartesian
    if you set Rscale=6371, it will give you the point and compute further use LonLatToAngleDistance
    http://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line
    """
    points = [point, point1, point2]
    points = np.array(points, 'f')
    p0 = points[0, :]
    p1 = points[1, :]
    p2 = points[2, :]
    v = p2 - p1
    w0 = p0 - p1
    w1 = p0 - p2
    dotprod = sum(w0 * v)
    if dotprod <= 0:
        #print 'beyond point1'
        Ploc = p1
        if Rscale != 1.0:
            a, hD, vD, az = LonLatToAngleDistance(Ploc, p0,
                                                  CalcRadius=False,
                                                  CalcDist=True,
                                                  CalcAzimuth=False,
                                                  Fast=True)
            dist = np.sqrt(hD**2 + vD**2)
        else:
            w0[:2] *= Rscale
            dist = np.sqrt(sum(w0**2))
    else:
        dotprod = sum(w1 * v)
        if dotprod >= 0:
            #   print 'beyond point2'
            Ploc = p2
            if Rscale != 1.0:
                a, hD, vD, az = LonLatToAngleDistance(Ploc, p0,
                                                      CalcRadius=False,
                                                      CalcDist=True,
                                                      CalcAzimuth=False,
                                                      Fast=True)
                dist = np.sqrt(hD**2 + vD**2)
            else:
                w1[:2] *= Rscale
                dist = np.sqrt(sum(w1**2))
        else:
            # get the projected point on the 3D line
            n = v / np.sqrt(sum(v * v)) # unit vector to show the line direction
            w10 = p1 - p0
            w10n = sum(w10 * n) * n
            Ploc = w10 - w10n + p0

            if Rscale != 1.0:
                a, hD, vD, az = LonLatToAngleDistance(Ploc, p0,
                                                      CalcRadius=False,
                                                      CalcDist=True,
                                                      CalcAzimuth=False,
                                                      Fast=True)
                dist = np.sqrt(hD**2 + vD**2)
            else:
                vn = np.cross(v, w0)
                u = p2 - p1
                u[:2] *= Rscale
                vn[:2] *= Rscale
                dist = np.sqrt(sum(vn * vn) / sum(u * u))
    return dist, Ploc

def minDistToLineSeg3D(loc, segs, Rscale=1.0, Debug=False):
    """
    Compute minimum distance between loc and a line made of segments
    Segments in line are contrained by two points loc1, loc2
    """
    Npoints = len(segs)
    minDist = 1000.
    if Debug:
        print('minDistToLineSeg3D Debug')
    for iseg in range(1, Npoints):
        p1 = segs[iseg - 1]
        p2 = segs[iseg]
        dist_tmp, Ploc = ptToLineSeg3D(loc, p1, p2, Rscale=Rscale)
        a, hD, vD, az = LonLatToAngleDistance(loc, Ploc, CalcRadius=False,
                                              CalcDist=True,
                                              CalcAzimuth=False,
                                              Fast=True)
        dist = np.sqrt(hD**2 + vD**2)
        if Debug:
            print('Line seg %s of %s:' % (iseg, Npoints - 1), Ploc, hD, vD)
        if dist <= minDist:
            minDist = dist

    return minDist

def ptToSurf3D(point0, point1, point2, point3):
    """
    point and plane relationship
    distance and the coordinate on a plane
    if points are on spherical earth, then you need do the correction before using this function
    reference:
    http://jtaylor1142001.net/calcjat/Solutions/VPlanes/VP3Pts.htm
    http://www.9math.com/book/projection-point-plane

    Input:
        point0: a point outside or within the plane
        points (3 row, 3 col LIST contains three know points on the plane)
    Output:
        D: distance between point0 and the plane
        point1: projection point of point0 on the plane
    """

    points = [point1, point2, point3]
    points = np.array(points)
    vn = np.cross(points[1, :] - points[0, :],
                  points[2, :] - points[0, :])  # plane normal
    vn = vn / np.sqrt(np.dot(vn, vn))   # get the normal vector of the plane

    a = vn[0]; b = vn[1]; c = vn[2]
    d = -a * points[0, 0] - b * points[0, 1] - c * points[0, 2]

    u, v, w = point0

    L1 = a * u + b * v + c * w + d
    L2 = a**2 + b**2 + c**2
    x = u - a * L1 / L2
    y = v - b * L1 / L2
    z = w - c * L1 / L2

    point1 = x, y, z
    dist = abs(L1) / np.sqrt(L2)

    return dist, point1

def minDistToSurfSeg(loc, segs, Rscale=1.0, Debug=False):
    """
    Compute minimum distance between loc and a line made of segments
    Segments in line are contrained by two points loc1, loc2
    """
    Nseg = len(segs)
    minDist = 1000.
    if Debug:
        print('minDistToSurfSeg Debug')
        Ppoints = []

    for iseg in range(Nseg):
        points = segs[iseg]  # there are fout points
        point1, point2, point3, point4 = points
        tmp_dist, Ppoint = ptToSurf3D(loc, point1, point2, point3)
        if Debug:
            Ppoints.append(Ppoint)

        a, hD, vD, az = LonLatToAngleDistance(loc, Ppoint, CalcRadius=False,
                                              CalcDist=True,
                                              CalcAzimuth=False, Fast=True)
        check = CheckPointInPolygon(Ppoint, points)
        if check:
            dist = np.sqrt(hD**2 + vD**2)
        else:
            pointsClosed = points + [points[0]]
            dist = minDistToLineSeg3D(loc, pointsClosed,
                                      Rscale=Rscale, Debug=Debug)
        if dist <= minDist:
            minDist = dist

        if Debug:
            print('Segment %s: ' % iseg, dist, check, Ppoint)

    if Debug:
        return minDist, Ppoints
    else:
        return minDist

def FaultTraceGen(Origin, Dims, Mech):
    """
    Generate top FaultTrace and required parameter for function
    SimpleFaultSurface
    """
    # Fault dimension and focal mechanism
    lon0, lat0 = Origin # top center lon/lat along the strike with depth ztor
    Fl, dfl, Fw, dfw, ztor = Dims # fault length along strike, fault width down dip
    strike, dip, rake = Mech   # fault mechanism

    # extend to get fault surface
    loc = lon0, lat0, ztor
    hD = Fl / 2.
    vD = 0.0
    strike *= np.pi / 180.   # strike is within 0,2pi
    dipRad = dip * np.pi / 180.   # strike is within 0,2pi

    # along strike extension
    vector = [strike, hD, vD]
    loc2 = EndLocation(loc, vector)
    vector = [strike + np.pi, hD * 2, vD]
    loc1 = EndLocation(loc2, vector)
    FaultTrace1 = [loc1, loc2]

    AveDip = dip  # in degree
    UpperSeisDepth = ztor
    LowerSeisDepth = Fw * np.sin(dipRad) + ztor
    GridSpaceAlongStrike = dfl
    GridSpaceDownDip = dfw

    return (FaultTrace1, UpperSeisDepth,
            LowerSeisDepth, AveDip,
            GridSpaceAlongStrike, GridSpaceDownDip)

# Surface extension (given FaultModel Dict)
def SimpleFaultSurface(FaultTrace, UpperSeisDepth, LowerSeisDepth,
                       AveDip, GridSpaceAlongStrike=None,
                       GridSpaceDownDip=None):
    """
    Extend fault surface based on
        a.UCERF-type fault model (fault trace, upper seis depth, lower seis depth, and average dip)
        b.BBP-type (including Fling)
    Using Stirling method (downdip extension perpendicular to the average strike when dealing with
    this could deal with multiple segments
    FaultTrace has to be a list which show the surface trace the rupture top
    """

    AveDip = AveDip * np.pi / 180. # To Radius
    vD = LowerSeisDepth - UpperSeisDepth
    hD = vD / np.tan(AveDip)

    FaultGeom = []

    Npoints = len(FaultTrace)
    Nsegs = Npoints - 1

    if Nsegs == 1:
        iseg = 0
        loc1 = FaultTrace[iseg]
        loc2 = FaultTrace[iseg + 1]
        a, hD0, vD0, az = LonLatToAngleDistance(loc1, loc2, CalcRadius=False,
                                                CalcAzimuth=True,
                                                CalcDist=True)
        AveStrike = az

        # extension if given grid space
        if GridSpaceAlongStrike and GridSpaceDownDip:
            daa = GridSpaceAlongStrike
            ddd = GridSpaceDownDip
        elif GridSpaceAlongStrike and not GridSpaceDownDip:
            daa = ddd = GridSpaceAlongStrike
        elif not GridSpaceAlongStrike and GridSpaceDownDip:
            daa = ddd = GridSpaceDownDip
        else:
            # no grid generation (just get four corner points)
            FaultTrace1 = []    # downdip extension points
            vector = [az + np.pi / 2, hD, vD]    # downdip extension

            loc3 = EndLocation(loc2, vector)
            FaultTrace1.append(loc3)

            loc4 = EndLocation(loc1, vector)
            FaultTrace1.append(loc4)

            FaultTraceSeg = []
            FaultTraceSeg.append([loc1, loc2, loc3, loc4])

            FaultTrace = FaultTrace + FaultTrace1
            return FaultTrace, FaultTraceSeg, AveStrike

        # Extended fault surface grid
        Ncol = int(hD0 / daa + 1)
        Nrow = int(vD / np.sin(AveDip) / ddd + 1)

        for irow in range(Nrow):
            vector0 = [az + np.pi / 2,
                       ddd * np.cos(AveDip),
                       ddd * np.sin(AveDip)]
            if irow == 0:
                FaultDD = FaultTrace[0]
            else:
                FaultDD = EndLocation(FaultDD, vector0)
            FaultAA = [FaultDD,]
            for icol in range(1, Ncol):
                vector1 = [az, daa, 0.0]
                FaultAA.append(EndLocation(FaultAA[icol - 1], vector1))
            FaultGeom.append(FaultAA)

        return FaultGeom   # should have shape like: (Nrow,Ncol,3)B

    else:
        # multiple segments
        azs = []; hDs = []; vDs = []      # or strikes
        AveStrike = 0
        for iseg in range(Nsegs):
            loc1 = np.array(FaultTrace[iseg])
            loc2 = np.array(FaultTrace[iseg+1])
            a, hD0, vD0, az = LonLatToAngleDistance(loc1, loc2,
                                                    CalcRadius=False,
                                                    CalcAzimuth=True,
                                                    CalcDist=True)

            az = (az + 2 * np.pi) % (2 * np.pi)
            hDs.append(hD0)
            AveStrike += az
            azs.append(az)

        AveStrike /= Nsegs    # average strike used by Stirling fault extension

        # extension if given grid space
        if GridSpaceAlongStrike and GridSpaceDownDip:
            daa = GridSpaceAlongStrike
            ddd = GridSpaceDownDip
        elif GridSpaceAlongStrike and not GridSpaceDownDip:
            daa = ddd = GridSpaceAlongStrike
        elif not GridSpaceAlongStrike and GridSpaceDownDip:
            daa = ddd = GridSpaceDownDip
        else:
            # no grid generation (just get four corner points)
            FaultTraceSeg = []
            for iseg in range(Nsegs):
                loc1 = FaultTrace[iseg]
                loc2 = FaultTrace[iseg + 1]
                vector = [AveStrike + np.pi / 2., hD, vD]
                loc3 = EndLocation(loc2, vector)
                loc4 = EndLocation(loc1, vector)
                FaultTraceSeg.append([loc1, loc2, loc3, loc4])
            FaultTrace1 = []
            for ipoint in range(Npoints):
                loc0 = FaultTrace[ipoint]
                vector = [AveStrike + np.pi / 2., hD, vD]
                loc00 = EndLocation(loc0, vector)
                FaultTrace1.append(loc00)
            FaultTrace1.reverse()
            FaultTrace = FaultTrace + FaultTrace1
            return FaultTrace, FaultTraceSeg, AveStrike

        # loop over segments
        Nrow = int(vD / np.sin(AveDip) / ddd + 1)
        Ncol = []
        for iseg in range(Nsegs):
            Ncol.append(int(hDs[iseg] / daa + 1))
        for irow in range(Nrow):
            vector0 = [AveStrike + np.pi / 2,
                       ddd * np.cos(AveDip),
                       ddd * np.sin(AveDip)]
            if irow == 0:
                FaultDD = FaultTrace[0]
            else:
                FaultDD = EndLocation(FaultDD, vector0)
            FaultAA = [FaultDD,]
            icount = 0
            for iseg in range(Nsegs):
                for icol in range(1, Ncol[iseg]):
                    vector1 = [azs[iseg], daa, 0.0]
                    loc0 = EndLocation(FaultAA[icount], vector1)
                    FaultAA.append(loc0)
                    icount += 1
            FaultGeom.append(FaultAA)

        return FaultGeom   # should have shape like: (Nrow,Ncol,3)

def srfFaultSurfaceExtract(SRFfile):
    """
    Generate fault surface from SRF file convention
    Following the Graves' SRF convention used in BBP and CyberShake
    """

    lines = open(SRFfile, 'r').readlines()
    Nseg = int(lines[1].strip().split()[1])

    # loop over segments to get (Nrow,Ncol) of each segments
    # fault surface for each segment will be read latter
    srfFaultSurface = {}
    srfFaultSurface['segments'] = {}

    dims = []
    dips = []
    ztors = []
    for iseg in range(Nseg):
        il0 = 2 * iseg + 2  # fault geometry info
        spl = lines[il0].strip().split()
        lon0, lat0, L, W, Ncol, Nrow = np.array(spl, 'f')
        Ncol, Nrow = int(Ncol), int(Nrow)
        dims.append([Ncol, Nrow])

        il1 = il0 + 1     # focal mechanism and hypocenter info
        spl = lines[il1].strip().split()
        strike, dip, ztor, hypoAS, hypoDD = np.array(spl, 'f')
        dips.append(dip)    # will be used to get the average dip angle (over segments)
        ztors.append(ztor)

    srfFaultSurface['segments']['dims'] = dims
    srfFaultSurface['segments']['dips'] = dips
    srfFaultSurface['segments']['ztors'] = ztors

    il0 = 2 * (Nseg + 1)
    Npoints = int(lines[il0].strip().split()[1])

    # jump to the data block (for each segments, there are a data block)
    il0 = il0 + 1

    locs = []; rakes = []
    while il0 < len(lines):
        spl = lines[il0].strip().split()
        lon, lat, dep, strike, dip, Area, Tinit, dt = np.array(spl, 'f')
        locs.append([lon, lat, dep])

        il0 = il0 + 1
        spl = lines[il0].strip().split()
        rake, slipA_AlongRake, Nt = np.array(spl[:3], 'f')
        rakes.append(rake) # will be used to get average rake (over points)
        dl = int(Nt / 6) + (Nt % 6 != 0) * 1
        il0 = il0 + dl + 1   # import (similar to the segments jump) ...

    Nrow1 = 0; Ncol1 = 0
    for iseg in range(Nseg):
        Nrow1 += dims[iseg][1]
        Ncol1 += dims[iseg][0]

    FaultGeom = np.array(locs).reshape((Nrow1, Ncol1, 3))
    srfFaultSurface['FaultGeom'] = FaultGeom
    srfFaultSurface['rakes'] = rakes

    return srfFaultSurface

def srfFaultSurfaceTest(SRFfile):
    """
    test function srfFaultSurfaceExtract
    """
    srfFaultSurface = srfFaultSurfaceExtract(SRFfile)
    FaultGeom = srfFaultSurface['FaultGeom']

    # fault surface geometry in 3D
    slon3d = FaultGeom[:, :, 0]
    slat3d = FaultGeom[:, :, 1]
    sdep3d = FaultGeom[:, :, 2]

    # fault surface projection in 2D (depth = 0 )
    sdep2d = np.zeros(len(slon3d)).tolist()

    # plot (depends on matplotlib)
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    fig = plt.figure(1)
    ax = Axes3D(fig)
    ax.plot(slon3d, slat3d, sdep2d, 'ro')
    linec = ax.plot_wireframe(slon3d, slat3d, -sdep3d)
    linec.set_color('b')

    ax.set_zlim3d(-50, 0)
    ax.set_xlabel('lon')
    ax.set_ylabel('lat')
    ax.set_zlabel('depth (km)')
    fig.savefig('./tmp/srfFaultSurfaceTest.png', format='png')

def DistanceX(SiteGeom, FaultTrace1, AveStrike=None, Fast=True,
              DealRx='Extension', Debug=False):
    """
    Compute Rx by extending the fault trace to infinity in both ends
    """
    Npoints = len(FaultTrace1)

    ps = FaultTrace1[0]
    pe = FaultTrace1[-1]

    if AveStrike == None:
        Radius, hD, vD, Azimuth1to2 = LonLatToAngleDistance(ps, pe,
                                                            CalcRadius=False,
                                                            CalcDist=True,
                                                            Fast=Fast,
                                                            CalcAzimuth=True,
                                                            Azimuth0to2PI=True)
        direction = Azimuth1to2
    else:
        direction = AveStrike

    vector = [direction + np.pi / 2., 1000, 0.0]
    Loc1 = ps; Loc2 = pe
    Loc3 = EndLocation(Loc1, vector)
    Loc4 = EndLocation(Loc2, vector)
    verts = []; segs = []
    for ipoint in range(Npoints):
        verts.append(FaultTrace1[ipoint])
        segs.append(FaultTrace1[ipoint])
    verts.append(Loc4)
    verts.append(Loc3)

    check = CheckPointInPolygon(SiteGeom[:2], np.array(verts)[:, :2])
    if check:
        # one special case: Site is within the surface projection of
        # the source (hanging wall)
        distToFaultTrace = minDistToLine2D(SiteGeom, segs,
                                           Fast=Fast, Debug=Debug)
        Rx = distToFaultTrace
    else:
        # when site is outside the surface projection of the fault:
        # 1. extend the fault to both end (very large area)
        # 2. use azimuth to judge consider hanning all or not
        if DealRx == 'Extension':
            # extend to both end (regardless azimuth)
            vector = [direction, 1000.0, 0.0]    # horizontal extension vector
            Loc2 = EndLocation(pe, vector)

            vector[0] = direction + np.pi   # flip over trace direction
            Loc1 = EndLocation(ps, vector)

            # down-dip extension
            vector[0] = direction + np.pi / 2.
            Loc3 = EndLocation(Loc1, vector)
            Loc4 = EndLocation(Loc2, vector)

            verts = []   # to form the polygon
            segs = []    # to form the extend fault trace segments
            verts.append(Loc1)
            segs.append(Loc1)
            for ipoint in range(Npoints):
                verts.append(FaultTrace1[ipoint])
                segs.append(FaultTrace1[ipoint])
            verts.append(Loc2)
            segs.append(Loc2)
            verts.append(Loc4)
            verts.append(Loc3)
            check = CheckPointInPolygon(SiteGeom[:2], np.array(verts)[:, :2])
            if Debug:
                # Test extended fault trace and fault surface projection (plot)
                import matplotlib.pyplot as plt
                print('Site is within the Extended Fault Surface projection:',
                      check)
                verts1 = np.array(verts)
                fig = plt.figure(10)
                ax = fig.add_subplot(111)
                ax.plot(verts1[:, 0], verts1[:, 1], 'ro')
                # plot the initial points (where to start from)
                ax.plot([ps[0], pe[0]], [ps[1], pe[1]], 'bx')
                ax.plot(SiteGeom[0], SiteGeom[1], 'rs')
                ax.set_title('Fault Extension for Rx Calculation')

            distToExtendedTrace = minDistToLineSeg2D(SiteGeom, segs,
                                                     Fast=Fast, Debug=Debug)
            if check or distToExtendedTrace == 0.0:
                Rx = distToExtendedTrace
            else:
                Rx = -distToExtendedTrace

        else:

            # azimuth
            # Rx check  (using azimuth when site is outside the fault surface) 90-azimuth >= 45
	    # ...

            pass

    return Rx

def DistanceToSimpleFaultSurface(SiteGeom, FaultTrace1, UpperSeisDepth,
                                 LowerSeisDepth, AveDip,
                                 GridSpaceAlongStrike=None,
                                 GridSpaceDownDip=None,
                                 Fast=True, Debug=False,
                                 RrupCalc=True, RxCalc=True):
    """
    Compute Rjb,Rrup,Rx for simple fault plane
    FaultTrace1 is just the top fault trace
    """

    Npoints = len(FaultTrace1)
    if GridSpaceAlongStrike == None and GridSpaceDownDip == None:
        FaultTrace, FaultSeg, AveStrike = SimpleFaultSurface(FaultTrace1,
                                                             UpperSeisDepth,
                                                             LowerSeisDepth,
                                                             AveDip)

        # Rjb: points to line (or seg)  (quite good!)
        verts = FaultTrace
        vertsClosed = FaultTrace + [FaultTrace[0]]
        minRjb = minDistToLineSeg2D(SiteGeom, vertsClosed,
                                    Fast=Fast, Debug=Debug)
        check = CheckPointInPolygon(SiteGeom, verts)

        if check:
            Rjb = 0.0
        else:
            Rjb = minRjb

        if RrupCalc:
            # Rrup: point to surface (like Rjb, but to fault surface)
            if Debug:
                Rrup, Ppoints = minDistToSurfSeg(SiteGeom, FaultSeg,
                                                 Rscale=DegToKm,
                                                 Debug=Debug)
            else:
                Rrup = minDistToSurfSeg(SiteGeom, FaultSeg,
                                        Rscale=DegToKm, Debug=Debug)
        else:
            Rrup = None

        if RxCalc:
            # Rx: points to line (consider the sign, relative to the strike)
            Rx = DistanceX(SiteGeom, FaultTrace1,
                           AveStrike=None, Fast=Fast, Debug=Debug)
        else:
            Rx = None

        if Debug:
            print('Site is within surface projection: ', check)
            import matplotlib.pyplot as plt
            from mpl_toolkits.mplot3d import Axes3D
            FaultTrace = np.array(FaultTrace)
            Ppoints = np.array(Ppoints)
            fig = plt.figure(1)
            ax = Axes3D(fig)
            verts1 = np.array(vertsClosed)
            ax.plot(verts1[:, 0], verts1[:, 1], -verts1[:, 2], 'b-')
            ax.plot(FaultTrace[:, 0], FaultTrace[:, 1],
                    FaultTrace[:, 2] * 0.0, 'ko')
            ax.plot([SiteGeom[0]], [SiteGeom[1]], [0.0], 'r^')
            ax.plot(Ppoints[:, 0], Ppoints[:, 1], Ppoints[:, 2], 'rs')
            plt.show()

        return Rjb, Rrup, Rx

    else:
        FaultGeom = SimpleFaultSurface(FaultTrace1, UpperSeisDepth,
                                       LowerSeisDepth, AveDip,
                                       GridSpaceAlongStrike=GridSpaceAlongStrike,
                                       GridSpaceDownDip=GridSpaceDownDip)
        FaultGeom = np.array(FaultGeom)
        Nrow, Ncol, Nelm = FaultGeom.shape

        Rjb, Rrup, Rx = DistanceToEvenlyGriddedSurface(SiteGeom, FaultGeom,
                                                       Fast=Fast,
                                                       RrupCalc=RrupCalc,
                                                       RxCalc=RxCalc)

        if Debug:
            Fault = FaultGeom.reshape((Nrow * Ncol, 3))
            import matplotlib.pyplot as plt
            fig = plt.figure(2)
            ax = Axes3D(fig)
            ax.plot(Fault[:, 0], Fault[:, 1], -Fault[:, 2], 'b.')
            ax.plot(Fault[:, 0], Fault[:, 1], Fault[:, 2] * 0.0, 'ko')
            plt.show()

        return Rjb, Rrup, Rx

# General distance calculation (before this, you need to generate FaultGeo)
# only requirement is the fault geometry (explicitly)
def DistanceToEvenlyGriddedSurface(SiteGeo, FaultGeo, Fast=True,
                                   RrupCalc=True, RxCalc=True):
    """
    Compute Rjb, Rrup, Rx explicitly given discretized fault (3D) surface
    geometry and site location (in lon/lat)

    Rx just use the fault trace along-strike and azimuth between
    the two points at the two ends of the fault and site location.

    Inputs:
        FaultGeo: faultGeo
                list faultGeo has dim: (Nrow,Ncol,3), fault1 surface discretization: (Nrow,Ncol)
                Nrow: down-dip direction grid points; Ncol: along-strike direction grid points
        AveDip: average dip (to determine the points to use)
        SiteGeo: siteGeo
                list siteGeo has elements: rlon,rlat,rdep to generally specify the location of site
    Outputs:
        Rjb, Rrup, Rx
        They all have the following shape (site-based):

    """
    FaultGeo = np.array(FaultGeo)
    SiteGeo = np.array(SiteGeo)
    loc1 = SiteGeo

    minRjb = 1000.    # in km
    minRrup = 1000.   # in km
    minRx = 1000.     # in km

    Nrow, Ncol, Nelm = FaultGeo.shape

    surf = FaultGeo.reshape((Nrow * Ncol, Nelm))
    for loc2 in surf:
        alpha, hD, vD, az = LonLatToAngleDistance(loc1, loc2,
                                                  CalcRadius=False,
                                                  CalcDist=True,
                                                  Fast=Fast,
                                                  CalcAzimuth=False)
        if hD <= minRjb:
            minRjb = hD
        totalDist = np.sqrt(hD**2 + vD**2)
        if totalDist <= minRrup:
            minRrup = totalDist
    if RrupCalc:
        Rrup = minRrup
    else:
        Rrup = None

    # check site is within the surface projection of the fault
    verts = []
    irow = 0
    for icol in range(Ncol):
        verts.append(FaultGeo[irow, icol][:2].tolist())
    icol = 0
    for irow in range(Nrow):
        verts.append(FaultGeo[irow, icol][:2].tolist())
    irow = Nrow - 1
    for icol in range(Ncol):
        verts.append(FaultGeo[irow, icol][:2].tolist())
    icol = Ncol - 1
    for irow in range(Nrow):
        verts.append(FaultGeo[irow, icol][:2].tolist())
    check = CheckPointInPolygon(SiteGeo[:2], verts)
    if check:
        Rjb = 0.0
    else:
        Rjb = minRjb

    if 0:
        DDtmp = 0
        for irow in range(Nrow - 1):
            loc1 = FaultGeo[irow, 0]
            loc2 = FaultGeo[irow + 1, 0]
            alpha, hD, vD, az = LonLatToAngleDistance(loc1, loc2,
                                                      CalcRadius=False,
                                                      CalcDist=True,
                                                      Fast=Fast)
            DDtmp += DDtmp + hD
        DownDipGridSpace = DDtmp / (Nrow - 1)
        AStmp = 0
        for icol in range(Ncol - 1):
            loc1 = FaultGeo[0, icol]
            loc2 = FaultGeo[0, icol + 1]
            alpha, hD, vD, az = LonLatToAngleDistance(loc1, loc2,
                                                      CalcRadius=False,
                                                      CalcDist=True,
                                                      Fast=True)
            AStmp += AStmp + hD
        AlongStrikeGridSpace = AStmp / (Ncol - 1)
        AveGridSpace = (DownDipGridSpace + AlongStrikeGridSpace) / 2.0
        if minRjb < AveGridSpace:
            Rjb = 0.0
        else:
            Rjb = minRjb

    # ============
    # compute Rx by extending fault trace and fault surface projection
    # ============
    if RxCalc:
        FaultTrace = FaultGeo[0] # first row to get the fault trace for Rx calculation
        Rx = DistanceX(SiteGeo, FaultTrace, AveStrike=None, Fast=Fast)
    else:
        Rx = None

    return Rjb, Rrup, Rx

# ==========================
# Site-Specific Parameters
# ==========================
def calc_Z1(Vs30, Z1model):
    """
    Compute Z1.0 parameter for AS and CY model if not specified
    Vs30 in m/s
    return in km
    """
    if Vs30 < 0:
        print(' Vs30 should be larger than 0')
        raise ValueError

    if Z1model == 'AS':
        if Vs30 < 180.:
            Z10 = 6.745
        elif 180 <= Vs30 <= 500.:
            Z10 = 6.745 - 1.35*np.log(Vs30 / 180.)
        elif Vs30 > 500:
            Z10 = 5.394 - 4.48 * np.log(Vs30 / 500.)
    elif Z1model == 'CY':
        Z10 = 28.5 - 3.82 / 8 * np.log(Vs30**8 + 378.7**8)
    return np.exp(Z10)   # in meter

def calc_Z25(Vs30, Z1=None, Z15=None, Z1model='CY'):
    """
    input Vs30 m/s; Z1 in m; Z15 in m
    """

    if Vs30 == None and Z1 == None and Z15 == None:
        print('Either Vs30, Z1.0 or Z1.5 should be specified')
        raise ValueError
    if Vs30 != None and Vs30 < 0:
        print('Vs30 should be larger than 0')
        raise ValueError
    if Z1 != None and Z1 < 0:
        print('Z1 should be larger than 0')
        raise ValueError
    if Z15 != None and Z15 < 0:
        print('Z15 should be larger than 0')
        raise ValueError
    if Z15 != None:
        Z25 = 636 + 1.549 * Z15
    elif Z1 != None:
        Z25 = 519 + 3.595 * Z1
    elif Vs30 != None:
        Z1 = calc_Z1(Vs30, Z1model)
        Z25 = 519 + 3.595 * Z1
    return Z25 / 1000.   # in km

# ====================
# Applications (other)
# ====================
def GetIntraInterResiduals(residualT, EQID, sigmaT, tau, sigma, AS=None):
    """
    Compute Normalized Total, Inter- and Intra residuals
    Input:
        residualT: un-normalized total residual
        EQID: earthquake id for each record (used for group bin)
        sigmaT: standard deviation (total residual) for each record
        tau: standard deviation (inter-event) for each record
        sigma: standard deviation (intra-event) for each record
        AS: None (just use average method): not None (Use AS 1992 method to compute)
    Output:
        epsilonT: normalized total residual
        eta: normalized inter-event residual
        epsilon: normalized intra-event residual
    """

    residualT = np.array(residualT)
    sigmaT = np.array(sigmaT)
    tau = np.array(tau)
    sigma = np.array(sigma)

    # group bin
    Events = group_list(EQID)
    Neq = len(Events)
    ID = []
    for ieq in range(Neq):
        ID.append(Events[ieq][2])

    # Compute residuals
    eta_EQinter = []
    for ieq in range(Neq):
        index = (np.array(RID) == ID[ieq]).nonzero()[0]
        N_EQ = len(index)

        if AS != None:
            tau0 = np.mean(tau[index])
            sigma0 = np.mean(sigma[index])
            residual_EQinter0 = (tau0**2 * sum(residualT[index]) /
                                 (N_EQ * tau0**2 + sigma0**2))
        else:
            residual_EQinter0 = np.mean(residualT[index])
        for isub in index:
            eta_EQinter.append(residual_EQinter0 / tau[isub]) # for each record in the event group

    eta = np.array(eta_EQinter)
    residual_EQintra = residual_total - eta * tau
    epsilon = residual_EQintra / sigma
    epsilonT = residual_total / sigmaT

    return epsilonT, eta, epsilon

if __name__ == '__main__':

    if 0:
        # test Vs30 and Z10, Z25 calculation
        Vs30 = 863
        Z10 = calc_Z1(Vs30, 'CY')
        print(Z10)
        Z10 = calc_Z1(Vs30, 'AS')
        print(Z10)
        Z25 = calc_Z25(Vs30, Z1model='CY')
        print(Z25)

    if 0:
        # Test srfFaultSurfaceExtract
        FilePath = './Validation/DistancesTestFiles/inputs/158_0'
        FileName = '158_0.txt.variation-s0000-h0000'
        SRFfile = os.path.join(FilePath, FileName)
        srfFaultSurfaceTest(SRFfile)

    if 0:
        x1, y1 = 0, 0
        x2, y2 = 1, 1
        px, py = 1, 0
        dist, PP = ptToLine2D(x1, y1, x2, y2, px, py)
        print(dist)
        print(PP)

    if 0:
        point0 = [1, 1, 1]
        points = [[0, 0, 1], [0, 1, 0], [-1.0, 0.5, 0.5]]
        dist, point1 = ptToSurf3D(point0, points[0],
                                  points[1], points[2])
        print(dist, point1)

    if 0:
        point0 = [0.5, 0.5, 1]
        point0 = [0, 2, 1]
        point1 = [0, 0, 0]
        point2 = [0, 1, 0]
        dist, Ploc = ptToLineSeg3D(point0, point1, point2)
        print(dist, Ploc)
    if 0:
        point0 = [-118.286, 34.0192, 0.0]
        point1 = [-118.23476100000001, 34.067455000000002, -3.0]
        point2 = [-118.29677, 34.112583000000001, -3.0]
        dist_tmp, Ploc = ptToLineSeg3D(point0, point1, point2,
                                       Rscale=111.12)
        a, hD, vD, az = LonLatToAngleDistance(point0, Ploc,
                                              CalcRadius=False,
                                              CalcDist=True,
                                              CalcAzimuth=False, Fast=True)
        dist = np.sqrt(hD**2 + vD**2)
        print(dist, Ploc)

    if 0:
        point0 = [-117.8, 33.5, 0.0]
        point1 = [-118.0, 33, 3.0]
        point2 = [-118.0, 34, 3.0]
        point3 = [-117.5, 33, 15.0]
        dist = ptToSurf3D(point0, point1, point2, point3)
        print(dist)
