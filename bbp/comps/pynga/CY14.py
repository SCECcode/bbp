#!/usr/bin/env python
from __future__ import division, print_function

import os
import numpy as np

from pynga.utils import (GetKey, mapfunc, calc_dip, calc_Ztor,
                         calc_Rx, calc_Rrup, calc_Zhypo, calc_W,
                         rake2ftype_CY)

class CY14_nga(object):
    """
    Class for Chiou and Youngs 2014 NGA model
    """
    def __init__(self):
        self.filepth = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                    'NGA_west2')
        self.CoefFile = os.path.join(self.filepth, 'CY14.csv')
        self.Coefs = {}
        self.read_model_coefs()
        self.countries = ['California', 'Japan']

        # period independent parameters
        self.c2 = 1.06
        self.c4 = -2.1
        self.c4a = -0.5
        self.cRB = 50
        self.c8 = 0.2153
        self.c8a = 0.2695

    def read_model_coefs(self):
        self.CoefKeys = open(self.CoefFile, 'r').readlines()[1].strip().split(',')[1:]
        inputs = np.loadtxt(self.CoefFile, skiprows=2, delimiter=',')
        self.periods = inputs[:, 0]
        coefs = inputs[:, 1:]
        for i in range(len(self.periods)):
            T1 = self.periods[i]
            Tkey = GetKey(T1)
            # periods list ( -2: PGV, -1: PGA ) (mapping between the
            # NGA models accordingly, -1: PGV, 0: PGA)
            if Tkey == '-1.000':
                Tkey = '-2.000'    # PGV
                self.periods[i] = -2
            if Tkey == '0.000':
                Tkey = '-1.000'    # PGA
                self.periods[i] = -1

            self.Coefs[Tkey] = {}
            for ikey in range(len(self.CoefKeys)):
                key = self.CoefKeys[ikey]
                cmd = "self.Coefs['%s']['%s'] = coefs[%i,%i]" % (Tkey, key,
                                                                 i, ikey)
                exec(cmd)

    # call the function
    def __call__(self, M, Rjb, Vs30, T, rake, Ftype=None,
                 Rrup=None, Rx=None, dip=None, Ztor=None, Z10=None,
                 W=None, Zhypo=None, azimuth=None, Fhw=None, D_DPP=0,
                 AS=0, VsFlag=1, country='California',
                 CoefTerms={'terms':(1, 1, 1, 1, 1, 1, 1), 'NewCoefs':None}):

        if T == -1:
            T = 0.01     # for CY model, PGA's coefficients share with SA(0.01)
        if T in self.periods:
            self.T = T
        else:
            print('T is not in periods list, try to interpolate')
            raise ValueError

        # required inputs
        self.M = M         # Moment Magnitude
        self.Rjb = Rjb     # Joyner-Boore distance (km)
        self.rake = rake   # rake angle
        self.Vs30 = Vs30   # site-condition (m/s)
        self.AS = AS       # Aftershock flag (0 or 1)  (depends on the earthquake itself)
        self.VsFlag = VsFlag # 0: inferred Vs30; 1: measured Vs30
        self.country = country

        terms = CoefTerms['terms']
        NewCoefs = CoefTerms['NewCoefs']

        # Obtain optional parameters
        if Ftype != None:
            self.Fnm = 1 * (Ftype == 'NM')
            self.Frv = 1 * (Ftype == 'RV')
        else:
            if rake == None or rake < -180 or rake > 180.:
                print('rake angle should be within [-180,180]')
                raise ValueError
            else:
                self.Frv, self.Fnm = rake2ftype_CY(self.rake)

        if W == None:
            if self.rake == None:
                print('you should give either the fault width W '
                      'or the rake angle')
                raise ValueError
            else:
                W = calc_W(self.M, self.rake)
        else:
            self.W = W

        if dip == None:
            if self.rake == None:
                print('you should give either the fault dip '
                      'angle or the rake angle')
                raise ValueError
            else:
                self.dip = calc_dip(self.rake)
        else:
            self.dip = dip

        if Zhypo == None:
            self.Zhypo = calc_Zhypo(self.M, self.rake)
        else:
            self.Zhypo = Zhypo

        if Ztor == None:
            if Zhypo == None:
                if self.rake == None:
                    print('you should give either the Ztor or the rake angle')
                    raise ValueError
                else:
                    Zhypo = calc_Zhypo(self.M, self.rake)
            self.Ztor = calc_Ztor(W, self.dip, Zhypo)
        else:
            self.Ztor = Ztor

        if Fhw == None:
            if azimuth == None and Rx == None:
                print('either one of azimuth angle, Rx and Fhw '
                      'has to be specified')
                raise ValueError

            if azimuth != None:
                if 0 <= azimuth <= 180. and dip != 90.:
                    Fhw = 1
                else:
                    Fhw = 0

            elif Rx != None:
                if Rx >= 0 and dip != 90.:
                    Fhw = 1
                else:
                    Fhw = 0

            if dip == 90:
                Fhw = 0
        else:
            if azimuth == None:
                if Fhw == 1:
                    azimuth = 50
                else:
                    azimuth = -50.
        self.Fhw = Fhw

        # Compute Rx and Rrup
        if azimuth is None:
            Rx = 0.0
        elif azimuth == 90.:
            Rx = (Rrup / np.sin(self.dip * np.pi / 180.) -
                  Ztor / np.tan(self.dip * np.pi / 180.))
        elif azimuth > 0.0:
            Rx = Rjb * np.tan(azimuth * np.pi / 180.)
        elif azimuth <= 0.0:
            Rx = 0.0
        if Rx == None:
            self.Rx = calc_Rx(self.Rjb, self.Ztor, W, self.dip, azimuth, Rrup)
        else:
            self.Rx = Rx
        if Rrup == None:
            self.Rrup = calc_Rrup(self.Rx, self.Ztor, W,
                                  self.dip, azimuth, self.Rjb)
        else:
            self.Rrup = Rrup

	# Z10 (empirical relationship depends on dataset used to
	# obtain the relationship)
        if Z10 == None:
            if country == 'Japan':
                self.Z10 = np.exp(-5.23 / 2. * np.log((Vs30**2 + 412.**2) /
                                                      (1360.**2 + 412.**2)))
            else:
                self.Z10 = np.exp(-7.15 / 4. * np.log((Vs30**4 + 571.**4) /
                                                      (1360.**4 + 571.**4)))
        else:
            self.Z10 = Z10   # Z10 should be in meter  (for CY14 model)
        self.Z10 = self.Z10 / 1000.

        # directivity parameter
        if D_DPP == None:
            # compute D_DPP  Chiou and Spudich 2013
            pass
        else:
            self.D_DPP = D_DPP

        # update coeficient
        if NewCoefs != None:
            NewCoefKeys = NewCoefs.keys()
            Tkey = GetKey(self.T)
            for key in NewCoefKeys:
                self.Coefs[Tkey][key] = NewCoefs[key]

        IM = self.compute_im() # in g
        sigma, tau, sigmaT = self.calc_sigma_tau() # in ln(g)

        return IM, sigmaT, tau, sigma

    def flt_function(self):
        Ti = GetKey(self.T)

        c1 = self.Coefs[Ti]['c1']
        c1a = self.Coefs[Ti]['c1a']
        c1b = self.Coefs[Ti]['c1b']
        c1c = self.Coefs[Ti]['c1c']
        c1d = self.Coefs[Ti]['c1d']
        c7 = self.Coefs[Ti]['c7']
        c7b = self.Coefs[Ti]['c7b']
        c11 = self.Coefs[Ti]['c11']
        c11b = self.Coefs[Ti]['c11b']
        tmp = np.cosh(2 * max([self.M - 4.5, 0]))
        term0 = (c1 + (c1a + c1c / tmp) * self.Frv +
                 (c1b + c1d / tmp) * self.Fnm) # faulting type

        MeanZtor = self.calc_MeanZtor()
        D_Ztor = self.Ztor - MeanZtor
        term1 = (c7 + c7b / tmp) * D_Ztor # Ztor

        # Dip related
        term2 = (c11 + c11b / tmp) * (np.cos(self.dip * np.pi / 180.))**2

        #print 'f_flt=',term0+term1+term2
        return term0 + term1 + term2

    def calc_MeanZtor(self, M=None, Frv=None):
        if M == None:
            M = self.M
        if Frv == None:
            Frv = self.Frv

        MeanZtor = (Frv * (max([2.704 - 1.226 * max(M - 5.849, 0), 0]))**2 +
                    (1 - Frv) * (max([2.673 - 1.136 * max(M - 4.970, 0),
                                      0]))**2)
        #print 'MeanZtor=',MeanZtor
        return MeanZtor

    def moment_function(self):
        Ti = GetKey(self.T)
        c3 = self.Coefs[Ti]['c3']
        cn = self.Coefs[Ti]['cn']
        cM = self.Coefs[Ti]['cM']
        term2 = (self.c2 * (self.M - 6) +
                 (self.c2 - c3) / cn * np.log(1 + np.exp(cn * (cM - self.M))))
        #print 'f_mag=', term2
        return term2

    def distance_function(self):
        Ti = GetKey(self.T)
        c5 = self.Coefs[Ti]['c5']
        c6 = self.Coefs[Ti]['c6']
        cHM = self.Coefs[Ti]['cHM']
        cg1 = self.Coefs[Ti]['cg1']
        cg2 = self.Coefs[Ti]['cg2']
        cg3 = self.Coefs[Ti]['cg3']
        term3 = (self.c4 *
                 np.log(self.Rrup + c5 * np.cosh(c6 * max(self.M-cHM, 0))))
        term4 = ((self.c4a - self.c4) *
                 np.log(np.sqrt(self.Rrup**2 + self.cRB**2)))
        term5 = (cg1 + cg2 / np.cosh(max(self.M - cg3, 0))) * self.Rrup
        #print 'f_dis=', term3+term4+term5
        return term3 + term4 + term5

    def directivity_function(self):
        Ti = GetKey(self.T)

        c8b = self.Coefs[Ti]['c8b']
        d_taper = max([1 - max([self.Rrup - 40, 0]) / 30., 0])
        m_taper = min([max([self.M - 5.5, 0]) / 0.8, 1])
        term6 = (self.c8 * d_taper * m_taper *
                 np.exp(-self.c8a * (self.M - c8b)**2) * self.D_DPP)
        #print 'f_dir=',term6
        return term6

    def hw_function(self):
        Ti = GetKey(self.T)
        c9 = self.Coefs[Ti]['c9']
        c9a = self.Coefs[Ti]['c9a']
        c9b = self.Coefs[Ti]['c9b']

        d = self.dip * np.pi / 180.
        term7 = (c9 * self.Fhw * np.cos(d) * (c9a + (1 - c9a) *
                                              np.tanh(self.Rx / c9b)) *
                 (1 - np.sqrt(self.Rjb**2 + self.Ztor**2) / (self.Rrup + 1)))
        # print 'ihw, f_hw: ',self.Fhw, term7
        return term7

    def lnYref(self):
        # assume the site effects and basin effects are zero here
        lnYref = (self.moment_function() + self.distance_function() +
                  self.flt_function() + self.directivity_function() +
                  self.hw_function())
        return lnYref

    def site_function(self):
        Ti = GetKey(self.T)

        f1 = self.Coefs[Ti]['phi1']
        f2 = self.Coefs[Ti]['phi2']
        f3 = self.Coefs[Ti]['phi3']
        f4 = self.Coefs[Ti]['phi4']

        lnY_ref = self.lnYref()
        # print "Yref = ", np.exp(lnY_ref)
        term8 = f1 * min([np.log(self.Vs30 / 1130.), 0])
        term9 = (f2 * (np.exp(f3 * (min(self.Vs30, 1130) - 360)) -
                       np.exp(f3 * (1130 - 360))) *
                 np.log((np.exp(lnY_ref) + f4) / f4))
        #  print "f_site = ", term8+term9
        return term8 + term9

    def calc_MeanZ10(self, Vs30=None, country='California'):
        if Vs30 == None:
            Vs30 = self.Vs30
        if country == 'California':
            MeanLnZ10 = (-7.15 / 4. * np.log((Vs30**4 + 571.**4) /
                                             (1360.**4 + 571.**4)))
        elif country == 'Japan':
            MeanLnZ10 = (-5.23 / 2. * np.log((Vs30**2 + 412.**2) /
                                             (1360.**2 + 412.**2)))
        else:
            # for other region, just use default California
            MeanLnZ10 = (-7.15 / 4. * np.log((Vs30**4 + 571.**4) /
                                             (1360.**4 + 571.**4)))
        return np.exp(MeanLnZ10) / 1000.      # in km

    def basin_function(self, Z10=None, Tother=None):
        if Tother != None:
            Ti = GetKey(Tother)
        else:
            Ti = GetKey(self.T)

        if Z10 != None:
            self.Z10 = Z10

        MeanZ10 = self.calc_MeanZ10(country=self.country)

        #D_Z10 = abs(self.Z10-MeanZ10)
        D_Z10 = self.Z10 - MeanZ10

        phi5 = self.Coefs[Ti]['phi5']
        phi6 = self.Coefs[Ti]['phi6']

        term10 = phi5 * (1 - np.exp(-D_Z10 / phi6))
        # print 'f_basin = ', term10
        return term10

    def compute_im(self, terms=(1, 1, 1, 1, 1, 1, 1)):
        # use this one
        IM = np.exp(terms[0] * self.moment_function() +
                    terms[1] * self.flt_function() +
                    terms[2] * self.hw_function() +
                    terms[3] * self.distance_function() +
                    terms[4] * self.directivity_function() +
                    terms[5] * self.basin_function() +
                    terms[6] * self.site_function())
        #print 'IM: ', IM
        return IM

    def calc_NL(self):

        Ti = GetKey(self.T)

        f2 = self.Coefs[Ti]['phi2']
        f3 = self.Coefs[Ti]['phi3']
        f4 = self.Coefs[Ti]['phi4']

        yref = np.exp(self.lnYref())

        b = (f2 * (np.exp(f3 * (min(self.Vs30, 1130) - 360)) -
                   np.exp(f3 * (1130 - 360))))  # Eqn10

        #print 'NL=', b*yref / (yref+f4)
        return b * yref / (yref + f4)

    def calc_sigma_tau(self):
        Ti = GetKey(self.T)
        if self.VsFlag == 0:
            Finfer = 1
            Fmeasure = 0
        else:
            Finfer = 0
            Fmeasure = 1

        sigma1 = self.Coefs[Ti]['sigma1']
        sigma2 = self.Coefs[Ti]['sigma2']
        sigma3 = self.Coefs[Ti]['sigma3']
        tau1 = self.Coefs[Ti]['tau1']
        tau2 = self.Coefs[Ti]['tau2']
        NL = self.calc_NL()

        tmp = min([max([self.M, 5]), 7.25]) - 5
        tau = tau1 + (tau2 - tau1) / 2.25 * tmp
        sigma = ((sigma1 + (sigma2 - sigma1) / 2.25 * tmp) *
                 np.sqrt(sigma3 * Finfer + 0.7 * Fmeasure + (1 + NL)**2))

        # correct tau
        tauNL = (1 + NL) * tau
        sigmaT = np.sqrt(sigma**2 + tauNL**2)

        #return (sigma, tau, sigmaT)
        return (sigma, tauNL, sigmaT)

def CY14nga_test(T, CoefTerms):
    """
    Test CY nga model
    """
    M = 7.0
    Zhypo = 20.
    Ztor = 3.0
    dip = 30
    Ftype = 'RV'
    rake = 0    # for specific rupture
    W = 50.0

    Rjb = 0.0
    Rrup = (W * np.sin(dip * np.pi / 180.) + Ztor) * np.cos(dip * np.pi / 180.)
    Rx = W * np.cos(dip * np.pi / 180.)

    #print "Rx", Rx
    #Rx = Rrup

    Vs30 = 748.0, 1200., 345., 160.
    Vs30 = 180.
    Z25 = 4.0
    Z10 = 1.0

    AS = 0
    VsFlag = 0

    CYnga = CY14_nga()

    kwds = {'Ftype':Ftype, 'Ztor':Ztor, 'dip':dip, 'Rrup':Rrup,
            'Rx':Rx, 'Z10':Z10, 'AS':AS, 'VsFlag':VsFlag,
            'CoefTerms':CoefTerms}
    values = mapfunc(CYnga, M, Rjb, Vs30, T, rake, **kwds)
    print('Median, SigmaT, Tau, Sigma')
    for i in range(len(values)):
        print(values[i])
    return CYnga

if __name__ == '__main__':

    NewCoefs = None
    CoefTerms = {'terms':(1, 1, 1, 1, 1, 1, 1), 'NewCoefs':NewCoefs}

    for T in [0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.1, 0.12,
              0.15, 0.17, 0.2, 0.25, 0.3, 0.4, 0.5, 0.75, 1.0,
              1.5, 2.0, 3.0, 4.0, 5.0, 7.5, 10.0]:
    #for T in [0.3, ]:
        print('CY GM at %s' % ('%3.2f' % T))
        CYnga = CY14nga_test(T, CoefTerms)


