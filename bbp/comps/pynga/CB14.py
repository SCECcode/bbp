#!/usr/bin/env python
"""
CB 2014 NGA model
"""
from __future__ import division, print_function

import os
import numpy as np

from pynga.utils import (GetKey, mapfunc, calc_Rx, calc_Rrup,
                         calc_Ztor, calc_Zhypo, calc_W, calc_dip,
                         rake2ftype_CB)

class CB14_nga(object):
    """
    Class of NGA model of Campbell and Bozorgnia 2008
    """
    def __init__(self):
        """
	Model initialization
	"""
        self.filepth = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'NGA_west2')
        self.CoefFile = os.path.join(self.filepth, 'CB14.csv')
        self.Coefs = {}
        self.read_model_coefs()
        # for distinguished anelastic attenuation for different countries
        self.regions = ['CA', 'JP', 'CH']
        self.c = 1.88
        self.n = 1.18

    def read_model_coefs(self):
        #print len(open(self.CoefFile,'r').readlines())
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

    def __call__(self, M, Rjb, Vs30, T, rake, Rrup=None, Ftype=None,
                 dip=None, Zhypo=None, Ztor=None, Z25=None,
                 W=None, Rx=None, azimuth=None, Fhw=0,
                 Z10=None, Z15=None, Arb=0, region='CA', SJ=0,
                 CoefTerms={'terms':(1, 1, 1, 1, 1, 1, 1, 1, 1),
                            'NewCoefs':None}):
        """
        Call the class to compute median ground-motion intensity
        You have to call the function here to make the class rich
        region: indicate the anelastic atteneuation
        SJ = 1: Japan basin effects; SJ = 0: other regions
        Note: Zhypo, dip, Rjb, and Rx now is required
        """
        # Those inputs have to be specified
        self.M = M # moment magnitude
        self.Rjb = float(Rjb) # Joyner-Boore distance (km)
        self.Vs30 = float(Vs30) # time-averaged shear wave velocity over 30m subsurface depth (m/s)
        self.T = T # select period (sec)
        self.rake = rake # rake could be None then you have to give the W and dip
        self.region = region
        self.SJ = SJ
        terms = CoefTerms['terms']
        NewCoefs = CoefTerms['NewCoefs']
        self.Fhw = Fhw

        # check inputs
        if T in self.periods:
            self.T = T
        else:
            print('T is not in periods list, try to interpolate')
            raise ValueError

        if self.M == None or self.M < 0:
            print('Moment magnitude must be a postive number')
            raise ValueError
        if self.Rjb == None or self.Rjb < 0:
            print('Joyner-Boore distance must be a non-negative number')
            raise ValueError
        if self.Vs30 == None or self.Vs30 < 0:
            print('Vs30 must be a positive number')
            raise ValueError

        # Determine the Fault-related parameters (if necessary)
        if Ftype != None:
            self.Fnm = 1 * (Ftype == 'NM')
            self.Frv = 1 * (Ftype == 'RV')
        else:
            if rake == None or rake < -180 or rake > 180.:
                print('rake angle should be within [-180,180]')
                raise ValueError
            else:
                self.Frv, self.Fnm = rake2ftype_CB(self.rake)

        if W == None:
            if self.rake == None:
                print('you should give either the fault width W '
                      'or the rake angle')
                raise ValueError
            else:
                self.W = calc_W(self.M, self.rake)
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

        # Determine Site-Source related parameters (if necessary)
        if Rrup == None:
            if azimuth == None:
                if Fhw != None:
                    if Fhw == 1:
                        azimuth = 50   # hanging wall site
                    else:
                        azimuth = -50. # footwall site
                else:
                    azimuth = -50.
                    Fhw = 0
            if self.Rjb == 0:
                Fhw = 1
                azimuth = 90
            Rx_tmp = calc_Rx(self.Rjb, self.Ztor, W,
                             self.dip, azimuth, Rrup=Rrup)
            self.Rrup = calc_Rrup(Rx_tmp, self.Ztor, W,
                                  self.dip, azimuth, Rjb=self.Rjb)
        else:
            self.Rrup = Rrup
        if azimuth is None:
            # For Python 3 compatibility
            Rx = 0.0
        elif azimuth == 90.:
            Rx = (Rrup / np.sin(self.dip * np.pi / 180.) -
                  Ztor / np.tan(self.dip * np.pi / 180.))
        elif azimuth > 0.0:
            Rx = Rjb * np.tan(azimuth * np.pi / 180.)
        elif azimuth <= 0.0:
            Rx = 0.0
        if Rx == None:
            if azimuth == None:
                if Fhw != None:
                    if Fhw == 1:
                        azimuth = 50   # hanging wall site
                    else:
                        azimuth = -50. # footwall site
                else:
                    azimuth = -50.
                    Fhw = 0
            if self.Rjb == 0:
                Fhw = 1
                azimuth = 90
            self.Rx = calc_Rx(self.Rjb, self.Ztor, W,
                              self.dip, azimuth, Rrup=Rrup)
        else:
            self.Rx = Rx
            Fhw = 0 * (Rx < 0) + 1*(Rx >= 0)
        self.Fhw = Fhw

        # Determine Site-Specific parameters (those empirical
        # relationships dependes on database)
        if Z25 == None:
            # if Z25 not provided, use the default values
            if region == 'CA':
                self.Z25 = np.exp(7.089 - 1.144 * np.log(Vs30))
            elif region == 'JP':
                self.Z25 = np.exp(5.359 - 1.102 * np.log(Vs30))
            else:
                self.Z25 = np.exp(6.510 - 1.181 * np.log(Vs30))
        else:
            self.Z25 = Z25  # input Z25 should be in km

        # update coeficient (use updated coefficients)
        if NewCoefs != None:
            NewCoefKeys = NewCoefs.keys()
            Tkey = GetKey(self.T)
            for key in NewCoefKeys:
                self.Coefs[Tkey][key] = NewCoefs[key]

        # Compute IM and Standard deviation
        IM = self.compute_im(terms=terms)
        sigma, tau, sigmaT, sigmaArb = self.sd_calc()
        if Arb == 0:
            return IM, sigmaT, tau, sigma
        else:
            return IM, sigmaArb, tau, sigma

    # ============================
    # Function used in this class
    # ============================
    def moment_function(self, Tother=None):
        """
        Moment term
        """
        if Tother != None:
            Ti = GetKey(Tother)
        else:
            Ti = GetKey(self.T)

        c0 = self.Coefs[Ti]['c0']
        c1 = self.Coefs[Ti]['c1']
        c2 = self.Coefs[Ti]['c2']
        c3 = self.Coefs[Ti]['c3']
        c4 = self.Coefs[Ti]['c4']

        if self.M <= 4.5:
            f_mag = c0 + c1 * self.M
        elif 4.5 < self.M <= 5.5:
            f_mag = c0 + c1 * self.M + c2 * (self.M - 4.5)
        elif 5.5 < self.M <= 6.5:
            f_mag = (c0 + c1 * self.M + c2 * (self.M - 4.5) +
                     c3 * (self.M - 5.5))
        else:
            f_mag = (c0 + c1 * self.M + c2 * (self.M - 4.5) +
                     c3 * (self.M - 5.5) + c4 * (self.M - 6.5))
        #print 'f_mag:',f_mag
        return f_mag

    def distance_function(self, Tother=None):
        """
        Geometrical attenuation
        """
        if Tother != None:
            Ti = GetKey(Tother)
        else:
            Ti = GetKey(self.T)

        c5 = self.Coefs[Ti]['c5']
        c6 = self.Coefs[Ti]['c6']
        c7 = self.Coefs[Ti]['c7']

        Rtmp = np.sqrt(self.Rrup**2 + c7**2)
        f_dis = (c5 + c6 * self.M) * np.log(Rtmp)
        #print 'f_dis:', f_dis
        return f_dis

    def attenuation_function(self, Tother=None):
        """
        Anelastic Attneuation
        """
        if Tother != None:
            Ti = GetKey(Tother)
        else:
            Ti = GetKey(self.T)
        c20 = self.Coefs[Ti]['c20']
        D_c20 = self.Coefs[Ti]['D_c20_%s' % self.region]
        if self.Rrup > 80:
            f_atn = (c20 + D_c20) * (self.Rrup - 80)
        else:
            f_atn = 0.0
        #print 'f_atn:', f_atn
        return f_atn

    def fault_function(self, Tother=None):
        """
        Fault mechanism term
        or style of the fault
        """
        if Tother != None:
            Ti = GetKey(Tother)
        else:
            Ti = GetKey(self.T)
        c8 = self.Coefs[Ti]['c8']
        c9 = self.Coefs[Ti]['c9']
        f_fltF = c8 * self.Frv + c9 * self.Fnm
        f_fltM = (0 * (self.M <= 4.5) +
                  (self.M - 4.5) * (4.5 < self.M <= 5.5) + 1 * (self.M > 5.5))
        f_flt = f_fltF * f_fltM
        #print 'f_flt:', f_flt
        return f_flt

    def hw_function(self, Tother=None):
        """
        Hanging Wall term
        """
        if Tother != None:
            Ti = GetKey(Tother)
        else:
            Ti = GetKey(self.T)

        c10 = self.Coefs[Ti]['c10']
        a2 = self.Coefs[Ti]['a2']
        h1 = self.Coefs[Ti]['h1']
        h2 = self.Coefs[Ti]['h2']
        h3 = self.Coefs[Ti]['h3']
        h4 = self.Coefs[Ti]['h4']
        h5 = self.Coefs[Ti]['h5']
        h6 = self.Coefs[Ti]['h6']

        #print 'Rx:', self.Rx
        R1 = self.W * np.cos(self.dip * np.pi / 180.)
        R2 = 62 * self.M - 350
        tmp1 = self.Rx / R1
        tmp2 = (self.Rx - R1) / (R2 - R1)
        f1_Rx = h1 + h2 * tmp1 + h3 * (tmp1)**2
        f2_Rx = h4 + h5 * tmp2 + h6 * (tmp2)**2
        f_hngRx = (0 * (self.Rx < 0) +
                   f1_Rx * (0.0 <= self.Rx < R1) +
                   max([f2_Rx, 0.0]) * (self.Rx >= R1))

        f_hngRrup = (1 * (self.Rrup == 0) +
                     (self.Rrup - self.Rjb) / self.Rrup * (self.Rrup > 0))
        f_hngM = (0 * (self.M <= 5.5) +
                  ((self.M - 5.5) * (1 + a2 * (self.M - 6.5))) *
                  (5.5 < self.M <= 6.5) +
                  (1 + a2 * (self.M - 6.5)) * (self.M > 6.5))
        f_hngZ = ((1 - 0.06 * self.Ztor) * (self.Ztor <= 16.66) +
                  0 * (self.Ztor > 16.66))
        f_hngD = (90 - self.dip) / 45.
        f_hng = c10 * f_hngRx * f_hngRrup * f_hngM * f_hngZ * f_hngD * self.Fhw
        #print 'f_hng: ', f_hng
        return f_hng

    def rup_dip_function(self, Tother=None):
        if Tother != None:
            Ti = GetKey(Tother)
        else:
            Ti = GetKey(self.T)
        c19 = self.Coefs[Ti]['c19']
        f_dip = (c19 * self.dip * (self.M <= 4.5) +
                 c19 * (5.5 - self.M) * self.dip *
                 (4.5 < self.M <= 5.5) + 0 * (self.M > 5.5))
        #print 'f_dip:', f_dip
        return f_dip

    def hypo_depth_function(self, Tother=None):
        if Tother != None:
            Ti = GetKey(Tother)
        else:
            Ti = GetKey(self.T)

        c17 = self.Coefs[Ti]['c17']
        c18 = self.Coefs[Ti]['c18']
        f_hypoH = (0 * (self.Zhypo <= 7) +
                   (self.Zhypo - 7) * (7 < self.Zhypo <= 20) +
                   13 * (self.Zhypo > 20))
        f_hypoM = (c17 * (self.M <= 5.5) +
                   (c17 + (c18 - c17) * (self.M - 5.5)) *
                   (5.5 < self.M <= 6.5) + c18 * (self.M > 6.5))
        f_hyp = f_hypoH * f_hypoM
        #print 'f_hyp:',f_hyp
        return f_hyp

    def basin_function(self, Tother=None, Z25=None):
        """
        Basin-effect term
        """
        if Tother != None:
            Ti = GetKey(Tother)
        else:
            Ti = GetKey(self.T)
        if Z25 == None:
            Z25 = self.Z25

        c14 = self.Coefs[Ti]['c14']
        c15 = self.Coefs[Ti]['c15']
        c16 = self.Coefs[Ti]['c16']
        k3 = self.Coefs[Ti]['k3']
        if Z25 <= 1.0:
            f_sed = (c14 + c15 * self.SJ) * (Z25 - 1.0)
        elif 1.0 < Z25 <= 3.0:
            f_sed = 0.0
        else:
            f_sed = (c16 * k3 * np.exp(-0.75) *
                     (1 - np.exp(-0.25 * (Z25 - 3.0))))
        #print 'f_sed:', f_sed
        return f_sed

    def A1100_calc(self):
        Tother = -1.0
        A1100 = np.exp(self.moment_function(Tother=Tother) +
                       self.distance_function(Tother=Tother) +
                       self.attenuation_function(Tother=Tother) +
                       self.fault_function(Tother=Tother) +
                       self.hw_function(Tother=Tother) +
                       self.rup_dip_function(Tother=Tother) +
                       self.hypo_depth_function(Tother=Tother) +
                       self.basin_function(Tother=Tother) +
                       self.site_function(A1100=0, Vs30=1100., Tother=Tother))
        #print 'PGA1100:', A1100
        #print '=================================='
        return A1100

    def site_function(self, A1100=None, Vs30=None, Tother=None):
        """
        Shallow site effect term
        Be careful to the input variables (they are keys, not arguments)
        """

        # PGA at reference rock that has Vs30 = 1100 (unit: m/s)
        if A1100 == None:
            A1100 = self.A1100_calc()

        if Vs30 == None:
            Vs30 = self.Vs30

        if Tother != None:
            Ti = GetKey(Tother)
        else:
            Ti = GetKey(self.T)

        c11 = self.Coefs[Ti]['c11']
        c12 = self.Coefs[Ti]['c12']
        c13 = self.Coefs[Ti]['c13']
        k1 = self.Coefs[Ti]['k1']
        k2 = self.Coefs[Ti]['k2']

        if Vs30 <= k1:
            f_siteG = (c11 * np.log(Vs30 / k1) +
                       k2 * (np.log(A1100 + self.c * (Vs30 / k1)**self.n) -
                             np.log(A1100 + self.c)))
        else:
            f_siteG = (c11 + k2 * self.n) * np.log(Vs30 / k1)

        if Vs30 <= 200:
            f_siteJ = ((c12 + k2 * self.n) *
                       (np.log(Vs30 / k1) - np.log(200. / k1)))
        else:
            f_siteJ = (c13 + k2 * self.n) * np.log(Vs30 / k1)

        f_site = f_siteG + self.SJ*f_siteJ
        #print 'f_site:', f_site
        return f_site

    # Final function to compute Sa, PGA, PGV
    def compute_im(self, terms=(1, 1, 1, 1, 1, 1, 1, 1, 1)):
        """
        Compute IM based on functional form of CB08 model
        """
        IM = np.exp(terms[0] * self.moment_function() +
                    terms[1] * self.fault_function() +
                    terms[2] * self.hw_function() +
                    terms[3] * self.rup_dip_function() +
                    terms[4] * self.hypo_depth_function() +
                    terms[5] * self.distance_function() +
                    terms[6] * self.attenuation_function() +
                    terms[7] * self.basin_function() +
                    terms[8] * self.site_function())
        if self.T <= 0.25 and self.T != -2: # and self.T != -1.0:
            Tother = -1.0
            # This is PGA itself
            IM1 = np.exp(terms[0] * self.moment_function(Tother=Tother) +
                         terms[1]*self.fault_function(Tother=Tother) +
                         terms[2]*self.hw_function(Tother=Tother) +
                         terms[3]*self.rup_dip_function(Tother=Tother) +
                         terms[4]*self.hypo_depth_function(Tother=Tother) +
                         terms[5]*self.distance_function(Tother=Tother) +
                         terms[6]*self.attenuation_function(Tother=Tother) +
                         terms[7]*self.basin_function(Tother=Tother) +
                         terms[8]*self.site_function(Tother=Tother))
            if IM < IM1:
                # This is for SA (not for PGA and PGV, since they are
                # computed above)
                IM = IM1
        #print 'IM: ', IM
        #print '=============================================='
        return IM

    # function used to compute standard deviation terms
    def alpha_calc(self, Vs30=None, Tother=None):
        if Vs30 == None:
            Vs30 = self.Vs30
        if Tother == None:
            Ti = GetKey(self.T)
        else:
            Ti = GetKey(Tother)

        k1 = self.Coefs[Ti]['k1']
        k2 = self.Coefs[Ti]['k2']

        A1100 = self.A1100_calc()

        # compute alpha
        if Vs30 < k1:
            alpha = (k2 * A1100 *
                     (1. / (A1100 + self.c * (Vs30 / k1)**self.n) -
                      1. / (A1100 + self.c)))
        else:
            alpha = 0

        return alpha

    def sigma_tau_lnY(self, Tother=None):
        if Tother != None:
            Ti = GetKey(Tother)
        else:
            Ti = GetKey(self.T)
        phi1 = self.Coefs[Ti]['phi1']
        phi2 = self.Coefs[Ti]['phi2']
        tau1 = self.Coefs[Ti]['tau1']
        tau2 = self.Coefs[Ti]['tau2']
        phi_lnY = (phi1 * (self.M <= 4.5) +
                   (phi2 + (phi1 - phi2) * (5.5 - self.M)) *
                   (4.5 < self.M <= 5.5) + phi2 * (self.M > 5.5))
        tau_lnY = (tau1 * (self.M <= 4.5) +
                   (tau2 + (tau1 - tau2) * (5.5 - self.M)) *
                   (4.5 < self.M <= 5.5) + tau2 * (self.M > 5.5))
        return phi_lnY, tau_lnY

    def sigma_tau_calc(self, Vs30=None, Tother=None):
        """
        Intra-event and inter-vent residual standard deviation
        """
        if Tother != None:
            Ti = GetKey(Tother)
        else:
            Ti = GetKey(self.T)
        if Vs30 == None:
            Vs30 = self.Vs30

        rho = self.Coefs[Ti]['rho']
        phi_lnAF = self.Coefs[Ti]['phi_lnAF']
        alpha = self.alpha_calc(Vs30=Vs30, Tother=Tother)

        phi_lnPGA, tau_lnPGA = self.sigma_tau_lnY(Tother=-1)
        phi_lnY, tau_lnY = self.sigma_tau_lnY()

        phi_lnYb = np.sqrt(phi_lnY**2 - phi_lnAF**2)
        phi_lnAb = np.sqrt(phi_lnPGA**2 - phi_lnAF**2) # Ab = PGA_b

        tau_lnYb = tau_lnY
        tau_lnAb = tau_lnPGA

        # Calculate sigma
        sigma = np.sqrt(phi_lnY**2 + (alpha * phi_lnAb)**2 +
                        2 * alpha * rho * phi_lnYb * phi_lnAb)
        # calculate tau
        tau = np.sqrt(tau_lnYb**2 + (alpha * tau_lnAb)**2 +
                      2 * alpha * rho * tau_lnYb * tau_lnAb)

        return sigma, tau

    def sd_calc(self, Vs30=None, Tother=None):

        if Tother != None:
            Ti = GetKey(Tother)
        else:
            Ti = GetKey(self.T)
        if Vs30 == None:
            Vs30 = self.Vs30

        # for RotD50
        sigma, tau = self.sigma_tau_calc(Vs30=Vs30, Tother=Tother)
        sigmaT = np.sqrt(sigma**2 + tau**2)

        # For arbitrary horizontal components
        phi_c = self.Coefs[Ti]['phi_c']
        sigmaArb = np.sqrt(sigmaT**2 + phi_c**2)

        # standard deviations are in logarithm scale !!!
        return (sigma, tau, sigmaT, sigmaArb)

def CB14nga_test(T, CoefTerms):
    """
    Test CB nga model
    """
    M = 8.8
    Zhypo = 35.
    Ztor = 4.098
    dip = 18
    Ftype = 'RV'
    rake = 107    # for specific rupture
    W = 200

    Rjb = 0.0
    Rrup = [5.23399,]
    azimuth = 90.
    Rx = -2
    Fhw = None
    #Fhw = 0
    Vs30 = 760
    Z25 = 1.0

    Arb = 0

    # How to use it
    CBnga = CB14_nga()
    kwds = {'Ftype':Ftype, 'Z25':Z25, 'Rrup':Rrup, 'Zhypo':Zhypo,
            'W':W, 'Fhw':Fhw, 'Ztor':Ztor, 'azimuth':azimuth,
            'dip':dip, 'Arb':Arb, 'CoefTerms':CoefTerms}
    values = mapfunc(CBnga, M, Rjb, Vs30, T, rake, **kwds)

    for i in range(len(values)):
        print(values[i])

    return CBnga

if __name__ == '__main__':
    CoefTerms = {'terms':(1, 1, 1, 1, 1, 1, 1, 1, 1), 'NewCoefs':None}
    #for T in [0.01, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25,
    #0.3, 0.4, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 7.5, 10.0,
    #-1,-2]:
    for T in [-1, 0.3, 1.0, 3.0]:
        print('CB SA at %s' % ('%3.2f' % T))
        CBnga = CB14nga_test(T, CoefTerms)
