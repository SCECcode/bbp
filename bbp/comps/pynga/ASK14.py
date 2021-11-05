#!/usr/bin/env python
from __future__ import division, print_function

import os
import numpy as np

from pynga.utils import (GetKey, calc_W, calc_dip, calc_Zhypo,
                         calc_Ztor, mapfunc, calc_Rx, calc_Rrup,
                         calc_Z1, rake2ftype_AS)

class ASK14_nga(object):
    """
    ASK14 NGA model class
    """
    def __init__(self):

        self.filepth = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                    'NGA_west2')
        self.CoefFile = os.path.join(self.filepth, 'ASK14.csv')
        self.Coefs = {}
        self.read_model_coefs()
        self.countries = ['California', 'Japan']

        # region for regional corrections (Vs30 and Rrup)
        self.regions = ['CA', 'TW', 'CN', 'JP',]

        # New period independent parameters:
        self.M2 = 5
        self.a4 = -0.1
        self.a5 = -0.49
        self.a7 = 0
        self.a2HW = 0.2
        self.n = 1.5
        self.a39 = 0.0

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

    def __call__(self, M, Rjb, Vs30, T, rake, Ftype=None,
                 Rrup=None, Rx=None, Ry0=None, dip=None,
                 Ztor=None, Z10=None,
                 W=None, Zhypo=None, azimuth=None, Fhw=None,
                 Fas=0, CRjb=15, VsFlag=1,
                 region='CA', country='California',
                 CoefTerms={'terms':(1, 1, 1, 1, 1, 1, 1), 'NewCoefs':None}):

        # input Z10 should be in meter !!!
        # CRjb is defined for aftershock, if you set Fas = 0, then
        # this will not be used

        self.M = M    # moment magnitude
        self.Rjb = Rjb     # JB distance (km)
        self.Vs30 = Vs30    # site condition (m/s)
        self.rake = rake   # rake andgle
        self.country = country
        self.region = region

        if T in self.periods:
            self.T = T
        else:
            print('T is not in periods list, try to interpolate')
            raise ValueError

        self.c = 2.4 * (self.T != -2) + 2400 * (self.T == -2)

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
                self.Frv, self.Fnm = rake2ftype_AS(self.rake)

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

        self.azimuth = azimuth   # use the original one if available

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

        if azimuth == None:
            if Fhw == 1:
                azimuth = 50
            else:
                azimuth = -50.

        self.Fhw = Fhw

        # Compute Rrup and Rx
        if azimuth == 90.:
            Rx = (Rrup / np.sin(self.dip * np.pi / 180.) -
                  Ztor / np.tan(self.dip * np.pi / 180.))
        elif azimuth > 0.0:
            Rx = Rjb * np.tan(azimuth * np.pi / 180.)
        elif azimuth <= 0.0:
            Rx = 0.0
        if Rx == None:
            self.Rx = calc_Rx(self.Rjb, self.Ztor, self.W,
                              self.dip, azimuth, Rrup)
        else:
            self.Rx = Rx
        if Rrup == None:
            self.Rrup = calc_Rrup(self.Rx, self.Ztor, self.W,
                                  self.dip, azimuth, self.Rjb)
        else:
            self.Rrup = Rrup

        if Ry0 == None:
            if self.azimuth != None and self.Rx != None:
                self.Ry0 = self.Rx * np.tan(self.azimuth * np.pi / 180.)
            else:
                self.Ry0 = None
        else:
            self.Ry0 = Ry0    # attention here (Ry0)

        # Z10
        if Z10 == None:
            if country == 'Japan':
                self.Z10 = (np.exp(-5.23 / 2. *
                                   np.log((Vs30**2 + 412.**2) /
                                          (1360.**2 + 412.**2))))
            else:
                self.Z10 = (np.exp(-7.67 / 4. *
                                   np.log((Vs30**4 + 610.**4) /
                                          (1360.**4 + 610.**4))))
        else:
            self.Z10 = Z10

        # for ASK14, the Z10 used in calculation is in km
        self.Z10 = self.Z10 / 1000.
        self.Fas = Fas    # aftershock flag (0 or 1)
        self.CRjb = CRjb
        self.VsFlag = VsFlag    # 0: estimated Vs30; 1: measured Vs30

        # update coeficient
        if NewCoefs != None:
            NewCoefKeys = NewCoefs.keys()
            Tkey = GetKey(self.T)
            for key in NewCoefKeys:
                self.Coefs[Tkey][key] = NewCoefs[key]

        # Compute IM and uncertainties
        IM = self.compute_im(terms=terms)
        sigma, tau, sigmaT = self.calc_sigma_tau()

        return IM, sigmaT, tau, sigma

    def base_model(self, Tother=None):
        # Basically, this is the distance-magnitude term
        if Tother != None:
            Ti = GetKey(Tother)
        else:
            Ti = GetKey(self.T)

        c4 = self.Coefs[str(Ti)]['c4']
        a1 = self.Coefs[str(Ti)]['a1']
        a2 = self.Coefs[str(Ti)]['a2']
        a3 = self.Coefs[str(Ti)]['a3']
        a6 = self.Coefs[str(Ti)]['a6']
        a8 = self.Coefs[str(Ti)]['a8']
        a17 = self.Coefs[str(Ti)]['a17']
        M1 = self.Coefs[str(Ti)]['M1']
        
        c4M = (c4 * (self.M > 5) +
               (c4 - (c4 - 1) * (5 - self.M)) * (4 < self.M <= 5) +
               1 * (self.M <= 4))
        Rtmp = np.sqrt(self.Rrup**2 + c4M**2)
        if self.M < self.M2:
            output = (a1 + self.a4 * (self.M2 - M1) +
                      a8 * (8.5 - self.M2)**2 +
                      a6 * (self.M - self.M2) +
                      self.a7 * (self.M - self.M2)**2 +
                      (a2 + a3 * (self.M2 - M1)) * np.log(Rtmp) +
                      a17 * self.Rrup)
        elif self.M2 <= self.M < M1:
            output = (a1 + self.a4 * (self.M - M1) +
                      a8 * (8.5 - self.M)**2 +
                      (a2 + a3 * (self.M - M1)) * np.log(Rtmp) +
                      a17 * self.Rrup)
        elif self.M >= M1:
            output = (a1 + self.a5 * (self.M - M1) +
                      a8 * (8.5 - self.M)**2 +
                      (a2 + a3 * (self.M - M1)) * np.log(Rtmp) +
                      a17 * self.Rrup)
        #print 'f_base=', output
        return output

    def flt_function(self, Tother=None):
        # fault type and aftershock flag
        if Tother != None:
            Ti = GetKey(Tother)
        else:
            Ti = GetKey(self.T)

        a11 = self.Coefs[Ti]['a11']
        a12 = self.Coefs[Ti]['a12']
        a14 = self.Coefs[Ti]['a14']

        f7 = (a11 * (self.M > 5) +
              a11 * (self.M - 4) * (4 < self.M <= 5) +
              0 * (self.M <= 4))
        f8 = (a12 * (self.M > 5) +
              a12 * (self.M - 4) * (4 < self.M <= 5) +
              0 * (self.M <= 4))
        f11 = (a14 * (self.CRjb <= 5) +
               a14 * (1 - (self.CRjb - 5) / 10.) * (5 < self.CRjb <= 15) +
               0 * (self.CRjb > 15))

        output = self.Frv * f7 + self.Fnm * f8 + self.Fas * f11
	#print 'f_flt=', output
        return output

    def ztor_function(self, Tother=None):
        # depth to top of rupture model

        if Tother != None:
            Ti = GetKey(Tother)
        else:
            Ti = GetKey(self.T)

        a15 = self.Coefs[Ti]['a15']

        if self.Ztor < 20:
            output = a15 * self.Ztor / 20.
        else:
            output = a15
        #print 'f_ztor=', output
        return output

    def hw_function(self, Tother=None):
        # hanging wall function
        if self.Rx < 0:
            output = 0.0
        else:
            if Tother != None:
                Ti = GetKey(Tother)
            else:
                Ti = GetKey(self.T)
            a13 = self.Coefs[Ti]['a13']

            # taper1
            if self.dip > 30:
                taper1 = (90 - self.dip) / 45.
            else:
                taper1 = 60. / 45.

            # taper2
            if self.M >= 6.5:
                taper2 = 1 + self.a2HW * (self.M - 6.5)
            elif 5.5 < self.M < 6.5:
                taper2 = (1 + self.a2HW * (self.M - 6.5) -
                          (1 - self.a2HW) * (self.M - 6.5)**2)
            else:
                taper2 = 0.0

            # taper3 (constrain the hanging wall effects, this should
            # decreasing the difference between simulation and NGA
            # GMPEs)
            R1 = self.W * np.cos(self.dip * np.pi / 180.)
            R2 = 3 * R1
            h1 = 0.25
            h2 = 1.5
            h3 = -0.75
            if self.Rx < R1:
                taper3 = h1 + h2 * (self.Rx / R1) + h3 * (self.Rx / R1)**2
            elif R1 <= self.Rx <= R2:
                taper3 = 1 - (self.Rx - R1) / (R2 - R1)
            else:
                taper3 = 0.0

            # taper4 (constrain based on ztor)
            if self.Ztor <= 10:
                taper4 = 1 - self.Ztor**2 / 100.
            else:
                taper4 = 0.0

            # taper5 (constrain azimuthally)
            Ry1 = self.Rx * np.tan(20 * np.pi / 180.)
            if self.Ry0 != None:
                if self.Ry0 < Ry1:
                    taper5 = 1.0
                elif self.Ry0 - Ry1 < 5:
                    taper5 = 1 - (self.Ry0 - Ry1) / 5.
                elif self.Ry0 - Ry1 > 5:
                    taper5 = 0.0
            else:
                taper5 = (1.0 * (self.Rjb == 0) +
                          (1 - self.Rjb / 30.) *
                          (self.Rjb < 30 and self.Rjb != 0) +
                          0.0 * (self.Rjb >= 30.))
            #print 'taper1, taper2, taper3, taper4, taper5, a13:',
            #taper1, taper2, taper3, taper4, taper5, a13
            output = a13 * taper1 * taper2 * taper3 * taper4 * taper5
        #print 'Rx, f_hng=', self.Rx, output
        return output

    def calc_MeanZ10(self, Vs30=None, country='California'):
        if Vs30 == None:
            Vs30 = self.Vs30

        if country == 'Japan':
            lnZ10 = (-5.23 / 2. * np.log((Vs30**2 + 412.**2) /
                                         (1360.**2 + 412.**2)))
        else:
            lnZ10 = (-7.67 / 4. * np.log((Vs30**4 + 610.**4) /
                                         (1360.**4 + 610.**4)))
        MeanZ10 = np.exp(lnZ10) / 1000. # convert to km for program to use
        return MeanZ10

    def soil_function(self, Z10=None, Vs30=None, Tother=None):
        # soil depth function
        if Tother != None:
            Ti = GetKey(Tother)
            T = Tother
        else:
            Ti = GetKey(self.T)
            T = self.T

        if Z10 == None:
            Z10 = self.Z10

        if Vs30 == None:
            Vs30 = self.Vs30

        a43 = self.Coefs[Ti]['a43']
        a44 = self.Coefs[Ti]['a44']
        a45 = self.Coefs[Ti]['a45']
        a46 = self.Coefs[Ti]['a46']

        MeanZ10 = self.calc_MeanZ10(Vs30=Vs30, country=self.country)

        term = np.log((Z10 + 0.01) / (MeanZ10 + 0.01))
        #  print 'Z10, Z10hat, a43, term', Z10, MeanZ10, a43, term
        output = (a43 * (Vs30 <= 200.) + a44 * (200 < Vs30 <= 300) +
                  a45 * (300 < Vs30 <= 500) + a46 * (Vs30 > 500))
        output = output * term
        #print 'f_soil=',output
        return output

    # compute Vs30* for soil and site effect functions
    def CalcVs30Star(self, Vs30, T):
        # compute V1 (used in soil-depth model)
        if T <= 0.5:
            V1 = 1500. # m/s
        elif 0.5 < T <= 3.0:
            V1 = np.exp(-0.35 * np.log(T / 0.5) + np.log(1500.))
        elif T >= 3.0:
            V1 = 800.

        # calculate Vs30*
        if Vs30 < V1:
            Vs30_1 = Vs30
        else:
            Vs30_1 = V1

        return V1, Vs30_1

    def site_model(self, SA1100, Vs30=None, Tother=None):
        # Site-response model

        if Tother != None:
            Ti = GetKey(Tother)
            T = Tother
        else:
            Ti = GetKey(self.T)
            T = self.T

        a10 = self.Coefs[Ti]['a10']
        b = self.Coefs[Ti]['b']
        Vlin = self.Coefs[Ti]['VLIN']

        if Vs30 == None:
            Vs30 = self.Vs30

        V1, Vs30_1 = self.CalcVs30Star(Vs30, T)

        if Vs30 < Vlin:
            output = (a10 * np.log(Vs30_1 / Vlin) -
                      b * np.log(SA1100 + self.c) +
                      b * np.log(SA1100 + self.c * (Vs30_1 / Vlin)**self.n))
        else:
            output = (a10 + b * self.n) * np.log(Vs30_1 / Vlin)
        #print 'f_site=',output
        return output

    def SA1100_calc(self):
        # compute SA1100 (different from AS08)
        SA1100Rock = 0.0
        Vs30Rock = 1100.
        Z10Rock = calc_Z1(Vs30Rock, 'AS') / 1000.   # attention here
        Tother = self.T # SA at current period !  cout << "Z10Rock: "
        #<< Z10Rock << ", Z10hat: " << Z10hat << ", a46: " <<
        #s_a46[iT] << ", term: " << tmp << endl;
        SA1100 = (self.base_model(Tother=Tother) +
                  self.flt_function(Tother=Tother) +
                  self.site_model(SA1100Rock, Vs30=Vs30Rock, Tother=Tother) +
                  self.Fhw * self.hw_function(Tother=Tother) +
                  self.ztor_function(Tother=Tother) +
                  self.soil_function(Z10=Z10Rock,
                                     Vs30=Vs30Rock,
                                     Tother=Tother))
        output = np.exp(SA1100)
        return output

    def RegionalCorrection(self, Vs30=None, Rrup=None, Tother=None):
        if Tother != None:
            Ti = GetKey(Tother)
            T = Tother
        else:
            Ti = GetKey(self.T)
            T = self.T
        if Vs30 == None:
            Vs30 = self.Vs30
        if Rrup == None:
            Rrup = self.Rrup

        for key in ['VLIN', 'a31', 'a28', 'a29', 'a36', 'a37',
                    'a38', 'a40', 'a41', 'a42']:
            cmd = "%s = self.Coefs['%s']['%s']" % (key, Ti, key)
            exec(cmd)

        if self.region == 'CA':
            return 0.0
        elif self.region == 'TW':
            f11 = a31 * np.log(Vs30 / Vlin)
            return f11 + a25 * Rrup
        elif self.region == 'CN':
            return a28 * Rrup
        elif self.region == 'JP':
            f12 = (a36 * (Vs30 < 200) +
                   a37 * (200 <= Vs30 < 300) +
                   a38 * (300 <= Vs30 < 400) +
                   self.a39 * (400 <= Vs30 < 500) +
                   a40 * (500 <= Vs30 < 700) +
                   a41 * (700 <= Vs30 < 1000) +
                   a42 * (Vs30 >= 1000))
            return f12 + a29 * Rrup

    # function to compute the intensity
    def logline(self, x1, x2, y1, y2, x):
        # linear interpolation
        k = (y2 - y1) / (x2 - x1)
        C = y1 - k * x1
        y = k * x + C
        return y

    def compute_im(self, terms=(1, 1, 1, 1, 1, 1, 1)):

        # print 'Compute SA1100_calc'
        SA1100 = self.SA1100_calc()
        #print 'SA1100=',SA1100
        # print '===================================='
        LnSa = (terms[0] * self.base_model() +
                terms[1] * self.flt_function() +
                terms[2] * (self.Fhw*self.hw_function() +
                            terms[3] * self.ztor_function()) +
                terms[4] * self.site_model(SA1100) +
                terms[5] * self.soil_function() +
                terms[6] * self.RegionalCorrection())
        IM = np.exp(LnSa)
        # print 'IM = ',IM
        return IM

    # compute standard deviations
    def calc_alpha(self):
        Ti = GetKey(self.T)
        Vlin = self.Coefs[Ti]['VLIN']
        b = self.Coefs[Ti]['b']

        SA1100 = self.SA1100_calc()
        if self.Vs30 >= Vlin:
            alpha = 0.0
        else:
            alpha = (-b * SA1100 / (SA1100 + self.c) +
                     b * SA1100 / (SA1100 + self.c *
                                   (self.Vs30 / Vlin)**self.n))
       # print 'Alpha=',alpha
        return alpha

    def calc_sigma_tau(self):
        Ti = GetKey(self.T)
        if self.country != 'Japan':
            if self.VsFlag == 0:
                s1 = self.Coefs[Ti]['s01']
                s2 = self.Coefs[Ti]['s02']
            if self.VsFlag == 1:
                s1 = self.Coefs[Ti]['s11']
                s2 = self.Coefs[Ti]['s12']
            phi_AL = (s1 * (self.M < 4) +
                      (s1 + (s2 - s1) / 2 * (self.M - 4)) *
                      (4 <= self.M <= 6) + s2 * (self.M > 6))
        else:
            s5 = self.Coefs[Ti]['s5']
            s6 = self.Coefs[Ti]['s6']
            phi_AL = (s5 * (self.Rrup < 30) +
                      (s5 + (s6 - s5) / 50. * (self.Rrup - 30)) *
                      (30 <= self.Rrup <= 80) + s6 * (self.Rrup > 80))

        s3 = self.Coefs[Ti]['s3']
        s4 = self.Coefs[Ti]['s4']
        tau_AL = (s3 * (self.M < 5) +
                  (s3 + (s4 - s3) / 2. * (self.M - 5)) *
                  (5 <= self.M < 7) + s4 * (self.M >= 7))

        phi_Amp = 0.4
        alpha = self.calc_alpha()
        phi_B = np.sqrt(phi_AL**2 - phi_Amp**2)
        sigma = np.sqrt((phi_B * (1 + alpha))**2 + phi_Amp**2)

        tau_B = tau_AL
        tau = tau_B * (1 + alpha)

        sigmaT = np.sqrt(sigma**2 + tau**2)
        return (sigma, tau, sigmaT)

def ASK14nga_test(T, CoefTerms):
    """
    Test AS nga model
    """
    Mw = 8.0
    Zhypo = 8.0
    Ztor = 0.0
    dip = 90
    Ftype = 'SS'
    rake = 0    # for specific rupture
    W = 100

    Rjb = 3.0
    Rrup = Rjb
    Rx = Rrup
    #Rrup = (W*np.sin(dip*np.pi/180.)+Ztor) * np.cos(dip*np.pi/180.)
    #Rx = W*np.cos(dip*np.pi/180.)

    #print "Rx", Rx
    #Rx = Rrup

    Vs30 = 748.0, 1200., 345., 160.
    Vs30 = 760.
    Z25 = Z10 = None

    Fas = 0
    VsFlag = 0

    ASKnga = ASK14_nga()

    kwds = {'Ftype':Ftype, 'Rrup':Rrup, 'Rx':Rx, 'dip':dip, 'Ztor':Ztor,
            'W':W, 'Z10':Z10, 'Fas':Fas, 'VsFlag':VsFlag,
            'CoefTerms':CoefTerms}
    values = mapfunc(ASKnga, Mw, Rjb, Vs30, T, rake, **kwds)
    print('Median, SigmaT, Tau, Sigma')
    for i in range(len(values)):
        print(values[i])

    return ASKnga

if __name__ == '__main__':
    if 0:
        # SA test
        #T = 0.1; NewCoefs = {'Vlin':500,'b':-1.024}
        #T = 0.1; NewCoefs = None
        #T = 0.1; NewCoefs = {'Vlin':1032.5,'b':-1.624}
        NewCoefs = None
        T = 1.0
        CoefTerms = {'terms':(1, 1, 1, 1, 1, 1, 1), 'NewCoefs':NewCoefs}

        print('AS SA at %s' % ('%3.2f' % T))
        AS14nga = AS14nga_test(T, CoefTerms)

        T = -1.0
        print('AS PGA:')
        CoefTerms = {'terms':(1, 1, 1, 1, 1, 1, 1), 'NewCoefs':None}
        AS14nga = ASK14nga_test(T, CoefTerms)
    else:
        # Notes: PGV for Vs30 = 760, Z10 = 24, the soil-depth function
        #should be the same as T=1.0
        # for T in [-1.0,-2.0, 0.01, 0.02,
        #0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.75,
        #1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 7.5, 10.0]:
        for T in [-1.0, ]:
            print('AS SA at %s' % ('%3.2f' % T))
            CoefTerms = {'terms':(1, 1, 1, 1, 1, 1, 1), 'NewCoefs':None}
            ASKnga = ASK14nga_test(T, CoefTerms)
