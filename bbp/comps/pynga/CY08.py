#!/usr/bin/env python

from __future__ import division, print_function

import numpy as np

from pynga.utils import (mapfunc, GetKey, calc_W, calc_dip, calc_Z1,
                         calc_Ztor, calc_Zhypo, calc_Rx, calc_Rrup,
                         rake2ftype_CY)

class CY08_nga(object):
    """
    Class for Chiou and Youngs 2008 NGA model
    """
    def __init__(self):

        # Model initialization

        # 1. MODEL COEFFICIENTS
        # periods list ( -1: PGV, -2: PGA )
        self.periods = [-1, 0.01, 0.02, 0.03, 0.04, 0.05, 0.075,
                        0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50,
                        0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 7.5,
                        10.0, -2.0]

        # model coefficients
        c1s = [-1.2687, -1.2687, -1.2515, -1.1744, -1.0671, -0.9464,
               -0.7051, -0.5747, -0.5309, -0.6352, -0.7766, -0.9278,
               -1.2176, -1.4695, -1.9278, -2.2453, -2.7307, -3.1413,
               -3.7413, -4.1814, -4.5187, -5.1224, -5.5872, 2.2884]

        c1as = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,
                0.0999, 0.0997, 0.0991, 0.0936, 0.0766, 0.0022,
                -0.0591, -0.0931, -0.0982, -0.0994, -0.0999, -0.1,
                0.1094]

        c1bs = [-0.255, -0.255, -0.255, -0.255, -0.255, -0.255,
                -0.254, -0.253, -0.25, -0.2449, -0.2382, -0.2313,
                -0.2146, -0.1972, -0.162, -0.14, -0.1184, -0.11,
                -0.104, -0.102, -0.101, -0.101, -0.1, -0.0626]

        c2s = [1.06, 1.06, 1.06, 1.06, 1.06, 1.06, 1.06, 1.06, 1.06,
               1.06, 1.06, 1.06, 1.06, 1.06, 1.06, 1.06, 1.06, 1.06,
               1.06, 1.06, 1.06, 1.06, 1.06, 1.06]

        c3s = [3.45, 3.45, 3.45, 3.45, 3.45, 3.45, 3.45, 3.45, 3.45,
               3.45, 3.45, 3.45, 3.45, 3.45, 3.45, 3.45, 3.45, 3.45,
               3.45, 3.45, 3.45, 3.45, 3.45, 3.45]

        cns = [2.996, 2.996, 3.292, 3.514, 3.563, 3.547, 3.448, 3.312,
               3.044, 2.831, 2.658, 2.505, 2.261, 2.087, 1.812, 1.648,
               1.511, 1.47, 1.456, 1.465, 1.478, 1.498, 1.502, 1.648]

        cMs = [4.1840, 4.1840, 4.1879, 4.1556, 4.1226, 4.1011, 4.0860,
               4.1030, 4.1717, 4.2476, 4.3184, 4.3844, 4.4979, 4.5881,
               4.7571, 4.8820, 5.0697, 5.2173, 5.4385, 5.5977, 5.7276,
               5.9891, 6.1930, 4.2979]

        c4s = [-2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1,
               -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1, -2.1,
               -2.1, -2.1, -2.1, -2.1, -2.1, -2.1]

        c4as = [-0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5,
                -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5,
                -0.5, -0.5, -0.5, -0.5, -0.5, -0.5]

        cRBs = [50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50,
                50, 50, 50, 50, 50, 50, 50, 50, 50, 50]

        c5s = [6.16, 6.16, 6.158, 6.155, 6.1508, 6.1441, 6.12, 6.085,
               5.9871, 5.8699, 5.7547, 5.6527, 5.4997, 5.4029, 5.29,
               5.248, 5.2194, 5.2099, 5.204, 5.202, 5.201, 5.2,
               5.2, 5.17]

        c6s = [0.4893, 0.4893, 0.4892, 0.489, 0.4888, 0.4884, 0.4872,
               0.4854, 0.4808, 0.4755, 0.4706, 0.4665, 0.4607, 0.4571,
               0.4531, 0.4517, 0.4507, 0.4504, 0.4501, 0.4501, 0.45,
               0.45, 0.45, 0.4407]

        cHMs = [3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
                3, 3, 3, 3, 3, 3, 3, 3]

        c7s = [0.0512, 0.0512, 0.0512, 0.0511, 0.0508, 0.0504,
               0.0495, 0.0489, 0.0479, 0.0471, 0.0464, 0.0458,
               0.0445, 0.0429, 0.0387, 0.035, 0.028, 0.0213,
               0.0106, 0.0041, 0.001, 0, 0, 0.0207]

        c7as = [0.0860, 0.0860, 0.0860, 0.0860, 0.0860, 0.0860, 0.0860,
                0.0860, 0.0860, 0.0860, 0.0860, 0.0860, 0.0850, 0.0830,
                0.0690, 0.0450, 0.0134, 0.0040, 0.0010, 0, 0, 0,
                0, 0.0437]

        c9s = [0.79, 0.79, 0.8129, 0.8439, 0.874, 0.8996, 0.9442,
               0.9677, 0.966, 0.9334, 0.8946, 0.859, 0.8019, 0.7578,
               0.6788, 0.6196, 0.5101, 0.3917, 0.1244, 0.0086,
               0, 0, 0, 0.3079]

        c9as = [1.5005, 1.5005, 1.5028, 1.5071, 1.5138, 1.523, 1.5597,
                1.6104, 1.7549, 1.9157, 2.0709, 2.2005, 2.3886, 2.5,
                2.6224, 2.669, 2.6985, 2.7085, 2.7145, 2.7164, 2.7172,
                2.7177, 2.718, 2.669]

        c10s = [-0.3218, -0.3218, -0.3323, -0.3394, -0.3453, -0.3502,
                -0.3579, -0.3604, -0.3565, -0.3470, -0.3379, -0.3314,
                -0.3256, -0.3189, -0.2702, -0.2059, -0.0852, 0.0160,
                0.1876, 0.3378, 0.4579, 0.7514, 1.1856, -0.1166]

        # NOTE:  The "g" refers to "gamma."
        cg1s = [-0.00804, -0.00804, -0.00811, -0.00839, -0.00875,
                -0.00912, -0.00973, -0.00975, -0.00883, -0.00778,
                -0.00688, -0.00612, -0.00498, -0.00420, -0.00308,
                -0.00246, -0.00180, -0.00147, -0.00117, -0.00107,
                -0.00102, -0.00096, -0.00094, -0.00275]

        cg2s = [-0.00785, -0.00785, -0.00792, -0.00819, -0.00855,
                -0.00891, -0.00950, -0.00952, -0.00862, -0.00759,
                -0.00671, -0.00598, -0.00486, -0.00410, -0.00301,
                -0.00241, -0.00176, -0.00143, -0.00115, -0.00104,
                -0.00099, -0.00094, -0.00091, -0.00625]

        cg3s = [4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
                4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4]

        # NOTE:  The coefficient "f" means "phi."
        f1s = [-0.4417, -0.4417, -0.434, -0.4177, -0.4, -0.3903,
               -0.404, -0.4423, -0.5162, -0.5697, -0.6109, -0.6444,
               -0.6931, -0.7246, -0.7708, -0.799, -0.8382, -0.8663,
               -0.9032, -0.9231, -0.9222, -0.8346, -0.7332, -0.7861]

        f2s = [-0.1417, -0.1417, -0.1364, -0.1403, -0.1591, -0.1862,
               -0.2538, -0.2943, -0.3113, -0.2927, -0.2662, -0.2405,
               -0.1975, -0.1633, -0.1028, -0.0699, -0.0425, -0.0302,
               -0.0129, -0.0016, 0, 0, 0, -0.0699]

        f3s = [-0.00701, -0.00701, -0.007279, -0.007354, -0.006977,
               -0.006467, -0.005734, -0.005604, -0.005845, -0.006141,
               -0.006439, -0.006704, -0.007125, -0.007435, -0.00812,
               -0.008444, -0.007707, -0.004792, -0.001828, -0.001523,
               -0.00144, -0.001369, -0.001361, -0.008444]

        f4s = [0.102151, 0.102151, 0.10836, 0.119888, 0.133641,
               0.148927, 0.190596, 0.230662, 0.266468, 0.255253,
               0.231541, 0.207277, 0.165464, 0.133828, 0.085153,
               0.058595, 0.031787, 0.019716, 0.009643, 0.005379,
               0.003223, 0.001134, 0.000515, 5.41]

        f5s = [0.2289, 0.2289, 0.2289, 0.2289, 0.2289, 0.229, 0.2292,
               0.2297, 0.2326, 0.2386, 0.2497, 0.2674, 0.312, 0.361,
               0.4353, 0.4629, 0.4756, 0.4785, 0.4796, 0.4799, 0.4799,
               0.48, 0.48, 0.2899]

        f6s = [0.014996, 0.014996, 0.014996, 0.014996, 0.014996,
               0.014996, 0.014996, 0.014996, 0.014988, 0.014964,
               0.014881, 0.014639, 0.013493, 0.011133, 0.006739,
               0.005749, 0.005544, 0.005521, 0.005517, 0.005517,
               0.005517, 0.005517, 0.005517, 0.006718]

        f7s = [580, 580, 580, 580, 579.9, 579.9, 579.6, 579.2,
               577.2, 573.9, 568.5, 560.5, 540, 512.9, 441.9,
               391.8, 348.1, 332.5, 324.1, 321.7, 320.9,
               320.3, 320.1, 459]

        f8s = [0.07, 0.07, 0.0699, 0.0701, 0.0702, 0.0701, 0.0686,
               0.0646, 0.0494, -0.0019, -0.0479, -0.0756, -0.096,
               -0.0998, -0.0765, -0.0412, 0.014, 0.0544,
               0.1232, 0.1859, 0.2295, 0.266, 0.2682, 0.1138]

        # uncertainties ( standard deviation ) (table 4 in CY 08 ES)
        self.tau1 = np.array([0.3437, 0.3437, 0.3471, 0.3603, 0.3718,
                              0.3848, 0.3878, 0.3835, 0.3719, 0.3601,
                              0.3522, 0.3438, 0.3351, 0.3353, 0.3429,
                              0.3577, 0.3769, 0.4023, 0.4406, 0.4784,
                              0.5074, 0.5328, 0.5542, 0.2539])

        self.tau2 = np.array([0.2637, 0.2637, 0.2671, 0.2803, 0.2918,
                              0.3048, 0.3129, 0.3152, 0.3128, 0.3076,
                              0.3047, 0.3005, 0.2984, 0.3036, 0.3205,
                              0.3419, 0.3703, 0.4023, 0.4406, 0.4784,
                              0.5074, 0.5328, 0.5542, 0.2381])

        self.sigma1 = np.array([0.4458, 0.4458, 0.4458, 0.4535, 0.4589,
                                0.4630, 0.4702, 0.4747, 0.4798, 0.4816,
                                0.4815, 0.4801, 0.4758, 0.4710, 0.4621,
                                0.4581, 0.4493, 0.4459, 0.4433, 0.4424,
                                0.4420, 0.4416, 0.4414, 0.4496])

        self.sigma2 = np.array([0.3459, 0.3459, 0.3459, 0.3537, 0.3592,
                                0.3635, 0.3713, 0.3769, 0.3847, 0.3902,
                                0.3946, 0.3981, 0.4036, 0.4079, 0.4157,
                                0.4213, 0.4213, 0.4213, 0.4213, 0.4213,
                                0.4213, 0.4213, 0.4213, 0.3554])

        self.sigma3 = np.array([0.8000, 0.8000, 0.8000, 0.8000, 0.8000,
                                0.8000, 0.8000, 0.8000, 0.8000, 0.8000,
                                0.7999, 0.7997, 0.7988, 0.7966, 0.7792,
                                0.7504, 0.7136, 0.7035, 0.7006, 0.7001,
                                0.7000, 0.7000, 0.7000, 0.7504])

        self.sigma4 = np.array([0.0663, 0.0663, 0.0663, 0.0663, 0.0663,
                                0.0663, 0.0663, 0.0663, 0.0612, 0.0530,
                                0.0457, 0.0398, 0.0312, 0.0255, 0.0175,
                                0.0133, 0.0090, 0.0068, 0.0045, 0.0034,
                                0.0027, 0.0018, 0.0014, 0.0133])

        # period match up (Old Coefs)
        self.Coefs = {}
        for i in range(len(self.periods)):
            T1 = self.periods[i]
            Tkey = GetKey(T1)
            self.Coefs[Tkey] = {}
            self.Coefs[Tkey]['c1'] = c1s[i]
            self.Coefs[Tkey]['c1a'] = c1as[i]
            self.Coefs[Tkey]['c1b'] = c1bs[i]
            self.Coefs[Tkey]['c2'] = c2s[i]
            self.Coefs[Tkey]['c3'] = c3s[i]
            self.Coefs[Tkey]['cn'] = cns[i]
            self.Coefs[Tkey]['cM'] = cMs[i]
            self.Coefs[Tkey]['c4'] = c4s[i]
            self.Coefs[Tkey]['c4a'] = c4as[i]
            self.Coefs[Tkey]['cRB'] = cRBs[i]
            self.Coefs[Tkey]['c5'] = c5s[i]
            self.Coefs[Tkey]['c6'] = c6s[i]
            self.Coefs[Tkey]['cHM'] = cHMs[i]
            self.Coefs[Tkey]['c7'] = c7s[i]
            self.Coefs[Tkey]['c7a'] = c7as[i]
            self.Coefs[Tkey]['c9'] = c9s[i]
            self.Coefs[Tkey]['c9a'] = c9as[i]
            self.Coefs[Tkey]['c10'] = c10s[i]
            self.Coefs[Tkey]['cg1'] = cg1s[i]
            self.Coefs[Tkey]['cg2'] = cg2s[i]
            self.Coefs[Tkey]['cg3'] = cg3s[i]
            self.Coefs[Tkey]['f1'] = f1s[i]
            self.Coefs[Tkey]['f2'] = f2s[i]
            self.Coefs[Tkey]['f3'] = f3s[i]
            self.Coefs[Tkey]['f4'] = f4s[i]
            self.Coefs[Tkey]['f5'] = f5s[i]
            self.Coefs[Tkey]['f6'] = f6s[i]
            self.Coefs[Tkey]['f7'] = f7s[i]
            self.Coefs[Tkey]['f8'] = f8s[i]
        self.CoefKeys = self.Coefs[list(self.Coefs.keys())[0]].keys()

    # call the function
    def __call__(self, M, Rjb, Vs30, T, rake, Ftype=None,
                 Rrup=None, Rx=None, dip=None, Ztor=None, Z10=None,
                 W=None, Zhypo=None, azimuth=None, Fhw=None,
                 AS=0, VsFlag=1,
                 CoefTerms={'terms':(1, 1, 1, 1, 1, 1), 'NewCoefs':None}):

        # call the function
        self.M = M         # Moment Magnitude
        self.Rjb = Rjb     # Joyner-Boore distance (km)
        self.rake = rake   # rake angle
        self.Vs30 = Vs30   # site-condition (m/s)

        if T in self.periods:
            self.T = T
        else:
            print('T is not in periods list, try to interpolate')
            raise ValueError

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
                print('you should give either the fault width W or the rake angle')
                raise ValueError
            else:
                W = calc_W(self.M, self.rake)
        else:
            self.W = W

        if dip == None:
            if self.rake == None:
                print('you should give either the fault dip angle or the rake angle')
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

        if Fhw == None:
            if azimuth == None and Rx == None:
                print('either one of azimuth angle, Rx and Fhw has to be specified')
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

        if self.Rjb == 0:
            azimuth = 90.
            Fhw = 1
        self.Fhw = Fhw

        # Compute Rx and Rrup
        if Rx == None:
            self.Rx = calc_Rx(self.Rjb, self.Ztor, W, self.dip, azimuth, Rrup)
        else:
            self.Rx = Rx
        if Rrup == None:
            self.Rrup = calc_Rrup(self.Rx, self.Ztor, W,
                                  self.dip, azimuth, self.Rjb)
        else:
            self.Rrup = Rrup

        # Z10
        if Z10 == None:
            self.Z10 = calc_Z1(self.Vs30, 'CY') # in meters
        else:
            self.Z10 = Z10 # Z10 should be in meter

        self.AS = AS # Aftershock flag (0 or 1) (depends on the earthquake itself)
        self.VsFlag = VsFlag # 0: inferred Vs30; 1: measured Vs30

        # update coeficient
        if NewCoefs != None:
            NewCoefKeys = NewCoefs.keys()
            Tkey = GetKey(self.T)
            for key in NewCoefKeys:
                self.Coefs[Tkey][key] = NewCoefs[key]

        IM = self.compute_im()
        sigma, tau, sigmaT = self.calc_sigma_tau()

        return IM, sigmaT, tau, sigma

    def flt_function(self):
        Ti = GetKey(self.T)

        c1 = self.Coefs[Ti]['c1']
        c1a = self.Coefs[Ti]['c1a']
        c1b = self.Coefs[Ti]['c1b']
        c7 = self.Coefs[Ti]['c7']
        c7a = self.Coefs[Ti]['c7a']
        c10 = self.Coefs[Ti]['c10']

        term0 = (c1 +
                 (c1a * self.Frv + c1b * self.Fnm + c7 * (self.Ztor - 4)) *
                 (1 - self.AS))
        term1 = (c10 + c7a * (self.Ztor - 4)) * self.AS
        return term0 + term1

    def moment_function(self):
        Ti = GetKey(self.T)

        c2 = self.Coefs[Ti]['c2']
        c3 = self.Coefs[Ti]['c3']
        cn = self.Coefs[Ti]['cn']
        cM = self.Coefs[Ti]['cM']

        term2 = (c2 * (self.M - 6) +
                 (c2 - c3) / cn * np.log(1 + np.exp(cn * (cM - self.M))))
        return term2

    def distance_function(self):
        Ti = GetKey(self.T)

        c4 = self.Coefs[Ti]['c4']
        c5 = self.Coefs[Ti]['c5']
        c6 = self.Coefs[Ti]['c6']
        cHM = self.Coefs[Ti]['cHM']
        c4a = self.Coefs[Ti]['c4a']
        cRB = self.Coefs[Ti]['cRB']
        cg1 = self.Coefs[Ti]['cg1']
        cg2 = self.Coefs[Ti]['cg2']
        cg3 = self.Coefs[Ti]['cg3']

        term3 = (c4 * np.log(self.Rrup +
                             c5 * np.cosh(c6 * max(self.M - cHM, 0))))
        term4 = (c4a - c4) * np.log(np.sqrt(self.Rrup**2 + cRB**2))
        term5 = (cg1 + cg2 / np.cosh(max(self.M - cg3, 0))) * self.Rrup
        return term3 + term4 + term5

    # make differences due the floating numbers computed using tanh
    def hw_function(self):
        Ti = GetKey(self.T)

        c9 = self.Coefs[Ti]['c9']
        c9a = self.Coefs[Ti]['c9a']

        d = self.dip * np.pi / 180.
        term6 = (c9 * self.Fhw *
                 np.tanh(self.Rx * np.cos(d) * np.cos(d) / c9a) *
                 (1 - np.sqrt(self.Rjb**2 + self.Ztor**2) /
                  (self.Rrup + 0.001)))
        return term6

    def lnYref(self):
        # equation 13a
        return (self.moment_function() +
                self.distance_function() +
                self.flt_function() +
                self.hw_function())

    def site_function(self):
        Ti = GetKey(self.T)

        f1 = self.Coefs[Ti]['f1']
        f2 = self.Coefs[Ti]['f2']
        f3 = self.Coefs[Ti]['f3']
        f4 = self.Coefs[Ti]['f4']

        lnY_ref = self.lnYref()
        term7 = f1 * min(np.log(self.Vs30 / 1130.), 0)
        term8 = (f2 * (np.exp(f3 * (min(self.Vs30, 1130) - 360)) -
                       np.exp(f3 * (1130 - 360))) *
                 np.log((np.exp(lnY_ref) + f4) / f4))
        return term7 + term8

    # make differences due the floating numbers computed using cosh
    def basin_function(self, Z10=None, Tother=None):

        if Tother != None:
            Ti = GetKey(Tother)
        else:
            Ti = GetKey(self.T)

        if Z10 != None:
            self.Z10 = Z10

        f5 = self.Coefs[Ti]['f5']
        f6 = self.Coefs[Ti]['f6']
        f7 = self.Coefs[Ti]['f7']
        f8 = self.Coefs[Ti]['f8']

        term9 = f5 * (1 - 1. / np.cosh(f6 * max(0, self.Z10 - f7)))
        # R use cosh(0.15*min(max(0, (Z1.0 - 15)), 300)) instead
        term10 = f8 / np.cosh(0.15 * max(0, self.Z10 - 15))
        # what R use
        #term10 = self.f8[Ti]/np.cosh(0.15*min(max(0, (selfZ10 - 15)), 300))
        return term9 + term10

    def compute_im(self, terms=(1, 1, 1, 1, 1, 1)):
        # use this one
        return np.exp(terms[0] * self.moment_function() +
                      terms[3] * self.distance_function() +
                      terms[1] * self.flt_function() +
                      terms[2] * self.hw_function() +
                      terms[4] * self.basin_function() +
                      terms[5] * self.site_function())

    def calc_NL(self):

        Ti = GetKey(self.T)

        f2 = self.Coefs[Ti]['f2']
        f3 = self.Coefs[Ti]['f3']
        f4 = self.Coefs[Ti]['f4']

        yref = np.exp(self.lnYref())

        b = f2 * (np.exp(f3 * (min(self.Vs30, 1130) - 360))
                  - np.exp(f3 * (1130 - 360))) # Eqn10
        c = f4

        return b * yref / (yref + c)

    def calc_sigma_tau(self):
        NL = self.calc_NL()
        if self.VsFlag == 0:
            Finfer = 1
            Fmeasure = 0
        else:
            Finfer = 0
            Fmeasure = 1

        indT = (np.array(self.periods) == self.T).nonzero()[0]

        # compute inter-event SD (Eqn 19)
        tau = (self.tau1[indT] +
               0.5 * (self.tau2[indT] - self.tau1[indT]) *
               (min(max(self.M, 5), 7) - 5))

        # compute intra-event SD (Eqn 20)
        sigma = ((self.sigma1[indT] +
                  0.5 * (self.sigma2[indT] - self.sigma1[indT]) *
                  (min(max(self.M, 5), 7) - 5) +
                  self.sigma4[indT] * self.AS) *
                 np.sqrt((self.sigma3[indT] * Finfer + 0.7 * Fmeasure) +
                         (1 + NL)**2))

        # correct tau
        tauNL = (1 + NL) * tau

        sigmaT = np.sqrt(sigma**2 + tauNL**2)

        #return (sigma, tau, sigmaT)
        return (sigma, tauNL, sigmaT)

def CY08nga_test(T, CoefTerms):
    """
    Test CY nga model
    """
    M = 7.75
    Rjb = 10.0

    #Vs30 = 748.0,1200.0, 356., 160.
    Vs30 = 865.0
    rake = 90    # for specific rupture

    W = 20

    Rrup = 21.0
    Rx = 20.0
    Ztor = 0.64274240
    dip = 45
    Z10 = None
    Z10 = 1000.0
    AS = 0
    VsFlag = 0

    CYnga = CY08_nga()

    kwds = {'Ztor':Ztor, 'dip':dip, 'Rrup':Rrup, 'Rx':Rx,
            'Z10':Z10, 'AS':AS, 'VsFlag':VsFlag, 'CoefTerms':CoefTerms}
    values = mapfunc(CYnga, M, Rjb, Vs30, T, rake, **kwds)
    print('Median, SigmaT, Tau, Sigma')
    for i in range(len(values)):
        print(values[i])
    return CYnga

if __name__ == '__main__':

    T = 0.1; NewCoefs = {'c1':-0.5747, 'c1a':0.1}
    T = 0.1; NewCoefs = {'c1':-0.6747, 'c1a':0.1}
    NewCoefs = None
    CoefTerms = {'terms':(1, 1, 1, 1, 1, 1), 'NewCoefs':NewCoefs}
    Z10 = 1000
    Ts = [3.0, 5.0, 10.0]
    for T in Ts:
        print('CY SA at %s' % ('%3.2f' % T))
        CYnga = CY08nga_test(T, CoefTerms)
        print(CYnga.basin_function(Z10=Z10, Tother=T))

    #T = -1.0
    #print 'CY PGA:'
    #CYnga = CY08nga_test(T)


