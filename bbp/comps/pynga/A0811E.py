#!/usr/bin/env python
"""
Atkinson 2008-2011 EPRI
"""

# Import Python modules
import math

class A0811E:
    """
    CENA Group 1 GMPE - Atkinson 2008-2011 EPRI Version
    """
    def __init__(self):
	"""
	Model initialization
	"""
	# 1. List of periods with defined coefficients
	self.periods = [ 0.01, 0.04, 0.10, 0.20, 0.40, 1.0, 2.0 ]

        # ===============================
	# period-dependent coefficients
	# ===============================
        self.c1s = [4, 8, 11, 13, 16, 19, 21]
        self.c2s = [0, 0, 0, 0, 0, 0, 0]
        self.c3s = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]
        self.c4s = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]
        self.c5s = [0, 0, 0, 0, 0, 0, 0]
        self.c6s = [0, 0, 0, 0, 0, 0, 0]
        self.c7s = [3.888, 3.888, 3.888, 3.888, 3.888, 3.888, 3.888]
        self.c8s = [-0.674, -0.674, -0.674, -0.674, -0.674, -0.674, -0.674]
        self.c9s = [2.933, 2.933, 2.933, 2.933, 2.933, 2.933, 2.933]
        self.c10s = [-0.51, -0.51, -0.51, -0.51, -0.51, -0.51, -0.51]
        self.c11s = [0.4190, 0.4170, 0.2450, 0.0420,
                     -0.1350, -0.2127, -0.1800]
        self.c12s = [0.00211, 0.00192, 0.00273, 0.00232,
                     0.00184, 0.00125, 0.00090]
        self.c13s = [1.600, 1.700, 0.960, 0.740, 0.710, 0.775, 0.900]
        self.c14s = [-0.30, -0.15, 0.00, 0.00, 0.00, 0.00, 0.00]
        
    def __call__(self, mag, Rjb, period):
	"""
	Compute IM for single period
	required inputs:
	mag, Rjb, period
	"""
        # Figure out period
        if period in self.periods:
            np = self.periods.index(period)
        else:
	    print 'period is not in periods list'
	    raise ValueError

        # Get BA 2008 B/C values
        kt = self.c1s[np]
        vs30 = 760.0
        fU = self.c2s[np]
        fSS = self.c3s[np]
        fRV = self.c4s[np]
        fNM = self.c5s[np]

        elogbza = self.BA2008(kt, mag, Rjb, vs30, fU, fSS, fRV, fNM)

        # Make BA08 adjustment
        BA2008ac = max(0.0, (self.c7s[np] + self.c8s[np] * mag))
        BA2008bc = max(0.0, (self.c9s[np] + self.c10s[np] * mag))
        BA2008af = BA2008ac - BA2008bc * math.log10(Rjb + 10.0)

        # Compute west to east from Atkinson 2011
        wtoe = self.c11s[np] + self.c12s[np] * Rjb
        # Compute AB11 BC to hard rock amplification (Table 2, factors inverted)
        logamp = (math.log(self.c13s[np] +
                           self.c14s[np] * math.log10(max(Rjb, 1))))

        im = elogbza + (BA2008af + wtoe) * math.log(10.0) + logamp
        
	return math.exp(im)

    def BA2008(self, kt, mag, rjb, vs30, fU, fSS, fRV, fNM):
        """
        Boore and Atkinson 2008 NGA period independent
        """
        a1 = 0.03
        pga_low = 0.06
        a2 = 0.09
        V1 = 180.0
        V2 = 300.0
        Vref = 760.0
        Mref = 4.5
        Rref = 1.0

        # Adjust kt from R's 1-based arrays to Python's 0-based arrays
        kt = kt - 1

        periods = [-10.,-1.,0.,0.01,0.02,0.025,0.03,0.04,0.05,0.075,0.1,
                    0.15,0.2,0.25,0.3,0.4,0.5,0.75,1.0,1.5,2.,3.,4.,5.,7.5,10.]

        blin = [0.,-0.6,-0.36,-0.36,-0.34,-0.334,-0.33,-0.307,-0.29,-0.23,-0.25,
                -0.28,-0.31,-0.39,-0.44,-0.5,-0.6,-0.69,-0.7,-0.72,-0.73,-0.74,
                -0.75,-0.75,-0.692,-0.65]

        b1 = [0.,-0.5,-0.64,-0.64,-0.63,-0.624,-0.62,-0.631,-0.64,-0.64,
              -0.6,-0.53,-0.52,-0.52,-0.52,-0.51,-0.5,-0.47,-0.44,-0.4,
              -0.38,-0.34,-0.31,-0.291,-0.247,-0.215]

        b2 = [0.,-0.06,-0.14,-0.14,-0.12,-0.114,-0.11,-0.11,-0.11,-0.11,
              -0.13,-0.18,-0.19,-0.16,-0.14,-0.1,-0.06,0.,0.,0.,0.,0.,0.,
              0.,0.,0.]

        c1 = [-0.55,-0.8737,-0.6605,-0.6622,-0.666,-0.6792,-0.6901,
               -0.7051,-0.717,-0.7205,-0.7081,-0.6961,-0.583,-0.5726,
               -0.5543,-0.6443,-0.6914,-0.7408,-0.8183,-0.8303,-0.8285,
               -0.7844,-0.6854,-0.5096,-0.3724,-0.09824]

        c2 = [0.,0.1006,0.1197,0.12,0.1228,0.1258,0.1283,0.1302,0.1317,
              0.1237,0.1117,0.09884,0.04273,0.02977,0.01955,0.04394,
              0.0608,0.07518,0.1027,0.09793,0.09432,0.07282,0.03758,
              -0.02391,-0.06568,-0.138]

        c3 = [-0.01151,-0.00334,-0.01151,-0.01151,-0.01151,-0.01151,
               -0.01151,-0.01151,-0.01151,-0.01151,-0.01151,-0.01113,
               -0.00952,-0.00837,-0.0075,-0.00626,-0.0054,-0.00409,
               -0.00334,-0.00255,-0.00217,-0.00191,-0.00191,-0.00191,
               -0.00191,-0.00191]

        h = [3.,2.54,1.35,1.35,1.35,1.35,1.35,1.35,1.35,1.55,1.68,1.86,
             1.98,2.07,2.14,2.24,2.32,2.46,2.54,2.66,2.73,2.83,2.89,2.93,
             3,3.04]

        e1 = [-0.03279,5.00121,-0.53804,-0.52883,-0.52192,-0.48270,
               -0.45285,-0.34873,-0.28476,0.00767,0.20109,0.46128,
               0.5718,0.51884,0.43825,0.3922,0.18957,-0.21338,
               -0.46896,-0.86271,-1.22652,-1.82979,-2.24656,
               -1.28408,-1.43145,-2.15446]

        e2 = [-0.03279,5.04727,-0.5035,-0.49429,-0.48508,-0.44711,
               -0.41831,-0.31319,-0.25022,0.04912,0.23102,0.48661,
               0.59253,0.53496,0.44516,0.40602,0.19878,-0.19496,
               -0.43443,-0.79593,-1.15514,-1.7469,
               -2.15906,-1.2127,-1.31632,-2.16137]

        e3 = [-0.03279,4.63188,-0.75472,-0.74551,-0.73906,-0.69862,
               -0.66722,-0.55727,-0.48462,-0.20578,0.03058,0.30185,
               0.4086,0.3388,0.25356,0.21398,0.00967,-0.49176,
               -0.78465,-1.20902,-1.57697,-2.22584,
               -2.58228,-1.50904,-1.81022,-2.53323]

        e4 = [-0.03279,5.0821,-0.5097,-0.49966,-0.48895,-0.45106,
               -0.42229,-0.32200,-0.26092,0.02706,0.22193,0.49328,
               0.61472,0.57747,0.5199,0.4608,0.26337,-0.10813,
               -0.3933,-0.88085,-1.27669,-1.91814,-2.38168,
               -1.41093,-1.59217,-2.14635]

        e5 = [0.29795,0.18322,0.28805,0.28897,0.25144,0.20904,0.17976,
              0.10021,0.06369,0.0117,0.04697,0.1799,0.52729,0.6088,
              0.64472,0.7861,0.76837,0.75179,0.6788,0.70689,0.77989,
              0.77966,1.24961,0.14271,0.52407,0.40387]

        e6 = [-0.20341,-0.12736,-0.10164,-0.10019,-0.11006,-0.11990,
               -0.12858,-0.14415,-0.15752,-0.17051,-0.15948,-0.14539,
               -0.12964,-0.13843,-0.15694,-0.07843,-0.09054,-0.14053,
               -0.18257,-0.2595,-0.29657,-0.45384,
               -0.35874,-0.39006,-0.37578,-0.48492]

        e7 = [0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.00102,0.08607,
              0.10601,0.02262,0.,0.10302,0.05393,0.19082,0.29888,
              0.67466,0.79508,0.,0.,0.]

        Mh = [7.,8.5,6.75,6.75,6.75,6.75,6.75,6.75,6.75,6.75,6.75,6.75,
              6.75,6.75,6.75,6.75,6.75,6.75,6.75,6.75,6.75,6.75,6.75,
              8.5,8.5,8.5]

        # Compute reference PGA4nl
        # Modified to use pga for pga4nl per
        # final version of BA 2008
        knl = 2 # (-1 from the R code as R arrays are 1-based)

        # Compute FMr
        if mag <= Mh[knl]:
            FMr = (e1[knl] * fU + e2[knl] * fSS + e3[knl] * fNM +
                   e4[knl] * fRV + e5[knl] * (mag - Mh[knl]) +
                   e6[knl] * (mag - Mh[knl]) * (mag - Mh[knl]))
        else:
            FMr = (e1[knl] * fU + e2[knl] * fSS + e3[knl] * fNM +
                   e4[knl] * fRV + e7[knl] * (mag - Mh[knl]))

        # Compute FDr
        Rr = math.sqrt(rjb * rjb + h[knl] * h[knl])
        FDr = ((c1[knl] + c2[knl] * (mag - Mref)) * math.log(Rr / Rref) +
               c3[knl] * (Rr - Rref))
        pga4nl = math.exp(FMr + FDr)

        # Compute PSA
        # Compute FM
        if mag <= Mh[kt]:
            FM = (e1[kt] * fU + e2[kt] * fSS + e3[kt] * fNM + e4[kt] * fRV +
                  e5[kt] * (mag - Mh[kt]) +
                  e6[kt] * (mag - Mh[kt]) * (mag - Mh[kt]))
        else:
            FM = (e1[kt] * fU + e2[kt] * fSS + e3[kt] * fNM +
                  e4[kt] * fRV + e7[kt] * (mag - Mh[kt]))

        # Compute FD
        R = math.sqrt(rjb * rjb + h[kt] * h[kt])
        FD = ((c1[kt] + c2[kt] * (mag - Mref)) * math.log(R / Rref) +
              c3[kt] * (R - Rref))

        # Compute Flin
        Flin = blin[kt] * math.log(vs30 / Vref)

        # Compute Fnl
        if vs30 <= V1:
            bnl = b1[kt]
        elif vs30 <= V2:
            bnl = ((b1[kt] - b2[kt]) * math.log(vs30 / V2) /
                   math.log(V1 / V2) + b2[kt])
        elif vs30 <= Vref:
            bnl = b2[kt] * math.log(vs30 / Vref) / math.log(V2 / Vref)
        else:
            bnl = 0.0

        if pga4nl <= a1:
            Fnl = bnl * math.log(pga_low / 0.1)
        elif pga4nl <= a2:
            delta_x = math.log(a2 / a1)
            delta_y = bnl * math.log(a2 / pga_low)
            c = (3.0 * delta_y - bnl * delta_x) / delta_x / delta_x
            d = -(2.0 * delta_y - bnl * delta_x) / delta_x / delta_x / delta_x
            pr = math.log(pga4nl / a1)
            Fnl = bnl * math.log(pga_low / 0.1) + c * pr * pr + d * pr * pr * pr
        else:
            Fnl = bnl * math.log(pga4nl / 0.1)

        value = FM + FD + Flin + Fnl
        return value            

if __name__ == '__main__':
    pass
