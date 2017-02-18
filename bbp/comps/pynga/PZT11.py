#!/usr/bin/env python
"""
Pezeshk et al. 2011 model
"""

# Import Python modules
import math

class PZT11:
    """ 
    CENA Group 1 GMPE - Pezeshk et al. 2011
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
        self.c1s = [1.5828, 1.6854, 0.9314, -0.4883, -2.3106, -5.4113, -6.9340]

        self.c2s = [0.2298, 0.2404, 0.3088, 0.6278, 1.0220, 1.6900, 1.9070]

        self.c3s = [-0.03847, -0.03578, -0.03844, -0.05654,
                     -0.07965, -0.11960, -0.12870]

        self.c4s = [-3.8325, -3.6129, -3.2926, -3.0304,
                     -2.9265, -2.8998, -3.0128]

        self.c5s = [0.3535, 0.3247, 0.3063, 0.2673, 0.2515, 0.2465, 0.2639]

        self.c6s = [0.3321, 0.2956, 0.7064, 0.5422, 0.4716, 0.3766, 0.3172]

        self.c7s = [-0.09165, -0.11800, -0.09521, -0.05347,
                     -0.04039, -0.02928, -0.02150]

        self.c8s = [-2.5517, -3.3320, -2.2090, -1.3516,
                     -1.0923, -0.9470, -0.8749]

        self.c9s = [0.18310, 0.19770, 0.14720, 0.08784,
                    0.06554, 0.05249, 0.04774]

        self.c10s = [-0.0004224, -0.0001113, -0.0009254, -0.0010450,
                      -0.0007853, -0.0004563, -0.0003025]

        self.c11s = [6.6521, 6.8113, 6.1621, 6.1905, 6.0263, 6.1234, 6.1355]

    def __call__(self, mag, Rrup, period):
	"""
	Compute IM for single period
	required inputs:
	mag, Rrup, period
	"""
        # Figure out period
        if period in self.periods:
            np = self.periods.index(period)
        else:
	    print 'period is not in periods list'
	    raise ValueError

        rr = math.sqrt(Rrup**2 + self.c11s[np]**2)
        R1 = math.log10(min(rr, 70))
        R2 = max(min(math.log10(rr/70.0), math.log10(140.0/70.0)), 0)
        R3 = max(math.log10(rr/140.0), 0)

        term = self.c1s[np] + self.c2s[np] * mag + self.c3s[np] * (mag**2)
        gterm = ((self.c4s[np] + self.c5s[np] * mag) * R1 +
                 (self.c6s[np] + self.c7s[np] * mag) * R2 +
                 (self.c8s[np] + self.c9s[np] * mag) * R3 + self.c10s[np] * rr)

        im = (term + gterm) * math.log(10)
        
	return math.exp(im)

if __name__ == '__main__':
    pass
