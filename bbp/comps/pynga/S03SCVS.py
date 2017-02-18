#!/usr/bin/env python
"""
Silva al. 2003 Single Corner Variable Stress Model
"""

# Import Python modules
import math

class S03SCVS:
    """ 
    CENA Group 1 GMPE - Silva et al. 2003's Single Corner Variable Stress
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
        self.c1s = [4.03930, 5.20890, 3.92885, 2.27495,
                    0.13162, -4.51914, -9.00041]
        self.c2s = [0.10412, 0.09698, 0.20331, 0.34400,
                    0.57890, 1.13220, 1.66899]
        self.c3s = [2.7, 2.8, 2.7, 2.6, 2.6, 2.6, 2.5]
        self.c4s = [-2.97465, -3.01742, -2.80630, -2.61448,
                     -2.45001, -2.16445, -1.86794]
        self.c5s = [0.19631, 0.19172, 0.17658, 0.16182,
                    0.14539, 0.11502, 0.08623]
        self.c6s = [-0.08874, -0.08150, -0.08961, -0.11211,
                     -0.16638, -0.29235, -0.37576]

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

        term = (self.c1s[np] + self.c2s[np] * mag +
                self.c6s[np] * ((mag - 6.0)**2))
        gterm = ((self.c4s[np] + self.c5s[np] * mag) *
                 math.log(Rjb + math.exp(self.c3s[np])))

        im = term + gterm
        
	return math.exp(im)

if __name__ == '__main__':
    pass
