#!/usr/bin/env python
"""
CB 08 NGA model
"""
from utils import *

class CB08_nga():
    """
    Class of NGA model of Campbell and Bozorgnia 2008
    """
    def __init__(self):
	"""
	Model initialization
	"""
        
	# ============
	# NGA models (parameters and coefficients)
	# ============
	# 0. period independent parameters
	self.c = 1.88
	self.n = 1.18

	# 1. List of periods with defined coefficients (PGA is -1; PGV is -2) 
	self.periods = [0.01, 0.02, 0.03, 0.05, 0.075, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 
		   0.50, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 7.5, 10.0, -1.0, -2.0 ] 

        # ===============================
	# period-dependent coefficients
	# ===============================
	c0s = [ -1.715, -1.68, -1.552, -1.209, -0.657, -0.314, -0.133, -0.486, -0.89, 
	    -1.171, -1.466, -2.569, -4.844, -6.406, -8.692, -9.701, -10.556, -11.212, 
	    -11.684, -12.505, -13.087, -1.715, 0.954]
	
	c1s = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.656, 0.972, 
	       1.196, 1.513, 1.6, 1.6, 1.6, 1.6, 1.6, 1.6, 0.5, 0.696] 

	c2s = [-0.53, -0.53, -0.53, -0.53, -0.53, -0.53, -0.53, -0.446, -0.362, -0.294, 
	   -0.186, -0.304, -0.578, -0.772, -1.046, -0.978, -0.638, -0.316, -0.07, 
	   -0.07, -0.07, -0.53, -0.309]

	c3s = [-0.262, -0.262, -0.262, -0.267, -0.302, -0.324, -0.339, -0.398, -0.458, 
	   -0.511, -0.592, -0.536, -0.406, -0.314, -0.185, -0.236, -0.491, -0.77, 
	   -0.986, -0.656, -0.422, -0.262, -0.019] 

	c4s = [-2.118, -2.123, -2.145, -2.199, -2.277, -2.318, -2.309, -2.22, -2.146, 
	       -2.095, -2.066, -2.041, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2.118, -2.016] 

	c5s = [0.17, 0.17, 0.17, 0.17, 0.17, 0.17, 0.17, 0.17, 0.17, 0.17, 0.17, 0.17, 
	   0.17, 0.17, 0.17, 0.17, 0.17, 0.17, 0.17, 0.17, 0.17, 0.17, 0.17] 

	c6s = [5.6, 5.6, 5.6, 5.74, 7.09, 8.05, 8.79, 7.6, 6.58, 6.04, 5.3, 4.73, 4, 
	   4, 4, 4, 4, 4, 4, 4, 4, 5.6, 4] 

	c7s = [0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 0.28, 
	   0.28, 0.255, 0.161, 0.094, 0, 0, 0, 0, 0, 0.28, 0.245] 

	c8s = [-0.12, -0.12, -0.12, -0.12, -0.12, -0.099, -0.048, -0.012, 0, 0, 0, 0, 
	   0, 0, 0, 0, 0, 0, 0, 0, 0, -0.12, 0] 

	c9s = [0.49, 0.49, 0.49, 0.49, 0.49, 0.49, 0.49, 0.49, 0.49, 0.49, 0.49, 0.49, 
	   0.49, 0.49, 0.49, 0.371, 0.154, 0, 0, 0, 0, 0.49, 0.358] 

	c10s = [1.058, 1.102, 1.174, 1.272, 1.438, 1.604, 1.928, 2.194, 2.351, 2.46, 
	   2.587, 2.544, 2.133, 1.571, 0.406, -0.456, -0.82, -0.82, -0.82, -0.82, 
	   -0.82, 1.058, 1.694] 

	c11s = [0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 
	  0.077, 0.15, 0.253, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.04, 0.092] 

	c12s = [0.61, 0.61, 0.61, 0.61, 0.61, 0.61, 0.61, 0.61, 0.7, 0.75, 0.85, 0.883, 
	  1, 1, 1, 1, 1, 1, 1, 1, 1, 0.61, 1] 

	k1s =  [865, 865, 908, 1054, 1086, 1032, 878, 748, 654, 587, 503, 457, 410, 400, 
	   400, 400, 400, 400, 400, 400, 400, 865, 400] 

	k2s = [-1.186, -1.219, -1.273, -1.346, -1.471, -1.624, -1.931, -2.188, -2.381, 
	   -2.518, -2.657, -2.669, -2.401, -1.955, -1.025, -0.299, 0, 0, 0, 0, 0, 
	   -1.186, -1.955] 

	k3s = [1.839, 1.84, 1.841, 1.843, 1.845, 1.847, 1.852, 1.856, 1.861, 1.865, 1.874, 
		1.883, 1.906, 1.929, 1.974, 2.019, 2.11, 2.2, 2.291, 2.517, 2.744, 1.839, 1.929]
	
	# aleatory uncertainty models
	# page 149 of CB 08 ES paper
	self.sigma_lnY  = [0.478, 0.480, 0.489, 0.510, 0.520, 0.531, 0.532, 0.534, 0.534,
		    0.544, 0.541, 0.550, 0.568, 0.568, 0.564, 0.571, 0.558, 0.576,
		    0.601, 0.628, 0.667, 0.478, 0.484, 0.667]

	# intra-event residual standard deviation
	self.tau_lnY =  [0.219, 0.219, 0.235, 0.258, 0.292, 0.286, 0.280, 0.249, 0.240,
		    0.215, 0.217, 0.214, 0.227, 0.255, 0.296, 0.296, 0.326, 0.297,
		    0.359, 0.428, 0.485, 0.219, 0.203 ]

	self.sigma_C = [0.166, 0.166, 0.165, 0.162, 0.158, 0.170, 0.180, 0.186, 0.191,
		    0.198, 0.206, 0.208, 0.221, 0.225, 0.222, 0.226, 0.229, 0.237,
		    0.237, 0.271, 0.290, 0.166, 0.190 ]

	self.rho =  [ 1.000, 0.999, 0.989, 0.963, 0.922, 0.898, 0.890, 0.871, 0.852,
		    0.831, 0.785, 0.735, 0.628, 0.534, 0.411, 0.331, 0.289, 0.261,
		    0.200, 0.174, 0.174, 1.000, 0.691 ]

	# Old Coefs (period match)
	self.Coefs = {}
	for i in xrange(len(self.periods)):
	    T1 = self.periods[i]
	    Tkey = GetKey(T1)
	    self.Coefs[Tkey] = {}
	    self.Coefs[Tkey]['c0']   = c0s[i]
	    self.Coefs[Tkey]['c1']   = c1s[i]
	    self.Coefs[Tkey]['c2']   = c2s[i]
	    self.Coefs[Tkey]['c3']   = c3s[i]
	    self.Coefs[Tkey]['c4']   = c4s[i]
	    self.Coefs[Tkey]['c5']   = c5s[i]
	    self.Coefs[Tkey]['c6']   = c6s[i]
	    self.Coefs[Tkey]['c7']   = c7s[i]
	    self.Coefs[Tkey]['c8']   = c8s[i]
	    self.Coefs[Tkey]['c9']   = c9s[i]
	    self.Coefs[Tkey]['c10']  = c10s[i]
	    self.Coefs[Tkey]['c11']  = c11s[i]
	    self.Coefs[Tkey]['c12']  = c12s[i]
	    self.Coefs[Tkey]['k1']   = k1s[i]
	    self.Coefs[Tkey]['k2']   = k2s[i]
	    self.Coefs[Tkey]['k3']   = k3s[i]
	self.CoefKeys = self.Coefs[self.Coefs.keys()[0]].keys()
	
    
    # Call to get the SA value
    def __call__(self,M,Rjb,Vs30,T, rake, Ftype=None, \
	    Rrup=None,dip=None,Ztor=None,Z25=None, \
	    W=None,Zhypo=None,azimuth=None,Fhw=0,\
	    Z10=None,Z15=None, Arb=0, \
	    CoefTerms={'terms':(1,1,1,1,1,1),'NewCoefs':None}):
	"""
	Call the class to compute median ground-motion intensity
	You have to call the function here to make the class rich
	"""
	# Those inputs have to be specified
        self.M = M    # moment magnitude
	self.Rjb = float(Rjb)   # Joyner-Boore distance (km)
	self.Vs30 = float(Vs30)  # time-averaged shear wave velocity over 30m subsurface depth (m/s)
	self.T = T   # select period (sec)
	
	self.rake = rake    # rake could be None then you have to give the W and dip

        terms = CoefTerms['terms']
	NewCoefs = CoefTerms['NewCoefs']

        # check inputs
	if T in self.periods:
	    self.T = T
	else:
	    print 'T is not in periods list, try to interpolate'
	    raise ValueError
	
	if self.M == None or self.M < 0:
	    print 'Moment magnitude must be a postive number'
	    raise ValueError
	if self.Rjb == None or self.Rjb < 0: 
	    print 'Joyner-Boore distance must be a non-negative number'
	    raise ValueError
	if self.Vs30 == None or self.Vs30 < 0: 
	    print 'Vs30 must be a positive number'
	    raise ValueError
	
	# Determine the Fault-related parameters (if necessary)
	if Ftype != None:
	    self.Fnm = 1*(Ftype == 'NM')
	    self.Frv = 1*(Ftype == 'RV')
	else: 
	    if rake == None or rake < -180 or rake > 180.:
		print 'rake angle should be within [-180,180]'
		raise ValueError
	    else: 
		self.Frv, self.Fnm = rake2ftype_CB( self.rake )
       
	if W == None:
	    if self.rake == None: 
		print 'you should give either the fault width W or the rake angle'
		raise ValueError
	    else:
		W = calc_W(self.M,self.rake)
	else: 
	    self.W = W 

	if dip == None:
	    if self.rake == None: 
		print 'you should give either the fault dip angle or the rake angle'
		raise ValueError
	    else:
		self.dip = calc_dip( self.rake )
	else:
	    self.dip = dip
	
	if Ztor == None:
	    if Zhypo == None:
		if self.rake == None: 
		    print 'you should give either the Ztor or the rake angle'
		    raise ValueError
		else:
		    Zhypo = calc_Zhypo( self.M, self.rake )
	    self.Ztor = calc_Ztor( W, self.dip, Zhypo )
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
	    if self.Rjb == 0:
		Fhw = 1
		azimuth = 90
            Rx = calc_Rx( self.Rjb,self.Ztor, W, self.dip, azimuth, Rrup=Rrup )
	    self.Rrup = calc_Rrup( Rx, self.Ztor, W, self.dip, azimuth, Rjb = self.Rjb )
	else:
	    self.Rrup = Rrup

        # Determine Site-Specific parameters
        if Z25 == None:
	    self.Z25 = calc_Z25(self.Vs30,Z1model='CY')

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
    def moment_function(self,M=None,Tother=None):
	"""
	Moment term
	"""
	if Tother != None:
	    Ti = GetKey( Tother )
	else:
	    Ti = GetKey( self.T )
	
	c0 = self.Coefs[Ti]['c0']
	c1 = self.Coefs[Ti]['c1']
	c2 = self.Coefs[Ti]['c2']
	c3 = self.Coefs[Ti]['c3']
	
        if M != None:
            self.M = M

	if self.M <= 5.5:
	    return c0 + c1 * self.M
	elif 5.5<self.M<=6.5:
	    return c0 + c1 * self.M + c2 * (self.M-5.5)
	else:
	    return c0 + c1 * self.M + c2 * (self.M-5.5) + c3*(self.M-6.5)


    def distance_function(self,M=None,Rrup=None,Tother=None):
	"""
	Distance term
	"""
	if Tother != None:
	    Ti = GetKey( Tother )
	else:
	    Ti = GetKey( self.T )
        if M != None:
            self.M = M
        if Rrup != None:
            self.Rrup = Rrup

	c4 = self.Coefs[Ti]['c4']
	c5 = self.Coefs[Ti]['c5']
	c6 = self.Coefs[Ti]['c6']
	
	Rtmp = np.sqrt( self.Rrup**2 + c6**2)
	return (c4+c5*self.M)*np.log(Rtmp)
    

    def fault_function(self,Tother=None):
	"""
	Fault mechanism term
	or style of the fault
	"""
	if Tother != None:
	    Ti = GetKey( Tother )
	else:
	    Ti = GetKey( self.T )

	c7 = self.Coefs[Ti]['c7']
	c8 = self.Coefs[Ti]['c8']
	
	if self.Ztor < 1:
	    f_fltz = self.Ztor
	else:
	    f_fltz = 1
	return c7*self.Frv*f_fltz+c8*self.Fnm


    def hw_function(self,Tother=None):
	"""
	Hanging Wall term
	"""
	if Tother != None:
	    Ti = GetKey( Tother )
	else:
	    Ti = GetKey( self.T )


	c9 = self.Coefs[Ti]['c9']
	
	if self.Rjb == 0:
	    f_hngr = 1
	elif self.Rjb > 0 and self.Ztor < 1:
	    f_hngr = (max(self.Rrup,np.sqrt( self.Rjb**2+1 ))-self.Rjb) / max(self.Rrup,np.sqrt(self.Rjb**2+1))
	elif self.Rjb >0 and self.Ztor >= 1:
	    f_hngr = (self.Rrup-self.Rjb)/self.Rrup
	else:
	    print 'Rjb should be larger or equal to 0'
	    raise ValueError
	
	if self.M <= 6.0:
	    f_hngm = 0
	elif 6.0<self.M<6.5:
	    f_hngm = 2*(self.M-6.0)
	else:
	    f_hngm = 1

	if self.Ztor >= 20:
	    f_hngz = 0
	elif 0 <= self.Ztor < 20:
	    f_hngz = (20-self.Ztor)/20
	else:
	    #print 'Ztor is less than 0'   # R code cannot handle this
	    f_hngz = 0

	if self.dip <= 70:
	    f_hngd = 1
	else:
	    f_hngd = (90-self.dip)/20
	
	return c9*f_hngr*f_hngm*f_hngz*f_hngd


    def basin_function(self,Tother=None,Z25=None):
	"""
	Basin-effect term
	"""
	if Tother != None:
	    Ti = GetKey( Tother )
	else:
	    Ti = GetKey( self.T )
	
	if Z25 != None: 
	    self.Z25 = Z25 

	c11 = self.Coefs[Ti]['c11']
	c12 = self.Coefs[Ti]['c12']
	k3  = self.Coefs[Ti]['k3']
	
	if self.Z25 < 1:
	    return c11 * (self.Z25-1)
	elif 1 <= self.Z25 <= 3:
	    return 0
	else:
	    return c12 * k3*np.exp(-0.75)*(1-np.exp(-0.25*(self.Z25-3)))


    def A1100_calc(self):
	Tother = -1.0
	A1100 =  np.exp( self.moment_function(Tother)+
			 self.distance_function(Tother)+
			 self.fault_function(Tother)+
			 self.hw_function(Tother)+
			 self.basin_function(Tother=Tother)+
	                 self.site_function(A1100=0,Vs30=1100.,Tother = Tother) )
        return A1100


    def site_function(self,A1100=None,Vs30=None,Tother=None):
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
	
	c10 = self.Coefs[Ti]['c10']
	k1 = self.Coefs[Ti]['k1']
	k2  = self.Coefs[Ti]['k2']

	if Vs30 < k1:
	    return c10 * np.log(Vs30/k1) + k2*(np.log(A1100+self.c*(Vs30/k1)**self.n)-np.log(A1100+self.c))
	elif k1 <= Vs30 < 1100.:
	    return (c10+k2*self.n)*np.log(Vs30 /k1)
	else:
	    return (c10+k2*self.n)*np.log(1100./k1)



    # Final function to compute Sa, PGA, PGV
    def compute_im(self,terms=(1,1,1,1,1,1)):
	"""
	Compute IM based on functional form of CB08 model
	"""
	IM =  np.exp(terms[0]*self.moment_function()+
		          terms[3]*self.distance_function()+
			  terms[1]*self.fault_function()+
			  terms[2]*self.hw_function()+
			  terms[5]*self.site_function()+
			  terms[4]*self.basin_function())

	if self.T <= 0.25: # and self.T != -1.0:
	    Tother=-1.0
	    # This is PGA itself
	    IM1 =  np.exp(terms[0]*self.moment_function(Tother)+
			  terms[3]*self.distance_function(Tother)+
			  terms[1]*self.fault_function(Tother)+
			  terms[2]*self.hw_function(Tother)+
			  terms[5]*self.site_function(Tother=Tother)+
			  terms[4]*self.basin_function(Tother=Tother))

	    if IM < IM1: 
		# This is for SA (not for PGA, since PGA is computed above)
		IM = IM1
	
	return IM

    # function used to compute standard deviation terms
    def alpha_calc( self, Vs30=None, Tother=None ):
	if Vs30 == None:
	    Vs30 = self.Vs30
	if Tother == None:
	    Ti = GetKey( self.T )
	else:
	    Ti = GetKey( Tother )
	
	k1 = self.Coefs[Ti]['k1']
	k2 = self.Coefs[Ti]['k2']

	A1100 = self.A1100_calc()

	# compute alpha
	if Vs30 < k1:
	    alpha = k2 * A1100 * (1./(A1100+self.c*(Vs30/k1)**self.n)-1./(A1100+self.c))
	else:
	    alpha = 0
	
        return alpha


    def sigma_calc( self, Vs30=None, Tother=None ):
	"""
	Intra-event residual standard deviation
	"""
	sigma_lnAF = 0.3
	sigma_lnYb = np.sqrt(self.sigma_lnY[(np.array(self.periods)==self.T).nonzero()[0]]**2-sigma_lnAF**2)
	sigma_lnAb = np.sqrt(self.sigma_lnY[(np.array(self.periods)==-1.0).nonzero()[0]]**2-sigma_lnAF**2)    # PGA
	alpha = self.alpha_calc()
	sigma = np.sqrt(sigma_lnYb**2+sigma_lnAF**2 + alpha**2*sigma_lnAb**2 + \
		2*alpha*self.rho[(np.array(self.periods)==self.T).nonzero()[0]]*sigma_lnYb*sigma_lnAb )    # Eqn (15) CB08 ES
	return sigma
    
    def sd_calc(self):
	# compute SD at interested period self.T 
	indT = (np.array(self.periods)==self.T).nonzero()[0]
	
	tau = self.tau_lnY[indT]
	sigma = self.sigma_calc()
	sigmaT = np.sqrt( sigma**2 + tau**2 )
	sigmaArb = np.sqrt( sigmaT**2 + self.sigma_C[indT]**2 )
        
	# standard deviations are in logarithm scale !!!
	return (sigma, tau, sigmaT, sigmaArb)

def CB08nga_test(T, CoefTerms):
    """
    Test CB nga model
    """
    M = 4.0
    Vs30 = 748.0,1200.,345.,
    Vs30 = 760.
    #Vs30 = c(748.0,1200.,345.,160.)

    Z25 = None
    Ztor = 3
    dip = 90
    
    Rjb = np.arange(1,200,5)
    Rrup = Rjb
    Ftype = 'SS'
    rake = 0
    Arb = 0
    W = 10.

    # How to use it
    CBnga = CB08_nga()
    kwds = {'Ftype':Ftype,'Z25':Z25,'Rrup':Rrup,'W':W,'Ztor':Ztor,'dip':dip,'Arb':Arb,'CoefTerms':CoefTerms}
    values = mapfunc( CBnga, M, Rjb, Vs30, T, rake,**kwds )

    for i in xrange( len(values) ):
	print Rrup[i], values[i]

    return CBnga

if __name__ == '__main__':
    T = 2.0; NewCoefs = {'c1':1.6,'c2':-0.978}
    T = 2.0; NewCoefs = {'c1':1.7,'c2':-0.648}
    T = 2.0; NewCoefs = None
    CoefTerms = {'terms':(1,1,1,1,1,1),'NewCoefs':NewCoefs}
    Ts = [2.0, 3.0, 4.0, 5.0, 7.5, 10.0]
    Ts = [0.3]
    for T in Ts:
        print 'CB SA at %s'%('%3.2f'%T)
        CBnga = CB08nga_test(T,CoefTerms)
    T = -1.0
    print 'CB PGA'
    CBnga = CB08nga_test(T,CoefTerms)

