#!/usr/bin/env python

from utils import *

# use can use the class to directly compute the IM and variances
# and you can use each individual functions in the class to do the regression ! 
#

class AS08_nga:
    """
    AS08 NGA model class
    """
    def __init__(self):
        
	# constant parameters
	self.c = 1.88
	self.c1 = 6.75
	self.c2 = 50
	self.c4 = 4.5
	self.n = 1.18

	# 1. MODEL COEFFICIENTS 
	self.periods =  [ -1.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.10, 0.15, 0.20, 0.25, 
		      0.30, 0.40, 0.50, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 7.5, 10.0, -2.0] 

	Vlins = [865.1, 865.1, 865.1, 907.8, 994.5, 1053.5, 1085.7, 1032.5, 877.6, 748.2, 
	   654.3, 587.1, 503.0, 456.6, 410.5, 400.0, 400.0, 400.0, 400.0, 400.0, 
	   400.0, 400.0, 400.0, 400.0] 

	bs = [ -1.186, -1.186, -1.219, -1.273, -1.308, -1.346, -1.471, -1.624, -1.931, 
	   -2.188, -2.381, -2.518, -2.657, -2.669, -2.401, -1.955, -1.025, -0.299, 
	   0.000, 0.000, 0.000, 0.000, 0.000, -1.955] 

	a1s = [ 0.804, 0.811, 0.855, 0.962, 1.037, 1.133, 1.375, 1.563, 1.716, 1.687, 
	   1.646, 1.601, 1.511, 1.397, 1.137, 0.915, 0.510, 0.192, -0.280, -0.639, 
	   -0.936, -1.527, -1.993, 5.7578] 

	a2s = [ -0.9679, -0.9679, -0.9774, -1.0024, -1.0289, -1.0508, -1.0810, -1.0833, 
	   -1.0357, -0.9700, -0.9202, -0.8974, -0.8677, -0.8475, -0.8206, -0.8088, 
	   -0.7995, -0.7960, -0.7960, -0.7960, -0.7960, -0.7960, -0.7960, -0.9046] 

	a3s = [ 0.265, 0.265, 0.265, 0.265, 0.265, 0.265, 0.265, 0.265, 0.265, 0.265, 
	   0.265, 0.265, 0.265, 0.265, 0.265, 0.265, 0.265, 0.265, 0.265, 0.265, 
	   0.265, 0.265, 0.265, 0.265] 

	a4s = [ -0.231, -0.231, -0.231, -0.231, -0.231, -0.231, -0.231, -0.231, -0.231, 
	   -0.231, -0.231, -0.231, -0.231, -0.231, -0.231, -0.231, -0.231, -0.231, 
	   -0.231, -0.231, -0.231, -0.231, -0.231, -0.231] 

	a5s = [-0.398, -0.398, -0.398, -0.398, -0.398, -0.398, -0.398, -0.398, -0.398, 
	   -0.398, -0.398, -0.398, -0.398, -0.398, -0.398, -0.398, -0.398, -0.398, 
	   -0.398, -0.398, -0.398, -0.398, -0.398, -0.398] 

	a8s = [ -0.0372, -0.0372, -0.0372, -0.0372, -0.0315, -0.0271, -0.0191, -0.0166, 
	   -0.0254, -0.0396, -0.0539, -0.0656, -0.0807, -0.0924, -0.1137, -0.1289, 
	   -0.1534, -0.1708, -0.1954, -0.2128, -0.2263, -0.2509, -0.2683, -0.1200] 

	a10s = [ 0.9445, 0.9445, 0.9834, 1.0471, 1.0884, 1.1333, 1.2808, 1.4613, 1.8071, 
	   2.0773, 2.2794, 2.4201, 2.5510, 2.5395, 2.1493, 1.5705, 0.3991, -0.6072, 
	   -0.9600, -0.9600, -0.9208, -0.7700, -0.6630, 1.5390]

	a12s =  [0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0181, 
	   0.0309, 0.0409, 0.0491, 0.0619, 0.0719, 0.0800, 0.0800, 0.0800, 0.0800, 
	   0.0800, 0.0800, 0.0800, 0.0800, 0.0800, 0.0800] 

	a13s = [ -0.0600, -0.0600, -0.0600, -0.0600, -0.0600, -0.0600, -0.0600, -0.0600, 
	   -0.0600, -0.0600, -0.0600, -0.0600, -0.0600, -0.0600, -0.0600, -0.0600, 
	   -0.0600, -0.0600, -0.0600, -0.0600, -0.0600, -0.0600, -0.0600, -0.0600] 

	a14s = [ 1.0800, 1.0800, 1.0800, 1.1331, 1.1708, 1.2000, 1.2000, 1.2000, 1.1683, 
	   1.1274, 1.0956, 1.0697, 1.0288, 0.9971, 0.9395, 0.8985, 0.8409, 0.8000, 
	   0.4793, 0.2518, 0.0754, 0.0000, 0.0000, 0.7000] 

	a15s = [ -0.3500, -0.3500, -0.3500, -0.3500, -0.3500, -0.3500, -0.3500, -0.3500, 
	   -0.3500, -0.3500, -0.3500, -0.3500, -0.3500, -0.3191, -0.2629, -0.2230, 
	   -0.1668, -0.1270, -0.0708, -0.0309, 0.0000, 0.0000, 0.0000, -0.3900] 

	a16s = [ 0.9000, 0.9000, 0.9000, 0.9000, 0.9000, 0.9000, 0.9000, 0.9000, 
	   0.9000, 0.9000, 0.9000, 0.9000, 0.8423, 0.7458, 0.5704, 0.4460, 0.2707, 
	   0.1463, -0.0291, -0.1535, -0.2500, -0.2500, -0.2500, 0.6300] 

	a18s = [-0.0067, -0.0067, -0.0067, -0.0067, -0.0067, -0.0076, -0.0093, -0.0093, 
		-0.0093, -0.0083, -0.0069, -0.0057, -0.0039, -0.0025, 0.0000, 0.0000, 
		0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000] 
	

	# SD (aleatory uncertainties)
        self.s1_est = [0.590, 0.590, 0.590, 0.605, 0.615, 0.623, 0.630, 0.630, 0.630, 0.630, 0.630,
		       0.630, 0.630, 0.630, 0.630, 0.630, 0.615, 0.604, 0.589, 0.578, 0.570,
		       0.611, 0.640, 0.590]                        

	self.s2_est = [0.470, 0.470, 0.470, 0.478, 0.483, 0.488, 0.495, 0.501, 0.509, 0.514, 0.518,
		       0.522, 0.527, 0.532, 0.539, 0.545, 0.552, 0.558, 0.565, 0.570, 0.587,
		       0.618, 0.640, 0.470]

	self.s1_meas = [0.576, 0.576, 0.576, 0.591, 0.602, 0.610, 0.617, 0.617, 0.616, 0.614, 0.612,
                       0.611, 0.608, 0.606, 0.602, 0.594, 0.566, 0.544, 0.527, 0.515, 0.510, 
		       0.572, 0.612, 0.576]

	self.s2_meas = [0.453, 0.453, 0.453, 0.461, 0.466, 0.471, 0.479, 0.485, 0.491, 0.495, 0.497,  
		        0.499, 0.501, 0.504, 0.506, 0.503, 0.497, 0.491, 0.500, 0.505, 0.529,0.579, 0.612, 0.453]
	
	self.s3 = [0.470, 0.420, 0.420, 0.462, 0.492, 0.515, 0.550, 0.550, 0.550, 0.520, 0.497,
		   0.479, 0.449, 0.426, 0.385, 0.350, 0.350, 0.350, 0.350, 0.350, 0.350,
	           0.350, 0.350, 0.420]

	self.s4 = [0.300, 0.300, 0.300, 0.305, 0.309, 0.312, 0.317, 0.321, 0.326, 0.329, 0.332,
		   0.335, 0.338, 0.341, 0.346, 0.350, 0.350, 0.350, 0.350, 0.350, 0.350,
		  0.350, 0.350, 0.300]

	self.rho = [1.000, 1.000, 1.000, 0.991, 0.982, 0.973, 0.952, 0.929, 0.896, 0.874, 0.856,
		    0.841, 0.818, 0.783, 0.680, 0.607, 0.504, 0.431, 0.328, 0.255, 0.200,
		    0.200, 0.200, 0.740]
        
	self.Coefs = {}
	for i in xrange( len(self.periods) ):
	    T1 = self.periods[i]
	    Tkey = GetKey( T1 )
	    self.Coefs[Tkey] = {}
		
	    self.Coefs[Tkey]['Vlin']  = Vlins[i] 
	    self.Coefs[Tkey]['b'] = bs[i]    
	    self.Coefs[Tkey]['a1'] = a1s[i]   
	    self.Coefs[Tkey]['a2']  = a2s[i]   
	    self.Coefs[Tkey]['a3']  = a3s[i]   
	    self.Coefs[Tkey]['a4']  = a4s[i]   
	    self.Coefs[Tkey]['a5']  = a5s[i]   
	    self.Coefs[Tkey]['a8']  = a8s[i]   
	    self.Coefs[Tkey]['a10'] = a10s[i]  
	    self.Coefs[Tkey]['a12'] = a12s[i]  
	    self.Coefs[Tkey]['a13']  = a13s[i]  
	    self.Coefs[Tkey]['a14']  = a14s[i]  
	    self.Coefs[Tkey]['a15'] = a15s[i]  
	    self.Coefs[Tkey]['a16']  = a16s[i]  
	    self.Coefs[Tkey]['a18'] = a18s[i]  
	self.CoefKeys = self.Coefs[self.Coefs.keys()[0]].keys()


    def __call__(self,M,Rjb,Vs30,T,rake, Ftype=None, \
	         Rrup=None,Rx=None,dip=None,Ztor=None,Z10=None,\
		 W=None,Zhypo=None,azimuth=None,Fhw=None, \
		 Fas=0,VsFlag=1, AS09=None,  \
		 CoefTerms={'terms':(1,1,1,1,1,1,1),'NewCoefs':None} \
		 ):
	
	self.M = M    # moment magnitude
	self.Rjb = Rjb     # JB distance (km)
	self.Vs30 = Vs30    # site condition (m/s)
	self.rake = rake   # rake andgle
	
	if T in self.periods:
	    self.T = T
	else:
	    print 'T is not in periods list, try to interpolate'
	    raise ValueError
	
	terms = CoefTerms['terms']
	NewCoefs = CoefTerms['NewCoefs']

        # Obtain optional parameters
	if Ftype != None:
	    self.Fnm = 1*(Ftype == 'NM')
	    self.Frv = 1*(Ftype == 'RV')
	else: 
	    if rake == None or rake < -180 or rake > 180.:
		print 'rake angle should be within [-180,180]'
		raise ValueError
	    else: 
		self.Frv, self.Fnm = rake2ftype_AS( self.rake )

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


	if Fhw == None:
	    if azimuth == None and Rx == None:
		print 'either one of azimuth angle, Rx and Fhw has to be specified'
		raise ValueError

	    if azimuth != None:
		if 0 <= azimuth <= 180. and dip != 90.:
		    Fhw = 1
		else:
		    Fhw = 0
	    
	    elif Rx != None:
		if Rx >=0 and dip != 90.:
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
	
	# Compute Rrup and Rx
	if Rx == None:                                                                
	    self.Rx = calc_Rx( self.Rjb, self.Ztor, self.W, self.dip, azimuth, Rrup )
	else:
	    self.Rx = Rx
	if Rrup == None:
	    self.Rrup = calc_Rrup( self.Rx, self.Ztor, self.W, self.dip, azimuth, self.Rjb )
	else:
	    self.Rrup = Rrup

	# Z10
	if Z10 == None:
	    self.Z10 = calc_Z1(self.Vs30,'AS')   # in meter
	else:
	    self.Z10 = Z10     # should be in meter

	self.Fas = Fas    # aftershock flag (0 or 1 )
        self.VsFlag = VsFlag    # 0: estimated Vs30; 1: measured Vs30
	self.AS09 = AS09    # 09 paper for the hanging wall

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

    
    def logline(self,x1,x2,y1,y2,x):
	# linear interpolation
	k = (y2-y1)/(x2-x1)
	C = y1-k*x1
	y = k*x+C
	return y


    def base_model(self,Tother=None):
	# Basically, this is the distance-magnitude term
	if Tother != None:
	    Ti = GetKey(Tother)
	else:
	    Ti = GetKey(self.T)

        a1 = self.Coefs[Ti]['a1']
        a2 = self.Coefs[Ti]['a2']
        a3 = self.Coefs[Ti]['a3']
        a4 = self.Coefs[Ti]['a4']
        a5 = self.Coefs[Ti]['a5']
        a8 = self.Coefs[Ti]['a8']

	Rtmp = np.sqrt(self.Rrup**2+self.c4**2)
	if self.M <= self.c1:
	    output = a1+ a4*(self.M-self.c1)+a8*(8.5-self.M)**2+ \
	          (a2+a3*(self.M-self.c1))*np.log(Rtmp)
	else:
	    output = a1+a5*(self.M-self.c1)+a8*(8.5-self.M)**2+ \
	          (a2+a3*(self.M-self.c1))*np.log(Rtmp)
	return output 


    def flt_function(self,Tother=None):
	# fault type and aftershock flag
	if Tother != None:
	    Ti = GetKey(Tother)
	else:
	    Ti = GetKey(self.T)
        
	a12 = self.Coefs[Ti]['a12']
	a13 = self.Coefs[Ti]['a13']
	a15 = self.Coefs[Ti]['a15']

        output = (a12*self.Frv + a13*self.Fnm + a15*self.Fas) 
	return output 

    def hw_function(self,Tother=None):
	# hanging wall function
	
	if self.Rx < 0:
	    output = 0.0

	else:

	    if Tother != None:
		Ti = GetKey( Tother )
	    else:
		Ti = GetKey( self.T )
            
	    a14 = self.Coefs[Ti]['a14']

	    if self.Rjb<30.:
		tape1 = 1-self.Rjb/30.
	    else:
		tape1 = 0

	    d = self.dip * np.pi/180.
	    if self.Rx <= self.W*np.cos(d):
		tape2 = 0.5+self.Rx/(2*self.W*np.cos(d))
	    if self.Rx > self.W*np.cos(d) or self.dip == 90:
		tape2 = 1.0
	    
	    if self.Rx >= self.Ztor:
		tape3 = 1
	    else:
		tape3 = self.Rx / self.Ztor
	    
	    if self.M <= 6:
		tape4 = 0.0
	    elif 6 < self.M < 7:
		tape4 = self.M - 6
	    else:
		tape4 = 1.0
	    
	    if self.AS09 == None:
		if self.dip >= 70:
		    tape5 = 1-(self.dip-70)/20   # AS 2008 Eq 12
		else:
		    tape5 = 1.0
	    else:
		# updated model
		if self.dip >=30: 
		    tape5 = 1-(self.dip-30)/60.  # AS 2009 Eq 5
		else:
		    tape5 = 1.0
	    output = a14*tape1*tape2*tape3*tape4*tape5	   
	return output 

    def ztor_function(self, Tother=None):
	# depth to top of rupture model

        if Tother != None:
	    Ti = GetKey( Tother )
	else:
	    Ti = GetKey( self.T )

	a16 = self.Coefs[Ti]['a16']

	if self.Ztor<10:
	    output = a16*self.Ztor/10.
	else:
	    output = a16
	return output


    def ld_function(self, Tother=None):
	# large distance function

        if Tother != None:
	    Ti = GetKey( Tother )
	else:
	    Ti = GetKey( self.T )
	
	a18 = self.Coefs[Ti]['a18']

	if self.M < 5.5:
	    tape6 = 1.0
	elif 5.5 <= self.M <= 6.5:
	    tape6 = 0.5*(6.5-self.M)+0.5
	else:
	    tape6 = 0.5
	
	if self.Rrup < 100:
	    output = 0.0 
	else:
	    output = a18*(self.Rrup-100)*tape6
        return output 


    # compute Vs30* for soil and site effect functions
    def CalcVs30Star( self, Vs30, T ):

	# compute V1 (used in soil-depth model)
	if T == -2.0:
	    V1 = 862     # For PGV (m/s)
	else: 
	    if T <= 0.5:
		V1 = 1500.    # m/s
	    elif 0.5 < T <= 1.0:
		V1 = np.exp(8.0-0.795*np.log(T/0.21))
	    elif 1 < T < 2:
		V1 = np.exp(6.76-0.297*np.log(T))
	    elif T >= 2:
		V1 = 700
	    else:
		pass
	
	# calculate Vs30*
	if Vs30 < V1:
	    Vs30_1 = Vs30
	else:
	    Vs30_1 = V1

	return V1, Vs30_1



    def site_model(self,PGA1100,Vs30=None,Tother=None):
	# Site-response model
	
	if Tother != None:
	    Ti = GetKey( Tother )
	    T = Tother
	else:
	    Ti = GetKey( self.T )
	    T = self.T
	
	a10 = self.Coefs[Ti]['a10']
	b = self.Coefs[Ti]['b']
	Vlin = self.Coefs[Ti]['Vlin']
	
	if Vs30 == None:
	    Vs30 = self.Vs30

	V1, Vs30_1 = self.CalcVs30Star(Vs30,T)
	
	if Vs30 < Vlin:
	    output = a10*np.log(Vs30_1/Vlin)-\
		    b*np.log(PGA1100+self.c) + \
		    b*np.log(PGA1100+self.c*(Vs30_1/Vlin)**self.n)
	else:
	    output = (a10+b*self.n)*np.log(Vs30_1/Vlin)
	return output 


    def soil_function(self,Z10=None,Vs30=None,Tother=None,coefs=None):
        # soil depth function
	if Tother != None:
	    Ti = GetKey( Tother )
	    T = Tother
	else:
	    Ti = GetKey( self.T )
	    T = self.T
	
	if coefs == None:
	    a10 = self.Coefs[Ti]['a10']
	    b = self.Coefs[Ti]['b']
	else: 
	    a10 = coefs['a10']
	    b = coefs['b']

	if Z10 == None:
	    Z10 = self.Z10

	if Vs30 == None:
	    Vs30 = self.Vs30
        
	V1, Vs30_1 = self.CalcVs30Star(Vs30,T)

	# compute \hat Z1.0
	if Vs30 < 180.:
	    Z11 = 6.745
	elif 180. <= Vs30 <= 500.:
	    Z11 = 6.745 - 1.35*np.log(Vs30/180.)
	else:
	    Z11 = 5.394 - 4.48*np.log(Vs30/500.)
	Z11 = np.exp(Z11)

        # compute e2
	if T == -2.0: 
	    # for PGV
	    e2 = -0.25 * np.log(Vs30/1000.)*np.log(1./0.35) 
	else: 
	    if T < 0.35 or Vs30 > 1000.:
		e2 = 0.0
	    elif 0.35 <= T <= 2.0:
		e2 = -0.25 * np.log(Vs30/1000.)*np.log(T/0.35)
	    else:
		e2 = -0.25 * np.log(Vs30/1000.)*np.log(2./0.35) 
	
	# compute a21
	tmp1 = (a10+b*self.n)*np.log(Vs30_1/min(V1,1000.))+e2*np.log((Z10+self.c2)/(Z11+self.c2))
	if Vs30 >= 1000.:
	    a21 = 0.0
	else:
	    if tmp1 < 0:    # what in R
		a21 = -(a10+b*self.n)*np.log(Vs30_1/min(V1,1000.)) / np.log((Z10+self.c2)/(Z11+self.c2))
	    else:
		a21 = e2

	# compute a22
	if T < 2:
	    # For PGV, a22 is computed using T=1.0
	    a22 = 0.0
	else:
	    a22 = 0.0625 * (T-2)
	
	term1 = a21*np.log((Z10+self.c2)/(Z11+self.c2))
	if Z10 >= 200:
	    output = term1 + a22*np.log(Z10/200.)
	else:
	    output = term1
        #print 'Soil Depth = ', Ti, output 

	return output 



    def PGA1100_calc(self):
	
	# compute PGA1100  
	# PGA when Vs30 = 1100
	PGA1100Rock = 0.0
	Vs30Rock = 1100.
	Z10Rock = calc_Z1( Vs30Rock, 'AS' )
	Tother = -1.0
	PGA1100 = self.base_model(Tother)+self.flt_function(Tother)+\
		  self.site_model(PGA1100Rock,Vs30=Vs30Rock,Tother=Tother)+\
		  self.Fhw*self.hw_function(Tother)+self.ztor_function(Tother)+self.ld_function(Tother)  + \
		  self.soil_function(Z10=Z10Rock,Vs30=Vs30Rock,Tother=Tother)
	output = np.exp( PGA1100 ) 
	return output 



    # function to compute the intensity
    def compute_im( self, terms=(1,1,1,1,1,1,1)):

	PGA1100 = self.PGA1100_calc()

        Td = -1.25+0.3*self.M
	Td = 10**(Td)
	Td = min(Td,10)

	if self.T <= Td:
	    LnSa = terms[0] * self.base_model() + \
		   terms[1] * self.flt_function() + \
		   terms[2] * (self.Fhw*self.hw_function() + self.ztor_function()) + \
		   terms[3] * self.ld_function() + \
		   terms[4] * self.site_model(PGA1100) + \
		   terms[5] * self.soil_function()
	else:

	    Vs30Rock = 1100
	    Z10Rock =  calc_Z1(Vs30Rock,'AS')
	    
	    periods1 = np.array( self.periods )
	    ind_low =  (periods1<Td).nonzero()[0]
	    ind_high = (periods1>Td).nonzero()[0]
	    Td1 = max( periods1[ind_low] ) 
	    Td2 = min( periods1[ind_high] ) 

	    Ti = GetKey( Td1 )
	    LnSa1 = terms[0] * self.base_model(Td1) + \
		    terms[1] * self.flt_function(Td1) + \
		    terms[2] * (self.Fhw*self.hw_function(Td1) + self.ztor_function(Td1)) + \
		    terms[3] * self.ld_function(Td1) + \
		    terms[4] * self.site_model(PGA1100,Vs30Rock,Td1)  + \
		    terms[5] * self.soil_function(Z10Rock,Vs30Rock,Td1)
	    
	    Ti = GetKey( Td2 )
	    LnSa2 = terms[0] * self.base_model(Td2) + \
		    terms[1] * self.flt_function(Td2) + \
		    terms[2] * (self.Fhw*self.hw_function(Td2) + self.ztor_function(Td2)) + \
		    terms[3] * self.ld_function(Td2) + \
		    terms[4] * self.site_model(PGA1100,Vs30Rock,Td2) + \
		    terms[5] * self.soil_function(Z10Rock,Vs30Rock,Td2)
	    
	    # The Rock Spectral acceleration is computed at T=Td
	    LnSaRockTd = self.logline(np.log(Td1),np.log(Td2),LnSa1,LnSa2,np.log(Td))    # log-linear interpolation
            
	    # the Rock Spectral acceleration at self.T is computed by scaling LnSaRockTd by (Td/T)**2
	    LnSaRockT = LnSaRockTd + np.log((Td/self.T)**2)  
	    
	    siteSoil = self.site_model(PGA1100) + self.soil_function()     # soil
	    siteRock = self.site_model(PGA1100,Vs30Rock) + self.soil_function(Z10Rock,Vs30Rock)    #  in R rock (For validation with R)
	    #siteRock = self.site_model(PGA1100,Vs30Rock)    
	    # what in Matlab and Java rock (for validation with Matlab)  (use this one)
	    # since when Vs30Rock = 1100, then soil_function = 0.0, so these two are equivalent !!!
	    soil_amp = siteSoil - siteRock  # both in R and Matlab
	    
	    LnSa = LnSaRockT + terms[6] * soil_amp
	    
	IM = np.exp(LnSa)
	return IM

    
    # compute standard deviations
    def alpha_calc(self):
	Ti = GetKey( self.T )
	Vlin = self.Coefs[Ti]['Vlin']
	b = self.Coefs[Ti]['b']

	PGA1100 = self.PGA1100_calc()
	if self.Vs30 >= Vlin:
	    return 0
	else:
	    return -b*PGA1100/(PGA1100+self.c)  + \
		    b*PGA1100/(PGA1100+self.c*(self.Vs30/Vlin)**self.n)    # as reported
     

    def calc_sigma0_tau0(self,Tother=None):
	if Tother == None:
	    indT = (np.array(self.periods)==self.T).nonzero()[0]
	else:
	    indT = (np.array(self.periods)==Tother).nonzero()[0]
	
	if self.VsFlag == 0:
	    s1 = self.s1_est[indT]
	    s2 = self.s2_est[indT]
	else:
	    s1 = self.s1_meas[indT]
	    s2 = self.s2_meas[indT]
        s3 = self.s3[indT]
	s4 = self.s4[indT]

	sigma0 = s1*(self.M<5) + (s1+0.5*(s2-s1)*(self.M-5))*(5<=self.M<=7) + s2*(self.M>7)   # Eqn 27
	tau0 = s3*(self.M<5) + (s3+0.5*(s4-s3)*(self.M-5))*(5<=self.M<=7) + s4*(self.M>7)     # Eqn 28
        return sigma0, tau0


    def calc_sigma_tau(self):
	sigma_amp = 0.3
	
	sigma0, tau0 = self.calc_sigma0_tau0()    # at T
	sigmaB_T = np.sqrt( sigma0**2-sigma_amp**2 )
	tauB_T = tau0
	
	Tother = -1.0      # for PGA
	sigma0, tau0 = self.calc_sigma0_tau0(Tother)
	sigmaB_PGA = np.sqrt( sigma0**2-sigma_amp**2 )
	tauB_PGA = tau0

	alpha = self.alpha_calc()
	indT = (np.array(self.periods)==self.T).nonzero()[0]
	
	# The paper of AS08 is incorrect, check the PEER report
	sigma = np.sqrt( sigmaB_T**2+sigma_amp**2+(alpha*sigmaB_PGA)**2 + \
	                 2.*alpha*sigmaB_T * sigmaB_PGA*self.rho[indT] )
	tau = np.sqrt( tauB_T**2 + (alpha*tauB_PGA)**2 + \
	                2.*alpha*tauB_T*tauB_PGA*self.rho[indT] )
	sigmaTot = np.sqrt( sigma**2 + tau**2 )
	
	return (sigma, tau, sigmaTot)


def AS08nga_test(T,CoefTerms):
    """
    Test AS nga model
    """
    if 0:
	# input list
	Mw = 7.75
	Rjb = 9.589290479428028 
	Vs30 = 354.0,768.,1200.,160. 
	rake = 180    # for specific rupture

	Rrup = 12.18438518124245 
	Rx = 10.099037378240684
	Ztor= 0.64274240
	dip = 79.39554
	W = 19.65842096
	Z10 = 0.0, 300., 500., 1000.

	Fas = 0
	VsFlag = 0
    else: 
	Mw = 6.93 
	Ztor = 3 
	Ftype = 'RV'; rake = None
	W = 3.85 
	dip = 70 
	Rrup = Rjb = Rx = 30 
	Fhw = 0 
	Vs30 = [860,300]
	Z10 = 0.024 * 1000   # in meter 
	Z25 = 2.974    # in km 
	VsFlag = 0 
	Fas = 0


    ASnga = AS08_nga()

    # Test of shallow soil depth site but different Vs30 (especially the small value in Vs30) to explain the smaller value in b(r) of AS-BA
    Vs30 = [150,200,300,400,500,600,800,1200]
    for vs in Vs30:
        print ASnga.soil_function(Z10=24., Vs30=vs, Tother=3.0)
    raw_input()

    kwds = {'Ftype':Ftype,'Rrup':Rrup,'Rx':Rx,'dip':dip,'Ztor':Ztor,'W':W,'Z10':Z10,'Fas':Fas,'VsFlag':VsFlag,'CoefTerms':CoefTerms}
    values = mapfunc( ASnga, Mw, Rjb, Vs30, T, rake, **kwds )
    print 'Median, SigmaT, Tau, Sigma'
    for i in xrange( len(values) ):
	print values[i]
    
    return ASnga



if __name__ == '__main__':
    if 0:
	# SA test 
	#T = 0.1; NewCoefs = {'Vlin':500,'b':-1.024}
	#T = 0.1; NewCoefs = None
	#T = 0.1; NewCoefs = {'Vlin':1032.5,'b':-1.624}
	NewCoefs = None
	T = 3.0
	CoefTerms = {'terms':(1,1,1,1,1,1,1), 'NewCoefs':NewCoefs}

	print 'AS SA at %s'%('%3.2f'%T)
	ASnga = AS08nga_test(T,CoefTerms)
	
	T = -1.0
	print 'AS PGA:'
	CoefTerms = {'terms':(1,1,1,1,1,1,1), 'NewCoefs':None}
	ASnga = AS08nga_test(T,CoefTerms)
    else: 
	# Notes: PGV for Vs30 = 760, Z10 = 24, the soil-depth function should be the same as T=1.0
	T = -2.0
	print 'AS PGV:'
	CoefTerms = {'terms':(1,1,1,1,1,1,1), 'NewCoefs':None}
	ASnga = AS08nga_test(T,CoefTerms)

