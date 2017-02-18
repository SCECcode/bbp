#!/usr/bin/env python
"""
BA 08 NGA model
"""
from utils import *
class BA08_nga:
    """ 
    Class of NGA model of Boore and Atkinson 2008
    """
    def __init__(self):
	"""
	Model initialization
	"""
	# 0. Given parameters (period independent parameters)
	self.a1 = 0.03     # in gravity (g)
	self.a2 = 0.09     # in gravity (g)
	self.pgalow = 0.06 # in gravity (g)
	self.V1 = 180.      # in m/s
	self.V2 = 300.      # in m/s
	self.Vref = 760.    # in m/s

	# 1. List of periods with defined coefficients (PGA is -1; PGV is -2) 
	self.periods = [ -2.0, -1.0, 0.01, 0.02, 0.03, 0.05, 0.075, 0.10, 0.15, 0.20, 0.25, 
		    0.30, 0.40, 0.50, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 7.5, 10.0 ]

        # ===============================
	# period-dependent coefficients
	# ===============================
	# 2. List of distance-scaling coefficients
	c1s = [ -0.87370, -0.66050, -0.66220, -0.66600, -0.69010, -0.71700, -0.72050, 
		-0.70810, -0.69610, -0.58300, -0.57260, -0.55430, -0.64430, -0.69140, 
		-0.74080, -0.81830, -0.83030, -0.82850, -0.78440, -0.68540, -0.50960, 
		-0.37240, -0.09824 ]
	 
	c2s = [ 0.10060, 0.11970, 0.12000, 0.12280, 0.12830, 0.13170, 0.12370, 0.11170, 
		0.09884, 0.04273, 0.02977, 0.01955, 0.04394, 0.06080, 0.07518, 0.10270, 
		0.09793, 0.09432, 0.07282, 0.03758, -0.02391, -0.06568, -0.13800 ] 
	 
	c3s = [ -0.00334, -0.01151, -0.01151, -0.01151, -0.01151, -0.01151, -0.01151, 
		-0.01151, -0.01113, -0.00952, -0.00837, -0.00750, -0.00626, -0.00540, 
		-0.00409, -0.00334, -0.00255, -0.00217, -0.00191, -0.00191, -0.00191, 
		-0.00191, -0.00191 ] 
	 
	hs = [ 2.54, 1.35, 1.35, 1.35, 1.35, 1.35, 1.55, 1.68, 1.86, 1.98, 2.07, 2.14, 
		2.24, 2.32, 2.46, 2.54, 2.66, 2.73, 2.83, 2.89, 2.93, 3.00, 3.04 ]     # in km

	e1s = [ 5.00121, -0.53804, -0.52883, -0.52192, -0.45285, -0.28476, 0.00767, 
		0.20109, 0.46128, 0.57180, 0.51884, 0.43825, 0.39220, 0.18957, -0.21338, 
		-0.46896, -0.86271, -1.22652, -1.82979, -2.24656, -1.28408, -1.43145, 
		-2.15446 ] 

	e2s = [ 5.04727, -0.50350, -0.49429, -0.48508, -0.41831, -0.25022, 0.04912, 
		0.23102, 0.48661, 0.59253, 0.53496, 0.44516, 0.40602, 0.19878, 
		-0.19496, -0.43443, -0.79593, -1.15514, -1.74690, -2.15906, -1.21270, 
		-1.31632, -2.16137 ]

	e3s = [ 4.63188, -0.75472, -0.74551, -0.73906, -0.66722, -0.48462, -0.20578, 
		0.03058, 0.30185, 0.4086, 0.3388, 0.25356, 0.21398, 0.00967, -0.49176, 
		-0.78465, -1.20902, -1.57697, -2.22584, -2.58228, -1.50904, -1.81022, 
		-2.53323 ] 

	e4s = [ 5.0821, -0.5097, -0.49966, -0.48895, -0.42229, -0.26092, 0.02706, 0.22193, 
		0.49328, 0.61472, 0.57747, 0.5199, 0.4608, 0.26337, -0.10813, -0.3933, 
		-0.88085, -1.27669, -1.91814, -2.38168, -1.41093, -1.59217, -2.14635 ]

	e5s = [ 0.18322, 0.28805, 0.28897, 0.25144, 0.17976, 0.06369, 0.0117, 0.04697, 
		0.1799, 0.52729, 0.6088, 0.64472, 0.7861, 0.76837, 0.75179, 0.6788, 
		0.70689, 0.77989, 0.77966, 1.24961, 0.14271, 0.52407, 0.40387 ] 

	e6s = [ -0.12736, -0.10164, -0.10019, -0.11006, -0.12858, -0.15752, -0.17051, 
		-0.15948, -0.14539, -0.12964, -0.13843, -0.15694, -0.07843, -0.09054, 
		-0.14053, -0.18257, -0.2595, -0.29657, -0.45384, -0.35874, -0.39006, 
		-0.37578, -0.48492 ]

	e7s = [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.00102, 0.08607, 0.10601, 0.02262, 0, 
		0.10302, 0.05393, 0.19082, 0.29888, 0.67466, 0.79508, 0, 0, 0 ]

	Mhs = [ 8.5, 6.75, 6.75, 6.75, 6.75, 6.75, 6.75, 6.75, 6.75, 6.75, 6.75, 6.75, 
		6.75, 6.75, 6.75, 6.75, 6.75, 6.75, 6.75, 6.75, 8.5, 8.5, 8.5 ] 

	self.Mref = 4.5     # moment magnitude
	self.Rref = 1.0     # in km 

	# 3. list of site-amplification coefficients (period-dependent)
	blins = [ -0.6, -0.36, -0.36, -0.34, -0.33, -0.29, -0.23, -0.25, -0.28, -0.31, 
		 -0.39, -0.44, -0.5, -0.6, -0.69, -0.7, -0.72, -0.73, -0.74, -0.75, 
		 -0.75, -0.692, -0.65 ] 
	 
	b1s = [ -0.5, -0.64, -0.64, -0.63, -0.62, -0.64, -0.64, -0.6, -0.53, -0.52, 
		     -0.52, -0.52, -0.51, -0.5, -0.47, -0.44, -0.4, -0.38, -0.34, -0.31, 
		      -0.291, -0.247, -0.215 ]
	 
	b2s = [ -0.06, -0.14, -0.14, -0.12, -0.11, -0.11, -0.11, -0.13, -0.18, 
		 -0.19, -0.16, -0.14, -0.1, -0.06, 0 , 0, 0, 0, 0, 0, 0, 0, 0 ]
        
	# 4. list of Aleatory uncertainties
	# intra-event residual standard deviation
	self.sigma0 = [ 0.5  ,  0.502,  0.502,  0.502,  0.507,  0.516,  0.513,  0.52 ,
		       0.518,  0.523,  0.527,  0.546,  0.541,  0.555,  0.571,  0.573,
		       0.566,  0.58 ,  0.566,  0.583,  0.601,  0.626,  0.645]

	# inter-event residual standard deviation (when fault type is not specified)
	self.tau_U = [ 0.286,  0.265,  0.267,  0.267,  0.276,  0.286,  0.322,  0.313,
	               0.288,  0.283,  0.267,  0.272,  0.267,  0.265,  0.311,  0.318,
		       0.382,  0.398,  0.41 ,  0.394,  0.414,  0.465,  0.355]
	
	# inter-event residual standard deviation (when fault type is specified)
	self.tau_M = [ 0.256,  0.26 ,  0.262,  0.262,  0.274,  0.286,  0.32 ,  0.318,
		       0.29 ,  0.288,  0.267,  0.269,  0.267,  0.265,  0.299,  0.302,
		       0.373,  0.389,  0.401,  0.385,  0.437,  0.477,  0.477]
        
	self.sigma_TU = [ 0.576,0.566,0.569,0.569,0.578,0.589,0.606,0.608,
		          0.592,0.596,0.592,0.608,0.603,0.615,0.649,0.654,
			  0.684,0.702,0.7,0.702,0.73,0.781,0.735 ]

	self.sigma_TM = [ 0.56,  0.564,  0.566,  0.566,  0.576,  0.589,  0.606,  0.608,  
		          0.594,  0.596,  0.592,  0.608,  0.603,  0.615,  0.645,  0.647,  
			  0.679,  0.7,  0.695,  0.698,  0.744,  0.787,  0.801 ]
	
	# Old Coefs (period match)
	self.Coefs = {}
	for i in xrange(len(self.periods)):
	    T1 = self.periods[i]
	    Tkey = GetKey(T1)
	    self.Coefs[Tkey] = {}
	    self.Coefs[Tkey]['c1']   = c1s[i]
	    self.Coefs[Tkey]['c2']   = c2s[i]
	    self.Coefs[Tkey]['c3']   = c3s[i]
	    self.Coefs[Tkey]['h']    = hs[i]
	    self.Coefs[Tkey]['e1']   = e1s[i]
	    self.Coefs[Tkey]['e2']   = e2s[i]
	    self.Coefs[Tkey]['e3']   = e3s[i]
	    self.Coefs[Tkey]['e4']   = e4s[i]
	    self.Coefs[Tkey]['e5']   = e5s[i]
	    self.Coefs[Tkey]['e6']   = e6s[i]
	    self.Coefs[Tkey]['e7']   = e7s[i]
	    self.Coefs[Tkey]['Mh']   = Mhs[i]
	    self.Coefs[Tkey]['blin'] = blins[i]
	    self.Coefs[Tkey]['b1']   = b1s[i]
	    self.Coefs[Tkey]['b2']   = b2s[i]
	self.CoefKeys = self.Coefs[self.Coefs.keys()[0]].keys()

        self.fault = ['unspecified','strike-slip','normal','reverse','U','NM','SS','RV']


    def __call__( self,M,Rjb,Vs30,T,rake, Mech=3, Ftype=None, AB11=None,CoefTerms={'terms':(1,1,1),'NewCoefs':None}):
	"""
	Compute IM for single period
	required inputs:
	M, Rjb, Vs30, T
	rake: rake angle (degree), default is None (Unspecified fault type)
	or give Mech instead of rake
	Mech: 
	     0: strike
	     1: normal
	     2: reverse
	     else: unspecified (U=1) (Default)
	Ftype = 'U', or 'SS', or 'RV', or 'NM'
	AB11: consider the recent correction to the median value
	"""
	# ==================
	# Input variables
	# ==================
	self.M = float(M)	     # moment magnitude
	self.Rjb = float(Rjb)	     # Joyner-Boore distance (km)
        self.Vs30 = float( Vs30 )    # 30 meter averaged S wave velocity (m/s)

        terms = CoefTerms['terms']
	NewCoefs = CoefTerms['NewCoefs']

	if T in self.periods:
	    self.T = T
	else:
	    print 'T is not in periods list, try to interpolate'
	    raise ValueError
	
	# check inputs
	if self.M == None or self.M < 0:
	    print 'Moment magnitude must be a postive number'
	    raise ValueError
	if self.Rjb == None or self.Rjb < 0:
	    print 'Joyner-Boore distance must be a non-negative number'
	    raise ValueError
	if self.Vs30 == None or self.Vs30 < 0:
	    print 'Vs30 must be a positive number'
	    raise ValueError

	self.rake = rake
	self.Mech = Mech
	
	if rake == None and Mech == None and Ftype == None:
	    print 'either rake or (U,SS,NM,RV) should be provided'
	    raise ValueError
	else: 
	    if Ftype != None: 
		self.U = 1*(Ftype == 'U')
		self.SS = 1*(Ftype == 'SS')
		self.NM = 1*(Ftype == 'NM')
		self.RV = 1*(Ftype == 'RV')
            else: 
		if Mech != None and rake != None:
		    # giveng Mech and rake at the same time, use Mech, not rake
		    rake = None

		if rake != None and Mech == None:
		    # Get ftype from rake
		    self.rake = rake
		    self.ftype()
		
		if rake == None and Mech != None:
		    self.U = 1*(Mech>2)
		    self.SS = 1*(Mech==0)
		    self.NM = 1*(Mech==1)
		    self.RV = 1*(Mech==2)

	self.AB11 = AB11

	# modify the coefficients
	if NewCoefs != None:
	    # only update Coefs given by NewCoefs (at self.T)
	    Tkey = GetKey( self.T )
	    NewCoefKeys = NewCoefs.keys()
	    for key in NewCoefKeys:
		self.Coefs[Tkey][key] = NewCoefs[key]
	
	# ======================
	# begin to compute IM
	# ======================
        IM = self.compute_im(terms=terms) 
        sigmaT, tau, sigma = self.compute_std()
        
	return IM, sigmaT, tau, sigma
                         
    # ============================
    # Functions used in the class
    # they could also be output for 
    # further regression analysis
    # ============================
    def ftype(self):
	"""
	Fault-Type 
	"""
	
	FT = rake2ftype_BA( self.rake )

	if FT not in self.fault:
	    print 'Invalid fault type!'
	    print 'It should be in one of the following list:'
	    print self.fault
	    raise ValueError
	else:
	    if FT == 'unspecified' or FT == 'U':
		self.U = 1
	    else:
		self.U = 0
	    
	    if FT == 'strike-slip' or FT == 'SS':
		self.SS = 1
	    else:
		self.SS = 0
	    
	    if FT == 'normal' or FT == 'NM':
		self.NM = 1
	    else:
		self.NM = 0

	    if FT == 'reverse' or FT == 'RV':
		self.RV = 1
	    else:
		self.RV = 0
        return FT   


    def moment_function(self, Tother=None):
	"""
	Magnitude-Moment scaling
	"""
	if Tother != None:
	    Ti = GetKey(Tother)
	else:
	    Ti = GetKey(self.T)

        e1 = self.Coefs[Ti]['e1']
        e2 = self.Coefs[Ti]['e2']
        e3 = self.Coefs[Ti]['e3']
        e4 = self.Coefs[Ti]['e4']
        e5 = self.Coefs[Ti]['e5']
        e6 = self.Coefs[Ti]['e6']
        e7 = self.Coefs[Ti]['e7']
        Mh = self.Coefs[Ti]['Mh']

	if self.M <= Mh:
	    return e1*self.U + e2*self.SS + e3*self.NM + e4*self.RV + \
		   e5*(self.M-Mh) + e6*(self.M-Mh)**2.
	else:
	    return e1*self.U + e2*self.SS + e3*self.NM + e4*self.RV + \
		   e7*(self.M-Mh)

    def distance_function(self,Tother=None):
	"""
	Distance function
	Geometrical spreading? (yes ~ ln(R))
	"""
	if Tother != None:
	    Ti = GetKey(Tother)
	else:
	    Ti = GetKey(self.T)
        
	h = self.Coefs[Ti]['h']
	c1 = self.Coefs[Ti]['c1']
	c2 = self.Coefs[Ti]['c2']
	c3 = self.Coefs[Ti]['c3']

	R = np.sqrt( self.Rjb**2 + h**2 )
	return (c1+c2*(self.M-self.Mref))*np.log(R/self.Rref)+c3*(R-self.Rref)


    def soil_function(self, Vs30=None, Tother=None):
	"""
	Site Amplification Function
	"""

	if Vs30 != None: 
	    self.Vs30 = Vs30 
	    
	if Tother != None: 
	    Ti = GetKey( Tother ) 
	else: 
	    Ti = GetKey(self.T )
	
	# linear term
	blin = self.Coefs[Ti]['blin']
	flin = blin * np.log(self.Vs30/self.Vref)

        # =================
	# non-linear term
	# =================
	# 1. compute pga4nl, which is defined as the media PGA when Vs30=Vref=760 m/s
	Tpga = -1.0    # compute PGA
	pga4nl = np.exp( self.moment_function(Tother=Tpga) + self.distance_function(Tother=Tpga) )
		        
	b1 = self.Coefs[Ti]['b1']
	b2 = self.Coefs[Ti]['b2']

	if self.Vs30 <= self.V1:
	    bnl = b1
	
	elif self.Vs30 > self.V1 and self.Vs30 <= self.V2:
	    bnl = (b1-b2)*np.log(self.Vs30/self.V2) / np.log(self.V1/self.V2) + b2 

	elif self.Vs30 > self.V2 and self.Vs30 < self.Vref:
	    bnl = b2*np.log( self.Vs30/self.Vref) / np.log(self.V2/self.Vref)

	else:
	    bnl = 0
	
	# 2. compute smoothing constants 
	dx = np.log( self.a2/self.a1 )
	dy = bnl*np.log(self.a2/self.pgalow)
	c = (3*dy-bnl*dx)/(dx**2)
	d = -(2*dy-bnl*dx)/(dx**3)

        # 3. final equation for nonlinear term
	if pga4nl <= self.a1:
	    fnl = bnl * np.log( self.pgalow/0.1 )
	elif pga4nl > self.a1 and pga4nl <= self.a2:
	    term = c*(np.log(pga4nl/self.a1))**2 + d * (np.log(pga4nl/self.a1))**3
	    fnl = bnl * np.log( self.pgalow/0.1) + term
	else:
	    fnl = bnl * np.log( pga4nl/0.1 )
	
	return flin+fnl


    def compute_im(self,terms=(1,1,1)):
        """
	Compute IM based on functional form of BA08 model
	"""
	IM =  np.exp(terms[0]*self.moment_function()+
	             terms[1]*self.distance_function()+
		     terms[2]*self.soil_function())
	if self.AB11 == None:
	    return IM

	else:
	    # BA 2011 correction for intermediate magnitudes
	    fba = max(0,3.888-0.674*self.M)-max(0,2.933-0.510*self.M)*np.log10(self.Rjb+10.)
	    fba = 10**fba
	    return fba * IM


    def compute_std(self):
	if self.rake == None:
	    if self.U == 1:
		FT = 'U'
	    if self.SS ==1:
		FT = 'SS'
	    if self.NM == 1:
		FT = 'NM'
	    if self.RV == 1:
		FT = 'RV'
	else:
	    FT = self.ftype()
	
	try: 
	    ind = (np.array( self.periods ) == self.T).nonzero()[0]
	    if FT == 'U':
		return (self.sigma_TU[ind], self.tau_U[ind], self.sigma0[ind])
	    else:
		return (self.sigma_TM[ind], self.tau_M[ind], self.sigma0[ind])
	except:
	    print 'inputed T not found in the available periods list, try to do interpolation'
	    raise ValueError


def BA08nga_test(T,CoefTerms):
    """
    Test BA features
    """
    # input parameter list
    Rjb = 200.
    Rjb = np.arange(1,200,5)
    Vs30 = 748.0,1200.,345.,160.
    Vs30 = 760.
    Mw = 4.0
    AB11 = None
    rake = 0
    Ftype = 'SS'
    kwds = {'Mech':None,'Ftype':Ftype,'AB11':AB11,'CoefTerms':CoefTerms}
    BAnga = BA08_nga()    # BA08nga instance
    values = mapfunc( BAnga, Mw, Rjb, Vs30, T, rake, **kwds )
    for ivalue in xrange( len(values) ):
	print Rjb[ivalue], values[ivalue]


if __name__ == '__main__':

    T = 10.0; NewCoefs = {'c1':-0.1,'c2':-0.14}  # use the updated one 
    T = 10.0; NewCoefs = {'c1':-0.09824,'c2':-0.13800} # use the updated one 
    T = 10.0; NewCoefs = {'c1':-0.1,'c2':-0.1000}  # use the updated one
    T = 0.3; NewCoefs = None     # pure one

    print 'BA SA at %s second'%('%3.2f'%T)
    CoefTerms={'terms':(1,1,1),'NewCoefs':NewCoefs}
    BAnga = BA08nga_test(T,CoefTerms)
    #BAnga = BA08nga_test(T,CoefTerms)

    T = -1.0
    CoefTerms={'terms':(1,1,1),'NewCoefs':None}
    print 'BA PGA at %s second'%('%3.2f'%T)
    BAnga = BA08nga_test(T,CoefTerms)

