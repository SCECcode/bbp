#!/usr/bin/env python
"""
SC08 directivity model
"""
from utils import *

# =================================
# Correction of directivity term
# If from 4 NGA model, we get median Y0 (exp), and from SC directicity model, we get fD
# then the resulting ground motion is: fD*Y0 or exp(ln(fD)+ln(Y0)) 
# =================================

# compute required parameters for SC08 model (refer to matlab utilites)
# ...

class SC08_model:
    """
    Directivity correction to four NGA models
    # There is a Matlab tool to compute fD at different periods for different NGA models
    # don't need to use this one here, just as reference
    """
    def __init__(self,NGA,\
	    cuteps={'vr/vs':0.8,'c_cut':2.45,'s_cut':75,'r_cut':0.2,'d_cut':[40,70],'m_cut':[5.6,6.0]}):

	# when initial the class to get a instance, use A = NGA08_directivity(NGA)
	# then A(M,Rrup,IDP,T)
	# then use A.fD to get correct to the median from each NGA model
        # residual = a0 + (a+b*IDP) + eps_ij + eta_i
	# data-NGA = sigma_o + 
	
	self.nga = NGA   # NGA model name (AS08,BA08,CB08,CY08)
        
	# In IDP computing
	self.ratio = cuteps['vr/vs'] # 0.8
	self.ccut = cuteps['c_cut'] # 2.45
        self.scut = cuteps['s_cut'] # 75 km
        self.rcut = cuteps['r_cut'] # 0.2 

        # in fD computing
	self.dcut = cuteps['d_cut'] # [40,70] km  the other option is [200,250] km for large distance
	self.mcut = cuteps['m_cut'] # [5.6,6.0]
	
	#print self.dcut

	# MODEL COEFFICIENTS 
	if self.nga == 'AS08':
	    self.periods = [ 0.50, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 7.5, 10.0 ] 
	    a1s = [ 0.0,-0.0447,-0765,-0.1213,-0.1531,-0.1979,-0.2296,-0.2542,-0.3636,-0.5755 ]
	    bs = [ 0.0,0.0298,0.051,0.0809,0.102,0.1319,0.153,0.1695,0.2411,0.3489 ]
	    
	    a0s = [ 0.0421,-0.02137,0.02305,0.0338,0.04692,0.03191,0.06307,0.08158,-0.05407,-0.1517 ]
	    self.sigma = [ 0.5414,0.5586,0.5553,0.5174,0.524,0.5178,0.5294,0.5285,0.5147,0.5566 ]
	    self.tau = [ 0.2247,0.246,0.3138,0.326,0.3711,0.3595,0.3922,0.4192,0.4332,0.3656 ]
	    self.sigma0 = [ 0.5414,0.5596,0.5598,0.5225,0.5341,0.5393,0.5506,0.551,0.556,0.6626 ]
	    self.sratio = [ 0,0.002,0.008,0.01,0.019,0.04,0.039,0.041,0.074,0.16 ]
	    
	if self.nga == 'BA08':
	    self.periods = [0.50, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 7.5, 10.0 ] 
	    a1s = [0.0,-0.0532,-0.091,-0.1443,-0.1821,-0.2353,-0.2731,-0.3021,-0.4627,-0.8285 ]
	    bs = [0.0,0.0355,0.0607,0.0962,0.1214,0.1569,0.1821,0.2015,0.2727,0.4141 ]
	    
	    self.sigma = [ 0.5212,0.5387,0.5278,0.5052,0.5191,0.5197,0.5247,0.5513,0.534,0.5171 ]
	    self.tau = [ 0.1945,0.2618,0.3068,0.3214,0.3815,0.3985,0.3637,0.3801,0.4514,0.3387 ]
	    a0s = [ -0.00887,-0.003463,0.05466,0.0836,0.09181,0.03557,0.03217,0.02285,0.004121,-0.121 ]
	    self.sigma0 = [ 0.5212,0.5394,0.5327,0.5091,0.5301,0.5484,0.5559,0.5973,0.6005,0.6503 ]
	    self.sratio = [ 0,0.001,0.009,0.008,0.021,0.052,0.056,0.077,0.111,0.205 ]
	
	if self.nga == 'CB08':
	    self.periods = [0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 7.5, 10.0] 
	    a1s = [0.0,-0.0329,-0.0795,-0.1125,-0.159,-0.1921,-0.2172,-0.3227,-0.6419]
	    bs = [0.0,0.022,0.053,0.075,0.106,0.128,0.145,0.2147,0.3522]
	    
	    self.sigma = [ 0.5298,0.5234,0.4889,0.4921,0.4964,0.5075,0.5206,0.5151,0.5365]
	    self.tau = [ 0.2247,0.2666,0.2699,0.2507,0.2556,0.2323,0.2677,0.358,0.4071]
	    a0s = [ 0.02525,0.0589,0.08493,0.07915,0.05891,0.008163,-0.04461,-0.08447,-0.2105]
	    self.sigma0 = [  0.5298,0.5243,0.4899,0.4972,0.5129,0.525,0.5481,0.5613,0.6497]
	    self.sratio = [ 0,0.002,0.002,0.01,0.032,0.033,0.05,0.082,0.174 ]
	
	if self.nga == 'CY08':
	    self.periods = [0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 7.5, 10.0] 
	    a1s = [0.0,-0.026,-0.0627,-0.0887,-0.1254,-0.1514,-0.1715,-0.2797,-0.4847]
	    bs = [0.0,0.02,0.0482,0.0682,0.0965,0.1165,0.132,0.1865,0.2933]
	    
	    self.sigma = [ 0.5428,0.5393,0.5097,0.5307,0.5311,0.5503,0.5527,0.5476,0.5819]
	    self.tau = [ 0.3615,0.4042,0.3998,0.4044,0.4335,0.4274,0.4895,0.4693,0.3077]
	    a0s = [ 0.08438,0.0866,0.09405,0.1002,0.1064,0.1389,0.05773,0.05731,0.01719]
	    self.sigma0 = [ 0.5428,0.5404,0.5113,0.5349,0.5431,0.5626,0.5655,0.5713,0.6454]
	    self.sratio = [ 0,0.002,0.003,0.008,0.022,0.022,0.023,0.041,0.098 ]
	

	# define dictionary for coefficients selection
	self.Coefs = {}
	for i in xrange( len(self.periods) ):
	    T1 = self.periods[i]
	    Tkey = GetKey(T1)
	    self.Coefs[Tkey] = {}
	    self.Coefs[Tkey]['a'] = a1s[i]
	    self.Coefs[Tkey]['a0'] = a0s[i]
	    self.Coefs[Tkey]['b'] = bs[i]
	self.CoefKeys = self.Coefs[self.Coefs.keys()[0]].keys()

	self.fD=[]    # directivity correction   ( remember to plus the average residual a0 term)
	self.IDP = []

    def __call__(self,M,Rrup,ctildepr,s,h,Rfn,Rfp,T,NewCoefs=None):

	# call the function
	self.M = M        # Moment magnitude
	self.Rrup = Rrup  # rupture distance
	self.ctildepr = ctildepr    # c_wiggle_prime
	self.s = s      # along-strike distance between hypocenter and the closest point on the fault
	self.h = h      # downdip depth of hypocenter
	self.Rfn = Rfn      # radiation pattern (fault-normal)
	self.Rfp = Rfp      # fault parallel radiation pattern
        
	if T in self.periods:
	    self.T = T
	else:
	    print 'T is not in periods list, try to interpolate'
	    raise ValueError
	
	# change the value of ctildepr due to the different velocity ratio
	if self.ratio != 0.8:
	    tmp = 1./0.8 - 1./ctildepr    # (Rhypo-Rseis)/D  (doesn't depend on model parameters, just fault and site pair)
	    self.ctildepr = 1./ (1./self.ratio-tmp)    # definition of the c'
	else:
	    self.ctildepr = ctildepr

        self.IDP.append( self.calc_IDP() )

	# modify the coefficients
	if NewCoefs != None:
	    # only update Coefs given by NewCoefs (at self.T)
	    Tkey = GetKey( self.T )
	    NewCoefKeys = NewCoefs.keys()
	    for key in NewCoefKeys:
		self.Coefs[Tkey][key] = NewCoefs[key]
	
	if self.T >=0.5:
	    self.fD.append( self.calc_fD() )  
	else:
	    # For other periods will be 0.0 for fD
	    self.fD.append( 0.0 )


    def calc_IDP(self):

	C = ( min(self.ctildepr, self.ccut) - self.ratio ) / (self.ccut - self.ratio)
	S = np.log( min(self.scut,max(self.s,self.h)) )
	Rri = max( np.sqrt( self.Rfn**2 + self.Rfp**2 ), self.rcut )
	IDP = C * S * Rri
	return IDP


    def calc_fD(self):
	Ti = GetKey(self.T)
	
	a = self.Coefs[Ti]['a']
	b = self.Coefs[Ti]['b']
	a0 = self.Coefs[Ti]['a0']

	# distance and magnitude cut (model limitation)
	#fr = max(0,1-max(0,self.Rrup-self.dcut[0])/(self.dcut[1]-self.dcut[0]))
	#fM = min(1,max(0,self.M-self.mcut[0])/(self.mcut[1]-self.mcut[0]))
	
	# attention here:
	fr = 1
	fM = 1
        
	fD = fr*fM*(a+b*self.calc_IDP() )   # attention to the unit here
	
	return np.exp(fD + a0)

