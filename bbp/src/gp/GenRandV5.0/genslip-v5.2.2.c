#include "include.h"
#include "structure.h"
#include "function.h"
#include "defs.h"
#include "getpar.h"

void get_rupt(struct velmodel *vm, float *h,
	      float *srcd, float *recd, float *srcr,
	      float *recr, double *p, double *rad, float *tt);

/*
************************************************************************************

   03/06/2000: version 2.1

   * Added option to read input fault description from Generic Slip Format (GSF)
     file.

   * Added option to specify "Brune" type STF.

************************************************************************************

   07/08/2009: version 3.0

   * Added flag to make SRF output optional.

   * Added tsfac depth and slip scaling for initiation time.

************************************************************************************

   2010/04/09: version 3.1

   * Added option to adjust rupture speed at segment boundaries, see seg_delay

************************************************************************************

   2010/07/21: version 3.2

   * Added option to adjust rise time with sqrt(slip) [previously done using generic_slip2srf]
     This is done using the parameter stfparams.rt_scalefac=1 (=0 gives no slip scaling)

   * Set the following as defaults:

	 kmodel=2
         flip_at_surface=1
	 stretch_kcorner=0
	 circular_average=0
	 modified_corners=0
         truncate_zero_slip=1
         rand_rake_degs=60.0
	 slip_sigma=0.85

	 tsfac_coef=1.8
	 tsfac_factor=1
	 tsfac=-ts_coef*1.0e-09*exp(log(moment)/3.0)

	 rvfrac=0.8
	 shal_vrup=0.7

         shal_vrup_dep=6.5
         shal_vrup_deprange=1.5

         tsfac_dep=6.5
         tsfac_deprange=1.5

         stfparams.trise=1.6*1.0e-09*exp(log(moment)/3.0)
         stfparams.rt_scalefac=1
         stfparams.risetimedep=6.5
         stfparams.risetimedep_range=1.5
         stfparams.risetimefac=2.0

         side_taper=0.05
         bot_taper=0.1
         top_taper=0.0

************************************************************************************

   2013/10/01: version 3.2.1

   * Changed default/recommended usage to directly produce SRF file (instead of GSF
     and then using genericslip2srf to get SRF).  All options existed in previous
     versions. To directly produce SRF, use these steps (see README also):

          1) Modify output format by setting:

             write_srf=1
             write_gsf=0

          so that the code will directly produce the SRF file without having to
          make subsequent call to 'generic_slip2srf'

          2) Additional parameters needed for SRF output are:

             stype=$STYPE
             dt=$DT
             risetime=$RISETIME
             plane_header=1
             risetimefac=$RTFAC
             risetimedep=$RTDEP
             risetimedep_range=$RTDEP_RANGE
             rt_scalefac=$RT_SCALEFAC

          These should already be defined in the script.

          3) The default setting with "write_srf=1" will write to stdout.  So the
          call to the code should be something like:

          genslip-v3.2.1 read_erf=0 write_srf=1 ... rt_scalefac=$RT_SCALEFAC > $SRFDIR/$SRFFILE

          And then remove the subsequent call to 'generic_slip2srf'.

   * Additional modifications based on results of SCEC BBP testing to reduce
     the strength/coherency of longer period radiation.

   * Changed the default tapering parameters to:

         side_taper=0.02
         bot_taper=0.0
         top_taper=0.0

   * Added following parameters to adjust rise time along bottom portion of fault
     (similar to shallow adjustments).  Defaults are:

         stfparams.deep_risetimedep=16.5      # or hypo_depth if greater
         stfparams.deep_risetimedep_range=1.5
         stfparams.deep_risetimefac=1.5

   * Also added parameter to set a minimum background slip level given as a percentage
     of the average slip amount (basically fills-in very low/zero slip patches with
     long rise time low slip. Testing indicates that this has little impact on overall
     results, so for now the default setting is to not use this option:

         slip_water_level=-1

************************************************************************************

   2013/11/19: version 3.3

   * Changed default/recommended slip rate function to 'Mliu' which is modified version
     of UCSB function.  See gen_Mliu_stf() in srf_subs.c for details.

   * Added parameter 'risetime_coef' which is used to calculate average risetime
     from moment using:

         stfparams.trise = risetime_coef*1.0e-09*exp(log(mom)/3.0);

     This option only works if risetime=-1 which is now the default.

   * Changed default/recommended coefficient for average rise time relation from '1.6'
     to '1.45' =>

         risetime_coef = 1.45

     This can be overridden by specifying as a getper parameter, e.g., for use in
     CEUS, use risetime_coef=3.75

   * Added 'alphaT' parameter for scaling of rise time for dip & rake.  Update of alphaT
     parameter defined by Graves and Pitarka (2010) =>

          alphaT = 1.0/(1.0 + fD*fR*Calpha)

          where Calpha = 0.1   (max value when fD=fR=1.0)
         
	        fD = 1.0 - (dip - 45.0)/45.0      for 45.0 < dip <= 90.0
	        fD = 1.0                          for dip <= 45.0
	        fD = 0.0                          otherwise
         
	        fR = 1.0 - (|rake - 90.0|)/90.0   for 0.0 <= rake <= 180.0
	        fR = 0.0                          otherwise

          Note: should have 0 <= dip <= 90 and -180 <= rake <= 180

     Default for Calpha is 0.1, but can be passed to code as getpar() variable.

   * Added random perturbations to tsfac so that it is not 1:1 correlated with slip.
     Perturbations are log normal with ln(sigma)=tsfac_rand*tsfac.  Default for
     tsfac_rand=0.2.

   * Added random perturbations to risetime so that it is not 1:1 correlated with
     sqrt(slip).  Perturbations are log normal with ln(sigma)=rt_rand*trise.  Default
     for rt_rand=0.5.

************************************************************************************

   2014/12/31: version 4.0

   * Added option to compute rupture initiation times using FD approach (from Afnimar and
     Koketsu, 2000).  This is done by:

         1. Computing rupture 'speed' from scaling of local Vs (see get_rspeed())
         2. Converting to rupture 'slowness' with optional randomization (see get_rslow())
         3. Computing rupture times (see wfront2d_())

     In step 1, the rupture speed can also be further scaled with the underlying slip
     (fdrup_scale_slip=1); however, default is to not use this (fdrup_scale_slip=0).

     I found that scaling rup speed with slip and/or adding random perturbations to
     rup speed (rslow) and then computing FD times diminishes the impact of these
     perturbations.  The reason why is that the FD computation "heals" the wave front
     pretty quickly so that the perturations/scaling only affect things significantly
     if they are large scale and/or large magnitude.  I found it is more effective to
     obtain short length scale timing variations by applying the slip-scaling and timing
     perturbations directly after calculating the initial rupture time values (ie., doing
     it the old way with tsfac and tsfac_rand).

     New parameters with defaults are:

         fdrup_time=0,1         => use FD approach to compute rupture times (default=1)
         fdrup_scale_slip=0,1   => scale rslow with slip prior to computing FD times (default=0)

   *** NOTE ***

        When using this option, the along strike and down-dip subfault sizes must be the same.

************************************************************************************

   2015/03/16: version 4.0.1

   * Finally (I think!) sorted out the difference of the 2*pi factor for the Somerville
     and Mai wavenumber corners.  Here is my current understanding:

     Somerville included factor of 2*pi in his experssions for wavenumber corners KCx & KCy:
         
	 log10 (KCx) = 1.72 - 0.5*M
	 log10 (KCy) = 1.93 - 0.5*M

     To get corresponding correlation lengths, the following should be used:

         Ax_som = 2*pi/KCx
	 Ay_som = 2*pi/KCy

     or
         
	 log10 (Ax_som) = 0.5*M - 1.72 - log10 (2*pi) = 0.5*M - 2.52
	 log10 (Ay_som) = 0.5*M - 1.93 - log10 (2*pi) = 0.5*M - 2.73
         
     Mai's expressions for correlation lengths are:

	 log10 (Ax_mai) = 0.5*M - 2.5
	 log10 (Ay_mai) = 0.3333*M - 1.5

     In versions 4.0 and earlier, I had incorrectly applied the 2*pi factor to Mai,
     instead of taking it out of Somerville.  This means that I was using correlation
     lengths that were a factor of 2*pi too large.  But, I managed to make things work
     OK, by using the combination of "determinstic" low wave number and "stochastic"
     high wave number described in the Appendix of Graves and Pitarka (2010). By
     including the lowest 1-2 wavenumbers from the input uniform slip distribution,
     this effectively pushed the corner out to a higher wave number (smaller correlation
     length) that was actually in the right ball park for the given magnitude.

     With the correct accounting for the 2*pi factor, this deterministic+stochastic
     adjustment is not needed.  The only requirement is to make sure the mean of the
     clip distribution is positive, which is done by forcing the DC component of the
     complex wavenumber spectrum to be that for the uniform slip input.  Thus, the
     wtD and wtD factors are no longer needed, see kfilt_gaus() and kfilt_rphs().

   * Added option to compute the random wavenumber spectrum using "gaussian" random
     numbers as described in Liu et al (2006).  This is done in the routine kfilt_gaus().
     Previously, random phase was determined explicitly by setting the phase angle
     to a random number.  Then the amplitude spectrum was fit to the desired target.
     This old way is still an option and uses the routine kfilt_rphs().

     New input parameter to select these is "use_gaus":

         use_gaus=1	=> use gaussian approach with kfilt_gaus() [default]
                 =0	=> do it the old way using kfilt_rphs()

************************************************************************************

   2015/04/01: version 5.0 (no fooling!)

   * Changed variable names for local fault coordiante system from (x,y) to (stk,dip),
     e.g.,

        nx => nstk
	ny => ndip
	dx => dstk
	dy => ddip
	etc...

     (x,y,z) will come back later when envoking fault roughness.

   * Changed and made more consistent array variables for various rupture parameters.
     Declared separate arrays for complex (wavenumber) and real (space) variables, denoted
     by '_c' and '_r' suffixes, respectively. These include:

        slip_c, slip_r:		slip values
	rake_c, rake_r:		rake values
	tsfac1_c, tsfac1_r:	1st level rupture time perturbations, correlated with slip
	rtime1_c, rtime1_r:	1st level rise time perturbations, correlated with slip

	rough_c, rough_r:	roughness perturbations
	tsfac2_c, tsfac2_r:	2nd level rupture time perturbations, correlated with roughness
				(currently set to zero via tsfac_rand=-1.0)
	rtime2_c, rtime2_r:	2nd level rise time perturbations, correlated with roughness
				(currently set to zero via rtime_rand=-1.0)

	stk_r, dip_r:		strike & dip variations derived from roughness, not yet set-up

   * Made complex array dimensions 25% larger than real array dimensions to suppress wrap
     arround effects.  This replaces the use of "flip_at_surface".

   * Allow specification of correlation levels between slip and tsfac1, rtime1.
     Correlation can be between 0 (uncorrelated) and 1.0 (1:1 correlation).  Parameters
     and defaults are:

        tsfac1_scor=0.8
        rtime1_scor=0.8

   * Changed tsfac implementation so that perturbations are given by gaussian random
     numbers (via tsfac1_r array) with zero mean and sigma set to tsfac1_sigma*tsfac.
     "tsfac" is determined as before:

        tsfac = -tsfac_coef*1.0e-09*exp(log(mom)/3.0);
     
     except the default coefficient is changed to 1.1 from 1.8:

        tsfac_coef = 1.1;
     
     These are implemented using

        psrc[ip].rupt = rt + tsfac*tsfac1_r[ip];

     where rt is the background rupture time (from FD computation).

     The sigma can be changed via tsfac1_sigma, with the default set as

        tsfac1_sigma = 1.0;

     This parameterization is similar to previous versions, except it allows for larger
     values instead of being capped at a maximum of "tsfac".

   * Changed rise time perturbation implementation so that perturbations are given by
     gaussian random numbers (via rtime_r array) with mean of 1.0 and sigma of rtime1_sigma.
     The risetime perturbations scale with sqrt(rtime1_r). The default is

        rtime1_sigma = 0.85;

     and the proportionality constant (k) is determined similar to before using rt_scalefac
     except here we explicity take into account the depth scaling of risetime.

     Note that rtime_r is kind of like a surrogate for slip in this implemention. Such that
     when rtime1_scor=1.0, this implementation is essentailly the same as the rise time
     scaling with sqrt(slip) used in previous versions.

************************************************************************************

   2015/12/18: version 5.1

   * Added function GGG to apply rupture delays at segment boundaries using FD rupture
     time calculator

************************************************************************************
*/

extern void wfront2d_(int *, int *, int *, int *, double *, int *, double *, double *, int *, double *, int *);

int main(int ac,char **av)
{
FILE *fpr, *fpw;
struct complex *slip_c, *rake_c, *tsfac1_c, *rtime1_c;
struct complex *tsfac2_c, *rtime2_c, *rough_c;
float *slip_r, *rake_r, *tsfac1_r, *rtime1_r;
float *tsfac2_r, *rtime2_r, *rough_r, *stk_r, *dip_r;
float flen, fwid, dstk, ddip, sval;
float dks2, dkd2;
float dstk3, ddip3, dks3, dkd3;
float xx, yy, zz;
float sum, neg_sum;
int nstk, ndip, nstk2, ndip2, nstk3, ndip3;
int i, j, k, ip, ip2, it;
char infile[1024], str[1024];
char init_slip_file[1024], outfile[1024];

float extend_fac = -1.0;

/* RWG 2016-06-28:
   Added flen_max & fwid_max for case of generating 1 segment of multi-segment fault,
   flen_max = total along strike length of fault, if faults abut along strike; otherwise flen_max=flen
   fwid_max = total down-dip width, if faults abut along down-dip direction; otherwise fwid_max=fwid
*/

float flen_max = -1.0;
float fwid_max = -1.0;

float r0, r1, r2, y0;

float bigM;
float clen_s, clen_d;

int stretch_kcorner = DEFAULT_STRETCH_KCORNER;

float pi = 3.14159265;
double rperd = 0.017453293;

float mag;
float side_taper = DEFAULT_SIDE_TAP;
float bot_taper = DEFAULT_BOT_TAP;
float top_taper = DEFAULT_TOP_TAP;

int truncate_zero_slip = DEFAULT_TRUNCATE_ZERO_SLIP;

int generate_seed = 0;
long seed = 0;
long starting_seed;
int dump_last_seed = 0;
char seedfile[1024];

int kmodel = MAI_FLAG;   /* default is mai */
int circular_average = DEFAULT_CIRCULAR_AVERAGE;
int modified_corners = DEFAULT_MODIFIED_CORNERS;
int use_gaus = 1;

float kx_corner, ky_corner, xmag_exponent, ymag_exponent;

float mag_area_Acoef = -1.0;
float mag_area_Bcoef = -1.0;
int use_median_mag = 0;

int random_hypo = 0;
int uniform_prob4hypo = 0;
struct hypo_distr_params hpar_as, hpar_dd;
int calc_shypo = 1;

float shypo = -1.0e+15;
float dhypo = -1.0e+15;

int nrup_min = 10;
float nrup_scale_fac = 0.5;

float avgstk, rake, rt, tsmin;
float xhypo, xdep;
struct velmodel rvmod;
double rayp, rupt_rad;

int seg_delay = 0;
int nseg, nseg_bounds, ig;
float delh, hx, hg, gwid2, *gbnd, *gwid, shypo_mseg, *rvfac_seg;

float rupture_delay = 0.0;

float sd_rand = -1.0;
float s_min, s_max, s_avg, s_fac;
float d_min, d_max, d_avg, d_fac;

float rk_min, rk_max, rk_avg, rk_sig;
float *psrc_rake;

float set_rake = -999.0;
float rand_rake_degs = DEFAULT_RAND_RAKE_DEGS;
int test;

float shypo_step = DEFAULT_SHYPO_STEP;
float shypo_min_off = DEFAULT_SHYPO_MIN_OFF;
float dhypo_frac = DEFAULT_DHYPO_FRAC;
int slips_to_hypos = DEFAULT_SLIPS_TO_HYPOS;

float sh0;
int nh = -1;
int ns = -1;
int ih, js;

float rvfrac = DEFAULT_VR_TO_VS_FRAC;
float shal_vrup = 0.6;
float htol = 0.1;

float rf, sf, tf, sabs, sden, snum;
float slp_max, slp_avg, slp_min, slp_sig;
float target_savg = -1.0;

/*
RWG 2016-06-17

   if moment_fraction < 0, scale slip to give total seismic moment determined from
                           input magnitude (default)
   if 0 < moment_fraction < 1.0, scale slip to give a moment equal to (moment_fraction*Mo)
                                 where Mo is total seismic moment determined from
				 input magnitude; useful for multi-segment ruptures
*/
float moment_fraction = -1.0;

char roughnessfile[1024];
float hcoef = 1.0;
float lambda_min = -1.0;
float lambda_min_default = 0.08; /* based on Shi and Day (2014) */
float lambda_max = -1.0;	/* default for lambda_max is 'flen' */
float alpha_rough = -1.0;
float alpha_rough_default = 0.0050119;   /* 10^(-2.3) Shi & Day 2014 */
float hrms, rgh_max, rgh_avg, rgh_min, rgh_sig;
float ss, dd, hh, ssfac, ddfac, hhfac;
float cosS, sinS, cosD, sinD;

float erad = 6378.139;
float mrot = -90.0;
double g0 = 0.0;
double b0 = 0.0;
double amat[9], ainv[9];

float shal_vrdep1, shal_vrdep2, deep_vrdep1, deep_vrdep2, vrfac;
float shal_vrup_dep = DEFAULT_DEPTH_SCALING_LEVEL;
float shal_vrup_deprange = DEFAULT_DEPTH_SCALING_RANGE;

/* 20140422 RWG: Deep rupture speed decrease */
float deep_vrup = 0.6;
float deep_vrup_dep = 17.5;
float deep_vrup_deprange = 2.5;

float tsfac = -1.0e+15;;
float tsfac_coef = 1.1;
float tsf1_max, tsf1_avg, tsf1_min, tsf1_sig;
float tsfac1_sigma = 1.0;
float tsfac1_scor = 0.8;
float tsfac_rand = -1.0;
float tsf2_max, tsf2_avg, tsf2_min, tsf2_sig;
float tsfac2_scor = 0.5;

/* 20131119 RWG:

   If risetime=-1 as input (default) code will calculate average risetime
   from moment using:

      risetime = risetime_coef*1.0e-09*exp(log(mom)/3.0);

   Default rise time relation coefficient modified from somerville, 1.6 -> 1.45
   risetime_coef = 1.45
*/

/* 20150325 RWG:

   Changed default rise time relation coefficient modified back to 1.6 (testing?)
   risetime_coef = 1.6
*/

float risetime_coef = 1.6;
float deep_risetimedep_saved, deep_vrup_dep_saved;
float dmin1, dmax1, rtfac1, dmin2, dmax2, rtfac2;

float fzero = 0.0;
float fone = 1.0;

float sigfac;
float slip_sigma = DEFAULT_SLIP_SIGMA;
float rake_sigma = 15.0;
float rtime1_sigma = DEFAULT_SLIP_SIGMA;
float rt1_max, rt1_avg, rt1_min, rt1_sig;
float rtime1_scor = 0.8;
float rtime_rand = -1.0;
float rt2_max, rt2_avg, rt2_min, rt2_sig;
float rtime2_scor = 0.5;

int svr_wt = 0;

float dtop, avgdip, mom, mag_med;
struct velmodel vmod;
char velfile[1024];

int read_erf = 1;

struct standrupformat srf;
struct srf_planerectangle *prect_ptr;
struct srf_prectsegments *prseg_ptr;
struct srf_allpoints *apnts_ptr;
struct srf_apointvalues *apval_ptr;

int read_gsf = 0;
int write_gsf = 0;

int write_srf = 1;

struct generic_slip gslip;

float *stf, elon, elat;

struct pointsource *psrc;
struct stfpar2 stfparams;

int outbin = 0;

float slip_water_level = -1;
float slipmin;

float Calpha = 0.1;
float alphaT, fD, fR, avgrak;

int fdrup_time = 0;
int fdrup_scale_slip = 0;
double dh, *fdrt, *rslw, *fspace;
float *rspd;
int nstk_fd, ndip_fd, istk_off, idip_off;
int ixs, iys, nsring, ntot, *ispace;
float rvfmin = 0.25;
float rvfmax = 1.414;
float rvel_rand = 0.0;

velfile[0] = '\0';
init_slip_file[0] = '\0';
roughnessfile[0] = '\0';

sprintf(infile,"stdin");
sprintf(outfile,"stdout");

sprintf(srf.version,"1.0");

stfparams.dt = DEFAULT_DT;
stfparams.nt = NTMAX;
stfparams.trise = -1.0;
stfparams.rt_scalefac = DEFAULT_RT_SCALEFAC;
stfparams.risetimedep = DEFAULT_DEPTH_SCALING_LEVEL;
stfparams.risetimedep_range = DEFAULT_DEPTH_SCALING_RANGE;
stfparams.risetimefac = 2.0;

/* 20140422 RWG: Changed default deep risetime parameters to following: */
stfparams.deep_risetimedep = 17.5;
stfparams.deep_risetimedep_range = 2.5;
stfparams.deep_risetimefac = 2.0;

stfparams.rt_rand = 0.0;          /* ln(sigma) */
sprintf(stfparams.stype,"Mliu");  /* default is modified UCSB sincos */

/* RWG 2014-02-20 randomized hypocenter
   default probability tapering for randomized hypocenter */

/* along strike -> */
hpar_as.x0 = 0.2;	/* default tapering starts at 20% of fault length at end end */
hpar_as.x1 = 0.8;
hpar_as.f0 = 0.1;	/* default probability at edge is 10% of probability in middle of fault */
hpar_as.f1 = 0.1;

/* down dip -> */
hpar_dd.x0 = 0.4;	/* default tapering starts at 40% along top edge of fault */
hpar_dd.x1 = 0.8;	/* default tapering starts at 20% from fault bottom */
hpar_dd.f0 = 0.01;	/* default probability at top edge is 1% of probability in middle of fault */
hpar_dd.f1 = 0.1;	/* default probability at bottom edge is 10% of probability in middle of fault */

sprintf(seedfile,"dump_last_seed.txt");

setpar(ac,av);

getpar("infile","s",infile);
getpar("outfile","s",outfile);
getpar("outbin","d",&outbin);

getpar("velfile","s",velfile);

getpar("read_erf","d",&read_erf);
getpar("read_gsf","d",&read_gsf);

getpar("write_srf","d",&write_srf);

getpar("srf_version","s",srf.version);
if(strcmp(srf.version,"1.0") != 0 && strcmp(srf.version,"2.0") != 0)
   sprintf(srf.version,"1.0");

if(read_gsf == 1)
   {
   read_erf = 0;
   mstpar("mag","f",&mag);

   mstpar("nstk","d",&nstk);
   mstpar("ndip","d",&ndip);
   getpar("write_gsf","d",&write_gsf);
   }

if(read_erf == 0 && read_gsf == 0)
   {
   mstpar("mag","f",&mag);

   mstpar("nstk","d",&nstk);
   mstpar("ndip","d",&ndip);
   mstpar("dstk","f",&dstk);
   mstpar("ddip","f",&ddip);

   mstpar("dtop","f",&dtop);
   mstpar("strike","f",&avgstk);
   mstpar("dip","f",&avgdip);
   mstpar("rake","f",&rake);
   mstpar("elon","f",&elon);
   mstpar("elat","f",&elat);
   }

getpar("dt","f",&stfparams.dt);
getpar("nt","d",&stfparams.nt);
getpar("risetime","f",&stfparams.trise);
getpar("risetime_coef","f",&risetime_coef);
getpar("risetimefac","f",&stfparams.risetimefac);
getpar("risetimedep","f",&stfparams.risetimedep);
getpar("risetimedep_range","f",&stfparams.risetimedep_range);
getpar("rt_scalefac","f",&stfparams.rt_scalefac);
getpar("rt_rand","f",&stfparams.rt_rand);
getpar("stype","s",stfparams.stype);

getpar("deep_risetimefac","f",&stfparams.deep_risetimefac);
getpar("deep_risetimedep","f",&stfparams.deep_risetimedep);
getpar("deep_risetimedep_range","f",&stfparams.deep_risetimedep_range);

deep_risetimedep_saved = stfparams.deep_risetimedep;

getpar("shypo_step","f",&shypo_step);
getpar("shypo_min_off","f",&shypo_min_off);
getpar("dhypo_frac","f",&dhypo_frac);
getpar("slips_to_hypos","d",&slips_to_hypos);
getpar("ns","d",&ns);
getpar("nh","d",&nh);

getpar("shypo","f",&shypo);
getpar("dhypo","f",&dhypo);

/*XXXX*/
/* RWG 2014-02-20 randomized hypocenter */
getpar("random_hypo","d",&random_hypo);
getpar("generate_seed","d",&generate_seed);
getpar("nrup_min","d",&nrup_min);
getpar("nrup_scale_fac","f",&nrup_scale_fac);

getpar("uniform_prob4hypo","d",&uniform_prob4hypo);

getpar("hypo_taperperc_left","f",&hpar_as.x0);
getpar("hypo_taperval_left","f",&hpar_as.f0);
getpar("hypo_taperperc_right","f",&hpar_as.x1);
getpar("hypo_taperval_right","f",&hpar_as.f1);

getpar("hypo_taperperc_top","f",&hpar_dd.x0);
getpar("hypo_taperval_top","f",&hpar_dd.f0);
getpar("hypo_taperperc_bot","f",&hpar_dd.x1);
getpar("hypo_taperval_bot","f",&hpar_dd.f1);
/*XXXX*/

getpar("dump_last_seed","d",&dump_last_seed);
getpar("seedfile","s",seedfile);

getpar("rvfrac","f",&rvfrac);
getpar("shal_vrup","f",&shal_vrup);
getpar("shal_vrup_dep","f",&shal_vrup_dep);
getpar("shal_vrup_deprange","f",&shal_vrup_deprange);

/* 20140422 RWG: Deep rupture speed decrease */
getpar("deep_vrup","f",&deep_vrup);
getpar("deep_vrup_dep","f",&deep_vrup_dep);
getpar("deep_vrup_deprange","f",&deep_vrup_deprange);

deep_vrup_dep_saved = deep_vrup_dep;

getpar("rupture_delay","f",&rupture_delay);

getpar("fdrup_time","d",&fdrup_time);
getpar("fdrup_scale_slip","d",&fdrup_scale_slip);
getpar("rvfmin","f",&rvfmin);
getpar("rvfmax","f",&rvfmax);

getpar("tsfac","f",&tsfac);
getpar("tsfac_coef","f",&tsfac_coef);
getpar("tsfac_rand","f",&tsfac_rand);
getpar("rtime_rand","f",&rtime_rand);
getpar("sd_rand","f",&sd_rand);

getpar("rand_rake_degs","f",&rand_rake_degs);
getpar("set_rake","f",&set_rake);

getpar("target_savg","f",&target_savg);
getpar("moment_fraction","f",&moment_fraction);

getpar("kmodel","d",&kmodel);
getpar("use_gaus","d",&use_gaus);

if(kmodel < 0)
   kmodel = INPUT_CORNERS_FLAG;
if(kmodel != MAI_FLAG && kmodel != INPUT_CORNERS_FLAG && kmodel != FRANKEL_FLAG && kmodel < 100)
   kmodel = SOMERVILLE_FLAG;

fprintf(stderr,"kmodel= %d\n",kmodel);

if(kmodel == INPUT_CORNERS_FLAG)
   {
   mstpar("kx_corner","f",&kx_corner);
   mstpar("ky_corner","f",&ky_corner);

   xmag_exponent = 0.5;
   ymag_exponent = 0.5;
   getpar("xmag_exponent","f",&xmag_exponent);
   getpar("ymag_exponent","f",&ymag_exponent);
   }

getpar("mag_area_Acoef","f",&mag_area_Acoef);
getpar("mag_area_Bcoef","f",&mag_area_Bcoef);
getpar("use_median_mag","d",&use_median_mag);

getpar("modified_corners","d",&modified_corners);
getpar("circular_average","d",&circular_average);
getpar("stretch_kcorner","d",&stretch_kcorner);
getpar("seed","d",&seed);
getpar("side_taper","f",&side_taper);
getpar("bot_taper","f",&bot_taper);
getpar("top_taper","f",&top_taper);

getpar("truncate_zero_slip","d",&truncate_zero_slip);
getpar("slip_water_level","f",&slip_water_level);

getpar("slip_sigma","f",&slip_sigma);
getpar("rake_sigma","f",&rake_sigma);
getpar("tsfac1_sigma","f",&tsfac1_sigma);
getpar("tsfac1_scor","f",&tsfac1_scor);
getpar("rtime1_sigma","f",&rtime1_sigma);
getpar("rtime1_scor","f",&rtime1_scor);
getpar("svr_wt","d",&svr_wt);

getpar("flen_max","f",&flen_max);
getpar("fwid_max","f",&fwid_max);
getpar("extend_fac","f",&extend_fac);
if(extend_fac < 0)
   extend_fac = 1.10;

getpar("roughnessfile","s",roughnessfile);
getpar("alpha_rough","f",&alpha_rough);
getpar("lambda_min","f",&lambda_min);
getpar("lambda_max","f",&lambda_max);

if(alpha_rough < 0.0)
   alpha_rough = alpha_rough_default;
if(lambda_min < 0.0)
   lambda_min = lambda_min_default;

getpar("tsfac2_scor","f",&tsfac2_scor);
getpar("rtime2_scor","f",&rtime2_scor);

getpar("init_slip_file","s",init_slip_file);

getpar("seg_delay","d",&seg_delay);
if(seg_delay == 1)
   {
   mstpar("nseg","d",&nseg);

   if(nseg <= 1)
      seg_delay = 0;
   else
      {
      nseg_bounds = nseg - 1;

      rvfac_seg = (float *)check_malloc(nseg_bounds*sizeof(float));
      gbnd = (float *)check_malloc(nseg_bounds*sizeof(float));
      gwid = (float *)check_malloc(nseg_bounds*sizeof(float));

      for(ig=0;ig<nseg_bounds;ig++)
         {
         rvfac_seg[ig] = 0.5;    /* default is 50% reduction of rupture speed at seg boundaries */
         gwid[ig] = 6.0;	/* default width of delay zone is 6 km (3 km on each seg) */
         }

      getpar("rvfac_seg","vf",rvfac_seg);
      getpar("gwid","vf",gwid);

      if(fdrup_time == 0) /* change rvfac_seg to equivalent distance adjustment */
         {
         for(ig=0;ig<nseg_bounds;ig++)
            rvfac_seg[ig] = 1.0/rvfac_seg[ig] - 1.0;
         }
      }
   }

endpar();

psrc = (struct pointsource *)NULL;
gslip.np = -1;
gslip.spar = (struct slippars *)NULL;

if(read_erf == 1)
   psrc = read_ruppars(infile,psrc,&mag,&nstk,&ndip,&dstk,&ddip,&dtop,&avgstk,&avgdip,&elon,&elat);
else if(read_gsf == 1)
   psrc = read_gsfpars(infile,psrc,&gslip,&dstk,&ddip,&dtop,&avgdip);
else
   psrc = set_ruppars(psrc,&mag,&nstk,&ndip,&dstk,&ddip,&dtop,&avgstk,&avgdip,&rake,&elon,&elat);

flen = nstk*dstk;
fwid = ndip*ddip;

if(flen_max < flen)
   flen_max = flen;
if(fwid_max < fwid)
   fwid_max = fwid;

if(lambda_max < 0.0)
   lambda_max = flen_max;

bigM = log(10.0);
mom = exp(bigM*1.5*(mag + 10.7));
mag_med = 3.98 + log(flen*fwid)/bigM;
/* update 12/2005 */
mag_med = 3.87 + 1.05*log(flen*fwid)/bigM;

if(tsfac < -1.0e+10)
   tsfac = -tsfac_coef*1.0e-09*exp(log(mom)/3.0);

/* update 12/2007 */
if(mag_area_Acoef < 0.0)
   {
   mag_area_Acoef = 4.0;
   mag_area_Bcoef = 1.0;
   }
mag_med = mag_area_Acoef + mag_area_Bcoef*log(flen*fwid)/bigM;

if(mag < 0.0 || use_median_mag)
   mag = mag_med;

init_plane_srf(&srf,&gslip,&elon,&elat,nstk,ndip,&flen,&fwid,&dstk,&ddip,&avgstk,&avgdip,&dtop,&shypo,&dhypo);

if(seg_delay == 1)
   {
   nseg_bounds = get_seg_bounds(&srf,gbnd);

   if(nseg_bounds != (nseg-1))
      {
      fprintf(stderr,"Expecting %d segments, found %d; exiting ...\n",nseg,nseg_bounds+1);
      exit(-1);
      }
   }

if(velfile[0] != '\0')
   read_Fvelmodel(velfile,&vmod);
   /*
   read_velmodel(velfile,&vmod);
   */
else
   default_velmodel(&vmod);

/* generic fault length/width scaling */
clen_s = flen;
clen_d = fwid;
clen_s = flen_max;
clen_d = fwid_max;

if(clen_s > 2.0*clen_d)
   clen_s = 2.0*clen_d;
if(clen_d > 2.0*clen_s)
   clen_d = 2.0*clen_s;

if(kmodel == MAI_FLAG || kmodel == FRANKEL_FLAG) /* mai scaling */
   {
   /* add factor of log10(2*pi)
   xl = exp(bigM*(0.5*mag - 2.50 + 0.79818));
   yl = exp(bigM*(0.3333*mag - 1.50 + 0.79818));
   */

/* 2015-03-13 RWG
   I Believe the above is wrong, but seems to work OK with the wtD/wtS Deterministic+Stochastic
   combination used prior to this date.  This combination artificially extends the corners to
   somewhat higher wavenumbers and things look better.  

   But I think Somerville is incorrect in that the wavenumber corners actually contain the extra
   factor of 2*pi, whereas Mai and Beroza don't.  Somerville worked OK but required the ad-hoc
   wtD + wtS adjustment to push the corners out a bit.

   Subtracting a 2*pi factor from Somerville now appears to work OK.

   Still retain the uniform slip input DC component but no wtD/wtS adjustment is needed.
   See kfilt_gaus() and kfilt_rphs().
   
*/

   clen_s = exp(bigM*(0.5*mag - 2.50));
   clen_d = exp(bigM*(0.3333*mag - 1.50));

   if(circular_average)
      clen_d = clen_s;
   }

if(kmodel == SOMERVILLE_FLAG) /* somerville scaling */
   {
   /*
   xl = exp(bigM*(0.5*mag - 1.72));
   yl = exp(bigM*(0.5*mag - 1.93));
   */
   /* subtract factor of log10(2*pi), see explanation above */
   clen_s = exp(bigM*(0.5*mag - 1.72 - 0.79818));
   clen_d = exp(bigM*(0.5*mag - 1.93 - 0.79818));

   if(circular_average)
      {
      clen_s = exp(bigM*(0.5*mag - 1.825 - 0.79818));
      clen_d = clen_s;
      }
   }

if(modified_corners)
   {
   clen_s = exp(bigM*(0.5*mag - 2.00));
   clen_d = exp(bigM*(0.5*mag - 2.00));
   }

if(kmodel == INPUT_CORNERS_FLAG)
   {
   clen_s = exp(bigM*(xmag_exponent*mag - kx_corner));
   clen_d = exp(bigM*(ymag_exponent*mag - ky_corner));
   }

if(stfparams.trise < 0.0)
   stfparams.trise = risetime_coef*1.0e-09*exp(log(mom)/3.0);

/* 20131118: rise time modification for rake and dip => update of GP2010 alphaT parameter */

Calpha = 0.1;
         
fD = 0.0;
if(avgdip <= 90.0 && avgdip > 45.0)
   fD = 1.0 - (avgdip - 45.0)/45.0;
else if(avgdip <= 45.0 && avgdip >= 0.0)
   fD = 1.0;

avgrak = 0.0;
for(j=0;j<ndip*nstk;j++)
   avgrak = avgrak + psrc[j].rak;

avgrak = avgrak/((float)(nstk*ndip));
while(avgrak < -180.0)
   avgrak = avgrak + 360.0;
while(avgrak > 180.0)
   avgrak = avgrak - 360.0;
         
fR = 0.0;
if(avgrak <= 180.0 && avgrak >= 0.0)
   fR = 1.0 - sqrt((avgrak - 90.0)*(avgrak - 90.0))/90.0;
         
alphaT = 1.0/(1.0 + fD*fR*Calpha);
stfparams.trise = alphaT*stfparams.trise;

/* alphaT done */

/*
  Prior to 2015-03-25: Doubled fault length and width for v5.0.
  Works well in eliminating wrap-around, but creates problems when correlating
  slip with other fields.

  2015-03-25: Try extending length & width by just 25%
  2015-04-22: Try extending length & width by just 10%

nstk2 = 2*nstk;
ndip2 = 2*ndip;
*/

/*
  RWG 2016-06-28: Add flen_max/flen and fwid_max/fwid adjustment for scaling lengths/widths
  This means the full fault length/width will be used for determining waenumber spectra,
  as opposed to limiting wavenumbers to that given by only individual segment

  Only makes sense when generating individual segment that is part of larger multi-segment fault.
*/

nstk2 = (int)((flen_max*extend_fac/flen)*nstk);
if(nstk2%2)
   nstk2++;

ndip2 = (int)((fwid_max*extend_fac/fwid)*ndip);
if(ndip2%2)
   ndip2++;

dstk3 = dstk/3.0;
ddip3 = ddip/3.0;

nstk3 = 3*nstk;
nstk3 = (int)((flen_max*extend_fac/flen)*nstk3);
if(nstk3%2)
   nstk3++;

ndip3 = 3*ndip;
ndip3 = (int)((fwid_max*extend_fac/fwid)*ndip3);
if(ndip3%2)
   ndip3++;

/*XXXXX
nstk3 = nstk;
ndip3 = ndip;
dstk3 = dstk;
ddip3 = ddip;
*/

if(nstk3 < ndip3)
   nstk3 = ndip3;
else if(ndip3 < nstk3)
   ndip3 = nstk3;

dks2 = 1.0/(nstk2*dstk);
dkd2 = 1.0/(ndip2*ddip);

dks3 = 1.0/(nstk3*dstk3);
dkd3 = 1.0/(ndip3*ddip3);

slip_c = (struct complex *) check_malloc (nstk2*ndip2*sizeof(struct complex));
rake_c = (struct complex *) check_malloc (nstk2*ndip2*sizeof(struct complex));
tsfac1_c = (struct complex *) check_malloc (nstk2*ndip2*sizeof(struct complex));
rtime1_c = (struct complex *) check_malloc (nstk2*ndip2*sizeof(struct complex));

rough_c = (struct complex *) check_malloc (nstk3*ndip3*sizeof(struct complex));
tsfac2_c = (struct complex *) check_malloc (nstk3*ndip3*sizeof(struct complex));
rtime2_c = (struct complex *) check_malloc (nstk3*ndip3*sizeof(struct complex));

slip_r = (float *) check_malloc (nstk*ndip*sizeof(float));
rake_r = (float *) check_malloc (nstk*ndip*sizeof(float));
tsfac1_r = (float *) check_malloc (nstk*ndip*sizeof(float));
rtime1_r = (float *) check_malloc (nstk*ndip*sizeof(float));

rough_r = (float *) check_malloc (nstk*ndip*sizeof(float));
tsfac2_r = (float *) check_malloc (nstk*ndip*sizeof(float));
rtime2_r = (float *) check_malloc (nstk*ndip*sizeof(float));
stk_r = (float *) check_malloc (nstk*ndip*sizeof(float));
dip_r = (float *) check_malloc (nstk*ndip*sizeof(float));

psrc_rake = (float *) check_malloc (nstk*ndip*sizeof(float));

for(j=0;j<ndip*nstk;j++)
   psrc_rake[j] = psrc[j].rak;

/* XXXX */
/* RWG 2014-02-20 randomized hypocenter */
if(generate_seed == 1)
   gseed(&seed,psrc,&flen,&fwid,&dtop,&mag);

/* RWG 2014-03-21 set starting seed */
starting_seed = seed;

if(random_hypo == 1)
   {
   hpar_as.x0 = hpar_as.x0*flen;
   hpar_as.x1 = hpar_as.x1*flen;
   hpar_as.xlen = flen;
   hpar_as.xshift = -0.5*flen;

   hpar_dd.x0 = hpar_dd.x0*fwid;
   hpar_dd.x1 = hpar_dd.x1*fwid;
   hpar_dd.xlen = fwid;
   hpar_dd.xshift = 0.0;
   }

if(read_erf == 1)
   {
   if(random_hypo == 1)
      {
      calc_shypo = 0;

      ns = (int)(0.1*nrup_scale_fac*flen*fwid + 0.5);
      if(ns < nrup_min)
         ns = nrup_min;

      nh = 1;
      }
   else
      {
      calc_shypo = 1;
      shypo = -1.0e+15;
      dhypo = -1.0e+15;

      if(nh < 0)
         nh = (int)((flen-2.0*shypo_min_off)/shypo_step) + 1;
      if(ns < 0)
         {
         ns = slips_to_hypos*nh;
         if(nstk == 1 && ndip == 1)
            ns = 1;
         }

      sh0 = 0.5*(flen - (nh-1)*shypo_step);
      }
   }

if(shypo > -1.0e+14)
   calc_shypo = 0;

if(dhypo < -1.0e+14)
   dhypo = dhypo_frac*fwid;
/* XXXX */

fprintf(stderr,"mag= %.2f median mag= %.2f nslip= %d nhypo= %d\n",mag,mag_med,ns,nh);

fprintf(stderr,"nstk=  %d ndip=  %d dstk=  %12.8f ddip=  %12.8f\n",nstk,ndip,dstk,ddip);
fprintf(stderr,"nstk2= %d ndip2= %d\n",nstk2,ndip2);
fprintf(stderr,"nstk3= %d ndip3= %d dstk3= %12.8f ddip3= %12.8f\n",nstk3,ndip3,dstk3,ddip3);

if(fdrup_time == 1)
   {
   if(dstk/ddip > 1.01 || dstk/ddip < 0.99)
      fprintf(stderr,"*** dstk= %13.5e != ddip= %13.5e; possible problems with fdrup_time\n",dstk,ddip);

   dh = dstk;
   nsring = 2;

   rspd = (float *)check_malloc(nstk*ndip*sizeof(float));

   rslw = NULL;
   fdrt = NULL;

   fspace = NULL;
   ispace = NULL;
   }

sval = 0.0;
for(js=0;js<ns;js++)    /* loop over slip/rupture realizations */
   {
/*
   RWG 2014-03-21 set initial seed using increments of starting seed.
   Allows reproducability of ruptures without having to generate entire set.
*/

   seed = starting_seed;
   for(k=0;k<10*js;k++)
      sval = sfrand(&seed);

   fprintf(stderr,"js= %d seed= %d ran= %10.6f\n",js,seed,sval);

/* RWG 2014-02-20 randomized hypocenter */
   if(random_hypo == 1)
      {
      if(uniform_prob4hypo == 1)
         rhypo_uniform(&seed,&shypo,&dhypo,&flen,&fwid);
      else
         {
         shypo = rhypo1_lintaper(&seed,&hpar_as);
         dhypo = rhypo1_lintaper(&seed,&hpar_dd);
	 }
      }

/* do slip */

   /*
   init_slip_IO(slip_c,nstk2,ndip2,&dstk,&ddip,0,init_slip_file);
   */

   for(j=0;j<ndip2*nstk2;j++)
      {
      slip_c[j].re = 1.0;
      slip_c[j].im = 0.0;
      }

   fft2d(slip_c,nstk2,ndip2,-1,&dstk,&ddip);

   if(kmodel < 100)
      {
      if(use_gaus)
         kfilt_gaus(slip_c,nstk2,ndip2,&dks2,&dkd2,&clen_s,&clen_d,&seed,kmodel);
      else
         kfilt_rphs(slip_c,nstk2,ndip2,&dks2,&dkd2,&clen_s,&clen_d,&seed,kmodel);
      }

   fft2d(slip_c,nstk2,ndip2,1,&dks2,&dkd2);

/*
   Copy in to real array.

   Compute average for initial field, if negative change polarity of slip_c
   to ensure positivity of later correlations
*/
   slp_avg = 0.0;
   for(j=0;j<ndip;j++)
      {
      for(i=0;i<nstk;i++)
         {
	 ip = i + j*nstk;
	 ip2 = i + j*nstk2;
         slip_r[ip] = slip_c[ip2].re;
         slp_avg = slp_avg + slip_r[ip];
         }
      }

   if(slp_avg < 0.0)
      {
      for(ip=0;ip<nstk2*ndip2;ip++)
         {
         slip_c[ip].re = -slip_c[ip].re;
         slip_c[ip].im = -slip_c[ip].im;
	 }

      for(ip=0;ip<nstk*ndip;ip++)
         slip_r[ip] = -slip_r[ip];

      slp_avg = -slp_avg;
      }

   if(kmodel == FRANKEL_FLAG)
      {
      slip_sigma = -1.0;
      slp_min = 1.0e+20;

      for(ip=0;ip<nstk*ndip;ip++)
         {
	 if(slip_r[ip] < slp_min)
	    slp_min = slip_r[ip];
	 }

      slp_avg = 0.0;
      for(ip=0;ip<nstk*ndip;ip++)
         {
         slip_r[ip] = slip_r[ip] - slp_min;
         slp_avg = slp_avg + slip_r[ip];
	 }
      }

   fprintf(stderr,"slp_avg= %13.5e, normalized to 1\n",slp_avg/(float)(nstk*ndip));

   slp_avg = slp_avg/(float)(nstk*ndip);
   for(ip=0;ip<ndip*nstk;ip++)
      slip_r[ip] = slip_r[ip]/slp_avg;

   slp_avg = 1.0;

/*
   scale to desired sigma
*/

   slp_sig = 0.0;
   for(j=0;j<ndip;j++)
      {
      for(i=0;i<nstk;i++)
         {
	 ip = i + j*nstk;
         slp_sig = slp_sig + (slip_r[ip]-slp_avg)*(slip_r[ip]-slp_avg);
         }
      }
   slp_sig = sqrt(slp_sig/(ndip*nstk))/slp_avg;

   if(slip_sigma > 0.0)
      {
      sigfac = slip_sigma/slp_sig;
      for(j=0;j<ndip;j++)
         {
         for(i=0;i<nstk;i++)
            {
	    ip = i + j*nstk;
            slip_r[ip] = sigfac*(slip_r[ip] - slp_avg) + slp_avg;
            }
         }
      }

/*
   truncate any negative slip values => should double check that spectra is
   not changed too much
*/

   neg_sum = 0.0;
   sum = 0.0;
   for(j=0;j<ndip;j++)
      {
      for(i=0;i<nstk;i++)
         {
	 ip = i + j*nstk;
         if(slip_r[ip] < 0.0 && truncate_zero_slip)
	    {
	    neg_sum = neg_sum - slip_r[ip];
            slip_r[ip] = 0.0;
	    /*
	    fprintf(stderr,"slip= %13.5e\n",slip_r[ip]);
	    */
	    }
	 sum = sum + slip_r[ip];
	 }
      }

   sum = sum/(float)(nstk*ndip);
   neg_sum = neg_sum/(float)(nstk*ndip);

   if(truncate_zero_slip)
      fprintf(stderr,"ratio (negative slip)/(positive slip)= %f\n",neg_sum/(sum+neg_sum));

   taper_slip_all_r(slip_r,nstk,ndip,&side_taper,&bot_taper,&top_taper);

   if(slip_water_level > 0.0)
      {
      slipmin = sum*slip_water_level;
      for(j=0;j<ndip;j++)
         {
         for(i=0;i<nstk;i++)
            {
            ip = i + j*nstk;
            if(slip_r[ip] < slipmin)
               slip_r[ip] = slipmin;
            }
         }
      }

/* check moment and scale slip */

   slp_avg = target_savg;
   if(moment_fraction > 0.0)
      mom = mom*moment_fraction;

   scale_slip_r(psrc,slip_r,nstk,ndip,0,&dstk,&ddip,&dtop,&avgdip,&mom,&vmod,&slp_avg,&slp_max);

   fprintf(stderr,"mom= %13.5e avgslip= %.0f maxslip= %.0f\n",mom,slp_avg,slp_max);
   fprintf(stderr,"orig_sigma= %f ... ",slp_sig);

/* recalculate just to check, plus normalize by savg */
   slp_sig = 0.0;
   for(j=0;j<ndip;j++)
      {
      for(i=0;i<nstk;i++)
         {
	 ip = i + j*nstk;
         slp_sig = slp_sig + (psrc[ip].slip-slp_avg)*(psrc[ip].slip-slp_avg);
         }
      }
   slp_sig = sqrt(slp_sig/(ndip*nstk))/slp_avg;

   fprintf(stderr,"new_sigma= %f\n",slp_sig);

/*
   transform slip_c back to wavenumber domain for later use in correlations
*/
   fft2d(slip_c,nstk2,ndip2,-1,&dstk,&ddip);

/* end slip */

/* now do rake */

   for(j=0;j<ndip2*nstk2;j++)
      {
      rake_c[j].re = 1.0;
      rake_c[j].im = 0.0;
      }

   fft2d(rake_c,nstk2,ndip2,-1,&dstk,&ddip);

   if(use_gaus)
      kfilt_gaus(rake_c,nstk2,ndip2,&dks2,&dkd2,&clen_s,&clen_d,&seed,kmodel);
   else
      kfilt_rphs(rake_c,nstk2,ndip2,&dks2,&dkd2,&clen_s,&clen_d,&seed,kmodel);

   fft2d(rake_c,nstk2,ndip2,1,&dks2,&dkd2);

   rk_avg = 0.0;
   for(j=0;j<ndip;j++)
      {
      for(i=0;i<nstk;i++)
         {
	 ip = i + j*nstk;
	 ip2 = i + j*nstk2;
         rake_r[ip] = rake_c[ip2].re;

         rk_avg = rk_avg + rake_r[ip];
         }
      }
   rk_avg = rk_avg/(float)(nstk*ndip);

   rk_sig = 0.0;
   for(ip=0;ip<ndip*nstk;ip++)
      {
      rake_r[ip] = rake_r[ip] - rk_avg;

      if(rake_r[ip] < rk_min)
         rk_min = rake_r[ip];
      if(rake_r[ip] > rk_max)
         rk_max = rake_r[ip];

      rk_sig = rk_sig + (rake_r[ip])*(rake_r[ip]);
      }
   rk_sig = sqrt(rk_sig/(ndip*nstk));

   if(rake_sigma > 0.0)
      {
      sigfac = rake_sigma/rk_sig;
      for(ip=0;ip<nstk*ndip;ip++)
         rake_r[ip] = sigfac*rake_r[ip];
      }

   rk_min = 1.0e+15;
   rk_max = -1.0e+15;
   rk_avg = 0.0;
   rk_sig = 0.0;
   for(ip=0;ip<ndip*nstk;ip++)
      {
      if(rake_r[ip] < rk_min)
         rk_min = rake_r[ip];
      if(rake_r[ip] > rk_max)
         rk_max = rake_r[ip];

      rk_sig = rk_sig + (rake_r[ip])*(rake_r[ip]);
      }
   rk_sig = sqrt(rk_sig/(ndip*nstk));

   fprintf(stderr,"rake: avg= %13.5e min= %13.5e max= %13.5e sig= %13.5e\n",rk_avg,rk_min,rk_max,rk_sig);

   if(set_rake > -900.0)
      {
      for(j=0;j<ndip*nstk;j++)
         psrc_rake[j] = set_rake;
      }

   for(j=0;j<ndip*nstk;j++)
      psrc[j].rak = rake_r[j] + psrc_rake[j];

/* end rake */

/* now do tsfac1 */

   for(j=0;j<ndip2*nstk2;j++)
      {
      tsfac1_c[j].re = 1.0;
      tsfac1_c[j].im = 0.0;
      }

   fft2d(tsfac1_c,nstk2,ndip2,-1,&dstk,&ddip);

   if(use_gaus)
      kfilt_gaus(tsfac1_c,nstk2,ndip2,&dks2,&dkd2,&clen_s,&clen_d,&seed,kmodel);
   else
      kfilt_rphs(tsfac1_c,nstk2,ndip2,&dks2,&dkd2,&clen_s,&clen_d,&seed,kmodel);

/*
   correlate with slip_c
*/

   sf = sqrt(1.0 - tsfac1_scor*tsfac1_scor);
   for(j=0;j<ndip2*nstk2;j++)
      {
      tsfac1_c[j].re = tsfac1_scor*slip_c[j].re + sf*tsfac1_c[j].re;
      tsfac1_c[j].im = tsfac1_scor*slip_c[j].im + sf*tsfac1_c[j].im;
      }

   fft2d(tsfac1_c,nstk2,ndip2,1,&dks2,&dkd2);

/*
   Copy in to real array.
   Compute average for initial field then remove mean.
*/
   tsf1_avg = 0.0;
   for(j=0;j<ndip;j++)
      {
      for(i=0;i<nstk;i++)
         {
	 ip = i + j*nstk;
	 ip2 = i + j*nstk2;
         tsfac1_r[ip] = tsfac1_c[ip2].re;

         tsf1_avg = tsf1_avg + tsfac1_r[ip];
         }
      }
   fprintf(stderr,"tsf1_avg= %13.5e, removed from field\n",tsf1_avg/(float)(nstk*ndip));

   tsf1_avg = tsf1_avg/(float)(nstk*ndip);
   for(ip=0;ip<ndip*nstk;ip++)
      tsfac1_r[ip] = tsfac1_r[ip] - tsf1_avg;

   tsf1_avg = 0.0;

/*
   scale to desired sigma
*/
   tsf1_sig = 0.0;
   for(ip=0;ip<ndip*nstk;ip++)
      tsf1_sig = tsf1_sig + (tsfac1_r[ip])*(tsfac1_r[ip]);

   tsf1_sig = sqrt(tsf1_sig/(ndip*nstk));

   if(tsfac1_sigma < 0.0)
      tsfac1_sigma = 0.0;

   sigfac = tsfac1_sigma/tsf1_sig;
   for(ip=0;ip<nstk*ndip;ip++)
      tsfac1_r[ip] = sigfac*tsfac1_r[ip];

   tsf1_min = 1.0e+15;
   tsf1_max = -1.0e+15;
   tsf1_sig = 0.0;
   for(ip=0;ip<nstk*ndip;ip++)
      {
      if(tsfac1_r[ip] < tsf1_min)
         tsf1_min = tsfac1_r[ip];
      if(tsfac1_r[ip] > tsf1_max)
         tsf1_max = tsfac1_r[ip];

      tsf1_sig = tsf1_sig + (tsfac1_r[ip])*(tsfac1_r[ip]);
      }

   tsf1_sig = sqrt(tsf1_sig/(ndip*nstk));

   fprintf(stderr,"tsfac1: avg= %13.5e min= %13.5e max= %13.5e sig= %13.5e\n",tsf1_avg,tsf1_min,tsf1_max,tsf1_sig);

/* end tsfac1 */

/* now do rtime1 */

   for(j=0;j<ndip2*nstk2;j++)
      {
      rtime1_c[j].re = 1.0;
      rtime1_c[j].im = 0.0;
      }

   fft2d(rtime1_c,nstk2,ndip2,-1,&dstk,&ddip);

   if(use_gaus)
      kfilt_gaus(rtime1_c,nstk2,ndip2,&dks2,&dkd2,&clen_s,&clen_d,&seed,kmodel);
   else
      kfilt_rphs(rtime1_c,nstk2,ndip2,&dks2,&dkd2,&clen_s,&clen_d,&seed,kmodel);

/*
   correlate with slip_c
*/
   sf = sqrt(1.0 - rtime1_scor*rtime1_scor);
   for(j=0;j<ndip2*nstk2;j++)
      {
      rtime1_c[j].re = rtime1_scor*slip_c[j].re + sf*rtime1_c[j].re;
      rtime1_c[j].im = rtime1_scor*slip_c[j].im + sf*rtime1_c[j].im;
      }

   fft2d(rtime1_c,nstk2,ndip2,1,&dks2,&dkd2);

/*
   Copy in to real array.
   Compute average for initial field then normalize by average,
   this will take care of sign change, if necessary
*/
   rt1_avg = 0.0;
   for(j=0;j<ndip;j++)
      {
      for(i=0;i<nstk;i++)
         {
         ip = i + j*nstk;
         ip2 = i + j*nstk2;
         rtime1_r[ip] = rtime1_c[ip2].re;

         rt1_avg = rt1_avg + rtime1_r[ip];
         }
      }
   fprintf(stderr,"rt1_avg= %13.5e, normalized to 1\n",rt1_avg/(float)(nstk*ndip));

   rt1_avg = rt1_avg/(float)(nstk*ndip);
   for(ip=0;ip<ndip*nstk;ip++)
      rtime1_r[ip] = rtime1_r[ip]/rt1_avg;

   rt1_avg = 1.0;

/*
   average should now be 1, but may change below when negative values are truncated
*/

/*
   scale to desired sigma
*/
   rt1_sig = 0.0;
   for(ip=0;ip<ndip*nstk;ip++)
      rt1_sig = rt1_sig + (rtime1_r[ip]-rt1_avg)*(rtime1_r[ip]-rt1_avg);

   rt1_sig = sqrt(rt1_sig/(ndip*nstk))/rt1_avg;

   if(rtime1_sigma < 0.0)
      rtime1_sigma = 0.0;

   sigfac = rtime1_sigma/rt1_sig;
   for(ip=0;ip<nstk*ndip;ip++)
      rtime1_r[ip] = sigfac*(rtime1_r[ip] - rt1_avg) + rt1_avg;

/*
   truncate any negative values, reset sig, min, max and avg
*/
   neg_sum = 0.0;
   rt1_avg = 0.0;
   for(ip=0;ip<nstk*ndip;ip++)
      {
      if(rtime1_r[ip] < 0.0)
         {
         neg_sum = neg_sum - rtime1_r[ip];
         rtime1_r[ip] = 0.0;
         }

      rt1_avg = rt1_avg + rtime1_r[ip];
      }
fprintf(stderr,"ratio (negative rt)/(positive rt)= %f\n",neg_sum/(rt1_avg+neg_sum));

   rt1_avg = rt1_avg/(float)(nstk*ndip);
   for(ip=0;ip<ndip*nstk;ip++)
      rtime1_r[ip] = rtime1_r[ip]/rt1_avg;

   rt1_avg = 1.0;

   rt1_min = 1.0e+15;
   rt1_max = -1.0e+15;
   rt1_sig = 0.0;
   for(ip=0;ip<nstk*ndip;ip++)
      {
      if(rtime1_r[ip] < rt1_min)
         rt1_min = rtime1_r[ip];
      if(rtime1_r[ip] > rt1_max)
         rt1_max = rtime1_r[ip];

      rt1_sig = rt1_sig + (rtime1_r[ip]-rt1_avg)*(rtime1_r[ip]-rt1_avg);
      }

   rt1_sig = sqrt(rt1_sig/(ndip*nstk))/rt1_avg;

   fprintf(stderr,"rtime1: avg= %13.5e min= %13.5e max= %13.5e sig= %13.5e\n",rt1_avg,rt1_min,rt1_max,rt1_sig);

   stfparams.deep_risetimedep = deep_risetimedep_saved;
   xhypo = dhypo*sin(avgdip*rperd) + dtop + stfparams.deep_risetimedep_range;
   if(xhypo > stfparams.deep_risetimedep)
      stfparams.deep_risetimedep = xhypo;

/*
   Compute rt_scalefac

   Here we treat rtime1_r as a surrogate for slip.  Since rtime1_r is correlated with slip,
   having risetime scale with sqrt(rtime1_r) is analogous to having it scale with sqrt(slip),
   like was done in GP2010.

   constant "k" = 1/stfparams.rt_scalefac

   if svr_wt == 1 then
      use (slip*vrup)-weighted averaging, including depth variation to determine "k"
      Note: increase of rise time is roughly balanced by decrease in vrup
   else
      simple averaging of risetime across fault is used to determine "k"

   sf contains scaling of rise time with depth
   rf contains scaling of rup. speed with depth
*/

   if(stfparams.rt_scalefac > 0)
      {
      dmin1 = stfparams.risetimedep - stfparams.risetimedep_range;
      dmax1 = stfparams.risetimedep + stfparams.risetimedep_range;
      rtfac1 = stfparams.risetimefac - 1.0;

      dmin2 = stfparams.deep_risetimedep - stfparams.deep_risetimedep_range;
      dmax2 = stfparams.deep_risetimedep + stfparams.deep_risetimedep_range;
      rtfac2 = stfparams.deep_risetimefac - 1.0;

      snum = 0.0;
      sden = 0.0;
      for(j=0;j<ndip;j++)
         {
	 /* scaling by actual vrup */

         k = 0;
         xdep = vmod.th[0];
         while(xdep < psrc[j*nstk].dep && k < vmod.nlay)
            {
            k++;
            xdep = xdep + vmod.th[k];
            }
	 rf = vmod.vs[k];

         if(psrc[j*nstk].dep <= dmin1)
            {
	    rf = rf*rvfrac*shal_vrup;
            sf = 1.0 + rtfac1;
	    }
         else if(psrc[j*nstk].dep < dmax1 && psrc[j*nstk].dep > dmin1)
            {
	    rf = rf*rvfrac*shal_vrup*(dmax1-(psrc[j*nstk].dep))/(dmax1-dmin1);
            sf = 1.0 + rtfac1*(dmax1-(psrc[j*nstk].dep))/(dmax1-dmin1);
	    }
         else if(psrc[j*nstk].dep >= dmax1 && psrc[j*nstk].dep <= dmin2)
            {
	    rf = rf*rvfrac;
            sf = 1.0;
	    }
         else if(psrc[j*nstk].dep < dmax2 && psrc[j*nstk].dep > dmin2)
            {
	    rf = rf*rvfrac*deep_vrup*((psrc[j*nstk].dep)-dmin2)/(dmax2-dmin2);
            sf = 1.0 + rtfac2*((psrc[j*nstk].dep)-dmin2)/(dmax2-dmin2);
	    }
         else
            {
	    rf = rf*rvfrac*deep_vrup;
            sf = 1.0 + rtfac2;
	    }

         for(i=0;i<nstk;i++)
            {
            ip = i + j*nstk;

	    if(svr_wt == 1)
               sabs = rf*sqrt(psrc[ip].slip*psrc[ip].slip);
	    else if(psrc[ip].slip > (float)(0.0))
	       sabs = 1.0;
	    else
	       sabs = 0.0;

            if(sabs > (float)(0.0))
               {
               snum = snum + sf*sabs*sqrt(rtime1_r[ip]/rt1_avg);
               sden = sden + sabs;
               }
            }
         }

      stfparams.rt_scalefac = snum/sden;
      fprintf(stderr,"rt_scalefac= %f\n",stfparams.rt_scalefac);
      }

/* end rtime1 */

/* now do roughness */

   for(j=0;j<ndip3*nstk3;j++)
      {
      rough_c[j].re = 0.0;
      rough_c[j].im = 0.0;
      }

   fft2d(rough_c,nstk3,ndip3,-1,&dstk3,&ddip3);

   hcoef = 1.0;
   kfilt_beta2(rough_c,nstk3,ndip3,&dks3,&dkd3,&hcoef,&lambda_max,&lambda_min,&seed);

   fft2d(rough_c,nstk3,ndip3,1,&dks3,&dkd3);

/*
   scale roughness field to match desired fractal samplitude-to-wavelength ratio 'alpha_rough'
   Shi & Day 2014
*/

   rgh_avg = 0.0;
   for(ip=0;ip<ndip3*nstk3;ip++)
      rgh_avg = rgh_avg + rough_c[ip].re;

   rgh_avg = rgh_avg/(float)(nstk3*ndip3);

/* RWG 2015-09-21
   Prior to 2015-09-21 the scaling used the 'extended' fault length given in the first
   line below.  However, it really should be done using the actual fault length given
   in the 2nd line. Change made to crrect formulation 2015-09-21.
*/
   hrms = alpha_rough*dstk3*nstk3;
   hrms = alpha_rough*flen;
   hrms = alpha_rough*flen_max;

   sf = 0.0;
   for(ip=0;ip<ndip3*nstk3;ip++)
      sf = sf + (rough_c[ip].re - rgh_avg)*(rough_c[ip].re - rgh_avg);

   sf = hrms/sqrt(sf/(float)(nstk3*ndip3));
   ssfac = 0.5*sf/dstk3;
   ddfac = 0.5*sf/ddip3;

   rgh_avg = 0.0;
   for(j=0;j<ndip;j++)
      {
      for(i=0;i<nstk;i++)
         {
	 ip = i + j*nstk;
	 ip2 = (3*i + 1) + (3*j + 1)*nstk3;

         rough_r[ip] = sf*rough_c[ip2].re;
         rgh_avg = rgh_avg + rough_r[ip];

/* compute change in subfault location */

         cosS = cos(rperd*psrc[ip].stk);
         sinS = sin(rperd*psrc[ip].stk);
         cosD = cos(rperd*psrc[ip].dip);
         sinD = sin(rperd*psrc[ip].dip);

	 xx = rough_r[ip]*sinD*sinS;
	 yy = -rough_r[ip]*sinD*cosS;

	 gen_matrices(amat,ainv,&mrot,&(psrc[ip].lon),&(psrc[ip].lat));
	 gcproj(&xx,&yy,&(psrc[ip].lon),&(psrc[ip].lat),&erad,&g0,&b0,amat,ainv,0);

	 psrc[ip].dep = psrc[ip].dep + rough_r[ip]*cosD;

/* compute unit normal to roughness surface in fault plane coordinates-
      ss = along strike (pos. in strike direction)
      dd = down-dip (pos. toward bottom fault edge)
      hh = perpendicular to fault surface (pos. toward footwall)
*/

	 ss = ssfac*(rough_c[ip2+1].re - rough_c[ip2-1].re);            /* df/ds */
	 dd = ddfac*(rough_c[ip2+nstk3].re - rough_c[ip2-nstk3].re);    /* df/dd */
	 hh = -1.0;

	 hhfac = 1.0/sqrt(ss*ss + dd*dd + hh*hh);

	 ss = hhfac*ss;
	 dd = hhfac*dd;
	 hh = hhfac*hh;

/* now, rotate unit normal into N,E,Z coordinates-
      xx = pos. toward North
      yy = pos. toward East
      zz = pos. down
*/

	 xx = ss*cosS - (dd*cosD - hh*sinD)*sinS;
	 yy = ss*sinS + (dd*cosD - hh*sinD)*cosS;
	 zz = dd*sinD + hh*cosD;

/* now, compute strike & dip using unit normal

      tan(stk+90) = yy/xx
      cos(180-dip) = zz

*/

	 psrc[ip].stk = atan2(yy,xx)/rperd - 90.0;
	 psrc[ip].dip = 180.0 - acos(zz)/rperd;

	 while(psrc[ip].stk < 0.0)
	    psrc[ip].stk = psrc[ip].stk + 360.0;
	 while(psrc[ip].stk > 360.0)
	    psrc[ip].stk = psrc[ip].stk - 360.0;
         }
      }
   rgh_avg = rgh_avg/(float)(nstk*ndip);

   rgh_min = 1.0e+15;
   rgh_max = -1.0e+15;
   rgh_sig = 0.0;
   for(ip=0;ip<ndip*nstk;ip++)
      {
      if(rough_r[ip] < rgh_min)
         rgh_min = rough_r[ip];
      if(rough_r[ip] > rgh_max)
         rgh_max = rough_r[ip];

      rgh_sig = rgh_sig + (rough_r[ip]-rgh_avg)*(rough_r[ip]-rgh_avg);
      }
   rgh_sig = sqrt(rgh_sig/(ndip*nstk));

   fprintf(stderr,"rough: avg= %13.5e min= %13.5e max= %13.5e sig= %13.5e\n",rgh_avg,rgh_min,rgh_max,rgh_sig);

   if(roughnessfile[0] != '\0')
      {
      fpr = fopfile(roughnessfile,"w");

      for(j=0;j<ndip;j++)
         {
	 dd = (j+0.5)*ddip;
         for(i=0;i<nstk;i++)
            {
	    ip = i + j*nstk;

	    ss = (i+0.5)*dstk;

            fprintf(fpr,"%13.5e %13.5e %13.5e\n",ss,dd,rough_r[ip]);
            }
         }
      fclose(fpr);
      }

/*
   transform rough_c back to wavenumber domain for later use in correlations
*/
   fft2d(rough_c,nstk3,ndip3,-1,&dstk3,&ddip3);

/* end roughness */

/* now do tsfac2 */

   if(tsfac_rand <= 0.0)
      {
      for(ip=0;ip<nstk*ndip;ip++)
         tsfac2_r[ip] = 0.0;
      }
   else
      {
      for(j=0;j<ndip3*nstk3;j++)
         {
         tsfac2_c[j].re = 0.0;
         tsfac2_c[j].im = 0.0;
         }

      fft2d(tsfac2_c,nstk3,ndip3,-1,&dstk3,&ddip3);

      hcoef = 1.0;
      kfilt_beta2(tsfac2_c,nstk3,ndip3,&dks3,&dkd3,&hcoef,&lambda_max,&lambda_min,&seed);

/*
   correlate with rough_c
*/

      sf = sqrt(1.0 - tsfac2_scor*tsfac2_scor);
      for(j=0;j<ndip3*nstk3;j++)
         {
         tsfac2_c[j].re = tsfac2_scor*rough_c[j].re + sf*tsfac2_c[j].re;
         tsfac2_c[j].im = tsfac2_scor*rough_c[j].im + sf*tsfac2_c[j].im;
         }

      fft2d(tsfac2_c,nstk3,ndip3,1,&dks3,&dkd3);

/*
   compute average for initial field then remove it (should be zero anyway)
*/
      tsf2_avg = 0.0;
      for(j=0;j<ndip;j++)
         {
         for(i=0;i<nstk;i++)
            {
	    ip = i + j*nstk;
	    /*XXXXX*/
	    ip2 = i + j*nstk;
	    ip2 = (3*i + 1) + (3*j + 1)*nstk3;

            tsfac2_r[ip] = tsfac2_c[ip2].re;

            tsf2_avg = tsf2_avg + tsfac2_r[ip];
            }
         }
      tsf2_avg = tsf2_avg/(float)(nstk*ndip);
      fprintf(stderr,"tsf2_avg= %13.5e\n",tsf2_avg);

      for(ip=0;ip<ndip*nstk;ip++)
         tsfac2_r[ip] = tsfac2_r[ip] - tsf2_avg;

      tsf2_avg = 0.0;

/*
   scale to desired sigma
*/
      tsf2_sig = 0.0;
      for(ip=0;ip<ndip*nstk;ip++)
         tsf2_sig = tsf2_sig + (tsfac2_r[ip]-tsf2_avg)*(tsfac2_r[ip]-tsf2_avg);

      tsf2_sig = sqrt(tsf2_sig/(ndip*nstk));

      sigfac = tsfac_rand/tsf2_sig;
      for(ip=0;ip<nstk*ndip;ip++)
         tsfac2_r[ip] = sigfac*(tsfac2_r[ip] - tsf2_avg) + tsf2_avg;

      tsf2_min = 1.0e+15;
      tsf2_max = -1.0e+15;
      tsf2_sig = 0.0;
      for(ip=0;ip<nstk*ndip;ip++)
         {
         if(tsfac2_r[ip] < tsf2_min)
            tsf2_min = tsfac2_r[ip];
         if(tsfac2_r[ip] > tsf2_max)
            tsf2_max = tsfac2_r[ip];

         tsf2_sig = tsf2_sig + (tsfac2_r[ip]-tsf2_avg)*(tsfac2_r[ip]-tsf2_avg);
         }

      tsf2_sig = sqrt(tsf2_sig/(ndip*nstk));

      fprintf(stderr,"tsfac2: avg= %13.5e min= %13.5e max= %13.5e sig= %13.5e\n",tsf2_avg,tsf2_min,tsf2_max,tsf2_sig);
      }
/* end tsfac2 */

/* now do rtime2 */

   if(rtime_rand <= 0.0)
      {
      for(ip=0;ip<nstk*ndip;ip++)
         rtime2_r[ip] = 0.0;
      }
   else
      {
      for(j=0;j<ndip3*nstk3;j++)
         {
         rtime2_c[j].re = 0.0;
         rtime2_c[j].im = 0.0;
         }

      fft2d(rtime2_c,nstk3,ndip3,-1,&dstk3,&ddip3);

      hcoef = 1.0;
      kfilt_beta2(rtime2_c,nstk3,ndip3,&dks3,&dkd3,&hcoef,&lambda_max,&lambda_min,&seed);

/*
   correlate with rough_c
*/

      sf = sqrt(1.0 - rtime2_scor*rtime2_scor);
      for(j=0;j<ndip3*nstk3;j++)
         {
         rtime2_c[j].re = rtime2_scor*rough_c[j].re + sf*rtime2_c[j].re;
         rtime2_c[j].im = rtime2_scor*rough_c[j].im + sf*rtime2_c[j].im;
         }

      fft2d(rtime2_c,nstk3,ndip3,1,&dks3,&dkd3);

/*
   compute average for initial field then remove it (should be zero anyway)
*/
      rt2_avg = 0.0;
      for(j=0;j<ndip;j++)
         {
         for(i=0;i<nstk;i++)
            {
	    ip = i + j*nstk;
	    /*XXXXX*/
	    ip2 = i + j*nstk;
	    ip2 = (3*i + 1) + (3*j + 1)*nstk3;

            rtime2_r[ip] = rtime2_c[ip2].re;

            rt2_avg = rt2_avg + rtime2_r[ip];
            }
         }
      rt2_avg = rt2_avg/(float)(nstk*ndip);
      fprintf(stderr,"rt2_avg= %13.5e\n",rt2_avg);

      for(ip=0;ip<ndip*nstk;ip++)
         rtime2_r[ip] = rtime2_r[ip] - rt2_avg;

      rt2_avg = 0.0;

/*
   scale to desired sigma
*/
      rt2_sig = 0.0;
      for(ip=0;ip<ndip*nstk;ip++)
         rt2_sig = rt2_sig + (rtime2_r[ip]-rt2_avg)*(rtime2_r[ip]-rt2_avg);

      rt2_sig = sqrt(rt2_sig/(ndip*nstk));

      sigfac = rtime_rand/rt2_sig;
      for(ip=0;ip<nstk*ndip;ip++)
         rtime2_r[ip] = sigfac*(rtime2_r[ip] - rt2_avg) + rt2_avg;

      rt2_min = 1.0e+15;
      rt2_max = -1.0e+15;
      rt2_sig = 0.0;
      for(ip=0;ip<nstk*ndip;ip++)
         {
         if(rtime2_r[ip] < rt2_min)
            rt2_min = rtime2_r[ip];
         if(rtime2_r[ip] > rt2_max)
            rt2_max = rtime2_r[ip];

         rt2_sig = rt2_sig + (rtime2_r[ip]-rt2_avg)*(rtime2_r[ip]-rt2_avg);
         }

      rt2_sig = sqrt(rt2_sig/(ndip*nstk));

      fprintf(stderr,"rtime2: avg= %13.5e min= %13.5e max= %13.5e sig= %13.5e\n",rt2_avg,rt2_min,rt2_max,rt2_sig);
      }
/* end rtime2 */

   if(write_srf)
      load_slip_srf_dd4(&srf,&stfparams,psrc,rtime1_r,rtime2_r,&vmod);
      /*
      load_slip_srf_dd3(&srf,&stfparams,psrc,rtime1_r,&vmod);
      */

   for(ih=0;ih<nh;ih++)    /* loop over hypocenter realizations */
      {
      if(random_hypo != 1 && calc_shypo == 1)
         shypo = sh0 + ih*shypo_step - 0.5*flen;

/* calculate rupture time */

/* 20140422 RWG: Deep rupture speed decrease */
      deep_vrup_dep = deep_vrup_dep_saved;
      xhypo = dhypo*sin(avgdip*rperd) + dtop + deep_vrup_deprange;
      if(xhypo > deep_vrup_dep)
         deep_vrup_dep = xhypo;

      shal_vrdep1 = shal_vrup_dep - shal_vrup_deprange;
      shal_vrdep2 = shal_vrup_dep + shal_vrup_deprange;
      deep_vrdep1 = deep_vrup_dep - deep_vrup_deprange;
      deep_vrdep2 = deep_vrup_dep + deep_vrup_deprange;

      if(fdrup_time == 1)
         {
         get_rspeed(&vmod,rspd,psrc,nstk,ndip,&slp_max,&slp_avg,&rvfrac,&shal_vrup,&shal_vrdep1,&shal_vrdep2,&deep_vrup,&deep_vrdep1,&deep_vrdep2,&rvfmin,&rvfmax,fdrup_scale_slip);

	 if(seg_delay == 1)
            get_rsegdelay(rspd,nstk,ndip,&dh,nseg_bounds,gbnd,gwid,rvfac_seg);

         istk_off = 0;
         nstk_fd = nstk;
	 ixs = (int)((shypo+0.5*flen)/dstk + 0.5);

	 if(ixs < (nsring+1))
	    {
	    istk_off = (nsring+1) - ixs;
	    nstk_fd = nstk + istk_off;
            ixs = nsring + 1;
	    }
         else if(ixs > nstk-(nsring+2))
	    {
	    istk_off = 0;
	    nstk_fd = ixs + (nsring + 2);
            ixs = nstk_fd - (nsring + 2);
	    }

         idip_off = 0;
         ndip_fd = ndip;
         iys = (int)(dhypo/ddip + 0.5);

	 if(iys < (nsring+1))
	    {
	    idip_off = (nsring+1) - iys;
	    ndip_fd = ndip + idip_off;
            iys = nsring + 1;
	    }
         else if(iys > ndip-(nsring+2))
	    {
	    idip_off = 0;
	    ndip_fd = iys + (nsring + 2);
            iys = ndip_fd - (nsring + 2);
	    }

fprintf(stderr,"nstk_fd= %d ndip_fd= %d\n",nstk_fd,ndip_fd);
fprintf(stderr,"ixs= %d iys= %d\n",ixs,iys);

         ntot = nstk_fd*ndip_fd;

         rslw = (double *)check_realloc(rslw,ntot*sizeof(double));
         fdrt = (double *)check_realloc(fdrt,ntot*sizeof(double));
         fspace = (double *)check_realloc(fspace,(nstk_fd+ndip_fd)*sizeof(double));
         ispace = (int *)check_realloc(ispace,(nstk_fd+ndip_fd)*sizeof(int));

         get_rslow_stretch(rspd,nstk,ndip,rslw,nstk_fd,ndip_fd,istk_off,idip_off,&rvel_rand,&seed);
	 wfront2d_(&nstk_fd,&ndip_fd,&ixs,&iys,&dh,&nsring,fdrt,rslw,&ntot,fspace,ispace);
	 }
      else
         conv2vrup_dd(&vmod,&rvmod,&avgdip,&dtop,&fwid,&rvfrac,&shal_vrup,&shal_vrdep1,&shal_vrdep2);

      tsmin = 1.0e+15;
      for(j=0;j<ndip;j++)
         {
         yy = (j + 0.5)*ddip;
         for(i=0;i<nstk;i++)
            {
	    ip = i + j*nstk;

	    if(fdrup_time == 0)
	       {
               xx = (i+0.5)*dstk - 0.5*flen;

	       shypo_mseg = shypo;
	       if(seg_delay == 1)
	          {
	          delh = 0.0;
	          hx = xx - shypo;

                  for(ig=0;ig<nseg_bounds;ig++)
                     {
                     hg = (gbnd[ig]-0.5*flen)-shypo;
	             gwid2 = 0.5*gwid[ig];

                     if(hx >= 0.0 && hg >= 0.0 && shypo <= ((gbnd[ig]-0.5*flen)-gwid2))
                        {
                        if(hx > (hg - gwid2) && hx < (hg + gwid2))
                           delh = delh + rvfac_seg[ig]*(hx - (hg - gwid2));
                        else if(hx >= (hg + gwid2))
                           delh = delh + rvfac_seg[ig]*gwid[ig];
                        }

                     else if(hx < 0.0 && hg < 0.0 && shypo >= ((gbnd[ig]-0.5*flen)+gwid2))
                        {
                        if(hx < (hg + gwid2) && hx > (hg - gwid2))
                           delh = delh + rvfac_seg[ig]*(hx - (hg + gwid2));
                        else if(hx <= (hg - gwid2))
                           delh = delh - rvfac_seg[ig]*gwid[ig];
                        }
		     }

	          shypo_mseg = xx - (hx + delh);
	          }

               get_rupt(&rvmod,&htol,&dhypo,&yy,&shypo_mseg,&xx,&rayp,&rupt_rad,&rt);

/* 20140422 RWG:

   Slow rupture for deeper portion of fault, below (deep_vrup_dep-deep_vrup_deprange)
   or dhypo, whichever is deeper.  Rupture time adjustment is approximation for reduced
   rupture velocity.  Should implement travel time calculator in future(?).

*/

               if(psrc[ip].dep >= deep_vrdep1)
                  {
                  y0 = (deep_vrdep1-dtop)/sin(avgdip*rperd);
            
                  r0 = sqrt((xx-shypo_mseg)*(xx-shypo_mseg) + (yy-dhypo)*(yy-dhypo));
                  r1 = r0*(y0-dhypo)/(yy-dhypo);
                  r2 = r0*(yy-y0)/(yy-dhypo);

                  if(psrc[ip].dep < deep_vrdep2)
                     vrfac = 1.0 + (deep_vrup - 1.0)*((psrc[ip].dep)-deep_vrdep1)/(deep_vrdep2-deep_vrdep1);
                  else
                     vrfac = deep_vrup;

                  rt = rt*(r1 + r2/vrfac)/r0;
                  }
               }
	    else
	       rt = fdrt[(i+istk_off) + (j+idip_off)*nstk_fd];

/* 20150320 RWG:

   New way with array of gaussian values

   Add random perturbations to tsfac so that it is not 1:1 correlated with slip.
   Perturbations are log normal with ln(sigma)=tsfac_rand*tsfac.  Default for
   tsfac_rand=0.2.
   
*/

	    tf = tsfac*exp(tsfac2_r[ip]);
            psrc[ip].rupt = rt + tf*tsfac1_r[ip];

	    if(psrc[ip].rupt < tsmin)
	       tsmin = psrc[ip].rupt;
            }
         }

      /* adjust to start at rt=0.0, but only if hypo is in this segment */

      if(shypo < -0.5*flen || shypo > 0.5*flen || dhypo < 0.0 || dhypo > fwid)
         tsmin = 0.0;

      for(j=0;j<nstk*ndip;j++)
         psrc[j].rupt = psrc[j].rupt - tsmin + rupture_delay;

/* RWG 2014-02-20 randomized hypocenter */
      if(write_srf)
	 {
         load_rupt_srf(&srf,psrc,&shypo,&dhypo);

         if(strcmp(outfile,"stdout") == 0)
            sprintf(str,"stdout");

         else if(read_erf == 1)
	    {
            if(random_hypo == 1)
               sprintf(str,"%s-r%.6d.srf",outfile,js);

            else             /* old way */
               sprintf(str,"%s-s%.4d-h%.4d",outfile,js,ih);
	    }

         else
            sprintf(str,"%s.srf",outfile);

         write2srf(&srf,str,outbin);
	 }

      else if(gslip.np > 0 && write_gsf)
         {
         if(strcmp(outfile,"stdout") == 0)
            sprintf(str,"stdout");
         else
            sprintf(str,"%s-s%.4d-h%.4d.gsf",outfile,js,ih);

         write2gsf(&gslip,psrc,infile,str);
	 }
      }

   free_srf_stf(&srf);
   }

if(dump_last_seed == 1)
   {
   fpr = fopfile(seedfile,"w");
   fprintf(fpr,"%lld\n",seed);
   fclose(fpr);
   }
}
