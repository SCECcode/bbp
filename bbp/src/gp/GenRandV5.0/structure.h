#include "../StandRupFormat/structure.h"

struct complex
   {
   float re;
   float im;
   };

#ifndef STRUCT_VELMODEL
#define STRUCT_VELMODEL

struct velmodel
   {
   int nlay;
   float *vp;
   double *vs;    /* need double for ray tracing to get ruptime */
   float *den;
   float *th;
   float *dep;
   float *mu;     /* in CMS units */
   double *invb2; /* need double for ray tracing to get ruptime */
   };

#endif

struct stfpar
   {
   int nt;
   float dt;
   float trise;
   float risetimedep;
   float risetimedep_range;
   float risetimefac;
   float rt_scalefac;
   char stype[32];
   };

struct stfpar2
   {
   int nt;
   float dt;
   float trise;
   float trise_orig;		/* start with v5.6.3 store original 'risetime' (rarely used) */
   float alphaT;		/* start with v5.6.3 store alphaT here */
   float risetime_coef_mom;	/* start with v5.6.3 replaces 'risetime_coef' */
   float risetime_coef_savg;	/* start with v5.6.3 get avg rtime from avg slip, now default */
   float risetimedep;
   float risetimedep_range;
   float risetimefac;
   float deep_risetimedep;
   float deep_risetimedep_range;
   float deep_risetimefac;
   float beta_depth;		/* starting with v5.6.1 replaced by 'beta_shal_depth' */
   float beta_depth_range;	/* starting with v5.6.1 replaced by 'beta_shal_depth_range' */
   float beta_shal_depth;
   float beta_shal_depth_range;
   float beta_mid_depth;
   float beta_mid_depth_range;
   float beta_shal;
   float beta_mid;
   float beta_deep;
   float beta_asp;
   float beta_subevt;
   float rtime2slip_exp;
   float rt_scalefac;
   float rt_rand;
   char stype[32];
   };

struct hybrid_correlation_length
   {
   int flag;
   int kmodel;
   float fac;
   float dep;
   float dep_range;
   float side_taper;
   float shal_wt_start;
   float deep_wt_start;
   float shal_wt_end;
   float deep_wt_end;
   };

struct hypo_distr_params
   {
   float x0;
   float x1;
   float f0;
   float f1;
   float xlen;
   float xshift;
   };

struct coord3
   {
   float x1;
   float x2;
   float x3;
   };
