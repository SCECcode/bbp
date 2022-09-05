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
   float risetimedep;
   float risetimedep_range;
   float risetimefac;
   float deep_risetimedep;
   float deep_risetimedep_range;
   float deep_risetimefac;
   float beta_depth;
   float beta_depth_range;
   float beta_shal;
   float beta_deep;
   float rt_scalefac;
   float rt_rand;
   char stype[32];
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
