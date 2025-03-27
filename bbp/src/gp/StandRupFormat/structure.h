//for compatibility with params.h
#define LV 10
#define NQ 2000
#define NP 500

//This structure contains the output from srf2stoch
struct slipfile
{
  int nseg;
  float elon[LV];
  float elat[LV];
  int nx[LV];
  int ny[LV];
  float dx[LV];
  float dy[LV];
  float strike[LV];
  float dip[LV];
  float ravg[LV];
  float dtop[LV];
  float dhypo[LV];
  float shypo[LV];
  float* sp;
  float* tr;
  float* ti;
  float qfexp;
};

struct stfparOLD
   {
   int nt;
   float dt;
   float trise;
   };

struct pointsource20171024
   {
   float lon;
   float lat;
   float dep;
   float stk;
   float dip;
   float rak;
   float area;
   float slip;
   float rupt;
   };

struct pointsource
   {
   float lon;
   float lat;
   float dep;
   float stk;
   float dip;
   float rak;
   float area;
   float slip;
   float rupt;
   float vs;
   float den;
   float mu;
   float beta;
   float asp;
   int asp_mask;
   float subevt;
   int subevt_mask;
   float aseis;
   float rvf;
   };

struct srf_apointvalues20141109
   {
   float lon;
   float lat;
   float dep;
   float stk;
   float dip;
   float area;
   float tinit;
   float dt;
   float rake;
   float slip1;
   int nt1;
   float slip2;
   int nt2;
   float slip3;
   int nt3;
   float *stf1;
   float *stf2;
   float *stf3;
   };

struct srf_apointvalues20241010
   {
   float lon;
   float lat;
   float dep;
   float stk;
   float dip;
   float area;
   float tinit;
   float dt;
   float vs;			/* added for V2.0 */
   float den;			/* added for V2.0 */
   float rake;
   float slip1;
   int nt1;
   float slip2;
   int nt2;
   float slip3;
   int nt3;
   float *stf1;
   float *stf2;
   float *stf3;
   };

struct srf_apointvalues20250127
   {
   float lon;
   float lat;
   float dep;
   float stk;
   float dip;
   float area;
   float tinit;
   float dt;
   float vp;			/* added for V3.0 */
   float vs;
   float den;
   float rake;
   float slip1;
   int nt1;
   float slip2;
   int nt2;
   float slip3;
   int nt3;
   float *stf1;
   float *stf2;
   float *stf3;
   double mnn;			/* moment-rate stuff added for V3.0 */
   double mee;
   double mdd;
   double mne;
   double mnd;
   double med;
   int ntmr;
   float *mrf;
   };

struct srf_apointvalues
   {
   float lon;
   float lat;
   float dep;
   float stk;
   float dip;
   float area;
   float tinit;
   float dt;
   float vp;			/* added for V3.0 */
   float vs;
   float den;
   float rake;
   float slip1;
   int nt1;
   float slip2;
   int nt2;
   float slip3;
   int nt3;
   float *stf1;
   float *stf2;
   float *stf3;
   double mnn;				/* moment-rate stuff added for V3.0 */
   double mee;
   double mdd;
   double mne;
   double mnd;
   double med;
   int ntmr;
   float *mrf;
   float *mr_nn;
   float *mr_ee;
   float *mr_dd;
   float *mr_ne;
   float *mr_nd;
   float *mr_ed;
   };

struct srf_allpoints
   {
   int np;
   struct srf_apointvalues *apntvals;
   };

struct srf_prectsegments
   {
   float elon;
   float elat;
   int nstk;
   int ndip;
   float flen;
   float fwid;
   float dlen;
   float dwid;
   float stk;
   float dip;
   float dtop;
   float shyp;
   float dhyp;
   };

struct srf_planerectangle
   {
   int nseg;
   struct srf_prectsegments *prectseg;
   };

struct srf_headercomment
   {
   int nline;
   char *cbuf;
   };

struct standrupformat20141109
   {
   char version[32];
   char type[32];
   struct srf_planerectangle srf_prect;
   struct srf_allpoints srf_apnts;
   };

struct standrupformat20241010
   {
   char version[32];
   char type[32];
   int nseg;			/* added for V2.0 */
   int *np_seg;			/* added for V2.0 */
   struct srf_headercomment srf_hcmnt;		/* added for V2.0 */
   struct srf_planerectangle srf_prect;
   struct srf_allpoints srf_apnts;
   };

struct standrupformat
   {
   char version[32];
   char type[32];

/*
   'src_format' added for V3.0, options are -
         "SLIP" : traditional slip on fault surface [default]
         "MOMENT" (or "MOMENT-1MECH") : moment-tensor with single time function
         "MOMENT-6MECH" : moment-tensor with separate time function for each component
*/
   char src_format[64];

   int nseg;
   int *np_seg;
   struct srf_headercomment srf_hcmnt;
   struct srf_planerectangle srf_prect;
   struct srf_allpoints srf_apnts;
   };

struct slippars20171024
   {
   float lon;
   float lat;
   float dep;
   float ds;
   float dw;
   float stk;
   float dip;
   float rake;
   float slip;
   float tinit;
   int segno;
   float trise;
   float tau1_ratio;
   };

struct slippars
   {
   float lon;
   float lat;
   float dep;
   float ds;
   float dw;
   float stk;
   float dip;
   float rake;
   float slip;
   float tinit;
   int segno;
   float trise;
   float tau1_ratio;
   float vs;
   float den;
   float rt_s;
   float rt_e;
   float aseis;
   float rvf;
   };

struct segpars
   {
   int nstk;
   int ndip;
   float dstk;
   float ddip;
   };

struct generic_slip
   {
   int np;
   struct slippars *spar;
   struct segpars *gpar;
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

struct vruppars
   {
   float rvfrac;
   float shal_vrup;
   float shal_vrup_dep;
   float shal_vrup_deprange;
   float deep_vrup;
   float deep_vrup_dep;
   float deep_vrup_deprange;
   };
