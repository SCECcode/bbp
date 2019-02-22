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

struct standrupformat
   {
   char version[32];
   char type[32];
   int nseg;			/* added for V2.0 */
   int *np_seg;			/* added for V2.0 */
   struct srf_headercomment srf_hcmnt;		/* added for V2.0 */
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
   };

struct generic_slip
   {
   int np;
   struct slippars *spar;
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
