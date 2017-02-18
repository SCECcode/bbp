struct complex
   {
   float re;
   float im;
   };

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

struct mechparam
   {
   int nmech;
   int flag[3];
   float stk;
   float dip;
   float rak;
   };

struct gfheader
   {
   int read_flag;
   int nt;
   int id;
   int ir;
   float olon;
   float olat;
   float slon;
   float slat;
   float dt;
   float north;
   float east;
   float rng;
   float dep;
   float xazim;
   float tst;
   float lam;
   float mu;
   float rho;
   float gft;
   float mom;
   float xmom;
   float ymom;
   float zmom;
   float wt;
   };

struct gfparam
   {
   int flag3d;
   int swap_flag;
   int nc;
   char gftype[16];
   char gflocs[128];
   char gftimes[128];
   float rtol;
   int ngfr;
   int ngfd;
   float *gfn;
   float *gfe;
   float *gfr;
   float *gfd;
   float *gft;
   };

struct beroza
   {
   int npstk;
   int npdip;
   int inc_stk;
   int inc_dip;
   int robstf;
   float generic_risetime;
   float generic_pulsedur;
   float generic_t2;
   float *as;
   float *dd;
   float *slip;
   float *sv;
   float *rupt;
   float *rist;
   float *tdur;
   float *t2;
   };

struct okumura
   {
   int nstk;
   int ndip;
   float flen;
   float fwid;
   float dlen;
   float dwid;
   float shypo;
   float dhypo;
   float vrup;
   float *as;
   float *dd;
   float *slip;
   float *sv;
   float *rist;
   float *rupt;
   };

struct gene
   {
   int nstk;
   int ndip;
   float flen;
   float fwid;
   float dlen;
   float dwid;
   float shypo;
   float dhypo;
   float trise;
   float tdel;
   int maxtw;
   float *as;
   float *dd;
   float *slip;
   float *rake;
   float *vrup;
   int *nt;
   float *swgt;
   };

struct rob
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
   float *as;
   float *dd;
   float *slip;
   float *rake;
   float *trise;
   float *vrup;
   float *tsfac;
   };
