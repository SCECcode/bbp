struct velmodel
   {
   int nlay;
   double *vs;
   double *invb2;
   float *th;
   };

struct complex
   {
   float re;
   float im;
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

struct srfseg
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
   float *tinit;
   int *nt;
   float *dt;
   float *stfbuf;
   float **stf;
   };

struct standrupform
   {
   char version[32];
   int nseg;
   struct srfseg *fseg;
   };

struct srf_apointvaluesOLD
   {
   float as;
   float dd;
   float lon;
   float lat;
   float dep;
   float mu;
   float stk;
   float dip;
   float rake;
   float area;
   float slip;
   float tinit;
   int nt;
   float dt;
   float *stf;
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

struct standrupformat
   {
   char version[32];
   char type[32];
   struct srf_planerectangle srf_prect;
   struct srf_allpoints srf_apnts;
   };
