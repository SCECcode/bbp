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

struct sgtfileparams
   {
   int xfdr;
   int yfdr;
   int zfdr;
   char xfile[256];
   char yfile[256];
   char zfile[256];
   off_t head_off;
   off_t xcur_off;
   off_t ycur_off;
   off_t zcur_off;
   };

struct sgtparams
   {
   int nsgt;
   long long indx[4];
   float wt[4];
   int master_ip[4];
   };

struct sgtmaster
   {
   int geoproj;     /* =0: RWG local flat earth; =1: RWG great circle arcs; =2: UTM */
   float modellon;  /* longitude of geographic origin */
   float modellat;  /* latitude of geographic origin */
   float modelrot;  /* rotation of y-axis from south (clockwise positive)   */
   float xshift;    /* xshift of cartesian origin from geographic origin */
   float yshift;    /* yshift of cartesian origin from geographic origin */
   int globnp;      /* total number of SGT locations (entire model) */
   int localnp;     /* local number of SGT locations (this file only) */
   int nt;          /* number of time points                                */
   };

struct sgtindex   /* indices for all 'globnp' SGT locations */
   {
   long long indx; /* indx= xsgt*10000000 + ysgt*1000 + zsgt */
   int xsgt;     /* x grid location */
   int ysgt;     /* y grid location */
   int zsgt;     /* z grid location */
   float h;         /* grid spacing                                         */
   };

struct sgtheader    /* sgt header for v1.14 */
   {
   long long indx;  /* index of this SGT */
   int geoproj;     /* =0: RWG local flat earth; =1: RWG great circle arcs; =2: UTM */
   float modellon;  /* longitude of geographic origin */
   float modellat;  /* latitude of geographic origin */
   float modelrot;  /* rotation of y-axis from south (clockwise positive)   */
   float xshift;    /* xshift of cartesian origin from geographic origin */
   float yshift;    /* yshift of cartesian origin from geographic origin */
   int nt;          /* number of time points                                */
   float xazim;     /* azimuth of X-axis in FD model (clockwise from north) */
   float dt;        /* time sampling                                        */
   float tst;       /* start time of 1st point in GF                        */
   float h;         /* grid spacing                                         */
   float src_lat;   /* site latitude */
   float src_lon;   /* site longitude */
   float src_dep;   /* site depth */
   int xsrc;        /* x grid location for source (station in recip. exp.)  */
   int ysrc;        /* y grid location for source (station in recip. exp.)  */
   int zsrc;        /* z grid location for source (station in recip. exp.)  */
   float sgt_lat;   /* SGT location latitude */
   float sgt_lon;   /* SGT location longitude */
   float sgt_dep;   /* SGT location depth */
   int xsgt;        /* x grid location for output (source in recip. exp.)   */
   int ysgt;        /* y grid location for output (source in recip. exp.)   */
   int zsgt;        /* z grid location for output (source in recip. exp.)   */
   float cdist;     /* straight-line distance btw site and SGT location */
   float lam;       /* lambda [in dyne/(cm*cm)] at output point             */
   float mu;        /* rigidity [in dyne/(cm*cm)] at output point           */
   float rho;       /* density [in gm/(cm*cm*cm)] at output point           */
   float xmom;      /* moment strength of x-oriented force in this run      */
   float ymom;      /* moment strength of y-oriented force in this run      */
   float zmom;      /* moment strength of z-oriented force in this run      */
   };

struct gfparam
   {
   int flag3d;
   int swap_flag;
   int use_depdir;
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

struct geoprojection
   {
   int geoproj;
   float modellon;
   float modellat;
   float modelrot;
   float xshift;
   float yshift;
   int center_origin;
   double rperd;
   float erad;
   float fc;
   float g2;
   float radc;
   float cosR;
   float sinR;
   float kmlon;
   float kmlat;
   double g0;
   double b0;
   double amat[9];
   double ainv[9];
   };
