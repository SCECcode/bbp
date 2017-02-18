#define FD_STATCHAR 8
#define STATCHAR 12
#define COMPCHAR 4
#define TITLCHAR 64

struct statdata
   {
   char stat[STATCHAR];
   char comp[COMPCHAR];
   char stitle[TITLCHAR];
   int nt;
   float dt;
   int hr;
   int min;
   float sec;
   float edist;
   float az;
   float baz;
   };

struct mtheader    /* header for moment tensor output information */
   {
   char title[128];
   int indx;
   int nt;
   float dt;
   float h;
   float modelrot;   /* rotation of y-axis from south (clockwise positive) */
   float modellat;
   float modellon;
   float mu;
   int xsrc;
   int ysrc;
   int zsrc;
   int xsta;
   int ysta;
   int zsta;
   float xmom;
   float ymom;
   float zmom;
   };

struct seisheader
   {
   int indx;        /* numerical index of this location in statcords file */
   int ix;          /* x grid location for output */
   int iy;          /* y grid location for output */
   int iz;          /* z grid location for output */
   int nt;          /* number of time points                                */
   float dt;        /* time sampling                                        */
   float h;         /* grid spacing                                         */
   float modelrot;  /* rotation of y-axis from south (clockwise positive)   */
   float modellat;  /* latitude of model origin                             */
   float modellon;  /* longitude of model origin                            */
   char name[FD_STATCHAR]; /* station name */
   };

struct tsheader    /* structure for time slice header information */
   {
   int ix0;          /* starting x grid location for output */
   int iy0;          /* starting y grid location for output */
   int iz0;          /* starting z grid location for output */
   int it0;          /* starting time step for output */
   int nx;          /* number of x points                                */
   int ny;          /* number of y points                                */
   int nz;          /* number of z points                                */
   int nt;          /* number of time points                                */
   float dx;         /* grid spacing                                         */
   float dy;         /* grid spacing                                         */
   float dz;         /* grid spacing                                         */
   float dt;        /* time sampling                                        */
   float modelrot;  /* rotation of y-axis from south (clockwise positive)   */
   float modellat;  /* latitude of model origin                             */
   float modellon;  /* longitude of model origin                            */
   };

struct tsheader_proc  /* structure for individual processor time slice header */
   {
   int ix0;          /* starting x grid location for output */
   int iy0;          /* starting y grid location for output */
   int iz0;          /* starting z grid location for output */
   int it0;          /* starting time step for output */
   int iyleft;          /* global iy of 1st plane in this file */
   int iyright;          /* global iy of last plane in this file */
   int localny;          /* local # of iy planes */
   int nx;          /* number of x points                                */
   int ny;          /* number of y points                                */
   int nz;          /* number of z points                                */
   int nt;          /* number of time points                                */
   float dx;         /* grid spacing                                         */
   float dy;         /* grid spacing                                         */
   float dz;         /* grid spacing                                         */
   float dt;        /* time sampling                                        */
   float modelrot;  /* rotation of y-axis from south (clockwise positive)   */
   float modellat;  /* latitude of model origin                             */
   float modellon;  /* longitude of model origin                            */
   };
