#define STATCHAR 8
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
   float modelrot;   /* rotation of y-axis from south (clockwise positive) */
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
