#include "include.h"
#include "structure.h"
#include "function.h"
#include "getpar.h"

int size_float = sizeof(float);
int size_int = sizeof(int);

int main(ac,av)
int ac;
char **av;
{
struct statdata *hptr, *shead1, *shead2;
float *fptr, *h1, *h2, *r1, *r2;
float ang1, ang2, ang, t1, t2, tt;
int nt, nts;
float rot = -999.9;

char datafile[128];
char filein1[128];
char filein2[128];
char fileout1[128];
char fileout2[128];

int inbin1 = 0;
int inbin2 = 0;
int outbin1 = 0;
int outbin2 = 0;

float pi = 3.14159265;

setpar(ac, av);
mstpar("filein1","s",filein1);
mstpar("filein2","s",filein2);
mstpar("fileout1","s",fileout1);
mstpar("fileout2","s",fileout2);
getpar("rot","f",&rot);
getpar("inbin1","d",&inbin1);
getpar("inbin2","d",&inbin2);
getpar("outbin1","d",&outbin1);
getpar("outbin2","d",&outbin2);
endpar();

shead1 = (struct statdata *) check_malloc (sizeof(struct statdata));
shead2 = (struct statdata *) check_malloc (sizeof(struct statdata));

h1 = NULL;
h1 = read_wccseis(filein1,shead1,h1,inbin1);
h2 = NULL;
h2 = read_wccseis(filein2,shead2,h2,inbin2);

if(shead1->comp[0] == 'n' || shead1->comp[0] == 'N')
   ang1 = 0.0;
else if(shead1->comp[0] == 'e'  || shead1->comp[0] == 'E' || strcmp(shead1->comp,"vx") == 0)
   ang1 = 90.0;
else if(shead1->comp[0] == 's'  || shead1->comp[0] == 'S' || strcmp(shead1->comp,"vy") == 0)
   ang1 = 180.0;
else if(shead1->comp[0] == 'w' || shead1->comp[0] == 'W')
   ang1 = 270.0;
else
   ang1 = atof(shead1->comp);

if(shead2->comp[0] == 'n' || shead2->comp[0] == 'N')
   ang2 = 0.0;
else if(shead2->comp[0] == 'e' || shead2->comp[0] == 'E' || strcmp(shead2->comp,"vx") == 0)
   ang2 = 90.0;
else if(shead2->comp[0] == 's' || shead2->comp[0] == 'S' || strcmp(shead2->comp,"vy") == 0)
   ang2 = 180.0;
else if(shead2->comp[0] == 'w' || shead2->comp[0] == 'W')
   ang2 = 270.0;
else
   ang2 = atof(shead2->comp);

t1 = shead1->hr*3600.0 + shead1->min*60.0 + shead1->sec;
t2 = shead2->hr*3600.0 + shead2->min*60.0 + shead2->sec;

if(t1 < t2)
   {
   nts = (t2-t1)/shead1->dt;
   h1 = h1 + nts;
   shead1->nt = shead1->nt - nts;

   shead1->hr = shead2->hr;
   shead1->min = shead2->min;
   shead1->sec = shead2->sec;
   }
else if(t1 > t2)
   {
   nts = (t1-t2)/shead1->dt;
   h2 = h2 + nts;
   shead2->nt = shead2->nt - nts;

   shead2->hr = shead1->hr;
   shead2->min = shead1->min;
   shead2->sec = shead1->sec;
   }

if(shead1->nt < shead2->nt)
   shead2->nt = shead1->nt;
if(shead1->nt > shead2->nt)
   shead1->nt = shead2->nt;

nt = shead1->nt;

r1 = (float *) check_malloc (nt*size_float);
r2 = (float *) check_malloc (nt*size_float);

while(ang1 > 360.0)
   ang1 = ang1 - 360.0;
while(ang2 > 360.0)
   ang2 = ang2 - 360.0;

if(ang2 < ang1 && ang2+270 != ang1)
   {
   fptr = h1; h1 = h2; h2 = fptr;
   hptr = shead1; shead1 = shead2; shead2 = hptr;
   ang = ang1; ang1 = ang2; ang2 = ang;
   }

if(ang2 > ang1 && ang1+90 != ang2)
   {
   fptr = h1; h1 = h2; h2 = fptr;
   hptr = shead1; shead1 = shead2; shead2 = hptr;
   ang = ang1; ang1 = ang2; ang2 = ang;
   }

if(ang2 != ang1+90.0 && ang2 != ang1-270.0)
   {
   fprintf(stderr,"ang1= %f ang2= %f\n",ang1,ang2);
   fprintf(stderr,"Input components differ by more than 90 deg., exiting...\n");
   exit(-1);
   }

if(rot < -999.0)
   {
   fprintf(stderr,"Using back-azimuth for rotation angle.\n");
   rot = shead1->baz - 180.0;
   }

ang = (rot-ang1)*pi/180.0;
rotate(nt,r1,r2,h1,h2,&ang);

ang = rot;
while(ang >= 360.0)
   ang -= 360.0;
while(ang < 0.0)
   ang += 360.0;

sprintf(shead1->comp,"%d",(int)(ang));
write_wccseis(fileout1,shead1,r1,outbin1);

ang = 90.0 + ang;
while(ang > 360.0)
   ang -= 360.0;
while(ang < 0.0)
   ang += 360.0;

sprintf(shead2->comp,"%d",(int)(ang));
write_wccseis(fileout2,shead2,r2,outbin2);
}

void rotate(n,r,t,north,east,a)
float *r, *t, *north, *east, *a;
int n;
{
float cosA, sinA;

cosA = cos(*a);
sinA = sin(*a);

while(n--)
   {
   t[0] = east[0]*cosA - north[0]*sinA;
   r[0] = east[0]*sinA + north[0]*cosA;

   r++; t++;
   north++; east++;
   }
}
