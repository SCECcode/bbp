#include "include.h"
#include "structure.h"
#include "function.h"
#include "getpar.h"

double frand(void);
double sfrand(long *);

void sum(s1,h1,f1,t1,s2,h2,f2,t2,p,hd)
struct statdata *h1, *h2, *hd;
float *s1, *s2, *p, *f1, *f2, *t1, *t2;
{
float time1, time2;
int it, n;

hd->dt = h1->dt;
hd->edist = h1->edist;
hd->az = h1->az;
hd->baz = h1->baz;

/* set output start time to be earliest start time of two input files */

time1 = h1->hr*3600.0 + h1->min*60.0 + h1->sec + *t1;
time2 = h2->hr*3600.0 + h2->min*60.0 + h2->sec + *t2;

n = 0;
if(time1 < time2)
   {
   hd->hr = h1->hr;
   hd->min = h1->min;
   hd->sec = h1->sec + *t1;

   n = ((time2-time1)/hd->dt + 0.5);

   for(it=0;it<n;it++) {
      p[it] = (*f1)*s1[it];
	}

   h1->nt = h1->nt - n;
   s1 = s1 + n;
   p = p + n;
   }
else
   {
   hd->hr = h2->hr;
   hd->min = h2->min;
   hd->sec = h2->sec + *t2;

   n = ((time1-time2)/hd->dt + 0.5);

   for(it=0;it<n;it++) {
      p[it] = (*f2)*s2[it];
   }

   h2->nt = h2->nt - n;;
   s2 = s2 + n;
   p = p + n;
   }

if(h1->nt > h2->nt)
   {
   hd->nt = h1->nt + n;

   for(it=0;it<h2->nt;it++)
      p[it] = (*f1)*s1[it] + (*f2)*s2[it];

   for(it=h2->nt;it<h1->nt;it++)
      p[it] = (*f1)*s1[it];
   }
else
   {
   hd->nt = h2->nt + n;

   for(it=0;it<h1->nt;it++)
      p[it] = (*f1)*s1[it] + (*f2)*s2[it];

   for(it=h1->nt;it<h2->nt;it++)
      p[it] = (*f2)*s2[it];
   }
}

static  long    frandx = 1;

/* frand() returns a uniform distribution of random numbers
 * in the range -1.0 -> 1.0.
 */
double frand(void)
{
frandx = (frandx * 1103515245 + 12345) & 0x7fffffff;
return((double)(frandx)/1073741824.0 - 1.0);
}

/* sfrand() returns a uniform distribution of random numbers
 * in the range -1.0 -> 1.0.
 */
double sfrand(long *seed)
{
*seed = ((*seed) * 1103515245 + 12345) & 0x7fffffff;
return((double)(*seed)/1073741824.0 - 1.0);
}

void wcc_add(int param_string_len, char** param_string, float* s1, struct statdata* shead1, float* s2, struct statdata* shead2, float* p, struct statdata* shead3) {
	float f1 = 1.0;
	float f2 = 1.0;

	float t1 = 0.0;
	float t2 = 0.0;

	int it;
	float add_rand = 0.0;

	setpar(param_string_len,param_string);
	getpar("f1","f",&f1);
	getpar("t1","f",&t1);
	getpar("f2","f",&f2);
	getpar("t2","f",&t2);
	getpar("add_rand","f",&add_rand);
	endpar();

	sum(s1,shead1,&f1,&t1,s2,shead2,&f2,&t2,p,shead3);

	if(add_rand > (float)(0.0))
	   {
	   for(it=0;it<shead3->nt;it++)
	      p[it] = p[it] + add_rand*frand();
	   }
}
