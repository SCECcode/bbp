/**********************************************************************/
/*                                                                    */
/*           wcc_tfilter - written by RWG 03/92                       */
/*                                                                    */
/*           N-th order high-pass, low-pass and band-pass             */
/*           Butterworth filters for time series.                     */
/*                                                                    */
/**********************************************************************/

#include "include.h"
#include "structure.h"
#include "function.h"
#include "getpar.h"

#define 	PI 		3.14159265
#define         MAXFILES        500


char stringbuf[256];

int size_float = sizeof(float);
int size_cx = sizeof(struct complex);

void czero(p,n)
struct complex *p;
int n;
{
float zap = 0.0;
int i;

for(i=0;i<n;i++)
   p[i].re = p[i].im = zap;
}

void hp_filter(q,p,n,alpha,beta,sgn)
struct complex *q, *p;
float *alpha, *beta;
int n, sgn;
{
float are, aim, bre, bim;
float tmpre, tmpim, denom;
float one = 1.0;
int i, n1, k;

tmpre = one - (*alpha);
tmpim = -(*beta);
denom = one/(tmpre*tmpre + tmpim*tmpim);
are = tmpre*denom;
aim = -tmpim*denom;

bre = one + (*alpha);
bim = (*beta);

n1 = n - 1;
if(sgn == 1)
   i = 1;
if(sgn == -1)
   i = n1 - 1;

k = i - sgn;

q[k].re = (are*p[k].re - aim*p[k].im);
q[k].im = (are*p[k].im + aim*p[k].re);

while(n1--)
   {
   tmpre = (p[i].re - p[k].re) + bre*q[k].re - bim*q[k].im;
   tmpim = (p[i].im - p[k].im) + bre*q[k].im + bim*q[k].re;

   q[i].re = are*tmpre - aim*tmpim;
   q[i].im = are*tmpim + aim*tmpre;

   i = i + sgn;
   k = i - sgn;
   }
}

void lp_filter(q,p,n,alpha,beta,sgn)
struct complex *q, *p;
float *alpha, *beta;
int n, sgn;
{
float are, aim, bre, bim, cre, cim;
float tmpre, tmpim, denom;
float one = 1.0;
int i, n1, k;

tmpre = one - (*alpha);
tmpim = -(*beta);
denom = one/(tmpre*tmpre + tmpim*tmpim);
are = tmpre*denom;
aim = -tmpim*denom;

bre = one + (*alpha);
bim = (*beta);

cre = -(*alpha);
cim = -(*beta);

n1 = n - 1;
if(sgn == 1)
   i = 1;
if(sgn == -1)
   i = n1 - 1;

k = i - sgn;

tmpre = cre*p[k].re - cim*p[k].im;
tmpim = cre*p[k].im + cim*p[k].re;

q[k].re = are*tmpre - aim*tmpim;
q[k].im = are*tmpim + aim*tmpre;

while(n1--)
   {
   tmpre = cre*(p[i].re + p[k].re) - cim*(p[i].im + p[k].im)
      + bre*q[k].re - bim*q[k].im;
   tmpim = cre*(p[i].im + p[k].im) + cim*(p[i].re + p[k].re)
      + bre*q[k].im + bim*q[k].re;

   q[i].re = are*tmpre - aim*tmpim;
   q[i].im = are*tmpim + aim*tmpre;

   i = i + sgn;
   k = i - sgn;
   }
}

void set_fullpath(fname,path,name)
char *path, *name, *fname;
{
int j;

j = 0;
while(path[j] != '\0')
   j++;
if(path[j-1] != '/')
   {  
   path[j] = '/';
   path[j+1] = '\0';
   }  

sprintf(fname,"%s%s",path,name);
}

void makedir(path)
char *path;
{
char stmp[256], str[1024];
int rtn, j; 
mode_t mode = 00777;

j = 0;   
while(path[j] != '\0')
   j++;  
 
j--;
while(path[j] == '/')
   j--;
path[j+1] = '\0';
path[j+2] = '\0';
         
j = 0;
while(path[j] != '\0')
   {
   while(path[j] != '/' && path[j] != '\0')
      j++;  
 
   if(j != 0)
      {
      strncpy(stmp,path,j);
      stmp[j] = '\0';
      rtn = mkdir(stmp,mode);
 
      if(rtn == -1)
         {
         if(errno != EEXIST)
            {
            sprintf(str,"makedir() cannot make directory %s, exiting",stmp);
            perror(str);
            exit(-1);
            }
         }
      }
   j++;
   }
}

char *readline(fp)
FILE *fp;
{
char *ptr;
signed char c;

ptr = stringbuf;
while((c = getc(fp)) != EOF)
   {
   if(c == '\n')
      {
      *ptr = '\0';
      return(stringbuf);
      }
   *ptr++ = c;
   }
return(NULL);
}

int getname(str,name,nshft)
char *str, *name;
int *nshft;
{
int nc, nc1;

nc1 = 0;
while(*str == ' ' || *str == '\t') /* remove blank space */
   {
   str++;
   nc1++;
   }
 
if(*str == '"')
   {
   str++;
   nc1++;
   nc = 0;
   while(str[nc] != '"')
      {
      if(str[nc] == '\0')
         break;
 
      nc++;
      }
 
   if(str[nc] == '"')
      nc1++;
   }
else  
   {
   nc = 0;
   while(str[nc] != '\0')
      {
      if(str[nc] == ' ' || str[nc] == '\t')
         break;
 
      nc++;
      }
   }

strncpy(name,str,nc);
name[nc] = '\0';
*nshft = *nshft + nc1 + nc;
nc++;
return(nc);
}

void wcc_tfilter (int param_string_len, char** param_string, float* s1, struct statdata* shead1) {
	struct complex *q, *p, *tmpptr;
	double trig_arg;
	int j, it, order, i;
	int rval;
	float wplo, wphi, cosA, sinA;
	float fhi, flo, fnyq;
	float are, aim;
	float one = 1.0;
	float two = 2.0;
	int phase = 0;

	order = 0;
	fhi = 0.0;
	flo = 1.0e+15;

	setpar(param_string_len, param_string);
	getpar("order","d",&order);
	getpar("fhi","f",&fhi);
	getpar("flo","f",&flo);
	getpar("phase","d",&phase);
	endpar();

	p = NULL;
	q = NULL;
	p = (struct complex *) check_realloc(p,shead1->nt*size_cx);
	q = (struct complex *) check_realloc(q,shead1->nt*size_cx);

	for(it=0;it<shead1->nt;it++)
	{
		p[it].re = s1[it];
		p[it].im = 0.0;
	}

	czero(q,shead1->nt);

	fnyq = one/(two*shead1->dt);

	if(fhi < fnyq)
		wphi = tan(PI*fhi*shead1->dt);
	else
	wphi = 1.0e+10;

	if(flo < fnyq)
		wplo = tan(PI*flo*shead1->dt);

        /*  forward pass  */

	for(j=0;j<order;j++)
	{
		trig_arg = 0.5*(two*j + one)*PI/order;
		cosA = cos(trig_arg);
		sinA = -sin(trig_arg);

		if(flo < fnyq)    /* low-pass filter */
		{
			are = wplo*sinA;
		        aim = -wplo*cosA;
			lp_filter(q,p,shead1->nt,&are,&aim,1);

		        tmpptr = p; p = q; q = tmpptr;
		}
		if(fhi >= 0.0)    /* high-pass filter */
         	{
         	are = wphi*sinA;
         	aim = -wphi*cosA;
         	hp_filter(q,p,shead1->nt,&are,&aim,1);

         	tmpptr = p; p = q; q = tmpptr;
         	}
	}

        /*  reverse pass to obtain zero-phase response  */

	if(phase == 0)
	{
   		for(j=0;j<order;j++)
        	{
        		trig_arg = 0.5*(two*j + one)*PI/order;
         		cosA = cos(trig_arg);
         		sinA = -sin(trig_arg);

         		if(flo < fnyq)    /* low-pass filter */
            		{
            			are = wplo*sinA;
            			aim = -wplo*cosA;
            			lp_filter(q,p,shead1->nt,&are,&aim,-1);

            			tmpptr = p; p = q; q = tmpptr;
            		}

         		if(fhi >= 0.0)    /* high-pass filter */
            		{
            			are = wphi*sinA;
            			aim = -wphi*cosA;
            			hp_filter(q,p,shead1->nt,&are,&aim,-1);

            			tmpptr = p; p = q; q = tmpptr;
            		}
         	}
      }
	for(it=0;it<shead1->nt;it++)
		s1[it] = p[it].re;

	free(p);
	free(q);
}
