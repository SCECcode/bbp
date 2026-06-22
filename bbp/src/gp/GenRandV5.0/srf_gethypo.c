#include "include.h"
#include "structure.h"
#include "function.h"
#include "defs.h"
#include "getpar.h"

void read_srf(struct standrupformat *srf,char *file,int bflag);

int main(int ac,char **av)
{
int i, ip, j;
char infile[1024];

struct standrupformat srf1;
struct srf_prectsegments *prseg_ptr1;
struct srf_apointvalues *apval_ptr1;

int inbin = 0;

float hlon, hlat, hdep;
float tmin = 1.0e+15;

/* 2025-12-04
   Added option to check slip-rate function to determine/modify rupture initation time.
   Parameters below are used for these options.  Default is now to use options but can
   be skipped by setting "check_tinit=0".
*/
int check_tinit = 1;
float *stfp, tt;
int it, nt, itmax, it_chk_tt;
double tmp_d, max_sliprate_d, min_sliprate_d;
double tol_d = 1.0e-20;

sprintf(infile,"stdin");

setpar(ac,av);
getpar("infile","s",infile);
getpar("check_tinit","d",&check_tinit);
endpar();

read_srf(&srf1,infile,inbin);

apval_ptr1 = srf1.srf_apnts.apntvals;

for(ip=0;ip<srf1.srf_apnts.np;ip++)
   {

/* 2025-12-04
   Added option to check slip-rate function to determine/modify rupture initation
   time. This is only important if the slip-rate is given with some number of leading zeros.
   Default is now to use this option, but can be skipped by setting "check_tinit=0".
*/

      if(check_tinit)
         {
         tt = 1.0e+15;
         if(apval_ptr1[ip].nt1 > 0 || apval_ptr1[ip].nt2 > 0)
            {
	    tt = apval_ptr1[ip].tinit;

            if(apval_ptr1[ip].slip1*apval_ptr1[ip].slip1 >= apval_ptr1[ip].slip2*apval_ptr1[ip].slip2)
               {
               stfp = apval_ptr1[ip].stf1;
               nt = apval_ptr1[ip].nt1;
               }
            else
               {
               stfp = apval_ptr1[ip].stf2;
               nt = apval_ptr1[ip].nt2;
               }

            max_sliprate_d = -1.0;
            for(it=0;it<nt;it++)
               {
               tmp_d = (double)(stfp[it])*(double)(stfp[it]);
               if(tmp_d > max_sliprate_d)
                  {
                  max_sliprate_d = tmp_d;
                  itmax = it;
                  }
               }

            min_sliprate_d = tol_d*max_sliprate_d;

            it_chk_tt = 1;
            tmp_d = (double)(stfp[it_chk_tt])*(double)(stfp[it_chk_tt]);
            while(tmp_d < min_sliprate_d && it_chk_tt < itmax)
               it_chk_tt++;

            if(tt < (it_chk_tt-1)*apval_ptr1[ip].dt)
               tt = tt + (it_chk_tt-1)*apval_ptr1[ip].dt;
            }
         } /* end "check_tinit" */
   else
      tt = apval_ptr1[ip].tinit;

   if(tt < tmin)
      {
      hlon = apval_ptr1[ip].lon;
      hlat = apval_ptr1[ip].lat;
      hdep = apval_ptr1[ip].dep;
      tmin = tt;
      }
   }

fprintf(stdout,"%.5f\t%.5f\t%.5f\t%.5f\n",hlon,hlat,hdep,tmin);
}
