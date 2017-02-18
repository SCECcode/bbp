#include "include.h"
#include "structure.h"
#include "function.h"

#define MAXFILE 100

char infilebuf[MAXFILE*256];

main(int ac,char **av)
{
FILE *fopfile(), *fpr, *fpw;
float dx, dy, *len, *wid, *strike, *dip;
float xlon, xlat, xdep, xrak, xslp, tdel;
int nfile, i, j, k, ip, nstk_tot, *nstk, *ndip, poff;
int nseg, ig;
char filelist[256], *infile[MAXFILE], outfile[256], segmentfile[256];
char string[1024];

struct generic_slip gslip;
struct slippars *spar, *spar2;

int reverse_strike = 0;

setpar(ac,av);
mstpar("filelist","s",filelist);
mstpar("outfile","s",outfile);
mstpar("segmentfile","s",segmentfile);
mstpar("tdel","f",&tdel);
getpar("reverse_strike","d",&reverse_strike);
endpar();

fpr = fopfile(filelist,"r");
nfile = 0;
while(fgets(string,1024,fpr) != NULL)
   {
   if(nfile == MAXFILE)
      {
      fprintf(stderr,"MAXFILE=%d exceeded, exiting...\n",MAXFILE);
      exit(-1);
      }

   infile[nfile] = infilebuf + nfile*256;
   sscanf(string,"%s",infile[nfile]);
   nfile++;
   }
fclose(fpr);

fpr = fopfile(segmentfile,"r");

fgets(string,1024,fpr);
sscanf(string,"%d",&nseg);

len = (float *)check_malloc(nseg*sizeof(float));
wid = (float *)check_malloc(nseg*sizeof(float));
strike = (float *)check_malloc(nseg*sizeof(float));
dip = (float *)check_malloc(nseg*sizeof(float));
nstk = (int *)check_malloc(nseg*sizeof(int));
ndip = (int *)check_malloc(nseg*sizeof(int));

nstk_tot = 0;
for(i=0;i<nseg;i++)
   {
   fgets(string,1024,fpr);
   sscanf(string,"%f %f %d %d %f %f",&len[i],&wid[i],&nstk[i],&ndip[i],&strike[i],&dip[i]);

   nstk_tot = nstk_tot + nstk[i];
   }

fclose(fpr);

gslip.np = nstk_tot*ndip[0];
gslip.spar = (struct slippars *)check_malloc(gslip.np*sizeof(struct slippars));
spar = gslip.spar;

for(ip=0;ip<gslip.np;ip++)
   {
   spar[ip].rake = 0.0;
   spar[ip].slip = -1;
   spar[ip].tinit = 0.0;
   }

for(k=0;k<nfile;k++)
   {
   fpr = fopfile(infile[k],"r");

   poff = 0;
   for(ig=0;ig<nseg;ig++)
      {
      dx = len[ig]/nstk[ig];
      dy = wid[ig]/ndip[ig];
      for(i=0;i<nstk[ig];i++)
         {
         for(j=0;j<ndip[ig];j++)
            {
            fgets(string,1024,fpr);
            sscanf(string,"%f %f %f %f %f",&xlon,&xlat,&xdep,&xrak,&xslp);

	    ip = poff + i + j*nstk_tot;

	    spar[ip].lon = xlon;
	    spar[ip].lat = xlat;
	    spar[ip].dep = xdep;
	    spar[ip].ds = dx;
	    spar[ip].dw = dy;
	    spar[ip].stk = strike[ig];
	    spar[ip].dip = dip[ig];
	    spar[ip].segno = ig;

	    if(xslp > 0.0)
	       {
	       if(spar[ip].slip < 0.0)
	          {
	          spar[ip].rake = xrak*xslp;
	          spar[ip].slip = xslp;
	          spar[ip].tinit = k*tdel;
	          }
	       else
	          {
	          spar[ip].rake = spar[ip].rake + xrak*xslp;
	          spar[ip].slip = spar[ip].slip + xslp;
	          }
	       }
            }
         }

      poff = poff + nstk[ig];
      }

   fclose(fpr);
   }

if(reverse_strike)
   {
   spar2 = (struct slippars *)check_malloc(gslip.np*sizeof(struct slippars));
   for(j=0;j<ndip[0];j++)
      {
      for(i=0;i<nstk_tot;i++)
         {
	 k = i + j*nstk_tot;
	 ip = ((nstk_tot-1) - i) + j*nstk_tot;

	 spar2[ip].lon = spar[k].lon;
	 spar2[ip].lat = spar[k].lat;
	 spar2[ip].dep = spar[k].dep;
	 spar2[ip].ds = spar[k].ds;
	 spar2[ip].dw = spar[k].dw;
	 spar2[ip].stk = spar[k].stk;
	 spar2[ip].dip = spar[k].dip;
	 spar2[ip].rake = spar[k].rake;
	 spar2[ip].slip = spar[k].slip;
	 spar2[ip].tinit = spar[k].tinit;
	 spar2[ip].segno = (nseg-1) - spar[k].segno;
	 }
      }
   spar = spar2;
   }

fpw = fopfile(outfile,"w");

fprintf(fpw,"%d\n",gslip.np);

for(ip=0;ip<gslip.np;ip++)
   {
   if(spar[ip].slip > 0.0)
      spar[ip].rake = spar[ip].rake/spar[ip].slip;

   fprintf(fpw,"%11.5f %11.5f %8.4f %8.4f %8.4f %6.1f %6.1f %6.1f %8.2f %8.3f %3d\n",
                                                        spar[ip].lon,
                                                        spar[ip].lat,
                                                        spar[ip].dep,
                                                        spar[ip].ds,
                                                        spar[ip].dw,
                                                        spar[ip].stk,
                                                        spar[ip].dip,
                                                        spar[ip].rake,
                                                        spar[ip].slip,
                                                        spar[ip].tinit,
                                                        spar[ip].segno);
   }

fclose(fpw);
}
