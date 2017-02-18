#include "include.h"
#include "structure.h"
#include "function.h"
#include "defs.h"
#include "getpar.h"

#define         RPERD           0.017453293
#define         DPERR           57.29577951

int main(int ac,char **av)
{
FILE *fopfile(), *fpr, *fpw;
float dx, dy, len, wid, xp, yp;
float dtop, strike, dip;
float rr, ss, tt, tl, ravg, savg, *slip, *trise, *tinit, *sp, *tr, *ti;
float s1, s2, xx, yy, rake, cosA, sinA, sum;
int i, is, id, ip, nstk, ndip, nseg;
int ix, iy, j, k, m, kp, noff;
int nx, ny, nxdiv, nydiv, nxsum, nysum;
float shypo, dhypo, elon, elat, fac;
char infile[256], outfile[256];
char string[1024];

struct standrupformat srf;
int inbin = 0;

float target_dx = -1.0;
float target_dy = -1.0;

float avgstk = -1.0e+15;
int revstk;

sprintf(infile,"stdin");
sprintf(outfile,"stdout");

setpar(ac,av);
getpar("infile","s",infile);
getpar("outfile","s",outfile);

getpar("target_dx","f",&target_dx);
if(target_dx < 0.0)
   mstpar("dx","f",&dx);
getpar("target_dy","f",&target_dy);
if(target_dy < 0.0)
   mstpar("dy","f",&dy);

getpar("inbin","d",&inbin);
getpar("avgstk","f",&avgstk);
endpar();

if(avgstk > -1.0e+14)
   {
   while(avgstk >= 360.0)
      avgstk = avgstk - 360.0;
   while(avgstk < 0.0)
      avgstk = avgstk + 360.0;
   }

read_srf(&srf,infile,inbin);

if(strcmp(outfile,"stdout") == 0)
   fpw = stdout;
else
   fpw = fopfile(outfile,"w");

fprintf(fpw,"%d\n",srf.srf_prect.nseg);

noff = 0;
for(i=0;i<srf.srf_prect.nseg;i++)
   {
   elon = srf.srf_prect.prectseg[i].elon;
   elat = srf.srf_prect.prectseg[i].elat;
   nstk = srf.srf_prect.prectseg[i].nstk;
   ndip = srf.srf_prect.prectseg[i].ndip;
   len = srf.srf_prect.prectseg[i].flen;
   wid = srf.srf_prect.prectseg[i].fwid;

   strike = srf.srf_prect.prectseg[i].stk;
   dip = srf.srf_prect.prectseg[i].dip;
   dtop = srf.srf_prect.prectseg[i].dtop;
   shypo = srf.srf_prect.prectseg[i].shyp;
   dhypo = srf.srf_prect.prectseg[i].dhyp;

   while(strike >= 360.0)
      strike = strike - 360.0;
   while(strike < 0.0)
      strike = strike + 360.0;

   if(target_dx > 0.0)
      {
      nx = (int)(len/target_dx + 0.5);
      dx = len/((float)(nx));
      fprintf(stderr,"target_dx= %f actual dx= %f\n",target_dx,dx);
      }

   nx = (len/dx + 0.5);
   if(nx > nstk)
      {
      nx = nstk;
      dx = len/nx;
      }

   nxdiv = 1;
   while((nstk*nxdiv)%nx)
      nxdiv++;

   nxsum = (nstk*nxdiv)/nx;

   if(target_dy > 0.0)
      {
      ny = (int)(wid/target_dy + 0.5);
      dy = wid/((float)(ny));
      fprintf(stderr,"target_dy= %f actual dy= %f\n",target_dy,dy);
      }

   ny = (wid/dy + 0.5);
   if(ny > ndip)
      {
      ny = ndip;
      dy = wid/ny;
      }

   nydiv = 1;
   while((ndip*nydiv)%ny)
      nydiv++;

   nysum = (ndip*nydiv)/ny;

   fprintf(stderr,"seg= %d\n",i);
   fprintf(stderr,"nstk= %d nx= %d nxdiv= %d nxsum= %d\n",nstk,nx,nxdiv,nxsum);
   fprintf(stderr,"ndip= %d ny= %d nydiv= %d nysum= %d\n",ndip,ny,nydiv,nysum);

   slip = (float *) check_malloc (nxdiv*nydiv*nstk*ndip*sizeof(float));
   trise = (float *) check_malloc (nxdiv*nydiv*nstk*ndip*sizeof(float));
   tinit = (float *) check_malloc (nxdiv*nydiv*nstk*ndip*sizeof(float));
   sp = (float *) check_malloc (nx*ny*sizeof(float));
   tr = (float *) check_malloc (nx*ny*sizeof(float));
   ti = (float *) check_malloc (nx*ny*sizeof(float));

   savg = 0.0;
   ravg = 0.0;
   for(id=0;id<ndip;id++)
      {
      for(is=0;is<nstk;is++)
         {
	 kp = noff + is + id*nstk;

	 tt = srf.srf_apnts.apntvals[kp].tinit;
	 s1 = srf.srf_apnts.apntvals[kp].slip1;
	 s2 = srf.srf_apnts.apntvals[kp].slip2;

	 rake = srf.srf_apnts.apntvals[kp].rake;

	 cosA = cos(rake*RPERD);
	 sinA = sin(rake*RPERD);

	 xx = -s2*sinA + s1*cosA;
	 yy =  s2*cosA + s1*sinA;

	 ss = sqrt(xx*xx + yy*yy);

	 rr = 90;
	 if(yy < 0.0)
	    rr = 270;

	 if(xx != 0.0)
	    {
	    rr = DPERR*atan(yy/xx);
	    if(xx < 0.0)
	       rr = rr + 180;
	    }

	 while(rr < 0.0)
	    rr = rr + 360.0;
	 while(rr > 360.0)
	    rr = rr - 360.0;

	 tl = (srf.srf_apnts.apntvals[kp].dt)*(srf.srf_apnts.apntvals[kp].nt1 - 1);
	 if(srf.srf_apnts.apntvals[kp].nt2 > srf.srf_apnts.apntvals[kp].nt1)
	    tl = (srf.srf_apnts.apntvals[kp].dt)*(srf.srf_apnts.apntvals[kp].nt2 - 1);
	 if(tl < 0.0)
	    tl = 0.0;
	 for(k=0;k<nydiv;k++)
	    {
	    for(j=0;j<nxdiv;j++)
	       {
	       ip = is*nxdiv + j + (id*nydiv + k)*nstk*nxdiv;

	       slip[ip] = ss;
	       trise[ip] = tl;
	       tinit[ip] = tt;
	       }
	    }

         ravg = ravg + rr*ss;
         savg = savg + ss;
         }
      }

   ravg = ravg/savg;
   savg = savg/(nstk*ndip);

   fac = 1.0/(float)(nxsum*nysum);

   for(iy=0;iy<ny;iy++)
      {
      for(ix=0;ix<nx;ix++)
         {
	 ip = ix + iy*nx;

         sp[ip] = 0.0;
         tr[ip] = 0.0;
         ti[ip] = 0.0;
	 sum = 0.0;
	 for(k=0;k<nysum;k++)
	    {
	    for(j=0;j<nxsum;j++)
	       {
	       m = ix*nxsum + j + (iy*nysum + k)*nx*nxsum;

	       sp[ip] = sp[ip] + slip[m];
	       ti[ip] = ti[ip] + tinit[m];
	       tr[ip] = tr[ip] + trise[m]*slip[m];
	       sum = sum + slip[m];
	       }
	    }

         sp[ip] = sp[ip]*fac;
         ti[ip] = ti[ip]*fac;

	 if(sum > 0.0)
            tr[ip] = tr[ip]/sum;
	 else
            tr[ip] = 1.0e-05;
         }
      }

   fprintf(fpw,"%10.4f %10.4f %5d %5d %8.2f %8.2f\n",elon,elat,nx,ny,dx,dy);
   fprintf(fpw,"%4.0f %4.0f %4.0f %8.2f %8.2f %8.2f\n",strike,dip,ravg,dtop,shypo,dhypo);

   revstk = 0;
   if(avgstk > -1.0e+14 && (strike > avgstk + 90.0 || strike < avgstk - 90.0))
      revstk = 1;

   for(iy=0;iy<ny;iy++)
      {
      if(revstk)
         {
         for(ix=nx-1;ix>=0;ix--)
            fprintf(fpw," %5.0f",sp[ix + iy*nx]);
	 }
      else
         {
         for(ix=0;ix<nx;ix++)
            fprintf(fpw," %5.0f",sp[ix + iy*nx]);
	 }

      fprintf(fpw,"\n");
      }

   for(iy=0;iy<ny;iy++)
      {
      if(revstk)
         {
         for(ix=nx-1;ix>=0;ix--)
            fprintf(fpw," %5.2f",tr[ix + iy*nx]);
	 }
      else
         {
         for(ix=0;ix<nx;ix++)
            fprintf(fpw," %5.2f",tr[ix + iy*nx]);
	 }

      fprintf(fpw,"\n");
      }

   for(iy=0;iy<ny;iy++)
      {
      if(revstk)
         {
         for(ix=nx-1;ix>=0;ix--)
            fprintf(fpw," %5.2f",ti[ix + iy*nx]);
	 }
      else
         {
         for(ix=0;ix<nx;ix++)
            fprintf(fpw," %5.2f",ti[ix + iy*nx]);
	 }

      fprintf(fpw,"\n");
      }

   free(slip);
   free(trise);
   free(tinit);
   free(sp);
   free(tr);
   free(ti);

   noff = noff + nstk*ndip;
   }

fclose(fpw);
}
