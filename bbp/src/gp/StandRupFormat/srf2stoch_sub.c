#include "include.h"
#include "structure.h"
#include "function.h"
#include "defs.h"
#include "getpar.h"

#define         RPERD           0.017453293
#define         DPERR           57.29577951

void srf2stoch(int param_string_len, char** param_string, struct standrupformat* srf, struct slipfile* sfile, int debug)
{
FILE *fopfile(), *fpr, *fpw;
float dx, dy, len, wid, xp, yp;
float dtop, strike, dip;
float rr, ss, tt, tl, ravg, savg, *slip, *trise, *tinit, *sp, *tr, *ti;
float s1, s2, xx, yy, rake, cosA, sinA, sum;
int i, is, id, nstk, ndip, nseg;
long long ip, m;
int ix, iy, j, k, kp, noff;
int nx, ny, nxdiv, nydiv, nxsum, nysum;
float shypo, dhypo, elon, elat, fac;
char string[1024];

float target_dx = -1.0;
float target_dy = -1.0;

float avgstk = -1.0e+15;
int revstk;

debug = 1;

setpar(param_string_len, param_string);

getpar("target_dx","f",&target_dx);
if(target_dx < 0.0)
   mstpar("dx","f",&dx);
getpar("target_dy","f",&target_dy);
if(target_dy < 0.0)
   mstpar("dy","f",&dy);

getpar("avgstk","f",&avgstk);
endpar();

if(avgstk > -1.0e+14)
   {
   while(avgstk >= 360.0)
      avgstk = avgstk - 360.0;
   while(avgstk < 0.0)
      avgstk = avgstk + 360.0;
   }


slip = NULL;
trise = NULL;
tinit = NULL;
sp = NULL;
tr = NULL;
ti = NULL;
sfile->sp = check_malloc(NP*NQ*LV*sizeof(float));
sfile->tr = check_malloc(NP*NQ*LV*sizeof(float));
sfile->ti = check_malloc(NP*NQ*LV*sizeof(float));

//Initialize arrays to 0
memset(sfile->sp, 0, NP*NQ*LV*sizeof(float));
memset(sfile->tr, 0, NP*NQ*LV*sizeof(float));
memset(sfile->ti, 0, NP*NQ*LV*sizeof(float));

sfile->nseg = srf->srf_prect.nseg;

noff = 0;
for(i=0;i<srf->srf_prect.nseg;i++)
   {
   elon = srf->srf_prect.prectseg[i].elon;
   elat = srf->srf_prect.prectseg[i].elat;
   nstk = srf->srf_prect.prectseg[i].nstk;
   ndip = srf->srf_prect.prectseg[i].ndip;
   len = srf->srf_prect.prectseg[i].flen;
   wid = srf->srf_prect.prectseg[i].fwid;

	printf("len=%f\n", len);

   strike = srf->srf_prect.prectseg[i].stk;
   dip = srf->srf_prect.prectseg[i].dip;
   dtop = srf->srf_prect.prectseg[i].dtop;
   shypo = srf->srf_prect.prectseg[i].shyp;
   dhypo = srf->srf_prect.prectseg[i].dhyp;

   while(strike >= 360.0)
      strike = strike - 360.0;
   while(strike < 0.0)
      strike = strike + 360.0;

   if(target_dx > 0.0)
      {
      nx = (int)(len/target_dx + 0.5);
      dx = len/((float)(nx));
      if(debug) {
	fprintf(stderr,"target_dx= %f actual dx= %f\n",target_dx,dx);
	}
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
      if(debug) {
        fprintf(stderr,"target_dy= %f actual dy= %f\n",target_dy,dy);
	}
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

   if(debug) {
     fprintf(stderr,"seg= %d\n",i);
     fprintf(stderr,"nstk= %d nx= %d nxdiv= %d nxsum= %d\n",nstk,nx,nxdiv,nxsum);
     fprintf(stderr,"ndip= %d ny= %d nydiv= %d nysum= %d\n",ndip,ny,nydiv,nysum);
   }

   printf("size of realloc: %ld\n", (size_t)nxdiv*nydiv*nstk*ndip*sizeof(float));
   fflush(stdout);

   slip = (float *) check_realloc (slip,(size_t)nxdiv*nydiv*nstk*ndip*sizeof(float));
   trise = (float *) check_realloc (trise,(size_t)nxdiv*nydiv*nstk*ndip*sizeof(float));
   tinit = (float *) check_realloc (tinit,(size_t)nxdiv*nydiv*nstk*ndip*sizeof(float));
   sp = (float *) check_realloc (sp,nx*ny*sizeof(float));
   tr = (float *) check_realloc (tr,nx*ny*sizeof(float));
   ti = (float *) check_realloc (ti,nx*ny*sizeof(float));

   savg = 0.0;
   ravg = 0.0;
   for(id=0;id<ndip;id++)
      {
      for(is=0;is<nstk;is++)
         {
	 kp = noff + is + id*nstk;

/*
fprintf(stderr,"id= %d is= %d kp= %d noff= %d\n",id,is,kp,noff);
*/

	 tt = srf->srf_apnts.apntvals[kp].tinit;
	 s1 = srf->srf_apnts.apntvals[kp].slip1;
	 s2 = srf->srf_apnts.apntvals[kp].slip2;

	 rake = srf->srf_apnts.apntvals[kp].rake;

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

	 tl = (srf->srf_apnts.apntvals[kp].dt)*(srf->srf_apnts.apntvals[kp].nt1);
	 if(srf->srf_apnts.apntvals[kp].nt2 > srf->srf_apnts.apntvals[kp].nt1)
	    tl = (srf->srf_apnts.apntvals[kp].dt)*(srf->srf_apnts.apntvals[kp].nt2);
	 if(tl < 0.0)
	    tl = 0.0;
	 for(k=0;k<nydiv;k++)
	    {
	    for(j=0;j<nxdiv;j++)
	       {
	       ip = (long long)is*nxdiv + j + (id*nydiv + k)*(long long)nstk*nxdiv;
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
	 ip = ix + (long long)iy*nx;

         sp[ip] = 0.0;
         tr[ip] = 0.0;
         ti[ip] = 0.0;
	 sum = 0.0;
	 for(k=0;k<nysum;k++)
	    {
	    for(j=0;j<nxsum;j++)
	       {
	       m = (long long)ix*nxsum + j + (iy*nysum + k)*(long long)nx*nxsum;

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

   sfile->elon[i] = elon;
   sfile->elat[i] = elat;
   sfile->nx[i] = nx;
   sfile->ny[i] = ny;
   sfile->dx[i] = dx;
   sfile->dy[i] = dy;
   sfile->strike[i] = strike;
   sfile->dip[i] = dip;
   sfile->ravg[i] = ravg;
   sfile->dtop[i] = dtop;
   sfile->shypo[i] = shypo;
   sfile->dhypo[i] = dhypo;

   revstk = 0;
   if(avgstk > -1.0e+14 && (strike > avgstk + 90.0 || strike < avgstk - 90.0))
      revstk = 1;

   for(iy=0;iy<ny;iy++)
      {
      if(revstk)
         {
         for(ix=nx-1;ix>=0;ix--)
		sfile->sp[(nx-1-ix) + iy*nx + i*NQ*NP] = sp[ix + iy*nx];
	 }
      else
         {
        for(ix=0;ix<nx;ix++) {
			sfile->sp[ix + iy*nx + i*NQ*NP] = sp[ix + iy*nx];
	//		printf("Storing sp index %d in sfile index %d\n", ix+iy*nx, ix+iy*nx+i*nx*ny);
		}
	 }

      }

   for(iy=0;iy<ny;iy++)
      {
      if(revstk)
         {
         for(ix=nx-1;ix>=0;ix--)
		sfile->tr[(nx-1-ix) + iy*nx + i*NQ*NP] = tr[ix + iy*nx];
	 }
      else
         {
         for(ix=0;ix<nx;ix++)
		sfile->tr[ix + iy*nx + i*NQ*NP] = tr[ix + iy*nx];
	 }

      }

   for(iy=0;iy<ny;iy++)
      {
      if(revstk)
         {
         for(ix=nx-1;ix>=0;ix--)
		sfile->ti[ix + iy*nx+i*NQ*NP] = ti[ix + iy*nx];	
	 }
      else
         {
         for(ix=0;ix<nx;ix++)
                sfile->ti[ix + iy*nx+i*NQ*NP] = ti[ix + iy*nx];
	 }

      }

   noff = noff + nstk*ndip;
   }

free(slip);
free(trise);
free(tinit);
free(sp);
free(tr);
free(ti);

//fclose(fpw);
}
