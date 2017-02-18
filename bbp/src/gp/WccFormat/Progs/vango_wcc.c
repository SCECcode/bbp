/********************************************************************/
/*                                                                  */
/*    vango_wcc: Reads fortran binary output from 'normod' or      */
/*                'norma_rob' and sums the responses using the      */
/*                the slip model in the file LISAOUT.DAT.  This     */
/*                code assumes that each line of the input binary   */
/*                files have been padded in the front and the       */
/*                end with four bytes [sizeof(int)] of blank        */
/*                space.  This is only true if the input files      */
/*                are binary AND are written from a Fortran         */
/*                code.                                             */
/*                                                                  */
/*                                  -RWG 06/23/93                   */
/*                                                                  */
/********************************************************************/

#include <sys/file.h>
#include <stdio.h>
#include <math.h>

#define		MAXF	2
#define		MAXPTS	1000000
#define		MAXSUB	100000
#define		MAXSTAT	50
#define		MAXTW	200

struct faultparam
   {
   int nx;
   int ny;
   int ndx;
   int ndy;
   float len;
   float wid;
   float dip;
   float dtot;
   };

struct modelparam
   {
   int nlay;
   float *th;
   float *vs;
   float *rho;
   float *dep;
   float *mu;
   };

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

void *check_malloc(int);
double frand(void);
double sfrand(long *);
void write_wccseis(char *, struct statdata *, float *, int);

int size_float = sizeof(float);
int size_int = sizeof(int);

main (ac, av)
int ac;
char **av;
{
struct statdata head1;
FILE *fp, *fopfile();
float space[50000];
struct faultparam fault;
struct modelparam model;
int fd[MAXF];
float gf[MAXPTS], *seis[3], *x, *aslip, tdel[MAXTW], *sptr, dt, fac;
float rng, gftst, tshft, avgslip, maxslip;
int itst, nseq;
int i, j, ij, ic, nc, nsum, nsumst, nsumend, nmech, nwin, nsub, ntmax;
float mompct, momtrgt, fmnt;
int imech, istat, smode, isub, nt, it, itshft, k, nn, nstat, iwin, nt6;
int indx[MAXSUB], nsum1[MAXSUB];
char statlist[512], mechname[256], str[512];
char *statname[MAXSTAT], statbuf[MAXSTAT*STATCHAR];
char *compname[3*MAXSTAT], compbuf[3*MAXSTAT*4];
char rtimesin[256], rtimesout[256];
char lisafile[512];
float sgn = 1.0;

float tstmin;
float tst = 0.0;
float tstart = -1.0e+15;

int do_tsfac = 0;
float tsfac = 0.0;
float avgrupv = 3.0;
float *tsf, xs[2];
float xmax, xavg, rtmin;
float *xp, *yp, *rt;

float trand = 0.0;
long seed = 1;

int normap = 0;
int bailey = 0;
int flip_vert = 1;

float targetslip = -1.0;
float epi = 0.0;
char title[64];

int outbin = 0;

title[0] = '\0';
rtimesin[0] = '\0';

setpar(ac,av);
getpar("normap","d",&normap);
getpar("bailey","d",&bailey);
if(bailey == 1)
   flip_vert = 0;
getpar("flip_vert","d",&flip_vert);
getpar("outbin","d",&outbin);
getpar("targetslip","f",&targetslip);
getpar("trand","f",&trand);
getpar("seed","d",&seed);
getpar("epi","f",&epi);
getpar("tstart","f",&tstart);
getpar("title","s",title);
getpar("do_tsfac","d",&do_tsfac);
getpar("tsfac","f",&tsfac);
getpar("avgrupv","f",&avgrupv);
getpar("rtimesin","s",rtimesin);
if(rtimesin[0] != '\0')
   mstpar("rtimesout","s",rtimesout);
endpar();

if(tstart > -1.0e+14)
   tst = tstart;

for(ic=0;ic<3;ic++)
   seis[ic] = (float *) check_malloc (size_float*MAXPTS);

printf("Enter the number of different mechanisms\n");
scanf("%d",&nmech);

for(i=0;i<nmech;i++)
   {
   printf("Enter the file containing subfault synthetics for mechanism #%1d\n",i+1);
   scanf("%s",mechname);

   fd[i] = opfile_ro(mechname);
   }

printf("Enter the number of subfault elements along strike\n");
scanf("%d",&fault.nx);

printf("Enter the number of intermediate point sources per subfault element along strike\n");
scanf("%d",&fault.ndx);

printf("Enter the number of subfault elements down dip\n");
scanf("%d",&fault.ny);

printf("Enter the number of intermediate point sources per subfault element down dip\n");
scanf("%d",&fault.ndy);

nsub = fault.nx*fault.ny;
if(bailey == 0)
   {
   fault.ndx++;
   fault.ndy++;
   }

printf("Enter the fault length, width, and dip\n");
scanf("%f %f %f",&fault.len,&fault.wid,&fault.dip);
fault.dip = fault.dip*3.14159/180.0;

printf("Enter the depth to the top of the fault\n");
scanf("%f",&fault.dtot);

printf("Enter the number of layers in the velocity model\n");
scanf("%d",&model.nlay);

model.th = (float *) check_malloc (model.nlay*size_float);
model.dep = (float *) check_malloc (model.nlay*size_float);
model.vs = (float *) check_malloc (model.nlay*size_float);
model.rho = (float *) check_malloc (model.nlay*size_float);
model.mu = (float *) check_malloc (model.nlay*size_float);
for(i=0;i<model.nlay;i++)
   {
   printf("Enter the thickness, shear velocity, and density for layer %d\n",i);
   scanf("%f %f %f",&model.th[i],&model.vs[i],&model.rho[i]);
   model.mu[i] = model.vs[i]*model.vs[i]*model.rho[i];
   if(i)
      model.dep[i] = model.dep[i-1] + model.th[i];
   else
      model.dep[i] = model.th[i];
   }

/*
for(i=0;i<model.nlay;i++)
   {
   printf("dep= %10.3f  vs= %10.3f rho= %10.3f\n",model.dep[i],model.vs[i],model.rho[i]);
   }
   exit(-99);
*/

printf("Enter the mode of subfault summation (0=all,1=select,2=sequence,3=specify rows)\n");
scanf("%d",&smode);
if(smode)
   {
   if(smode == 1)
      {
      printf("Enter the number of subfaults to sum\n");
      scanf("%d",&nsum);
      printf("Enter the subfault indices\n");
      for(i=0;i<nsum;i++)
         scanf("%d",&indx[i]);
      }
   else if(smode == 2)
      {
      nsum = 0;
      printf("Enter the number of sequences to read\n");
      scanf("%d",&nseq);
      while(nseq--)
	 {
         printf("Enter the first and last subfaults in sequence\n");
         scanf("%d %d",&nsumst,&nsumend);
	 k = nsum;
         for(i=nsumst;i<nsumend;i++)
	    {
            indx[k+i-nsumst] = i;
	    nsum++;
	    }
	 }
      }
   else if(smode == 3)
      {
      printf("Enter the number of rows to read\n");
      scanf("%d",&nsum);
      printf("Enter the row indices\n");
      for(i=0;i<nsum;i++)
	 {
         scanf("%d",&nseq);

	 for(j=0;j<fault.nx;j++)
	    indx[i*fault.nx + j] = j*fault.ny + nseq;
	 }
      nsum = nsum*fault.nx;
      }

   for(i=0;i<nsub;i++)
      nsum1[i] = 0;

   for(i=0;i<nsum;i++)
      nsum1[indx[i]-1] = 1;
   }
else
   {
   for(i=0;i<nsub;i++)
      nsum1[i] = 1;
   }

printf("Enter the number of time windows\n");
scanf("%d",&nwin);
printf("*** Number of time windows= %d\n",nwin);
if(nwin > MAXTW)
   {
   fprintf(stderr,"*** Maximum number of time windows exceeded, exiting...\n");
   exit(-1);
   }

for(i=0;i<nwin;i++)
   {
   printf("Enter the time delay for window #%1d\n",i+1);
   scanf("%f",&tdel[i]);
   }

nn = nsub*nwin*nmech;
x = (float *) check_malloc (size_float*nn);
aslip = (float *) check_malloc (size_float*nsub);

printf("Enter the number of stations\n");
scanf("%d",&nstat);

if(nstat > 0)
   {
   for(i=0;i<nstat;i++)
      {
      printf("Enter the name for station #%2d\n",i+1);
      statname[i] = statbuf + i*STATCHAR;
      scanf("%s",statname[i]);
      for(ic=0;ic<3;ic++)
         {
         printf("Enter the 3 letter code for %s comp #%1d\n",statname[i],ic+1);
         compname[3*i + ic] = compbuf + (3*i + ic)*4;
         scanf("%s",compname[3*i + ic]);
         }
      }
   }
else
   {
   scanf("%s",statlist);
   fp = fopfile(statlist,"r");

   fscanf(fp,"%d",&nstat);
   for(i=0;i<nstat;i++)
      {
      statname[i] = statbuf + i*STATCHAR;
      fscanf(fp,"%s",statname[i]);
      for(ic=0;ic<3;ic++)
         {
         compname[3*i + ic] = compbuf + (3*i + ic)*4;
         fscanf(fp,"%s",compname[3*i + ic]);
         }
      }
   fclose(fp);
   }

printf("Enter the target moment and the moment of the input Green`s functions\n");
scanf("%f %f",&momtrgt,&fmnt);

printf("Enter name of slipmodel file (LISAOUT format)-\n");
scanf("%s",lisafile);

fp = fopfile(lisafile,"r");
fgets(str,256,fp);
fgets(str,256,fp);
fgets(str,256,fp);
fgets(str,256,fp);
fgets(str,256,fp);

for(i=0;i<nsub;i++)
   {
   for(j=0;j<nmech*nwin;j++)
      {
      ij = j + i*nmech*nwin;
      fscanf(fp,"%f",&x[ij]);
      }
   }
fclose(fp);

if(momtrgt > 0.0)
   normslip(&fault,&model,x,aslip,&avgslip,&maxslip,nmech,nwin,&momtrgt,&fmnt,nsum1,&mompct,normap,&targetslip);
else
   getmom(&fault,&model,x,aslip,&avgslip,&maxslip,nmech,nwin,&momtrgt,nsum1,&mompct,&targetslip);

printf("***** all slip values in meters\n");
printf("avg. slip= %6.2f\n",avgslip);
printf("max. slip= %6.2f\n",maxslip);
printf("*****\n");

printf("\n");
printf("*** moment used in summation= %8.2e\n",mompct*momtrgt);

for(j=0;j<fault.ny;j++)
   {
   for(i=0;i<fault.nx;i++)
      printf("%6.2f",aslip[i + j*fault.nx]);

   printf("\n");
   }

/* calculate time shift adjustments based on slip amount */

if(do_tsfac != 0 && tsfac == (float)(0.0))
   tsfac = -0.5*((fault.len/fault.nx) + (fault.wid/fault.ny))/avgrupv;

fprintf(stderr,"tsfac= %f\n",tsfac);

xmax = 0.0;
xavg = 0.0;
for(isub=0;isub<nsub;isub++)
   {
   xs[0] = 0.0;
   xs[1] = 0.0;
   for(imech=0;imech<nmech;imech++)
      {
      for(iwin=0;iwin<nwin;iwin++)
         {
         j = isub*nwin*nmech + nwin*imech + iwin;

	 xs[imech] = xs[imech] + x[j];
	 }
      }

   xs[0] = sqrt(xs[0]*xs[0] + xs[1]*xs[1]);

   xavg = xavg + xs[0];

   if(xs[0] > xmax)
      xmax = xs[0];
   }

xavg = xavg/(float)(nsub);

if((xmax-xavg) != (float)(0.0))
   tsfac = tsfac/(xmax-xavg);

tsf = (float *) check_malloc (nsub*sizeof(float));

for(isub=0;isub<nsub;isub++)
   {
   xs[0] = 0.0;
   xs[1] = 0.0;
   for(imech=0;imech<nmech;imech++)
      {
      for(iwin=0;iwin<nwin;iwin++)
         {
         j = isub*nwin*nmech + nwin*imech + iwin;

	 xs[imech] = xs[imech] + x[j];
	 }
      }

   xs[0] = sqrt(xs[0]*xs[0] + xs[1]*xs[1]);
   tsf[isub] = (xs[0]-xavg)*tsfac;

/*
   fprintf(stderr,"%13.5e %9.3f %9.3f %9.3f\n",tsf[isub],xs[0],xmax,xavg);
*/
   }

if(rtimesin[0] != '\0')
   {
   xp = (float *) check_malloc (nsub*fault.ndx*fault.ndy*sizeof(float));
   yp = (float *) check_malloc (nsub*fault.ndx*fault.ndy*sizeof(float));
   rt = (float *) check_malloc (nsub*fault.ndx*fault.ndy*sizeof(float));

   fp = fopfile(rtimesin,"r");

   rtmin = 1.0e+10;
   for(i=0;i<nsub;i++)
      {
      for(j=0;j<fault.ndx*fault.ndy;j++)
         {
	 k = j + i*fault.ndx*fault.ndy;
         fscanf(fp,"%f %f %f",&xp[k],&yp[k],&rt[k]);
	 rt[k] = rt[k] + tsf[i];

	 if(rt[k] < rtmin)
	    rtmin = rt[k];
         }
      }
   fclose(fp);

   fp = fopfile(rtimesout,"w");
   for(i=0;i<nsub;i++)
      {
      tsf[i] = tsf[i] - rtmin;
      for(j=0;j<fault.ndx*fault.ndy;j++)
         {
	 k = j + i*fault.ndx*fault.ndy;
         fprintf(fp,"%13.5e %13.5e %13.5e\n",xp[k],yp[k],rt[k]-rtmin);
         }
      }
   fclose(fp);
   }

ntmax = 0;
for(istat=0;istat<nstat;istat++)
   {
   tstmin = 1.0e+15;
   for(ic=0;ic<3;ic++)
      {
      sptr = seis[ic];
      for(i=0;i<MAXPTS;i++)
         sptr[i] = 0.0;
      }

   for(imech=0;imech<nmech;imech++)
      {
      for(isub=0;isub<nsub;isub++)
         {
	 lseek(fd[imech],size_int,1); /* skip padding in Fortran binary */
	 reed(fd[imech],&nc,size_int);

	 lseek(fd[imech],size_int,1); /* skip padding in Fortran binary */
	 for(ic=0;ic<nc;ic++)
	    {
	    sptr = seis[ic];
	    lseek(fd[imech],size_int,1); /* skip padding in Fortran binary */
	    reed(fd[imech],&rng,size_float);
	    reed(fd[imech],&gftst,size_float);

	    if(gftst < tstmin)
	       tstmin = gftst;

	    lseek(fd[imech],size_int,1); /* skip padding in Fortran binary */

	    lseek(fd[imech],size_int,1); /* skip padding in Fortran binary */
	    reed(fd[imech],&nt,size_int);
	    reed(fd[imech],&dt,size_float);
	    lseek(fd[imech],size_int,1); /* skip padding in Fortran binary */

	    lseek(fd[imech],size_int,1); /* skip padding in Fortran binary */
	    reed(fd[imech],gf,nt*size_float);
	    lseek(fd[imech],size_int,1); /* skip padding in Fortran binary */

	    if(nsum1[isub])
	       {
               for(iwin=0;iwin<nwin;iwin++)
                  {
                  j = isub*nwin*nmech + nwin*imech + iwin;
 
                  tshft = gftst + tdel[iwin] + trand*sfrand(&seed)*dt - tst;
		  tshft = tshft + tsf[isub];

                  if(tshft < 0.0)
		     {
		     itshft = (int)(tshft/dt - 0.5);
                     itst = -itshft;
		     }
                  else
		     {
		     itshft = (int)(tshft/dt + 0.5);
                     itst = 0;
		     }

                  for(it=itst;it<nt;it++)
                     {
                     k = it + itshft;
                     sptr[k] = sptr[k] + x[j]*gf[it];
		     }

	          if(k > ntmax)
	             ntmax = k;
                  }
	       }
	    }
         }
      }

   if(tstart < -1.0e+14) /* use tst from GFs */
      {
      tst = tstmin;
      itshft = (int)(tst/dt + 0.5);

      ntmax = ntmax - itshft;

      for(ic=0;ic<nc;ic++)
         {
         sptr = seis[ic];
         if(itshft > 0)
	    {
            for(it=0;it<ntmax;it++)
               sptr[it] = sptr[it + itshft];
	    }

         if(itshft < 0)
	    {
            for(it=ntmax-1;it>=-itshft;it--)
               sptr[it] = sptr[it + itshft];

            for(it=0;it<-itshft;it++)
               sptr[it] = 0.0;
	    }
         }
      }

   strcpy(head1.stat,statname[istat]);
   strcpy(head1.stitle,title);

   head1.nt = ntmax;
   head1.dt = dt;

   head1.hr = (int)(tst/3600.0);
   head1.min = (int)((tst - 3600.0*head1.hr)/60.0);
   head1.sec = (float)(tst - 60.0*head1.min - 3600.0*head1.hr);

   head1.edist = epi;
   head1.az = 0.0;
   head1.baz = 0.0;

   printf("%s %f\n",head1.stat,tst);

   for(ic=0;ic<nc;ic++)
      {
      sptr = seis[ic];
      strcpy(head1.comp,compname[istat*3 + ic]);

      if(flip_vert == 1 && ic == 2) /*makes vert positive up for norma output*/
	 {
         for(i=0;i<ntmax;i++)
            sptr[i] = -1.0*sptr[i];
         }

      sprintf(str,"%s.%s",head1.stat,head1.comp);
      write_wccseis(str,&head1,sptr,outbin);
      }
   }
}

normslip(fault,model,x,as,avgs,mxs,nm,nw,mtrg,fmnt,nsum,momp,normap,targs)
struct faultparam *fault;
struct modelparam *model;
float *x, *mtrg, *fmnt, *as, *avgs, *mxs, *momp, *targs;
int nm, nw, *nsum, normap;
{
float h, c1, tslip, slip[2], sinD, dh, ddh, sum, r, sum1;
float h0, h1, mubar, th0, th1;
int l0, l1;
int ll, im, iw, i, j, k, l, npnt;
float acnt, acnt0, darea, scale;
float momt;

float conv = 0.01;  /* cm -> meters */

*avgs = 0.0;
*mxs = 0.0;

sum = 0.0;
sum1 = 0.0;
slip[0] = 0.0;
slip[1] = 0.0;

sinD = sin(fault->dip);
dh = sinD*fault->wid/(fault->ny);
ddh = dh/(fault->ndy);
npnt = (fault->ndx)*(fault->ndy);

darea = fault->len*fault->wid/(fault->nx*fault->ny);

acnt0 = 0.0;
acnt = 0.0;
for(j=0;j<fault->nx;j++)
   {
   for(i=0;i<fault->ny;i++)
      {
      /* use averge of rigidity over subfault */

      mubar = 0.0;
      for(k=0;k<(fault->ndy);k++)
	 {
	 h0 = fault->dtot + i*dh + (k+0.5)*ddh;

         l0 = 0;
         while(h0 > model->dep[l0])
	    {
	    l0++;
	    if(l0 == model->nlay)
	       exit(-1);
	    }
	 mubar = mubar + model->mu[l0];
	 }
      mubar = mubar/(fault->ndy);

      for(im=0;im<nm;im++)
	 {
	 slip[im] = 0.0;
	 for(iw=0;iw<nw;iw++)
	    {
	    k = iw + im*nw + (i + j*fault->ny)*nm*nw;
	    slip[im] = slip[im] + x[k];
	    }
	 }

      tslip = sqrt(slip[0]*slip[0] + slip[1]*slip[1]);

      if(normap)
         sum1 = sum1 + mubar*npnt;
      else
	 sum1 = 1.0;

      if(*targs > 0.0)
         sum = sum + darea*(*targs)*mubar*tslip;
      else
         sum = sum + mubar*tslip;

      acnt = acnt + 1.0;
      if(tslip != 0.0)
	 acnt0 = acnt0 + 1.0;
      }
   }

if(*targs > 0.0)
   {
   r = (*mtrg)/(sum*1.0e+20);
   scale = conv*r;
   }
else
   {
   r = *mtrg/(npnt*(*fmnt));
   r = r*sum1/sum;
   scale = conv*(*mtrg)/(darea*sum*1.0e+20);
   }

fprintf(stderr,"****scale=%13.5e\n",scale);

for(j=0;j<fault->nx;j++)
   {
   for(i=0;i<fault->ny;i++)
      {
      /* use averge of rigidity over subfault */

      mubar = 0.0;
      for(k=0;k<(fault->ndy);k++)
	 {
	 h0 = fault->dtot + i*dh + (k+0.5)*ddh;

         l0 = 0;
         while(h0 > model->dep[l0])
	    {
	    l0++;
	    if(l0 == model->nlay)
	       exit(-1);
	    }
	 mubar = mubar + model->mu[l0];
	 }
      mubar = mubar/(fault->ndy);

if(j==0)
fprintf(stderr,"h0=%f mu=%f\n",h0,mubar);
/*
*/

      for(im=0;im<nm;im++)
	 {
	 slip[im] = 0.0;
	 for(iw=0;iw<nw;iw++)
	    {
	    k = iw + im*nw + (i + j*fault->ny)*nm*nw;

	    slip[im] = slip[im] + x[k];

	    if(normap || *targs > 0.0)
	       x[k] = x[k]*r;
	    else
	       x[k] = x[k]*r*mubar;
	    }
	 }

      /* calculate absolute slip in meters, store as rows */

      tslip = sqrt(slip[0]*slip[0] + slip[1]*slip[1]);
      as[j + i*fault->nx] = scale*tslip;

      *avgs = *avgs + as[j + i*fault->nx];
      if(as[j + i*fault->nx] > *mxs)
	 *mxs = as[j + i*fault->nx];
      }
   }

/* if subfault area has zero slip, don't add contribution of sources */
*avgs = (*avgs/(fault->nx*fault->ny))*(acnt/acnt0);

momt = *momp = 0.0;
for(j=0;j<fault->nx;j++)
   {
   for(i=0;i<fault->ny;i++)
      {
      h = fault->dtot + (i+0.5)*dh;
      l = 0;
      while(h > model->dep[l])
	 {
	 l++;
	 if(l == model->nlay)
	    exit(-1);
	 }

      ll = i + j*fault->ny;
      for(im=0;im<nm;im++)
	 {
	 slip[im] = 0.0;
	 for(iw=0;iw<nw;iw++)
	    {
	    k = iw + im*nw + ll*nm*nw;

	    slip[im] = slip[im] + x[k];
	    }
	 }

      tslip = sqrt(slip[0]*slip[0] + slip[1]*slip[1]);
      momt = momt + tslip;
      if(nsum[ll])
	 *momp = *momp + tslip;
      }
   }
*momp = (*momp)/momt;
}

getmom(fault,model,x,as,avgs,mxs,nm,nw,mtrg,nsum,momp,targs)
struct faultparam *fault;
struct modelparam *model;
float *x, *mtrg, *as, *avgs, *mxs, *momp, *targs;
int nm, nw, *nsum;
{
float h, c1, tslip, slip[2], sinD, dh, ddh, sum, r, sum1;
float h0, h1, mubar, th0, th1;
int l0, l1;
int ll, im, iw, i, j, k, l, npnt;
float acnt, acnt0, darea, scale;
float momt;

float conv = 0.01;  /* cm -> meters */

*avgs = 0.0;
*mxs = 0.0;

sum = 0.0;
sum1 = 0.0;
slip[0] = 0.0;
slip[1] = 0.0;

sinD = sin(fault->dip);
dh = sinD*fault->wid/(fault->ny);
ddh = dh/(fault->ndy);
npnt = (fault->ndx)*(fault->ndy);

darea = fault->len*fault->wid/(fault->nx*fault->ny);

acnt0 = 0.0;
acnt = 0.0;
for(j=0;j<fault->nx;j++)
   {
   for(i=0;i<fault->ny;i++)
      {
      /* use averge of rigidity over subfault */

      mubar = 0.0;
      for(k=0;k<(fault->ndy);k++)
	 {
	 h0 = fault->dtot + i*dh + (k+0.5)*ddh;

         l0 = 0;
         while(h0 > model->dep[l0])
	    {
	    l0++;
	    if(l0 == model->nlay)
	       exit(-1);
	    }
	 mubar = mubar + model->mu[l0];
	 }
      mubar = mubar/(fault->ndy);

      for(im=0;im<nm;im++)
	 {
	 slip[im] = 0.0;
	 for(iw=0;iw<nw;iw++)
	    {
	    k = iw + im*nw + (i + j*fault->ny)*nm*nw;
	    slip[im] = slip[im] + x[k];
	    }
	 }

      tslip = sqrt(slip[0]*slip[0] + slip[1]*slip[1]);
      as[j + i*fault->nx] = (*targs)*tslip*conv;

      if(as[j + i*fault->nx] > *mxs)
	 *mxs = as[j + i*fault->nx];

      *avgs = *avgs + as[j + i*fault->nx];

      sum = sum + darea*mubar*as[j + i*fault->nx];
      }
   }

*momp = 1.0;
*mtrg = sum*1.0e+22;
*avgs = *avgs/(fault->nx*fault->ny);
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
