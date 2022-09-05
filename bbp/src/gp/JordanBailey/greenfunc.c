#include "include.h"
#include "structure.h"
#include "function.h"

void get_gfpars(struct gfparam *gfp)
{
FILE *fpr, *fopfile();
float rng, z0;
int i, j;

if(gfp->flag3d == 0)
   {
   fpr = fopfile(gfp->gflocs,"r");

   fscanf(fpr,"%d",&gfp->ngfr);
   gfp->gfr = (float *) check_malloc (gfp->ngfr*sizeof(float));
   for(i=0;i<gfp->ngfr;i++)
      fscanf(fpr,"%f",&gfp->gfr[i]);

   fscanf(fpr,"%d",&gfp->ngfd);
   gfp->gfd = (float *) check_malloc (gfp->ngfd*sizeof(float));
   for(i=0;i<gfp->ngfd;i++)
      fscanf(fpr,"%f",&gfp->gfd[i]);

   fclose(fpr);

   fpr = fopfile(gfp->gftimes,"r");

   gfp->gft = (float *) check_malloc (gfp->ngfd*gfp->ngfr*sizeof(float));
   for(j=0;j<gfp->ngfd;j++)
      {
      for(i=0;i<gfp->ngfr;i++)
         {
         fscanf(fpr,"%*f %*f %*f");
         fscanf(fpr,"%f %f %f",&z0,&rng,&gfp->gft[i+j*(gfp->ngfr)]);

         if(gfp->gfr[i] < rng*(1.0-0.01) || gfp->gfr[i] > rng*(1.0+0.01) ||
            gfp->gfd[j] < z0*(1.0-0.01) || gfp->gfd[j] > z0*(1.0+0.01))
            {
            fprintf(stderr,"*** mismatch between gflocs and gftimes, exiting...\n");
            fprintf(stderr,"*** i= %d j= %d\n",i,j);
            exit(-1);
            }
         }
      }

   fclose(fpr);
   }
else
   {
   fpr = fopfile(gfp->gflocs,"r");

   fscanf(fpr,"%d",&gfp->ngfr);
   gfp->gfn = (float *) check_malloc (gfp->ngfr*sizeof(float));
   gfp->gfe = (float *) check_malloc (gfp->ngfr*sizeof(float));
   gfp->gfr = (float *) check_malloc (gfp->ngfr*sizeof(float));
   for(i=0;i<gfp->ngfr;i++)
      fscanf(fpr,"%f %f %f",&gfp->gfn[i],&gfp->gfe[i],&gfp->gfr[i]);

   fscanf(fpr,"%d",&gfp->ngfd);
   gfp->gfd = (float *) check_malloc (gfp->ngfd*sizeof(float));
   for(i=0;i<gfp->ngfd;i++)
      fscanf(fpr,"%f",&gfp->gfd[i]);

   fclose(fpr);
   }
}

void find_4gf(struct gfparam gfp,struct gfheader *gfh,float *rng,float *dep,float *se,float *sn)
{
int i, ir, ir0, ir1, id0, id1;
float sum, rr0, rr1, dd0, dd1;
float dr, dr0, dr1, drng;

if(gfp.flag3d == 0)
   {
   ir1 = 0;
   while(gfp.gfr[ir1] < *rng && ir1 < (gfp.ngfr)-1)
      ir1++;

   ir0 = ir1 - 1;

   if(ir1 == 0)
      {
      ir0 = ir1;
      fprintf(stderr,"*** range= %f < gf range limits\n",*rng);
      }

   if(ir1 == (gfp.ngfr)-1 && gfp.gfr[ir1] < *rng)
      {
      ir0 = ir1;
      fprintf(stderr,"*** range= %f > gf range limits\n",*rng);
      }

   rr0 = gfp.gfr[ir0] - *rng;
   rr1 = gfp.gfr[ir1] - *rng;
   rr0 = rr0*rr0;
   rr1 = rr1*rr1;
   }
else
   {
   ir0 = -1;
   dr0 = 1.0e+10;
   ir1 = -1;
   dr1 = 1.0e+10;
   for(ir=0;ir<gfp.ngfr;ir++)
      {
      drng = (*rng) - gfp.gfr[ir];

      if(drng < 0.0 && -drng <= gfp.rtol)
         {
	 dr = (gfp.gfn[ir]-(*sn))*(gfp.gfn[ir]-(*sn))
	       + (gfp.gfe[ir]-(*se))*(gfp.gfe[ir]-(*se));

	 if(dr < dr0)
	    {
	    dr0 = dr;
	    ir0 = ir;
	    }
	 }
      if(drng >= 0.0 && drng <= gfp.rtol)
         {
	 dr = (gfp.gfn[ir]-(*sn))*(gfp.gfn[ir]-(*sn))
	       + (gfp.gfe[ir]-(*se))*(gfp.gfe[ir]-(*se));

	 if(dr < dr1)
	    {
	    dr1 = dr;
	    ir1 = ir;
	    }
	 }
      }

   if(ir0 == -1)
      ir0 = ir1;
   if(ir1 == -1)
      ir1 = ir0;

   rr0 = (gfp.gfn[ir0]-(*sn))*(gfp.gfn[ir0]-(*sn))
          + (gfp.gfe[ir0]-(*se))*(gfp.gfe[ir0]-(*se));
   rr1 = (gfp.gfn[ir1]-(*sn))*(gfp.gfn[ir1]-(*sn))
          + (gfp.gfe[ir1]-(*se))*(gfp.gfe[ir1]-(*se));
   }

id1 = 0;
while(gfp.gfd[id1] < *dep && id1 < (gfp.ngfd)-1)
   id1++;

id0 = id1 - 1;

if(id1 == 0)
   id0 = id1;

if(id1 == (gfp.ngfd)-1 && gfp.gfd[id1] < *dep)
   id0 = id1;

dd0 = gfp.gfd[id0] - *dep;
dd1 = gfp.gfd[id1] - *dep;
dd0 = dd0*dd0;
dd1 = dd1*dd1;

gfh[0].read_flag = 0;
if(gfh[0].id == id0 && gfh[0].ir == ir0) /* already read-in */
   gfh[0].read_flag = 1;

gfh[1].read_flag = 0;
if(gfh[1].id == id0 && gfh[1].ir == ir1) /* already read-in */
   gfh[1].read_flag = 1;

gfh[2].read_flag = 0;
if(gfh[2].id == id1 && gfh[2].ir == ir0) /* already read-in */
   gfh[2].read_flag = 1;

gfh[3].read_flag = 0;
if(gfh[3].id == id1 && gfh[3].ir == ir1) /* already read-in */
   gfh[3].read_flag = 1;

/* even if read_flag=1 still set wt's since dists. may be slightly different */
   
gfh[0].id = id0; gfh[0].ir = ir0; gfh[0].wt = sqrt(dd0 + rr0);
gfh[1].id = id0; gfh[1].ir = ir1; gfh[1].wt = sqrt(dd0 + rr1);
gfh[2].id = id1; gfh[2].ir = ir0; gfh[2].wt = sqrt(dd1 + rr0);
gfh[3].id = id1; gfh[3].ir = ir1; gfh[3].wt = sqrt(dd1 + rr1);

if(gfh[0].wt == 0.0)
   { gfh[0].wt = 1.0; gfh[1].wt = 0.0; gfh[2].wt = 0.0; gfh[3].wt = 0.0; }
else if(gfh[1].wt == 0.0)
   { gfh[0].wt = 0.0; gfh[1].wt = 1.0; gfh[2].wt = 0.0; gfh[3].wt = 0.0; }
else if(gfh[2].wt == 0.0)
   { gfh[0].wt = 0.0; gfh[1].wt = 0.0; gfh[2].wt = 1.0; gfh[3].wt = 0.0; }
else if(gfh[3].wt == 0.0)
   { gfh[0].wt = 0.0; gfh[1].wt = 0.0; gfh[2].wt = 0.0; gfh[3].wt = 1.0; }
else
   {
   sum = 1.0/gfh[0].wt + 1.0/gfh[1].wt + 1.0/gfh[2].wt + 1.0/gfh[3].wt;
   sum = 1.0/sum;
   for(i=0;i<4;i++)
      gfh[i].wt = sum/gfh[i].wt;
   }
}

void read_4gf(char *gfpath,char *gfname,float *gf,int nts,struct gfheader *gfh,struct gfparam gfpar,float *maxgft,int *maxnt,float *dtout,float *space)
{
float *gfptr, fnt, pbar;
int kg, ig, resamp, ntpad, ntrsmp;
int gnt;

/* RWG 20140314: added following to turn print statements on(=1)/off(=0,default) */
int print_comments = 0;

float tol = 1.0e-02;

float pratio_tol, mratio_tol;
float ratio_tol = 0.00001;

pratio_tol = 1.0 + ratio_tol;
mratio_tol = 1.0 - ratio_tol;

*maxgft = -1.0e+15;
*maxnt = -1;

for(kg=0;kg<4;kg++)
   {
   if(gfh[kg].read_flag == 1)
      {
      if(print_comments == 1)
         fprintf(stderr,"*** already read-in\n");
      }
   else      /* new location, need to read-in GF */
      {
      if(print_comments == 1)
         fprintf(stderr,"*** new GF to read-in\n");

      gfptr = gf + gfpar.nc*kg*nts;
      zapit(gfptr,gfpar.nc*nts);

      read_gf(gfpath,gfname,gfptr,nts,&gfh[kg],gfpar);
      }

   if(gfpar.flag3d == 0 && gfh[kg].gft > *maxgft)
      *maxgft = gfh[kg].gft;

   /* resample if needed */
   if((*dtout)/gfh[kg].dt > pratio_tol || (*dtout)/gfh[kg].dt < mratio_tol)
      {
      fprintf(stderr,"*** RESAMPLED diff= %13.5e ratio= %13.5e\n",(*dtout)-gfh[kg].dt,(*dtout)/gfh[kg].dt);

      ntpad = 2*gfh[kg].nt;
      fnt = ntpad*gfh[kg].dt/(*dtout);
      gnt = (int)(fnt + 0.5);

      while(nt_tol(fnt,gnt) > tol)
         {
         ntpad++;
         fnt = ntpad*gfh[kg].dt/(*dtout);
         gnt = (int)(fnt + 0.5);
         }

/*
      while(modff(fnt,&fnt) != (float)(0.0))
         {
         ntpad++;
         fnt = ntpad*gfh[kg].dt/(*dtout);
         }
*/

      ntrsmp = (int)(fnt);
      if(ntrsmp > nts)
         {
         fprintf(stderr,"*** resampled nt > ntsum, exiting...\n");
         exit(-1);
         }

      if((*dtout) < gfh[kg].dt)
         resamp = 1;
      else
         resamp = -1;

      for(ig=0;ig<gfpar.nc;ig++)
	 {
	 gfptr = gf + (ig + gfpar.nc*kg)*nts;
         resample(gfptr,gfh[kg].nt,&gfh[kg].dt,resamp,ntpad,ntrsmp,dtout,space);
	 }

      gfh[kg].nt = ntrsmp;
      gfh[kg].dt = (*dtout);
      }

   if(gfh[kg].nt > *maxnt)
      *maxnt = gfh[kg].nt;
   }

if(gfpar.flag3d == 1)  /* adjust stime arrivals using average Vs near GF locations */
   {
   pbar = 0.0;
   for(kg=0;kg<4;kg++)
      pbar = pbar + sqrt(gfh[kg].rho/gfh[kg].mu);  /* slowness */
   pbar = 0.25*pbar;

   for(kg=0;kg<4;kg++)
      {
      gfh[kg].gft = pbar*sqrt(gfh[kg].rng*gfh[kg].rng + gfh[kg].dep*gfh[kg].dep);
      if(gfh[kg].gft > *maxgft)
         *maxgft = gfh[kg].gft;
      }
   }
}

void read_gf(char *gfpath,char *gfname,float *gf,int nts,struct gfheader *gfh,struct gfparam gfpar)
{
float *gfp;
int fdr, it, ig, swap_flag;
char ndx[16], edx[16], ddx[16], str[2048];

if(gfpar.flag3d == 0)
   {
   sprintf(str,"%s/%s%.4d%.4d",gfpath,gfname,(1 + gfh->id),(1 + gfh->ir));

   fdr = opfile_ro(str);

   reed(fdr,&gfh->nt,sizeof(int));
   reed(fdr,&gfh->dt,sizeof(float));
   reed(fdr,&gfh->dep,sizeof(float));
   reed(fdr,&gfh->rng,sizeof(float));
   reed(fdr,&gfh->tst,sizeof(float));
   reed(fdr,&gfh->mom,sizeof(float));
   reed(fdr,&gfh->mu,sizeof(float));

   if(gfpar.swap_flag == 0)  /* byte-order is OK, no swapping needed */
      swap_flag = 0;
   else                        /* byte-order not OK, swapping needed */
      {
      swap_in_place(1,(char *)(&gfh->nt));
      swap_in_place(1,(char *)(&gfh->dt));
      swap_in_place(1,(char *)(&gfh->dep));
      swap_in_place(1,(char *)(&gfh->rng));
      swap_in_place(1,(char *)(&gfh->tst));
      swap_in_place(1,(char *)(&gfh->mom));
      swap_in_place(1,(char *)(&gfh->mu));

      swap_flag = 1;
      }

   gfh->gft = gfpar.gft[gfh->ir + (gfh->id)*(gfpar.ngfr)];

   for(ig=0;ig<gfpar.nc;ig++)
      {
      gfp = gf + ig*nts;
      reed(fdr,gfp,gfh->nt*sizeof(float));

      if(swap_flag == 1)
         swap_in_place(gfh->nt,(char *)(gfp));
      }
   }
else
   {
   if(gfpar.gfn[gfh->ir] < 0.0)
      sprintf(ndx,"-%.6d",(int)(-1000.0*gfpar.gfn[gfh->ir] + 0.5));
   else
      sprintf(ndx,"+%.6d",(int)(1000.0*gfpar.gfn[gfh->ir] + 0.5));

   if(gfpar.gfe[gfh->ir] < 0.0)
      sprintf(edx,"-%.6d",(int)(-1000.0*gfpar.gfe[gfh->ir] + 0.5));
   else
      sprintf(edx,"+%.6d",(int)(1000.0*gfpar.gfe[gfh->ir] + 0.5));

   sprintf(ddx,"%.5d",(int)(1000.0*gfpar.gfd[gfh->id] + 0.5));

   if(gfpar.use_depdir == 1)
      sprintf(str,"%s/D%s/%s_n%se%sd%s.sgt",gfpath,ddx,gfname,ndx,edx,ddx);
   else
      sprintf(str,"%s/%s_n%se%sd%s.sgt",gfpath,gfname,ndx,edx,ddx);

   fprintf(stderr,"GF= %s\n",str);

   fdr = opfile_ro(str);

   reed(fdr,&swap_flag,sizeof(int));
   reed(fdr,&gfh->olon,sizeof(float));
   reed(fdr,&gfh->olat,sizeof(float));
   reed(fdr,&gfh->slon,sizeof(float));
   reed(fdr,&gfh->slat,sizeof(float));
   reed(fdr,&gfh->north,sizeof(float));
   reed(fdr,&gfh->east,sizeof(float));
   reed(fdr,&gfh->dep,sizeof(float));
   reed(fdr,&gfh->xazim,sizeof(float));
   reed(fdr,&gfh->lam,sizeof(float));
   reed(fdr,&gfh->mu,sizeof(float));
   reed(fdr,&gfh->rho,sizeof(float));
   reed(fdr,&gfh->gft,sizeof(float));
   reed(fdr,&gfh->xmom,sizeof(float));
   reed(fdr,&gfh->ymom,sizeof(float));
   reed(fdr,&gfh->zmom,sizeof(float));
   reed(fdr,&gfh->tst,sizeof(float));
   reed(fdr,&gfh->nt,sizeof(int));
   reed(fdr,&gfh->dt,sizeof(float));

   if(swap_flag == SWAP_FLAG)  /* byte-order is OK, no swapping needed */
      swap_flag = 0;
   else                        /* byte-order not OK, swapping needed */
      {
      swap_in_place(1,(char *)(&swap_flag));
      swap_in_place(1,(char *)(&gfh->olon));
      swap_in_place(1,(char *)(&gfh->olat));
      swap_in_place(1,(char *)(&gfh->slon));
      swap_in_place(1,(char *)(&gfh->slat));
      swap_in_place(1,(char *)(&gfh->north));
      swap_in_place(1,(char *)(&gfh->east));
      swap_in_place(1,(char *)(&gfh->dep));
      swap_in_place(1,(char *)(&gfh->xazim));
      swap_in_place(1,(char *)(&gfh->lam));
      swap_in_place(1,(char *)(&gfh->mu));
      swap_in_place(1,(char *)(&gfh->rho));
      swap_in_place(1,(char *)(&gfh->gft));
      swap_in_place(1,(char *)(&gfh->xmom));
      swap_in_place(1,(char *)(&gfh->ymom));
      swap_in_place(1,(char *)(&gfh->zmom));
      swap_in_place(1,(char *)(&gfh->tst));
      swap_in_place(1,(char *)(&gfh->nt));
      swap_in_place(1,(char *)(&gfh->dt));

      if(swap_flag != SWAP_FLAG)  /* byte-order still not OK */
         {
	 fprintf(stderr,"Problem with swap_flag= %d, exiting...\n",swap_flag);
	 exit(-1);
         }
      else
         swap_flag = 1;
      }

   gfh->rng = sqrt((gfh->north)*(gfh->north) + (gfh->east)*(gfh->east));

fprintf(stderr,"xm= %13.5e ymom= %13.5e zmom= %13.5e\n",gfh->xmom,gfh->ymom,gfh->zmom);
   if(gfh->xmom >= 0.0)
      {
      for(ig=0;ig<6;ig++)
         {
         gfp = gf + ig*nts;
         reed(fdr,gfp,gfh->nt*sizeof(float));

         if(swap_flag == 1)
            swap_in_place(gfh->nt,(char *)(gfp));
         }
      }

   if(gfh->ymom >= 0.0)
      {
      for(ig=6;ig<12;ig++)
         {
         gfp = gf + ig*nts;
         reed(fdr,gfp,gfh->nt*sizeof(float));

         if(swap_flag == 1)
            swap_in_place(gfh->nt,(char *)(gfp));
         }
      }

   if(gfh->zmom >= 0.0)
      {
      for(ig=12;ig<18;ig++)
         {
         gfp = gf + ig*nts;
         reed(fdr,gfp,gfh->nt*sizeof(float));

         if(swap_flag == 1)
            swap_in_place(gfh->nt,(char *)(gfp));
         }
      }
   }

close(fdr);
}

void sum_4gf(float *seis,int ntout,float *gf,struct gfheader *gfh,int ntsum,int maxnt,float *rupt,float *maxgft,float *tstart,int tdflag,struct mechparam mp)
{
int ig, it, im;
float backt0, t0[4];
float *sptr, *gfptr;

backt0 = 0.0;
for(ig=0;ig<4;ig++)
   {
   t0[ig] = *maxgft - gfh[ig].gft;
   backt0 = backt0 - t0[ig]*gfh[ig].wt;
   t0[ig] = t0[ig] + gfh[ig].tst;
   }
backt0 = backt0 + *rupt - *tstart;

for(im=0;im<mp.nmech;im++)
   {
   sptr = seis + 3*im*ntout;
   gfptr = gf + 12*im*ntsum;
   timeshift_gf(sptr,ntout,gfptr,gfh,ntsum,t0,&backt0,tdflag);
   }
}

void timeshift_gf(float *seis,int ntout,float *gf,struct gfheader *gfh,int ntsum,float *t0,float *bt0,int tdf)
{
struct complex *gc0, *gc1, *gc2;
float *gf0, *gf1, *gf2, *sv, *sn, *se, *gfv, *gfn, *gfe;
float cosA, sinA, arg, fac, norm, tmpre, scale, tsh;
int i, ig, tapst, it, nts3, nts6, nts9;
int itshift, it0, maxnt, ntp2, nf2;

int taplen = 10;
float zap = 0.0;
float half = 0.5;
float one = 1.0;
float two = 2.0;
float pi = 3.141592654;

maxnt = 0;

sv = seis;
sn = seis + ntout;
se = seis + 2*ntout;

if(tdf)
   {
   for(ig=0;ig<4;ig++)
      {
      gf0 = gf + 3*ig*ntsum;
      gf1 = gf + 3*ig*ntsum + ntsum;
      gf2 = gf + 3*ig*ntsum + 2*ntsum;

   /* taper */

      tapst = gfh[ig].nt - taplen;
      arg = pi/(float)(taplen);
      for(it=tapst+1;it<gfh[ig].nt;it++)
         {
         fac = half*(one + cos(arg*(it-tapst)));
   
         gf0[it] = fac*gf0[it];
         gf1[it] = fac*gf1[it];
         gf2[it] = fac*gf2[it];
         }

/* apply time shift */

      tsh = t0[ig] + *bt0;
      if(tsh >= 0.0)
         itshift = (int)(tsh/gfh[ig].dt + 0.5);
      else
         itshift = (int)(tsh/gfh[ig].dt - 0.5);

      it0 = gfh[ig].nt + itshift;

      if(it0 > ntsum)
         it0 = ntsum;

      if(it0 > ntout)
         it0 = ntout;

      if(itshift < 0)
         {
         for(i=0;i<it0;i++)
            {
	    sv[i] = sv[i] + gf0[i-itshift];
	    sn[i] = sn[i] + gf1[i-itshift];
	    se[i] = se[i] + gf2[i-itshift];
	    }
         }
      else
         {
         for(i=it0-1;i>=itshift;i--)
            {
	    sv[i] = sv[i] + gf0[i-itshift];
	    sn[i] = sn[i] + gf1[i-itshift];
	    se[i] = se[i] + gf2[i-itshift];
	    }
         }

      if(it0 > maxnt)
         maxnt = it0;
      }
   }
else
   {
   ntp2 = 2;
   while(ntp2 < ntsum)
      ntp2 = ntp2*2;

   nf2 = ntp2/2;

   for(ig=0;ig<4;ig++)
      {
      gf0 = gf + 3*ig*ntsum;
      gf1 = gf + 3*ig*ntsum + ntsum;
      gf2 = gf + 3*ig*ntsum + 2*ntsum;

      /* taper and pad */

      tapst = gfh[ig].nt - taplen;
      arg = pi/(float)(taplen);
      for(it=tapst+1;it<gfh[ig].nt;it++)
         {
         fac = half*(one + cos(arg*(it-tapst)));
   
         gf0[it] = fac*gf0[it];
         gf1[it] = fac*gf1[it];
         gf2[it] = fac*gf2[it];
         }
   
      for(it=gfh[ig].nt;it<ntp2;it++)
         {
         gf0[it] = zap;
         gf1[it] = zap;
         gf2[it] = zap;
         }
   
   /* apply time shift */
   
      gc0 = (struct complex *) gf0;
      gc1 = (struct complex *) gf1;
      gc2 = (struct complex *) gf2;
   
      forfft(gc0,ntp2,-1);
      forfft(gc1,ntp2,-1);
      forfft(gc2,ntp2,-1);
   
      /* zero nyquist */
      gc0[0].im = zap;
      gc1[0].im = zap;
      gc2[0].im = zap;
   
      fac = -(t0[ig] + *bt0)*two*pi/(ntp2*gfh[ig].dt);
      for(i=1;i<nf2;i++)
         {
         arg = fac*i;
         cosA = cos(arg);
         sinA = sin(arg);
   
         tmpre     = gc0[i].re*cosA - gc0[i].im*sinA;
         gc0[i].im = gc0[i].re*sinA + gc0[i].im*cosA;
         gc0[i].re = tmpre;
   
         tmpre     = gc1[i].re*cosA - gc1[i].im*sinA;
         gc1[i].im = gc1[i].re*sinA + gc1[i].im*cosA;
         gc1[i].re = tmpre;
   
         tmpre     = gc2[i].re*cosA - gc2[i].im*sinA;
         gc2[i].im = gc2[i].re*sinA + gc2[i].im*cosA;
         gc2[i].re = tmpre;
         }
      }
   
   /* reset pointers to first group of 8 GFs */
   
   gf0 = gf;
   gf1 = gf + ntsum;
   gf2 = gf + 2*ntsum;
   
   gc0 = (struct complex *) gf0;
   gc1 = (struct complex *) gf1;
   gc2 = (struct complex *) gf2;
   
   /* sum over 4 GFs to interpolated response */
   
   nts3 = 3*nf2;
   nts6 = 6*nf2;
   nts9 = 9*nf2;
   for(i=0;i<nf2;i++)
      {
      gc0[i].re = gc0[i].re + gc0[i+nts3].re + gc0[i+nts6].re + gc0[i+nts9].re;
      gc0[i].im = gc0[i].im + gc0[i+nts3].im + gc0[i+nts6].im + gc0[i+nts9].im;
   
      gc1[i].re = gc1[i].re + gc1[i+nts3].re + gc1[i+nts6].re + gc1[i+nts9].re;
      gc1[i].im = gc1[i].im + gc1[i+nts3].im + gc1[i+nts6].im + gc1[i+nts9].im;
   
      gc2[i].re = gc2[i].re + gc2[i+nts3].re + gc2[i+nts6].re + gc2[i+nts9].re;
      gc2[i].im = gc2[i].im + gc2[i+nts3].im + gc2[i+nts6].im + gc2[i+nts9].im;
      }

   invfft(gc0,ntp2,1);
   invfft(gc1,ntp2,1);
   invfft(gc2,ntp2,1);

   maxnt = ntout;
   if(ntp2 < maxnt)
      maxnt = ntp2;

   gfv = gf;
   gfn = gf + ntsum;
   gfe = gf + 2*ntsum;

   scale = one/(float)(ntp2);
   for(it=0;it<maxnt;it++)
      {
      sv[it] = sv[it] + scale*gfv[it];
      sn[it] = sn[it] + scale*gfn[it];
      se[it] = se[it] + scale*gfe[it];
      }
   }
}

void mech_4gf(float *gfmech,float *gf,struct gfheader *gfh,struct gfparam gfp,int nts,struct mechparam mp,float *azi,float *scl)
{
float *zdd, *rdd, *zds, *rds, *tds, *zss, *rss, *tss;
float *axx, *ayy, *azz, *axy, *axz, *ayz;
float *bxx, *byy, *bzz, *bxy, *bxz, *byz;
float *cxx, *cyy, *czz, *cxy, *cxz, *cyz;
float *gfn, *gfe, *gfv, *gfmptr;
float f1, f2, f3, f4, f5;
float cxS, sxS, cxD, sxD, cx2D, sx2D, cxL, sxL;
float cxT, sxT, cx2T, sx2T;
float arg, cosA, sinA, rad, tan, scale;
float xamp, yamp, zamp, sx, sy, sz;
float mxx, myy, mzz, mxy, mxz, myz;
float sum, rake;
float u1, u2, u3, vx, vy, vz, l2m;
float us, ud, ux, uy, uz;
int it, ig, im;

float half = 0.5;
float two = 2.0;
float rperd = 0.017453293;

arg = (mp.dip)*rperd;
cxD = cos(arg);
sxD = sin(arg);

cx2D = cxD*cxD - sxD*sxD;
sx2D = two*sxD*cxD;

sum = 0.0;
for(im=0;im<mp.nmech;im++)
   {
   u1 = u2 = u3 = 0;
   if(mp.flag[im] == U1FLAG)
      u1 = 1;
   else if(mp.flag[im] == U2FLAG)
      u2 = 1;
   else if(mp.flag[im] == U3FLAG)
      u3 = 1;

   gfmptr = gfmech + im*12*nts;
   zapit(gfmptr,12*nts);

   if(gfp.flag3d == 0)
      {
      arg = *azi - (mp.stk)*rperd;
      cxT = cos(arg);
      sxT = sin(arg);

      cx2T = cxT*cxT - sxT*sxT;
      sx2T = two*sxT*cxT;

      if(u1)
         rake = mp.rak;
      else if(u2)
         rake = mp.rak + 90.0;

      arg = rake*rperd;
      cxL = cos(arg);
      sxL = sin(arg);

      cosA = cos(*azi);
      sinA = sin(*azi);

      for(ig=0;ig<4;ig++)
         {
         gfv = gfmptr + 3*ig*nts;
         gfn = gfmptr + 3*ig*nts + nts;
         gfe = gfmptr + 3*ig*nts + 2*nts;

         zdd = gf + gfp.nc*ig*nts;
         rdd = gf + gfp.nc*ig*nts + nts;
         zds = gf + gfp.nc*ig*nts + 2*nts;
         rds = gf + gfp.nc*ig*nts + 3*nts;
         tds = gf + gfp.nc*ig*nts + 4*nts;
         zss = gf + gfp.nc*ig*nts + 5*nts;
         rss = gf + gfp.nc*ig*nts + 6*nts;
         tss = gf + gfp.nc*ig*nts + 7*nts;

/*
fprintf(stderr,"%d: scl= %.5e mu= %.5e wt = %.5e\n",ig,(*scl),gfh[ig].mu,gfh[ig].wt,gfh[ig].mom);
*/

         sum = sum + (*scl)*(gfh[ig].mu)*(gfh[ig].wt);
         scale = (*scl)*(gfh[ig].mu)*(gfh[ig].wt)/(gfh[ig].mom);

         f1 = half*scale*sxL*sx2D;
         f2 = -scale*(sxT*sxL*cx2D - cxT*cxL*cxD);
         f3 = -scale*(sx2T*cxL*sxD + half*cx2T*sxL*sx2D);
         f4 = scale*(cxT*sxL*cx2D + sxT*cxL*cxD);
         f5 = scale*(cx2T*cxL*sxD - half*sx2T*sxL*sx2D);

         for(it=0;it<gfh[ig].nt;it++)
            {
            gfv[it] = f1*zdd[it] + f2*zds[it] + f3*zss[it];

            rad = f1*rdd[it] + f2*rds[it] + f3*rss[it];
            tan = f4*tds[it] + f5*tss[it];

            gfn[it] = -tan*sinA + rad*cosA;
            gfe[it] = tan*cosA + rad*sinA;
            }
         }
      }
   else
      {
      arg = (mp.rak)*rperd;
      cxL = cos(arg);
      sxL = sin(arg);

      sum = 0.0;
      for(ig=0;ig<4;ig++)
         {
         arg = (mp.stk - gfh[ig].xazim)*rperd;
         cxT = cos(arg);
         sxT = sin(arg);
/*
         cx2T = cxT*cxT - sxT*sxT;
         sx2T = two*sxT*cxT;

         mxx = -(sxD*cxL*sx2T + sx2D*sxL*sxT*sxT);
         myy = (sxD*cxL*sx2T - sx2D*sxL*cxT*cxT);
         mzz = sx2D*sxL;
         mxy = (sxD*cxL*cx2T + half*sx2D*sxL*sx2T);
         mxz = -(cxD*cxL*cxT + cx2D*sxL*sxT);
         myz = -(cxD*cxL*sxT - cx2D*sxL*cxT);
*/

         vx = -sxD*sxT;
         vy =  sxD*cxT;
         vz = -cxD;

         us = u1*cxL - u2*sxL;
         ud = u1*sxL + u2*cxL;

         ux = -(u3*sxD - ud*cxD)*sxT + us*cxT;
         uy =  (u3*sxD - ud*cxD)*cxT + us*sxT;
         uz = -(u3*cxD + ud*sxD);

	 l2m = gfh[ig].lam + two*gfh[ig].mu;

         mxx = l2m*vx*ux + (gfh[ig].lam)*vy*uy + (gfh[ig].lam)*vz*uz;
         myy = (gfh[ig].lam)*vx*ux + l2m*vy*uy + (gfh[ig].lam)*vz*uz;
         mzz = (gfh[ig].lam)*vx*ux + (gfh[ig].lam)*vy*uy + l2m*vz*uz;
         mxy = (gfh[ig].mu)*(vx*uy + vy*ux);
         mxz = (gfh[ig].mu)*(vx*uz + vz*ux);
         myz = (gfh[ig].mu)*(vy*uz + vz*uy);

         arg = gfh[ig].xazim*rperd;
         cosA = cos(arg);
         sinA = sin(arg);

         gfv = gfmptr + 3*ig*nts;
         gfn = gfmptr + 3*ig*nts + nts;
         gfe = gfmptr + 3*ig*nts + 2*nts;

         axx = gf + gfp.nc*ig*nts;
         ayy = gf + gfp.nc*ig*nts + nts;
         azz = gf + gfp.nc*ig*nts + 2*nts;
         axy = gf + gfp.nc*ig*nts + 3*nts;
         axz = gf + gfp.nc*ig*nts + 4*nts;
         ayz = gf + gfp.nc*ig*nts + 5*nts;
         bxx = gf + gfp.nc*ig*nts + 6*nts;
         byy = gf + gfp.nc*ig*nts + 7*nts;
         bzz = gf + gfp.nc*ig*nts + 8*nts;
         bxy = gf + gfp.nc*ig*nts + 9*nts;
         bxz = gf + gfp.nc*ig*nts + 10*nts;
         byz = gf + gfp.nc*ig*nts + 11*nts;
         cxx = gf + gfp.nc*ig*nts + 12*nts;
         cyy = gf + gfp.nc*ig*nts + 13*nts;
         czz = gf + gfp.nc*ig*nts + 14*nts;
         cxy = gf + gfp.nc*ig*nts + 15*nts;
         cxz = gf + gfp.nc*ig*nts + 16*nts;
         cyz = gf + gfp.nc*ig*nts + 17*nts;

         sum = sum + (*scl)*(gfh[ig].mu)*(gfh[ig].wt);

         xamp = (*scl)*(gfh[ig].wt)/(gfh[ig].xmom);
         yamp = (*scl)*(gfh[ig].wt)/(gfh[ig].ymom);
         zamp = (*scl)*(gfh[ig].wt)/(gfh[ig].zmom);

         for(it=0;it<gfh[ig].nt;it++)
            {
            sx = xamp*(axx[it]*mxx + ayy[it]*myy + azz[it]*mzz
               + axy[it]*mxy + axz[it]*mxz + ayz[it]*myz);

            sy = yamp*(bxx[it]*mxx + byy[it]*myy + bzz[it]*mzz
               + bxy[it]*mxy + bxz[it]*mxz + byz[it]*myz);

            sz = zamp*(cxx[it]*mxx + cyy[it]*myy + czz[it]*mzz
               + cxy[it]*mxy + cxz[it]*mxz + cyz[it]*myz);

            gfe[it] = sx*sinA + sy*cosA;
            gfn[it] = sx*cosA - sy*sinA;
            gfv[it] = -sz;
            }
         }
      }
   }

*scl = sum;  /* scl now contains the moment released for this point source */
}

void getname_gf(char *str,char *gfname,struct gfheader *gfh,struct gfparam gfpar)
{
char ndx[16], edx[16], ddx[16];

if(gfpar.flag3d == 0)
   {
   sprintf(str,"%s%.4d%.4d",gfname,(1 + gfh->id),(1 + gfh->ir));
   }
else
   {
   if(gfpar.gfn[gfh->ir] < 0.0)
      sprintf(ndx,"-%.6d",(int)(-1000.0*gfpar.gfn[gfh->ir] + 0.5));
   else
      sprintf(ndx,"+%.6d",(int)(1000.0*gfpar.gfn[gfh->ir] + 0.5));

   if(gfpar.gfe[gfh->ir] < 0.0)
      sprintf(edx,"-%.6d",(int)(-1000.0*gfpar.gfe[gfh->ir] + 0.5));
   else
      sprintf(edx,"+%.6d",(int)(1000.0*gfpar.gfe[gfh->ir] + 0.5));

   sprintf(ddx,"%.5d",(int)(1000.0*gfpar.gfd[gfh->id] + 0.5));

   sprintf(str,"D%s/%s_n%se%sd%s.sgt",ddx,gfname,ndx,edx,ddx);
   }
}

int check_name(char *str,char *list,int n,int blen)
{
int i, lf;
char *lptr;

lf = 0;
i = 0;
while(i < n)
   {
   lptr = list + i*blen;
   if(strcmp(str,lptr) == 0)
      {
      lf = 1;
      i = n;
      }
    else
      i++;
   }

return(lf);
}

void find_4gf_use_closest(struct gfparam gfp,struct gfheader *gfh,float *rng,float *dep,float *se,float *sn)
{
int i, ir, ir0, ir1, id0, id1;
float sum, rr0, rr1, dd0, dd1;
float dr, dr0, dr1, drng;
int imin;

if(gfp.flag3d == 0)
   {
   ir1 = 0;
   while(gfp.gfr[ir1] < *rng && ir1 < (gfp.ngfr)-1)
      ir1++;

   ir0 = ir1 - 1;

   if(ir1 == 0)
      {
      ir0 = ir1;
      fprintf(stderr,"*** range= %f < gf range limits\n",*rng);
      }

   if(ir1 == (gfp.ngfr)-1 && gfp.gfr[ir1] < *rng)
      {
      ir0 = ir1;
      fprintf(stderr,"*** range= %f > gf range limits\n",*rng);
      }

   rr0 = gfp.gfr[ir0] - *rng;
   rr1 = gfp.gfr[ir1] - *rng;
   rr0 = rr0*rr0;
   rr1 = rr1*rr1;
   }
else
   {
   ir0 = -1;
   dr0 = 1.0e+10;
   ir1 = -1;
   dr1 = 1.0e+10;
   for(ir=0;ir<gfp.ngfr;ir++)
      {
      drng = (*rng) - gfp.gfr[ir];

      if(drng < 0.0 && -drng <= gfp.rtol)
         {
	 dr = (gfp.gfn[ir]-(*sn))*(gfp.gfn[ir]-(*sn))
	       + (gfp.gfe[ir]-(*se))*(gfp.gfe[ir]-(*se));

	 if(dr < dr0)
	    {
	    dr0 = dr;
	    ir0 = ir;
	    }
	 }
      if(drng >= 0.0 && drng <= gfp.rtol)
         {
	 dr = (gfp.gfn[ir]-(*sn))*(gfp.gfn[ir]-(*sn))
	       + (gfp.gfe[ir]-(*se))*(gfp.gfe[ir]-(*se));

	 if(dr < dr1)
	    {
	    dr1 = dr;
	    ir1 = ir;
	    }
	 }
      }

   if(ir0 == -1)
      ir0 = ir1;
   if(ir1 == -1)
      ir1 = ir0;

   rr0 = (gfp.gfn[ir0]-(*sn))*(gfp.gfn[ir0]-(*sn))
          + (gfp.gfe[ir0]-(*se))*(gfp.gfe[ir0]-(*se));
   rr1 = (gfp.gfn[ir1]-(*sn))*(gfp.gfn[ir1]-(*sn))
          + (gfp.gfe[ir1]-(*se))*(gfp.gfe[ir1]-(*se));
   }

id1 = 0;
while(gfp.gfd[id1] < *dep && id1 < (gfp.ngfd)-1)
   id1++;

id0 = id1 - 1;

if(id1 == 0)
   id0 = id1;

if(id1 == (gfp.ngfd)-1 && gfp.gfd[id1] < *dep)
   id0 = id1;

dd0 = gfp.gfd[id0] - *dep;
dd1 = gfp.gfd[id1] - *dep;
dd0 = dd0*dd0;
dd1 = dd1*dd1;

gfh[0].read_flag = 0;
if(gfh[0].id == id0 && gfh[0].ir == ir0) /* already read-in */
   gfh[0].read_flag = 1;

gfh[1].read_flag = 0;
if(gfh[1].id == id0 && gfh[1].ir == ir1) /* already read-in */
   gfh[1].read_flag = 1;

gfh[2].read_flag = 0;
if(gfh[2].id == id1 && gfh[2].ir == ir0) /* already read-in */
   gfh[2].read_flag = 1;

gfh[3].read_flag = 0;
if(gfh[3].id == id1 && gfh[3].ir == ir1) /* already read-in */
   gfh[3].read_flag = 1;

/* even if read_flag=1 still set wt's since dists. may be slightly different */
   
gfh[0].id = id0; gfh[0].ir = ir0; gfh[0].wt = sqrt(dd0 + rr0);
gfh[1].id = id0; gfh[1].ir = ir1; gfh[1].wt = sqrt(dd0 + rr1);
gfh[2].id = id1; gfh[2].ir = ir0; gfh[2].wt = sqrt(dd1 + rr0);
gfh[3].id = id1; gfh[3].ir = ir1; gfh[3].wt = sqrt(dd1 + rr1);

imin = 0;
for(ir=1;ir<4;ir++)
   {
   if(gfh[ir].wt < gfh[imin].wt)
      imin = ir;
   }

for(ir=0;ir<4;ir++)
   gfh[ir].wt = 0.0;

gfh[imin].wt = 1.0;
}

void reset_wt_4gf(struct gfheader *gfh,float *rng,float *minr,float *maxr)
{
float del, twt[4];
int ir, imax;

if(*rng >= *minr)
   {
   imax = 0;
   for(ir=1;ir<4;ir++)
      {
      if(gfh[ir].wt > gfh[imax].wt)
         imax = ir;
      }

   if(*rng >= *maxr)
      {
      for(ir=0;ir<4;ir++)
         gfh[ir].wt = 0.0;

      gfh[imax].wt = 1.0;
      }
   else
      {
      del = (*rng - *minr)/(*maxr - *minr);

      for(ir=0;ir<4;ir++)
         twt[ir] = 0.0;

      twt[imax] = 1.0;

      for(ir=0;ir<4;ir++)
         gfh[ir].wt = gfh[ir].wt + del*(twt[ir] - gfh[ir].wt);
      }
   }
}

void read_4gf_v2(char *gfpath,char *gfname,float *gf,int nts,struct gfheader *gfh,struct gfparam gfpar,int *maxnt,float *dtout,float *space)
{
float *gfptr, fnt, pbar;
int kg, ig, resamp, ntpad, ntrsmp;
int gnt;

/* RWG 20140314: added following to turn print statements on(=1)/off(=0,default) */
int print_comments = 0;

float tol = 1.0e-02;

float pratio_tol, mratio_tol;
float ratio_tol = 0.00001;

pratio_tol = 1.0 + ratio_tol;
mratio_tol = 1.0 - ratio_tol;

*maxnt = -1;

for(kg=0;kg<4;kg++)
   {
   if(gfh[kg].read_flag == 0 && gfh[kg].wt > (float)(0.0))
      {
      if(print_comments == 1)
         fprintf(stderr,"*** new GF to read-in\n");

      gfptr = gf + gfpar.nc*kg*nts;
      zapit(gfptr,gfpar.nc*nts);

      read_gf(gfpath,gfname,gfptr,nts,&gfh[kg],gfpar);
      gfh[kg].read_flag = 1;

      /* resample if needed */
      if((*dtout)/gfh[kg].dt > pratio_tol || (*dtout)/gfh[kg].dt < mratio_tol)
         {
         fprintf(stderr,"*** RESAMPLED diff= %13.5e ratio= %13.5e\n",(*dtout)-gfh[kg].dt,(*dtout)/gfh[kg].dt);

         ntpad = 2*gfh[kg].nt;
         fnt = ntpad*gfh[kg].dt/(*dtout);
         gnt = (int)(fnt + 0.5);

         while(nt_tol(fnt,gnt) > tol)
            {
            ntpad++;
            fnt = ntpad*gfh[kg].dt/(*dtout);
            gnt = (int)(fnt + 0.5);
            }

         ntrsmp = (int)(fnt);
         if(ntrsmp > nts)
            {
            fprintf(stderr,"*** resampled nt > ntsum, exiting...\n");
            exit(-1);
            }

         if((*dtout) < gfh[kg].dt)
            resamp = 1;
         else
            resamp = -1;

         for(ig=0;ig<gfpar.nc;ig++)
	    {
	    gfptr = gf + (ig + gfpar.nc*kg)*nts;
            resample(gfptr,gfh[kg].nt,&gfh[kg].dt,resamp,ntpad,ntrsmp,dtout,space);
	    }

         gfh[kg].nt = ntrsmp;
         gfh[kg].dt = (*dtout);
         }
      }

   if(gfh[kg].read_flag == 1 && gfh[kg].wt > (float)(0.0))
      {
      if(gfh[kg].nt > *maxnt)
         *maxnt = gfh[kg].nt;
      }
   }
}

void mech_4gf_v2(float *gfmech,float *gf,struct gfheader *gfh,struct gfparam gfp,int nts,struct mechparam mp,float *azi,float *scl)
{
float *zdd, *rdd, *zds, *rds, *tds, *zss, *rss, *tss;
float *axx, *ayy, *azz, *axy, *axz, *ayz;
float *bxx, *byy, *bzz, *bxy, *bxz, *byz;
float *cxx, *cyy, *czz, *cxy, *cxz, *cyz;
float *gfn, *gfe, *gfv, *gfmptr;
float f1, f2, f3, f4, f5;
float cxS, sxS, cxD, sxD, cx2D, sx2D, cxL, sxL;
float cxT, sxT, cx2T, sx2T;
float arg, cosA, sinA, rad, tan, scale;
float xamp, yamp, zamp, sx, sy, sz;
float mxx, myy, mzz, mxy, mxz, myz;
float sum, rake;
float u1, u2, u3, vx, vy, vz, l2m;
float us, ud, ux, uy, uz;
int it, ig, im;

float half = 0.5;
float two = 2.0;
float rperd = 0.017453293;

arg = (mp.dip)*rperd;
cxD = cos(arg);
sxD = sin(arg);

cx2D = cxD*cxD - sxD*sxD;
sx2D = two*sxD*cxD;

sum = 0.0;
for(im=0;im<mp.nmech;im++)
   {
   u1 = u2 = u3 = 0;
   if(mp.flag[im] == U1FLAG)
      u1 = 1;
   else if(mp.flag[im] == U2FLAG)
      u2 = 1;
   else if(mp.flag[im] == U3FLAG)
      u3 = 1;

   gfmptr = gfmech + im*12*nts;
   zapit(gfmptr,12*nts);

   if(gfp.flag3d == 0)
      {
      arg = *azi - (mp.stk)*rperd;
      cxT = cos(arg);
      sxT = sin(arg);

      cx2T = cxT*cxT - sxT*sxT;
      sx2T = two*sxT*cxT;

      if(u1)
         rake = mp.rak;
      else if(u2)
         rake = mp.rak + 90.0;

      arg = rake*rperd;
      cxL = cos(arg);
      sxL = sin(arg);

      cosA = cos(*azi);
      sinA = sin(*azi);

      for(ig=0;ig<4;ig++)
         {
	 if(gfh[ig].wt > (float)(0.0))
	    {
            gfv = gfmptr + 3*ig*nts;
            gfn = gfmptr + 3*ig*nts + nts;
            gfe = gfmptr + 3*ig*nts + 2*nts;

            zdd = gf + gfp.nc*ig*nts;
            rdd = gf + gfp.nc*ig*nts + nts;
            zds = gf + gfp.nc*ig*nts + 2*nts;
            rds = gf + gfp.nc*ig*nts + 3*nts;
            tds = gf + gfp.nc*ig*nts + 4*nts;
            zss = gf + gfp.nc*ig*nts + 5*nts;
            rss = gf + gfp.nc*ig*nts + 6*nts;
            tss = gf + gfp.nc*ig*nts + 7*nts;

            sum = sum + (*scl)*(gfh[ig].mu)*(gfh[ig].wt);
            scale = (*scl)*(gfh[ig].mu)*(gfh[ig].wt)/(gfh[ig].mom);

            f1 = half*scale*sxL*sx2D;
            f2 = -scale*(sxT*sxL*cx2D - cxT*cxL*cxD);
            f3 = -scale*(sx2T*cxL*sxD + half*cx2T*sxL*sx2D);
            f4 = scale*(cxT*sxL*cx2D + sxT*cxL*cxD);
            f5 = scale*(cx2T*cxL*sxD - half*sx2T*sxL*sx2D);

            for(it=0;it<gfh[ig].nt;it++)
               {
               gfv[it] = f1*zdd[it] + f2*zds[it] + f3*zss[it];

               rad = f1*rdd[it] + f2*rds[it] + f3*rss[it];
               tan = f4*tds[it] + f5*tss[it];

               gfn[it] = -tan*sinA + rad*cosA;
               gfe[it] = tan*cosA + rad*sinA;
               }
	    }
         }
      }
   else
      {
      arg = (mp.rak)*rperd;
      cxL = cos(arg);
      sxL = sin(arg);

      sum = 0.0;
      for(ig=0;ig<4;ig++)
         {
         arg = (mp.stk - gfh[ig].xazim)*rperd;
         cxT = cos(arg);
         sxT = sin(arg);

         vx = -sxD*sxT;
         vy =  sxD*cxT;
         vz = -cxD;

         us = u1*cxL - u2*sxL;
         ud = u1*sxL + u2*cxL;

         ux = -(u3*sxD - ud*cxD)*sxT + us*cxT;
         uy =  (u3*sxD - ud*cxD)*cxT + us*sxT;
         uz = -(u3*cxD + ud*sxD);

	 l2m = gfh[ig].lam + two*gfh[ig].mu;

         mxx = l2m*vx*ux + (gfh[ig].lam)*vy*uy + (gfh[ig].lam)*vz*uz;
         myy = (gfh[ig].lam)*vx*ux + l2m*vy*uy + (gfh[ig].lam)*vz*uz;
         mzz = (gfh[ig].lam)*vx*ux + (gfh[ig].lam)*vy*uy + l2m*vz*uz;
         mxy = (gfh[ig].mu)*(vx*uy + vy*ux);
         mxz = (gfh[ig].mu)*(vx*uz + vz*ux);
         myz = (gfh[ig].mu)*(vy*uz + vz*uy);

         arg = gfh[ig].xazim*rperd;
         cosA = cos(arg);
         sinA = sin(arg);

         gfv = gfmptr + 3*ig*nts;
         gfn = gfmptr + 3*ig*nts + nts;
         gfe = gfmptr + 3*ig*nts + 2*nts;

         axx = gf + gfp.nc*ig*nts;
         ayy = gf + gfp.nc*ig*nts + nts;
         azz = gf + gfp.nc*ig*nts + 2*nts;
         axy = gf + gfp.nc*ig*nts + 3*nts;
         axz = gf + gfp.nc*ig*nts + 4*nts;
         ayz = gf + gfp.nc*ig*nts + 5*nts;
         bxx = gf + gfp.nc*ig*nts + 6*nts;
         byy = gf + gfp.nc*ig*nts + 7*nts;
         bzz = gf + gfp.nc*ig*nts + 8*nts;
         bxy = gf + gfp.nc*ig*nts + 9*nts;
         bxz = gf + gfp.nc*ig*nts + 10*nts;
         byz = gf + gfp.nc*ig*nts + 11*nts;
         cxx = gf + gfp.nc*ig*nts + 12*nts;
         cyy = gf + gfp.nc*ig*nts + 13*nts;
         czz = gf + gfp.nc*ig*nts + 14*nts;
         cxy = gf + gfp.nc*ig*nts + 15*nts;
         cxz = gf + gfp.nc*ig*nts + 16*nts;
         cyz = gf + gfp.nc*ig*nts + 17*nts;

         sum = sum + (*scl)*(gfh[ig].mu)*(gfh[ig].wt);

         xamp = (*scl)*(gfh[ig].wt)/(gfh[ig].xmom);
         yamp = (*scl)*(gfh[ig].wt)/(gfh[ig].ymom);
         zamp = (*scl)*(gfh[ig].wt)/(gfh[ig].zmom);

         for(it=0;it<gfh[ig].nt;it++)
            {
            sx = xamp*(axx[it]*mxx + ayy[it]*myy + azz[it]*mzz
               + axy[it]*mxy + axz[it]*mxz + ayz[it]*myz);

            sy = yamp*(bxx[it]*mxx + byy[it]*myy + bzz[it]*mzz
               + bxy[it]*mxy + bxz[it]*mxz + byz[it]*myz);

            sz = zamp*(cxx[it]*mxx + cyy[it]*myy + czz[it]*mzz
               + cxy[it]*mxy + cxz[it]*mxz + cyz[it]*myz);

            gfe[it] = sx*sinA + sy*cosA;
            gfn[it] = sx*cosA - sy*sinA;
            gfv[it] = -sz;
            }
         }
      }
   }

*scl = sum;  /* scl now contains the moment released for this point source */
}

void sum_4gf_v2(float *seis,int ntout,float *gf,struct gfheader *gfh,int ntsum,int maxnt,float *rupt,float *maxgft,float *tstart,int tdflag,struct mechparam mp, float *tadjust)
{
int ig, it, im;
float backt0, t0[4];
float *sptr, *gfptr;

backt0 = 0.0;
for(ig=0;ig<4;ig++)
   {
   if(gfh[ig].wt > (float)(0.0))
      {
      t0[ig] = *maxgft - gfh[ig].gft;
      backt0 = backt0 - t0[ig]*gfh[ig].wt;
      t0[ig] = t0[ig] + gfh[ig].tst + *tadjust;
      }
   }
backt0 = backt0 + *rupt - *tstart;

for(im=0;im<mp.nmech;im++)
   {
   sptr = seis + 3*im*ntout;
   gfptr = gf + 12*im*ntsum;
   timeshift_gf_v2(sptr,ntout,gfptr,gfh,ntsum,t0,&backt0,tdflag);
   }
}

void timeshift_gf_v2(float *seis,int ntout,float *gf,struct gfheader *gfh,int ntsum,float *t0,float *bt0,int tdf)
{
struct complex *gc0, *gc1, *gc2;
float *gf0, *gf1, *gf2, *sv, *sn, *se, *gfv, *gfn, *gfe;
float cosA, sinA, arg, fac, norm, tmpre, scale, tsh;
int i, ig, tapst, it, nts3, nts6, nts9;
int itshift, it0, maxnt, ntp2, nf2;

int taplen = 10;
float zap = 0.0;
float half = 0.5;
float one = 1.0;
float two = 2.0;
float pi = 3.141592654;

maxnt = 0;

sv = seis;
sn = seis + ntout;
se = seis + 2*ntout;

if(tdf)
   {
   for(ig=0;ig<4;ig++)
      {
      if(gfh[ig].wt > (float)(0.0))
         {
         gf0 = gf + 3*ig*ntsum;
         gf1 = gf + 3*ig*ntsum + ntsum;
         gf2 = gf + 3*ig*ntsum + 2*ntsum;

   /* taper */

         tapst = gfh[ig].nt - taplen;
         arg = pi/(float)(taplen);
         for(it=tapst+1;it<gfh[ig].nt;it++)
            {
            fac = half*(one + cos(arg*(it-tapst)));
   
            gf0[it] = fac*gf0[it];
            gf1[it] = fac*gf1[it];
            gf2[it] = fac*gf2[it];
            }

/* apply time shift */

         tsh = t0[ig] + *bt0;
         if(tsh >= 0.0)
            itshift = (int)(tsh/gfh[ig].dt + 0.5);
         else
            itshift = (int)(tsh/gfh[ig].dt - 0.5);

         it0 = gfh[ig].nt + itshift;

         if(it0 > ntsum)
            it0 = ntsum;

         if(it0 > ntout)
            it0 = ntout;

         if(itshift < 0)
            {
            for(i=0;i<it0;i++)
               {
	       sv[i] = sv[i] + gf0[i-itshift];
	       sn[i] = sn[i] + gf1[i-itshift];
	       se[i] = se[i] + gf2[i-itshift];
	       }
            }
         else
            {
            for(i=it0-1;i>=itshift;i--)
               {
	       sv[i] = sv[i] + gf0[i-itshift];
	       sn[i] = sn[i] + gf1[i-itshift];
	       se[i] = se[i] + gf2[i-itshift];
	       }
            }

         if(it0 > maxnt)
            maxnt = it0;
         }
      }
   }
else
   {
   ntp2 = 2;
   while(ntp2 < ntsum)
      ntp2 = ntp2*2;

   nf2 = ntp2/2;

   for(ig=0;ig<4;ig++)
      {
      if(gfh[ig].wt > (float)(0.0))
         {
         gf0 = gf + 3*ig*ntsum;
         gf1 = gf + 3*ig*ntsum + ntsum;
         gf2 = gf + 3*ig*ntsum + 2*ntsum;

      /* taper and pad */

         tapst = gfh[ig].nt - taplen;
         arg = pi/(float)(taplen);
         for(it=tapst+1;it<gfh[ig].nt;it++)
            {
            fac = half*(one + cos(arg*(it-tapst)));
   
            gf0[it] = fac*gf0[it];
            gf1[it] = fac*gf1[it];
            gf2[it] = fac*gf2[it];
            }
   
         for(it=gfh[ig].nt;it<ntp2;it++)
            {
            gf0[it] = zap;
            gf1[it] = zap;
            gf2[it] = zap;
            }
   
   /* apply time shift */
   
         gc0 = (struct complex *) gf0;
         gc1 = (struct complex *) gf1;
         gc2 = (struct complex *) gf2;
   
         forfft(gc0,ntp2,-1);
         forfft(gc1,ntp2,-1);
         forfft(gc2,ntp2,-1);
   
         /* zero nyquist */
         gc0[0].im = zap;
         gc1[0].im = zap;
         gc2[0].im = zap;
   
         fac = -(t0[ig] + *bt0)*two*pi/(ntp2*gfh[ig].dt);
         for(i=1;i<nf2;i++)
            {
            arg = fac*i;
            cosA = cos(arg);
            sinA = sin(arg);
   
            tmpre     = gc0[i].re*cosA - gc0[i].im*sinA;
            gc0[i].im = gc0[i].re*sinA + gc0[i].im*cosA;
            gc0[i].re = tmpre;
   
            tmpre     = gc1[i].re*cosA - gc1[i].im*sinA;
            gc1[i].im = gc1[i].re*sinA + gc1[i].im*cosA;
            gc1[i].re = tmpre;
   
            tmpre     = gc2[i].re*cosA - gc2[i].im*sinA;
            gc2[i].im = gc2[i].re*sinA + gc2[i].im*cosA;
            gc2[i].re = tmpre;
            }
         }
   
   /* reset pointers to first group of 8 GFs */
   
      gf0 = gf;
      gf1 = gf + ntsum;
      gf2 = gf + 2*ntsum;
   
      gc0 = (struct complex *) gf0;
      gc1 = (struct complex *) gf1;
      gc2 = (struct complex *) gf2;
   
   /* sum over 4 GFs to interpolated response */
   
      nts3 = 3*nf2;
      nts6 = 6*nf2;
      nts9 = 9*nf2;
      for(i=0;i<nf2;i++)
         {
         gc0[i].re = gc0[i].re + gc0[i+nts3].re + gc0[i+nts6].re + gc0[i+nts9].re;
         gc0[i].im = gc0[i].im + gc0[i+nts3].im + gc0[i+nts6].im + gc0[i+nts9].im;
   
         gc1[i].re = gc1[i].re + gc1[i+nts3].re + gc1[i+nts6].re + gc1[i+nts9].re;
         gc1[i].im = gc1[i].im + gc1[i+nts3].im + gc1[i+nts6].im + gc1[i+nts9].im;
   
         gc2[i].re = gc2[i].re + gc2[i+nts3].re + gc2[i+nts6].re + gc2[i+nts9].re;
         gc2[i].im = gc2[i].im + gc2[i+nts3].im + gc2[i+nts6].im + gc2[i+nts9].im;
         }

      invfft(gc0,ntp2,1);
      invfft(gc1,ntp2,1);
      invfft(gc2,ntp2,1);

      maxnt = ntout;
      if(ntp2 < maxnt)
         maxnt = ntp2;

      gfv = gf;
      gfn = gf + ntsum;
      gfe = gf + 2*ntsum;

      scale = one/(float)(ntp2);
      for(it=0;it<maxnt;it++)
         {
         sv[it] = sv[it] + scale*gfv[it];
         sn[it] = sn[it] + scale*gfn[it];
         se[it] = se[it] + scale*gfe[it];
         }
      }
   }
}

void find_4gf_v2(struct gfparam gfp,struct gfheader *gfh,float *rng,float *dep,float *se,float *sn,float *maxgft)
{
int i, ir, ir0, ir1, id0, id1;
float sum, rr0, rr1, dd0, dd1;
float dr, dr0, dr1, drng;

if(gfp.flag3d == 0)
   {
   ir1 = 0;
   while(gfp.gfr[ir1] < *rng && ir1 < (gfp.ngfr)-1)
      ir1++;

   ir0 = ir1 - 1;

   if(ir1 == 0)
      {
      ir0 = ir1;
      fprintf(stderr,"*** range= %f < gf range limits\n",*rng);
      }

   if(ir1 == (gfp.ngfr)-1 && gfp.gfr[ir1] < *rng)
      {
      ir0 = ir1;
      fprintf(stderr,"*** range= %f > gf range limits\n",*rng);
      }

   rr0 = gfp.gfr[ir0] - *rng;
   rr1 = gfp.gfr[ir1] - *rng;
   rr0 = rr0*rr0;
   rr1 = rr1*rr1;
   }
else
   {
   ir0 = -1;
   dr0 = 1.0e+10;
   ir1 = -1;
   dr1 = 1.0e+10;
   for(ir=0;ir<gfp.ngfr;ir++)
      {
      drng = (*rng) - gfp.gfr[ir];

      if(drng < 0.0 && -drng <= gfp.rtol)
         {
	 dr = (gfp.gfn[ir]-(*sn))*(gfp.gfn[ir]-(*sn))
	       + (gfp.gfe[ir]-(*se))*(gfp.gfe[ir]-(*se));

	 if(dr < dr0)
	    {
	    dr0 = dr;
	    ir0 = ir;
	    }
	 }
      if(drng >= 0.0 && drng <= gfp.rtol)
         {
	 dr = (gfp.gfn[ir]-(*sn))*(gfp.gfn[ir]-(*sn))
	       + (gfp.gfe[ir]-(*se))*(gfp.gfe[ir]-(*se));

	 if(dr < dr1)
	    {
	    dr1 = dr;
	    ir1 = ir;
	    }
	 }
      }

   if(ir0 == -1)
      ir0 = ir1;
   if(ir1 == -1)
      ir1 = ir0;

   rr0 = (gfp.gfn[ir0]-(*sn))*(gfp.gfn[ir0]-(*sn))
          + (gfp.gfe[ir0]-(*se))*(gfp.gfe[ir0]-(*se));
   rr1 = (gfp.gfn[ir1]-(*sn))*(gfp.gfn[ir1]-(*sn))
          + (gfp.gfe[ir1]-(*se))*(gfp.gfe[ir1]-(*se));
   }

id1 = 0;
while(gfp.gfd[id1] < *dep && id1 < (gfp.ngfd)-1)
   id1++;

id0 = id1 - 1;

if(id1 == 0)
   id0 = id1;

if(id1 == (gfp.ngfd)-1 && gfp.gfd[id1] < *dep)
   id0 = id1;

dd0 = gfp.gfd[id0] - *dep;
dd1 = gfp.gfd[id1] - *dep;
dd0 = dd0*dd0;
dd1 = dd1*dd1;

if(gfh[0].id != id0 || gfh[0].ir != ir0) /* need to read-in */
   gfh[0].read_flag = 0;

if(gfh[1].id != id0 || gfh[1].ir != ir1) /* need to read-in */
   gfh[1].read_flag = 0;

if(gfh[2].id != id1 || gfh[2].ir != ir0) /* need to read-in */
   gfh[2].read_flag = 0;

if(gfh[3].id != id1 || gfh[3].ir != ir1) /* need to read-in */
   gfh[3].read_flag = 0;

/* even if read_flag=1 still set wt's since dists. may be slightly different */
   
gfh[0].id = id0; gfh[0].ir = ir0; gfh[0].wt = sqrt(dd0 + rr0);
gfh[1].id = id0; gfh[1].ir = ir1; gfh[1].wt = sqrt(dd0 + rr1);
gfh[2].id = id1; gfh[2].ir = ir0; gfh[2].wt = sqrt(dd1 + rr0);
gfh[3].id = id1; gfh[3].ir = ir1; gfh[3].wt = sqrt(dd1 + rr1);

*maxgft = -1.0e+15;
for(i=0;i<4;i++)
   {
   gfh[i].gft = gfp.gft[gfh[i].ir + (gfh[i].id)*(gfp.ngfr)];

   if(gfh[i].gft > *maxgft)
      *maxgft = gfh[i].gft;
   }

if(gfh[0].wt == 0.0)
   { gfh[0].wt = 1.0; gfh[1].wt = 0.0; gfh[2].wt = 0.0; gfh[3].wt = 0.0; }
else if(gfh[1].wt == 0.0)
   { gfh[0].wt = 0.0; gfh[1].wt = 1.0; gfh[2].wt = 0.0; gfh[3].wt = 0.0; }
else if(gfh[2].wt == 0.0)
   { gfh[0].wt = 0.0; gfh[1].wt = 0.0; gfh[2].wt = 1.0; gfh[3].wt = 0.0; }
else if(gfh[3].wt == 0.0)
   { gfh[0].wt = 0.0; gfh[1].wt = 0.0; gfh[2].wt = 0.0; gfh[3].wt = 1.0; }
else
   {
   sum = 1.0/gfh[0].wt + 1.0/gfh[1].wt + 1.0/gfh[2].wt + 1.0/gfh[3].wt;
   sum = 1.0/sum;
   for(i=0;i<4;i++)
      gfh[i].wt = sum/gfh[i].wt;
   }
}

void reset_wt_4gf_v2(struct gfheader *gfh,float *rng,float *minr,float *maxr,float *tadjust)
{
float del, twt[4];
int ir, imax;

float tt;

*tadjust = 0.0;

tt = 0.0;
for(ir=0;ir<4;ir++)
   tt = tt + gfh[ir].gft*gfh[ir].wt;

if(*rng >= *minr)
   {
   imax = 0;
   for(ir=1;ir<4;ir++)
      {
      if(gfh[ir].wt > gfh[imax].wt)
         imax = ir;
      }

   if(*rng >= *maxr)
      {
      for(ir=0;ir<4;ir++)
         gfh[ir].wt = 0.0;

      gfh[imax].wt = 1.0;
      *tadjust = tt - gfh[imax].gft;
      }
   else
      {
      del = (*rng - *minr)/(*maxr - *minr);

      for(ir=0;ir<4;ir++)
         twt[ir] = 0.0;

      twt[imax] = 1.0;

      for(ir=0;ir<4;ir++)
         gfh[ir].wt = gfh[ir].wt + del*(twt[ir] - gfh[ir].wt);
      }
   }
}
