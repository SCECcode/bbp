#include "include.h"
#include "structure.h"
#include "function.h"

void timeshift_sgt(float *seis,int ntout,float *gf,struct sgtheader *gfh,int ntsum,float *t0,float *bt0,int nsgt);

void *sgt_subset(struct sgtfileparams *sgtfpar,struct sgtfileparams *sgtextract,struct sgtmaster *sgtmast,struct sgtindex *sgtindx,int nm,long long *mindx,char *dir)
{
struct sgtmaster exmast;
struct sgtindex *exindx;
struct sgtheader sgthead;
float *sgt;
int im, ip, nreed;
off_t blen, off;
char ofile[512];

exmast.geoproj = sgtmast->geoproj;
exmast.modellon = sgtmast->modellon;
exmast.modellat = sgtmast->modellat;
exmast.modelrot = sgtmast->modelrot;
exmast.xshift = sgtmast->xshift;
exmast.yshift = sgtmast->yshift;
exmast.globnp = nm;
exmast.localnp = nm;
exmast.nt = sgtmast->nt;

exindx = (struct sgtindex *) check_malloc ((exmast.globnp)*sizeof(struct sgtindex));

for(im=0;im<nm;im++)
   {
   ip = 0;
   while(mindx[im] > sgtindx[ip].indx && ip < (sgtmast->globnp)-1)
      ip++;

   exindx[im].indx = sgtindx[ip].indx;
   exindx[im].xsgt = sgtindx[ip].xsgt;
   exindx[im].ysgt = sgtindx[ip].ysgt;
   exindx[im].zsgt = sgtindx[ip].zsgt;
   exindx[im].h = sgtindx[ip].h;

   if(mindx[im] != exindx[im].indx)
      {
      fprintf(stderr,"houston, we have a problem...\n");
      exit(-1);
      }
   }

sgt = (float *) check_malloc (6*(sgtmast->nt)*sizeof(float));
blen = (off_t)(sizeof(struct sgtheader)) + (off_t)(6*(sgtmast->nt)*sizeof(float));

if(strcmp(dir,".") != 0)
   makedir(dir);

if(sgtfpar->xfile[0] != '\0')
   {
   if(sgtextract->xfile[0] == '/')
      sprintf(ofile,"%s",sgtextract->xfile);
   else
      sprintf(ofile,"%s/%s",dir,sgtextract->xfile);

   sgtextract->xfdr = croptrfile(ofile);

   rite(sgtextract->xfdr,&exmast,sizeof(struct sgtmaster));
   rite(sgtextract->xfdr,exindx,(exmast.globnp)*sizeof(struct sgtindex));

   for(im=0;im<nm;im++)
      {
      ip = 0;
      while(mindx[im] > sgtindx[ip].indx && ip < (sgtmast->globnp)-1)
         ip++;

      off = sgtfpar->head_off + (off_t)(ip)*blen - sgtfpar->xcur_off;
      lseek(sgtfpar->xfdr,off,SEEK_CUR);
      sgtfpar->xcur_off = sgtfpar->head_off + (off_t)(ip)*blen;

      nreed = reed(sgtfpar->xfdr,&sgthead,sizeof(struct sgtheader));
      sgtfpar->xcur_off = sgtfpar->xcur_off + (off_t)(nreed);

      nreed = reed(sgtfpar->xfdr,sgt,6*(sgtmast->nt)*sizeof(float));
      sgtfpar->xcur_off = sgtfpar->xcur_off + (off_t)(nreed);

      rite(sgtextract->xfdr,&sgthead,sizeof(struct sgtheader));
      rite(sgtextract->xfdr,sgt,6*(sgtmast->nt)*sizeof(float));
      }

   close(sgtfpar->xfdr);
   close(sgtextract->xfdr);
   }

if(sgtfpar->yfile[0] != '\0')
   {
   if(sgtextract->yfile[0] == '/')
      sprintf(ofile,"%s",sgtextract->yfile);
   else
      sprintf(ofile,"%s/%s",dir,sgtextract->yfile);

   sgtextract->yfdr = croptrfile(ofile);

   rite(sgtextract->yfdr,&exmast,sizeof(struct sgtmaster));
   rite(sgtextract->yfdr,exindx,(exmast.globnp)*sizeof(struct sgtindex));

   for(im=0;im<nm;im++)
      {
      ip = 0;
      while(mindx[im] > sgtindx[ip].indx && ip < (sgtmast->globnp)-1)
         ip++;

      off = sgtfpar->head_off + (off_t)(ip)*blen - sgtfpar->ycur_off;
      lseek(sgtfpar->yfdr,off,SEEK_CUR);
      sgtfpar->ycur_off = sgtfpar->head_off + (off_t)(ip)*blen;

      nreed = reed(sgtfpar->yfdr,&sgthead,sizeof(struct sgtheader));
      sgtfpar->ycur_off = sgtfpar->ycur_off + (off_t)(nreed);

      nreed = reed(sgtfpar->yfdr,sgt,6*(sgtmast->nt)*sizeof(float));
      sgtfpar->ycur_off = sgtfpar->ycur_off + (off_t)(nreed);

      rite(sgtextract->yfdr,&sgthead,sizeof(struct sgtheader));
      rite(sgtextract->yfdr,sgt,6*(sgtmast->nt)*sizeof(float));
      }

   close(sgtfpar->yfdr);
   close(sgtextract->yfdr);
   }

if(sgtfpar->zfile[0] != '\0')
   {
   if(sgtextract->zfile[0] == '/')
      sprintf(ofile,"%s",sgtextract->zfile);
   else
      sprintf(ofile,"%s/%s",dir,sgtextract->zfile);

   sgtextract->zfdr = croptrfile(ofile);

   rite(sgtextract->zfdr,&exmast,sizeof(struct sgtmaster));
   rite(sgtextract->zfdr,exindx,(exmast.globnp)*sizeof(struct sgtindex));

   for(im=0;im<nm;im++)
      {
      ip = 0;
      while(mindx[im] > sgtindx[ip].indx && ip < (sgtmast->globnp)-1)
         ip++;

      off = sgtfpar->head_off + (off_t)(ip)*blen - sgtfpar->zcur_off;
      lseek(sgtfpar->zfdr,off,SEEK_CUR);
      sgtfpar->zcur_off = sgtfpar->head_off + (off_t)(ip)*blen;

      nreed = reed(sgtfpar->zfdr,&sgthead,sizeof(struct sgtheader));
      sgtfpar->zcur_off = sgtfpar->zcur_off + (off_t)(nreed);

      nreed = reed(sgtfpar->zfdr,sgt,6*(sgtmast->nt)*sizeof(float));
      sgtfpar->zcur_off = sgtfpar->zcur_off + (off_t)(nreed);

      rite(sgtextract->zfdr,&sgthead,sizeof(struct sgtheader));
      rite(sgtextract->zfdr,sgt,6*(sgtmast->nt)*sizeof(float));
      }

   close(sgtfpar->zfdr);
   close(sgtextract->zfdr);
   }
}

struct sgtindex *get_sgtpars(struct sgtfileparams *sgtfpar,struct sgtmaster *sgtmast,struct sgtindex *sgtindx)
{
struct sgtmaster tmast;
struct sgtindex *tindx;
int ip;

int eflag = 0;
int rflag = 0;

if(sgtfpar->xfile[0] != '\0')
   {
   sgtfpar->xfdr = opfile_ro(sgtfpar->xfile);

   reed(sgtfpar->xfdr,sgtmast,sizeof(struct sgtmaster));
   sgtindx = (struct sgtindex *) check_malloc ((sgtmast->globnp)*sizeof(struct sgtindex));

   reed(sgtfpar->xfdr,sgtindx,(sgtmast->globnp)*sizeof(struct sgtindex));

   sgtfpar->head_off = (off_t)(sizeof(struct sgtmaster)) + (off_t)((sgtmast->globnp)*sizeof(struct sgtindex));
   sgtfpar->xcur_off = sgtfpar->head_off;

   rflag = 1;
   }

if(sgtfpar->yfile[0] != '\0')
   {
   sgtfpar->yfdr = opfile_ro(sgtfpar->yfile);

   if(rflag == 0)
      {
      reed(sgtfpar->yfdr,sgtmast,sizeof(struct sgtmaster));
      sgtindx = (struct sgtindex *) check_malloc ((sgtmast->globnp)*sizeof(struct sgtindex));
      reed(sgtfpar->yfdr,sgtindx,(sgtmast->globnp)*sizeof(struct sgtindex));

      sgtfpar->head_off = (off_t)(sizeof(struct sgtmaster)) + (off_t)((sgtmast->globnp)*sizeof(struct sgtindex));
      }
   else
      {
      reed(sgtfpar->yfdr,&tmast,sizeof(struct sgtmaster));

      eflag = 0;
      if(tmast.geoproj != sgtmast->geoproj)
         eflag = -1;
      if(tmast.modellon > sgtmast->modellon + 0.001 || tmast.modellon < sgtmast->modellon - 0.001)
         eflag = -1;
      if(tmast.modellat > sgtmast->modellat + 0.001 || tmast.modellat < sgtmast->modellat - 0.001)
         eflag = -1;
      if(tmast.modelrot > sgtmast->modelrot + 0.001 || tmast.modelrot < sgtmast->modelrot - 0.001)
         eflag = -1;
      if(tmast.xshift > sgtmast->xshift + 0.001 || tmast.xshift < sgtmast->xshift - 0.001)
         eflag = -1;
      if(tmast.yshift > sgtmast->yshift + 0.001 || tmast.yshift < sgtmast->yshift - 0.001)
         eflag = -1;
      if(tmast.globnp != sgtmast->globnp)
         eflag = -1;
      if(tmast.localnp != sgtmast->localnp)
         eflag = -1;
      if(tmast.nt != sgtmast->nt)
         eflag = -1;

      if(eflag != 0)
         {
	 fprintf(stderr,"sgtmaster inconsistency in yfile= %s, exiting ...\n",sgtfpar->yfile);
	 exit(-1);
	 }

      tindx = (struct sgtindex *) check_malloc ((tmast.globnp)*sizeof(struct sgtindex));
      reed(sgtfpar->yfdr,tindx,(tmast.globnp)*sizeof(struct sgtindex));

      for(ip=0;ip<sgtmast->globnp;ip++)
         {
         if(tindx[ip].indx != sgtindx[ip].indx)
            eflag = -1;
	 }

      if(eflag != 0)
         {
	 fprintf(stderr,"sgtindex inconsistency in yfile= %s, exiting ...\n",sgtfpar->yfile);
	 exit(-1);
	 }

      free(tindx);
      }

   sgtfpar->ycur_off = sgtfpar->head_off;
   rflag = 1;
   }

if(sgtfpar->zfile[0] != '\0')
   {
   sgtfpar->zfdr = opfile_ro(sgtfpar->zfile);

   if(rflag == 0)
      {
      reed(sgtfpar->zfdr,sgtmast,sizeof(struct sgtmaster));
      sgtindx = (struct sgtindex *) check_malloc ((sgtmast->globnp)*sizeof(struct sgtindex));
      reed(sgtfpar->zfdr,sgtindx,(sgtmast->globnp)*sizeof(struct sgtindex));

      sgtfpar->head_off = (off_t)(sizeof(struct sgtmaster)) + (off_t)((sgtmast->globnp)*sizeof(struct sgtindex));
      }
   else
      {
      reed(sgtfpar->zfdr,&tmast,sizeof(struct sgtmaster));

      eflag = 0;
      if(tmast.geoproj != sgtmast->geoproj)
         eflag = -1;
      if(tmast.modellon > sgtmast->modellon + 0.001 || tmast.modellon < sgtmast->modellon - 0.001)
         eflag = -1;
      if(tmast.modellat > sgtmast->modellat + 0.001 || tmast.modellat < sgtmast->modellat - 0.001)
         eflag = -1;
      if(tmast.modelrot > sgtmast->modelrot + 0.001 || tmast.modelrot < sgtmast->modelrot - 0.001)
         eflag = -1;
      if(tmast.xshift > sgtmast->xshift + 0.001 || tmast.xshift < sgtmast->xshift - 0.001)
         eflag = -1;
      if(tmast.yshift > sgtmast->yshift + 0.001 || tmast.yshift < sgtmast->yshift - 0.001)
         eflag = -1;
      if(tmast.globnp != sgtmast->globnp)
         eflag = -1;
      if(tmast.localnp != sgtmast->localnp)
         eflag = -1;
      if(tmast.nt != sgtmast->nt)
         eflag = -1;

      if(eflag != 0)
         {
	 fprintf(stderr,"sgtmaster inconsistency in zfile= %s, exiting ...\n",sgtfpar->zfile);
	 exit(-1);
	 }

      tindx = (struct sgtindex *) check_malloc ((tmast.globnp)*sizeof(struct sgtindex));
      reed(sgtfpar->zfdr,tindx,(tmast.globnp)*sizeof(struct sgtindex));

      for(ip=0;ip<sgtmast->globnp;ip++)
         {
         if(tindx[ip].indx != sgtindx[ip].indx)
            eflag = -1;
	 }

      if(eflag != 0)
         {
	 fprintf(stderr,"sgtindex inconsistency in zfile= %s, exiting ...\n",sgtfpar->zfile);
	 exit(-1);
	 }

      free(tindx);
      }

   sgtfpar->zcur_off = sgtfpar->head_off;
   rflag = 1;
   }

return(sgtindx);
}

void find_sgt(struct sgtparams *sgtpar,struct sgtmaster *sgtmast,struct sgtindex *sgtindx,struct sgtindex *eqindx,struct sgtindex *statindx,float *maxd,float *fwt)
{
float xx, yy, zz, rng, zdp, rexact, zexact, del, delta[4];
float sum, mind, xwt;
int ip, i;
int p0, p1, p2, p3;
int zflag;

/* first see if there is an exact match */

ip = 0;
while(eqindx->indx > sgtindx[ip].indx && ip < (sgtmast->globnp)-1)
   ip++;

if(eqindx->indx == sgtindx[ip].indx)   /* great, exact match, this will be easy */
   {
   sgtpar->nsgt = 1;
   sgtpar->indx[0] = sgtindx[ip].indx;
   sgtpar->wt[0] = 1.0;

   /*
   fprintf(stderr,"SGT EXACT: eqindx= %Ld ip= %d\n",eqindx->indx,ip);
   */
   }
else  /* more difficult, find up to 4 SGT that bracket point in range and depth */
   {
   sgtpar->nsgt = 0;
   p0 = -1; p1 = -1; p2 = -1; p3 = -1;

   delta[0] = 1.0e+15;
   delta[1] = 1.0e+15;
   delta[2] = 1.0e+15;
   delta[3] = 1.0e+15;

   xx = (eqindx->xsgt - statindx->xsgt);
   yy = (eqindx->ysgt - statindx->ysgt);
   rexact = xx*xx + yy*yy;
   zexact = (eqindx->zsgt - statindx->zsgt);

   for(ip=0;ip<sgtmast->globnp;ip++)
      {
      xx = (sgtindx[ip].xsgt - statindx->xsgt);
      yy = (sgtindx[ip].ysgt - statindx->ysgt);
      rng = xx*xx + yy*yy;
      zdp = (sgtindx[ip].zsgt - statindx->zsgt);

      xx = (sgtindx[ip].xsgt - eqindx->xsgt);
      yy = (sgtindx[ip].ysgt - eqindx->ysgt);
      zz = (sgtindx[ip].zsgt - eqindx->zsgt);
      del = xx*xx + yy*yy + zz*zz;

      if(rng <= rexact && zdp <= zexact)
         {
	 if(p0 < 0)   /* first time here */
	    {
	    p0 = sgtpar->nsgt;
	    sgtpar->nsgt = sgtpar->nsgt + 1;
	    }

	 if(del < delta[p0])
	    {
	    delta[p0] = del;
	    sgtpar->indx[p0] = sgtindx[ip].indx;
	    sgtpar->wt[p0] = sqrt(del);
	    }
	 }
      else if(rng > rexact && zdp <= zexact)
         {
	 if(p1 < 0)   /* first time here */
	    {
	    p1 = sgtpar->nsgt;
	    sgtpar->nsgt = sgtpar->nsgt + 1;
	    }

	 if(del < delta[p1])
	    {
	    delta[p1] = del;
	    sgtpar->indx[p1] = sgtindx[ip].indx;
	    sgtpar->wt[p1] = sqrt(del);
	    }
	 }
      else if(rng <= rexact && zdp > zexact)
         {
	 if(p2 < 0)   /* first time here */
	    {
	    p2 = sgtpar->nsgt;
	    sgtpar->nsgt = sgtpar->nsgt + 1;
	    }

	 if(del < delta[p2])
	    {
	    delta[p2] = del;
	    sgtpar->indx[p2] = sgtindx[ip].indx;
	    sgtpar->wt[p2] = sqrt(del);
	    }
	 }
      else /* should be (rng > rexact && zdp > zexact)  */
         {
	 if(p3 < 0)   /* first time here */
	    {
	    p3 = sgtpar->nsgt;
	    sgtpar->nsgt = sgtpar->nsgt + 1;
	    }

	 if(del < delta[p3])
	    {
	    delta[p3] = del;
	    sgtpar->indx[p3] = sgtindx[ip].indx;
	    sgtpar->wt[p3] = sqrt(del);
	    }
	 }
      }

   zflag = -1;
   for(i=0;i<sgtpar->nsgt;i++)
      {
      if(sgtpar->wt[i] == 0.0)
         zflag = i;
      }

   if(zflag >= 0)
      {
      for(i=0;i<sgtpar->nsgt;i++)
         sgtpar->wt[i] = 0.0;

      sgtpar->wt[zflag] = 1.0;
      }
   else
      {
      sum = 0.0;
      for(i=0;i<sgtpar->nsgt;i++)
         sum = sum + 1.0/sgtpar->wt[i];

      sum = 1.0/sum;
      for(i=0;i<sgtpar->nsgt;i++)
         sgtpar->wt[i] = sum/sgtpar->wt[i];
      }

   mind = 1.0e+15;
   for(i=0;i<sgtpar->nsgt;i++)
      {
      if(delta[i] < mind)
	 {
         mind = delta[i];
	 xwt = sgtpar->wt[i];
	 }
      }
   if(mind > *maxd)
      {
      *maxd = mind;
      *fwt = xwt;
      }

   /*
   for(i=0;i<sgtpar->nsgt;i++)
      {
      if(delta[i] > *maxd)
         *maxd = delta[i];
      }
      */

   if(sgtpar->nsgt < 4)
      fprintf(stderr,"*** tried to find 4 SGT, but only found %d: eq.zsgt= %d\n",sgtpar->nsgt,eqindx->zsgt);

      /*
      {
for(i=0;i<sgtpar->nsgt;i++)
   {
   fprintf(stderr,"%d) sgti=%Ld eqi=%Ld delta=%13.5e maxd=%13.5e\n",i,sgtpar->indx[i],eqindx->indx,delta[i],*maxd);
   }
   exit(-1);
   }
*/

   }
}

void read_sgt(struct sgtfileparams *sgtfpar,struct sgtmaster *sgtmast,struct sgtindex *sgtindx,struct sgtheader *sgthead,float *sgtbuf)
{
int ip;
float *sgtptr;

float xmom = 0.0;
float ymom = 0.0;

for(ip=0;ip<18*(sgtmast->globnp)*(sgtmast->nt);ip++)
   sgtbuf[ip] = 0.0;

if(sgtfpar->xfile[0] != '\0')
   {
   lseek(sgtfpar->xfdr,(sgtfpar->head_off),SEEK_SET);

   for(ip=0;ip<(sgtmast->globnp);ip++)
      {
      sgtptr = sgtbuf + ip*18*(sgtmast->nt);
      reed(sgtfpar->xfdr,&sgthead[ip],sizeof(struct sgtheader));
      reed(sgtfpar->xfdr,sgtptr,6*(sgtmast->nt)*sizeof(float));
      }
   close(sgtfpar->xfdr);

   xmom = sgthead[0].xmom;
   }

if(sgtfpar->yfile[0] != '\0')
   {
   lseek(sgtfpar->yfdr,(sgtfpar->head_off),SEEK_SET);

   for(ip=0;ip<(sgtmast->globnp);ip++)
      {
      sgtptr = sgtbuf + (ip*18 + 6)*(sgtmast->nt);
      reed(sgtfpar->yfdr,&sgthead[ip],sizeof(struct sgtheader));
      reed(sgtfpar->yfdr,sgtptr,6*(sgtmast->nt)*sizeof(float));

      sgthead[ip].xmom = xmom;
      }
   close(sgtfpar->yfdr);

   ymom = sgthead[0].ymom;
   }

if(sgtfpar->zfile[0] != '\0')
   {
   lseek(sgtfpar->zfdr,(sgtfpar->head_off),SEEK_SET);

   for(ip=0;ip<(sgtmast->globnp);ip++)
      {
      sgtptr = sgtbuf + (ip*18 + 12)*(sgtmast->nt);
      reed(sgtfpar->zfdr,&sgthead[ip],sizeof(struct sgtheader));
      reed(sgtfpar->zfdr,sgtptr,6*(sgtmast->nt)*sizeof(float));

      sgthead[ip].xmom = xmom;
      sgthead[ip].ymom = ymom;
      }
   close(sgtfpar->zfdr);
   }
}

void sum_sgt(float *seis,int ntout,float *gfmech,struct sgtparams *sgtpar,struct sgtheader *sgthead,int ntsum,float *rupt,float *tstart,struct mechparam mp)
{
int ig, ip, it, im;
float pbar, maxgft, backt0, t0[4], gft[4];
float *sptr, *gfptr;

pbar = 0.0;
for(ig=0;ig<sgtpar->nsgt;ig++)
   {
   ip = sgtpar->master_ip[ig];
   pbar = pbar + sqrt(sgthead[ip].rho/sgthead[ip].mu);  /* slowness */
   }
pbar = pbar/(float)(sgtpar->nsgt);

maxgft = -1.0e+15;
for(ig=0;ig<sgtpar->nsgt;ig++)
   {
   ip = sgtpar->master_ip[ig];
   gft[ig] = pbar*sgthead[ip].cdist;
   if(gft[ig] > maxgft)
      maxgft = gft[ig];
   }

backt0 = 0.0;
for(ig=0;ig<sgtpar->nsgt;ig++)
   {
   ip = sgtpar->master_ip[ig];
   t0[ig] = maxgft - gft[ig];
   backt0 = backt0 - t0[ig]*sgtpar->wt[ig];
   t0[ig] = t0[ig] + sgthead[ip].tst;
   }
backt0 = backt0 + *rupt - *tstart;

for(im=0;im<mp.nmech;im++)
   {
   sptr = seis + 3*im*ntout;
   gfptr = gfmech + 12*im*ntsum;
   timeshift_sgt(sptr,ntout,gfptr,sgthead,ntsum,t0,&backt0,sgtpar->nsgt);
   }
}

void timeshift_sgt(float *seis,int ntout,float *gf,struct sgtheader *gfh,int ntsum,float *t0,float *bt0,int nsgt)
{
struct complex *gc0, *gc1, *gc2;
float *gf0, *gf1, *gf2, *sv, *sn, *se, *gfv, *gfn, *gfe;
float cosA, sinA, arg, fac, norm, tmpre, scale, tsh;
int i, ig, tapst, it, nts3, nts6, nts9;
int itshift, it0, nf2;

int taplen = 10;
float zap = 0.0;
float half = 0.5;
float one = 1.0;
float two = 2.0;
float pi = 3.141592654;

sv = seis;
sn = seis + ntout;
se = seis + 2*ntout;

for(ig=0;ig<nsgt;ig++)
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
   }
}

void mech_sgt(float *gfmech,float *sgtbuf,struct sgtheader *sgthead,struct sgtparams *sgtpar,int nts,struct mechparam mp,float *scl)
{
struct sgtheader *sgtheadptr;
float *sgtbufptr;
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

   arg = (mp.rak)*rperd;
   cxL = cos(arg);
   sxL = sin(arg);

   sum = 0.0;
   for(ig=0;ig<(sgtpar->nsgt);ig++)
      {
      sgtheadptr = sgthead + sgtpar->master_ip[ig];
      sgtbufptr = sgtbuf + 18*sgtheadptr->nt*sgtpar->master_ip[ig];

      arg = (mp.stk - sgtheadptr->xazim)*rperd;
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

      l2m = sgtheadptr->lam + two*sgtheadptr->mu;

      mxx = l2m*vx*ux + (sgtheadptr->lam)*vy*uy + (sgtheadptr->lam)*vz*uz;
      myy = (sgtheadptr->lam)*vx*ux + l2m*vy*uy + (sgtheadptr->lam)*vz*uz;
      mzz = (sgtheadptr->lam)*vx*ux + (sgtheadptr->lam)*vy*uy + l2m*vz*uz;
      mxy = (sgtheadptr->mu)*(vx*uy + vy*ux);
      mxz = (sgtheadptr->mu)*(vx*uz + vz*ux);
      myz = (sgtheadptr->mu)*(vy*uz + vz*uy);

      arg = sgtheadptr->xazim*rperd;
      cosA = cos(arg);
      sinA = sin(arg);

      gfv = gfmptr + 3*ig*nts;
      gfn = gfmptr + 3*ig*nts + nts;
      gfe = gfmptr + 3*ig*nts + 2*nts;

      axx = sgtbufptr;
      ayy = sgtbufptr + (sgtheadptr->nt);
      azz = sgtbufptr + 2*(sgtheadptr->nt);
      axy = sgtbufptr + 3*(sgtheadptr->nt);
      axz = sgtbufptr + 4*(sgtheadptr->nt);
      ayz = sgtbufptr + 5*(sgtheadptr->nt);
      bxx = sgtbufptr + 6*(sgtheadptr->nt);
      byy = sgtbufptr + 7*(sgtheadptr->nt);
      bzz = sgtbufptr + 8*(sgtheadptr->nt);
      bxy = sgtbufptr + 9*(sgtheadptr->nt);
      bxz = sgtbufptr + 10*(sgtheadptr->nt);
      byz = sgtbufptr + 11*(sgtheadptr->nt);
      cxx = sgtbufptr + 12*(sgtheadptr->nt);
      cyy = sgtbufptr + 13*(sgtheadptr->nt);
      czz = sgtbufptr + 14*(sgtheadptr->nt);
      cxy = sgtbufptr + 15*(sgtheadptr->nt);
      cxz = sgtbufptr + 16*(sgtheadptr->nt);
      cyz = sgtbufptr + 17*(sgtheadptr->nt);

      sum = sum + (*scl)*(sgtheadptr->mu)*(sgtpar->wt[ig]);

      xamp = 0.0;
      if(sgtheadptr->xmom > 0.0)
         xamp = (*scl)*(sgtpar->wt[ig])/(sgtheadptr->xmom);
	 
      yamp = 0.0;
      if(sgtheadptr->ymom > 0.0)
         yamp = (*scl)*(sgtpar->wt[ig])/(sgtheadptr->ymom);
	 
      zamp = 0.0;
      if(sgtheadptr->zmom > 0.0)
         zamp = (*scl)*(sgtpar->wt[ig])/(sgtheadptr->zmom);

      for(it=0;it<sgtheadptr->nt;it++)
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

*scl = sum;  /* scl now contains the moment released for this point source */
}

void *get_sgt_subset(struct sgtfileparams *sgtfpar,struct sgtmaster *sgtmast,struct sgtindex *sgtindx,struct sgtheader *sgthead,int nm,long long *mindx,float *sgtbuf)
{
struct sgtmaster exmast;
struct sgtindex *exindx;
float *sgtptr;
int im, ip, nreed;
off_t blen, off;

float xmom = 0.0;
float ymom = 0.0;
float zmom = 0.0;

exmast.geoproj = sgtmast->geoproj;
exmast.modellon = sgtmast->modellon;
exmast.modellat = sgtmast->modellat;
exmast.modelrot = sgtmast->modelrot;
exmast.xshift = sgtmast->xshift;
exmast.yshift = sgtmast->yshift;
exmast.globnp = nm;
exmast.localnp = nm;
exmast.nt = sgtmast->nt;

for(ip=0;ip<18*nm*(sgtmast->nt);ip++)
   sgtbuf[ip] = 0.0;

exindx = (struct sgtindex *) check_malloc ((exmast.globnp)*sizeof(struct sgtindex));

for(im=0;im<nm;im++)
   {
   ip = 0;
   while(mindx[im] > sgtindx[ip].indx && ip < (sgtmast->globnp)-1)
      ip++;

   exindx[im].indx = sgtindx[ip].indx;
   exindx[im].xsgt = sgtindx[ip].xsgt;
   exindx[im].ysgt = sgtindx[ip].ysgt;
   exindx[im].zsgt = sgtindx[ip].zsgt;
   exindx[im].h = sgtindx[ip].h;

   if(mindx[im] != exindx[im].indx)
      {
      fprintf(stderr,"houston, we have a problem...\n");
      exit(-1);
      }
   }

blen = (off_t)(sizeof(struct sgtheader)) + (off_t)(6*(sgtmast->nt)*sizeof(float));

if(sgtfpar->xfile[0] != '\0')
   {
   for(im=0;im<nm;im++)
      {
      ip = 0;
      while(mindx[im] > sgtindx[ip].indx && ip < (sgtmast->globnp)-1)
         ip++;

      off = sgtfpar->head_off + (off_t)(ip)*blen - sgtfpar->xcur_off;
      lseek(sgtfpar->xfdr,off,SEEK_CUR);
      sgtfpar->xcur_off = sgtfpar->head_off + (off_t)(ip)*blen;

      nreed = reed(sgtfpar->xfdr,&sgthead[im],sizeof(struct sgtheader));
      sgtfpar->xcur_off = sgtfpar->xcur_off + (off_t)(nreed);

      sgtptr = sgtbuf + im*18*(sgtmast->nt);
      nreed = reed(sgtfpar->xfdr,sgtptr,6*(sgtmast->nt)*sizeof(float));
      sgtfpar->xcur_off = sgtfpar->xcur_off + (off_t)(nreed);
      }

   close(sgtfpar->xfdr);

   xmom = sgthead[0].xmom;
   }

if(sgtfpar->yfile[0] != '\0')
   {
   for(im=0;im<nm;im++)
      {
      ip = 0;
      while(mindx[im] > sgtindx[ip].indx && ip < (sgtmast->globnp)-1)
         ip++;

      off = sgtfpar->head_off + (off_t)(ip)*blen - sgtfpar->ycur_off;
      lseek(sgtfpar->yfdr,off,SEEK_CUR);
      sgtfpar->ycur_off = sgtfpar->head_off + (off_t)(ip)*blen;

      nreed = reed(sgtfpar->yfdr,&sgthead[im],sizeof(struct sgtheader));
      sgtfpar->ycur_off = sgtfpar->ycur_off + (off_t)(nreed);

      sgtptr = sgtbuf + (im*18 + 6)*(sgtmast->nt);
      nreed = reed(sgtfpar->yfdr,sgtptr,6*(sgtmast->nt)*sizeof(float));
      sgtfpar->ycur_off = sgtfpar->ycur_off + (off_t)(nreed);
      }

   close(sgtfpar->yfdr);

   ymom = sgthead[0].ymom;
   }

if(sgtfpar->zfile[0] != '\0')
   {
   for(im=0;im<nm;im++)
      {
      ip = 0;
      while(mindx[im] > sgtindx[ip].indx && ip < (sgtmast->globnp)-1)
         ip++;

      off = sgtfpar->head_off + (off_t)(ip)*blen - sgtfpar->zcur_off;
      lseek(sgtfpar->zfdr,off,SEEK_CUR);
      sgtfpar->zcur_off = sgtfpar->head_off + (off_t)(ip)*blen;

      nreed = reed(sgtfpar->zfdr,&sgthead[im],sizeof(struct sgtheader));
      sgtfpar->zcur_off = sgtfpar->zcur_off + (off_t)(nreed);

      sgtptr = sgtbuf + (im*18 + 12)*(sgtmast->nt);
      nreed = reed(sgtfpar->zfdr,sgtptr,6*(sgtmast->nt)*sizeof(float));
      sgtfpar->zcur_off = sgtfpar->zcur_off + (off_t)(nreed);
      }

   close(sgtfpar->zfdr);

   zmom = sgthead[0].zmom;
   }

sgtmast->globnp = exmast.globnp;
sgtmast->localnp = exmast.localnp;

for(im=0;im<nm;im++)
   {
   sgtindx[im].indx = exindx[im].indx;
   sgtindx[im].xsgt = exindx[im].xsgt;
   sgtindx[im].ysgt = exindx[im].ysgt;
   sgtindx[im].zsgt = exindx[im].zsgt;
   sgtindx[im].h = exindx[im].h;

/* replace just in case values were read-over */
   sgthead[im].xmom = xmom;
   sgthead[im].ymom = ymom;
   sgthead[im].zmom = zmom;
   }
}

void point_mech_sgt(float *gfmech,float *sgtbuf,struct sgtheader *sgthead,struct sgtparams *sgtpar,int nts,struct mechparam mp,float *mom)
{
struct sgtheader *sgtheadptr;
float *sgtbufptr;
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

for(im=0;im<mp.nmech;im++)
   {
   gfmptr = gfmech + im*12*nts;
   zapit(gfmptr,12*nts);

   arg = (mp.rak)*rperd;
   cxL = cos(arg);
   sxL = sin(arg);

   sum = 0.0;
   for(ig=0;ig<(sgtpar->nsgt);ig++)
      {
      sgtheadptr = sgthead + sgtpar->master_ip[ig];
      sgtbufptr = sgtbuf + 18*sgtheadptr->nt*sgtpar->master_ip[ig];

      arg = (mp.stk - sgtheadptr->xazim)*rperd;
      cxT = cos(arg);
      sxT = sin(arg);

      cx2T = cxT*cxT - sxT*sxT;
      sx2T = two*sxT*cxT;

      mxx = -(sxD*cxL*sx2T + sx2D*sxL*sxT*sxT);
      myy = (sxD*cxL*sx2T - sx2D*sxL*cxT*cxT);
      mzz = sx2D*sxL;
      mxy = (sxD*cxL*cx2T + half*sx2D*sxL*sx2T);
      mxz = -(cxD*cxL*cxT + cx2D*sxL*sxT);
      myz = -(cxD*cxL*sxT - cx2D*sxL*cxT);

/*
fprintf(stderr,"mxx= %13.5e\n",mxx);
fprintf(stderr,"myy= %13.5e\n",myy);
fprintf(stderr,"mzz= %13.5e\n",mzz);
fprintf(stderr,"mxy= %13.5e\n",mxy);
fprintf(stderr,"mxz= %13.5e\n",mxz);
fprintf(stderr,"myz= %13.5e\n",myz);
*/

      arg = sgtheadptr->xazim*rperd;
      cosA = cos(arg);
      sinA = sin(arg);

      gfv = gfmptr + 3*ig*nts;
      gfn = gfmptr + 3*ig*nts + nts;
      gfe = gfmptr + 3*ig*nts + 2*nts;

      axx = sgtbufptr;
      ayy = sgtbufptr + (sgtheadptr->nt);
      azz = sgtbufptr + 2*(sgtheadptr->nt);
      axy = sgtbufptr + 3*(sgtheadptr->nt);
      axz = sgtbufptr + 4*(sgtheadptr->nt);
      ayz = sgtbufptr + 5*(sgtheadptr->nt);
      bxx = sgtbufptr + 6*(sgtheadptr->nt);
      byy = sgtbufptr + 7*(sgtheadptr->nt);
      bzz = sgtbufptr + 8*(sgtheadptr->nt);
      bxy = sgtbufptr + 9*(sgtheadptr->nt);
      bxz = sgtbufptr + 10*(sgtheadptr->nt);
      byz = sgtbufptr + 11*(sgtheadptr->nt);
      cxx = sgtbufptr + 12*(sgtheadptr->nt);
      cyy = sgtbufptr + 13*(sgtheadptr->nt);
      czz = sgtbufptr + 14*(sgtheadptr->nt);
      cxy = sgtbufptr + 15*(sgtheadptr->nt);
      cxz = sgtbufptr + 16*(sgtheadptr->nt);
      cyz = sgtbufptr + 17*(sgtheadptr->nt);

      xamp = 0.0;
      if(sgtheadptr->xmom > 0.0)
         xamp = (*mom)*(sgtpar->wt[ig])/(sgtheadptr->xmom);
	 
      yamp = 0.0;
      if(sgtheadptr->ymom > 0.0)
         yamp = (*mom)*(sgtpar->wt[ig])/(sgtheadptr->ymom);
	 
      zamp = 0.0;
      if(sgtheadptr->zmom > 0.0)
         zamp = (*mom)*(sgtpar->wt[ig])/(sgtheadptr->zmom);

      for(it=0;it<sgtheadptr->nt;it++)
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
