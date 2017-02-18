#include "include.h"
#include "structure.h"
#include "function.h"
#include "defs.h"

void init_plane_srf(struct standrupformat *srf,struct generic_slip *gslip,float *elon,float *elat,int nx,int ny,float *fl,float *fw,float *dx,float *dy,float *stk,float *dip,float *dtop,float *sh,float *dh)
{
struct srf_prectsegments *prseg_ptr;
float avglon, avglat, dmin, dold, se, sn;
int ig, ip_seg, nseg, i, ip;
struct slippars *spar;

double dbllon, dbllat, dbldlen, dbldwid, dblstk, dbldip;

float rperd = 0.017453293;

sprintf(srf[0].type,"PLANE");

if(gslip->np > 0)
   {
   spar = gslip->spar;

   nseg = 0;
   for(i=0;i<gslip->np;i++)
      {
      if(spar[i].segno > nseg)
         nseg = spar[i].segno;
      }
   nseg++;

   srf[0].srf_prect.nseg = nseg;
   srf[0].srf_prect.prectseg = (struct srf_prectsegments *)check_malloc(srf[0].srf_prect.nseg*sizeof(struct srf_prectsegments));

   prseg_ptr = srf[0].srf_prect.prectseg;

   for(i=0;i<nseg;i++)
      {
      dbllon = 0.0; 
      dbllat = 0.0; 
      dblstk = 0.0;
      dbldip = 0.0;
      dbldlen = 0.0;
      dbldwid = 0.0;

      ip_seg = 0;
      dold = -1;

      prseg_ptr[i].nstk = 0;
      prseg_ptr[i].ndip = 0;
      for(ip=0;ip<gslip->np;ip++)
         {
         if(spar[ip].segno == i)
            {
            if(ip_seg == 0)
               dmin = spar[ip].dep;

            if(spar[ip].dep == dmin)
               {
               dbllon = dbllon + spar[ip].lon;
               dbllat = dbllat + spar[ip].lat;
               prseg_ptr[i].nstk++;
               }

            if(spar[ip].dep != dold)
               {
               dold = spar[ip].dep;
               prseg_ptr[i].ndip++;
               }

            ip_seg++;

            if(ip_seg == 1)
               {
               dblstk = spar[ip].ds;
               dbldip = spar[ip].dw;
               dblstk = spar[ip].stk;
               dbldip = spar[ip].dip;
               }
            else
               {
               dbldlen = dbldlen + spar[ip].ds;
               dbldwid = dbldwid + spar[ip].dw;

               if(spar[ip].stk > dblstk/(float)(ip_seg-1.0) + 90.0)
                  {
                  dblstk = dblstk + spar[ip].stk - 180.0;
                  dbldip = dbldip + 180.0 - spar[ip].dip;
                  }
               else if(spar[ip].stk < dblstk/(float)(ip_seg-1.0) - 90.0)
                  {
                  dblstk = dblstk + spar[ip].stk + 180.0;
                  dbldip = dbldip + 180.0 - spar[ip].dip;
                  }
               else
                  {
                  dblstk = dblstk + spar[ip].stk;
                  dbldip = dbldip + spar[ip].dip;
                  }
               }
            }
         }

      prseg_ptr[i].dlen = dbldlen/(float)(ip_seg);
      prseg_ptr[i].dwid = dbldwid/(float)(ip_seg);

      prseg_ptr[i].stk = dblstk/(float)(ip_seg);
      prseg_ptr[i].dip = dbldip/(float)(ip_seg);

      if(prseg_ptr[i].dip > 90.0)
         {
         prseg_ptr[i].dip = 180.0 - prseg_ptr[i].dip;
         prseg_ptr[i].stk = prseg_ptr[i].stk - 180.0;
         }

      while(prseg_ptr[i].stk < 0.0)
         prseg_ptr[i].stk = prseg_ptr[i].stk + 360.0;
      while(prseg_ptr[i].stk >= 360.0)
         prseg_ptr[i].stk = prseg_ptr[i].stk - 360.0;

      prseg_ptr[i].flen = prseg_ptr[i].nstk*prseg_ptr[i].dlen;
      prseg_ptr[i].fwid = prseg_ptr[i].ndip*prseg_ptr[i].dwid;

      prseg_ptr[i].dtop = dmin - 0.5*prseg_ptr[i].dwid*sin(rperd*prseg_ptr[i].dip);
      if(prseg_ptr[i].dtop < 0.0)
         prseg_ptr[i].dtop = 0.0;

      avglon = dbllon/(float)(prseg_ptr[i].nstk);
      avglat = dbllat/(float)(prseg_ptr[i].nstk);
      se = -prseg_ptr[i].dwid*sin(rperd*prseg_ptr[i].dip)*cos(rperd*prseg_ptr[i].stk);
      sn = prseg_ptr[i].dwid*sin(rperd*prseg_ptr[i].dip)*sin(rperd*prseg_ptr[i].stk);

      set_ll(&avglon,&avglat,&prseg_ptr[i].elon,&prseg_ptr[i].elat,&sn,&se);

      prseg_ptr[i].shyp = -999.9;
      prseg_ptr[i].dhyp = -999.9;

      fprintf(stderr,"%3d: %11.5f %11.5f %2d %2d %10.4f %10.4f %.1f %.1f %10.4f\n",i,prseg_ptr[i].elon,prseg_ptr[i].elat,prseg_ptr[i].nstk,prseg_ptr[i].ndip,prseg_ptr[i].flen,prseg_ptr[i].fwid,prseg_ptr[i].stk,prseg_ptr[i].dip,prseg_ptr[i].dtop);
      }
   }
else
   {
   srf[0].srf_prect.nseg = 1;
   srf[0].srf_prect.prectseg = (struct srf_prectsegments *)check_malloc(srf[0].srf_prect.nseg*sizeof(struct srf_prectsegments));
   prseg_ptr = srf[0].srf_prect.prectseg;

   prseg_ptr[0].elon = *elon;
   prseg_ptr[0].elat = *elat;
   prseg_ptr[0].nstk = nx;
   prseg_ptr[0].ndip = ny;
   prseg_ptr[0].flen = *fl;
   prseg_ptr[0].fwid = *fw;
   prseg_ptr[0].dlen = *dx;
   prseg_ptr[0].dwid = *dy;
   prseg_ptr[0].stk = *stk;
   prseg_ptr[0].dip = *dip;
   prseg_ptr[0].dtop = *dtop;
   prseg_ptr[0].shyp = *sh;
   prseg_ptr[0].dhyp = *dh;
   }

srf[0].srf_apnts.np = 0;
for(ig=0;ig<srf[0].srf_prect.nseg;ig++)
   srf[0].srf_apnts.np = srf[0].srf_apnts.np + (prseg_ptr[ig].nstk)*(prseg_ptr[ig].ndip);

srf[0].srf_apnts.apntvals = (struct srf_apointvalues *)check_malloc((srf[0].srf_apnts.np)*sizeof(struct srf_apointvalues));
}

void init_plane_srf_pre20150316(struct standrupformat *srf,struct generic_slip *gslip,float *elon,float *elat,int nx,int ny,float *fl,float *fw,float *dx,float *dy,float *stk,float *dip,float *dtop,float *sh,float *dh)
{
struct srf_prectsegments *prseg_ptr;
float avglon, avglat, dmin, dold, se, sn;
int ig, iseg, nseg, nc, i, ip;
struct slippars *spar, *spar2;

float rperd = 0.017453293;

sprintf(srf[0].type,"PLANE");

if(gslip->np > 0)
   {
   spar = gslip->spar;
   spar2 = (struct slippars *)check_malloc(gslip->np*sizeof(struct slippars));

   nseg = 0;
   for(i=0;i<gslip->np;i++)
      {
      if(spar[i].segno > nseg)
         nseg = spar[i].segno;
      }
   nseg++;

   srf[0].srf_prect.nseg = nseg;
   srf[0].srf_prect.prectseg = (struct srf_prectsegments *)check_malloc(srf[0].srf_prect.nseg*sizeof(struct srf_prectsegments));

   prseg_ptr = srf[0].srf_prect.prectseg;

   nc = 0; 
   for(i=0;i<nseg;i++)
      {
      avglon = 0.0; 
      avglat = 0.0; 
      prseg_ptr[i].stk = 0.0;
      prseg_ptr[i].dip = 0.0;
      prseg_ptr[i].dlen = 0.0;
      prseg_ptr[i].dwid = 0.0;
      iseg = 0;
      dold = -1;
      prseg_ptr[i].nstk = 0;
      prseg_ptr[i].ndip = 0;
      for(ip=0;ip<gslip->np;ip++)
         {
         if(spar[ip].segno == i)
            {
            spar2[nc].lon = spar[ip].lon;
            spar2[nc].lat = spar[ip].lat;
            spar2[nc].dep = spar[ip].dep;
            spar2[nc].ds = spar[ip].ds;
            spar2[nc].dw = spar[ip].dw;
            spar2[nc].stk = spar[ip].stk;
            spar2[nc].dip = spar[ip].dip;
            spar2[nc].rake = spar[ip].rake;
            spar2[nc].slip = spar[ip].slip;
            spar2[nc].tinit = spar[ip].tinit;
            spar2[nc].segno = spar[ip].segno;

            if(iseg == 0)
               dmin = spar[ip].dep;

            if(spar[ip].dep == dmin)
               {
               avglon = avglon + spar[ip].lon;
               avglat = avglat + spar[ip].lat;
               prseg_ptr[i].nstk++;
               }

            if(spar[ip].dep != dold)
               {
               dold = spar[ip].dep;
               prseg_ptr[i].ndip++;
               }

            prseg_ptr[i].dlen = prseg_ptr[i].dlen + spar[ip].ds;
            prseg_ptr[i].dwid = prseg_ptr[i].dwid + spar[ip].dw;

            nc++;
            iseg++;

            if(iseg == 1)
               {
               prseg_ptr[i].stk = spar[ip].stk;
               prseg_ptr[i].dip = spar[ip].dip;
               }
            else
               {
               if(spar[ip].stk > prseg_ptr[i].stk/(float)(iseg-1.0) + 90.0)
                  {
                  prseg_ptr[i].stk = prseg_ptr[i].stk + spar[ip].stk - 180.0;
                  prseg_ptr[i].dip = prseg_ptr[i].dip + 180.0 - spar[ip].dip;
                  }
               else if(spar[ip].stk < prseg_ptr[i].stk/(float)(iseg-1.0) - 90.0)
                  {
                  prseg_ptr[i].stk = prseg_ptr[i].stk + spar[ip].stk + 180.0;
                  prseg_ptr[i].dip = prseg_ptr[i].dip + 180.0 - spar[ip].dip;
                  }
               else
                  {
                  prseg_ptr[i].stk = prseg_ptr[i].stk + spar[ip].stk;
                  prseg_ptr[i].dip = prseg_ptr[i].dip + spar[ip].dip;
                  }
               }
            }
         }

      prseg_ptr[i].stk = prseg_ptr[i].stk/(float)(iseg);
      prseg_ptr[i].dip = prseg_ptr[i].dip/(float)(iseg);

      if(prseg_ptr[i].dip > 90.0)
         {
         prseg_ptr[i].dip = 180.0 - prseg_ptr[i].dip;
         prseg_ptr[i].stk = prseg_ptr[i].stk - 180.0;
         }

      while(prseg_ptr[i].stk < 0.0)
         prseg_ptr[i].stk = prseg_ptr[i].stk + 360.0;
      while(prseg_ptr[i].stk >= 360.0)
         prseg_ptr[i].stk = prseg_ptr[i].stk - 360.0;

      prseg_ptr[i].dlen = prseg_ptr[i].dlen/(float)(iseg);
      prseg_ptr[i].dwid = prseg_ptr[i].dwid/(float)(iseg);

      prseg_ptr[i].flen = prseg_ptr[i].nstk*prseg_ptr[i].dlen;
      prseg_ptr[i].fwid = prseg_ptr[i].ndip*prseg_ptr[i].dwid;

      prseg_ptr[i].dtop = dmin - 0.5*prseg_ptr[i].dwid*sin(rperd*prseg_ptr[i].dip);
      if(prseg_ptr[i].dtop < 0.0)
         prseg_ptr[i].dtop = 0.0;

      avglon = avglon/(float)(prseg_ptr[i].nstk);
      avglat = avglat/(float)(prseg_ptr[i].nstk);
      se = -prseg_ptr[i].dwid*sin(rperd*prseg_ptr[i].dip)*cos(rperd*prseg_ptr[i].stk);
      sn = prseg_ptr[i].dwid*sin(rperd*prseg_ptr[i].dip)*sin(rperd*prseg_ptr[i].stk);

      set_ll(&avglon,&avglat,&prseg_ptr[i].elon,&prseg_ptr[i].elat,&sn,&se);

      prseg_ptr[i].shyp = -999.9;
      prseg_ptr[i].dhyp = -999.9;

      fprintf(stderr,"%3d: %11.5f %11.5f %2d %2d %10.4f %10.4f %.1f %.1f %10.4f\n",i,prseg_ptr[i].elon,prseg_ptr[i].elat,prseg_ptr[i].nstk,prseg_ptr[i].ndip,prseg_ptr[i].flen,prseg_ptr[i].fwid,prseg_ptr[i].stk,prseg_ptr[i].dip,prseg_ptr[i].dtop);
      }

   free(spar2);
   }
else
   {
   srf[0].srf_prect.nseg = 1;
   srf[0].srf_prect.prectseg = (struct srf_prectsegments *)check_malloc(srf[0].srf_prect.nseg*sizeof(struct srf_prectsegments));
   prseg_ptr = srf[0].srf_prect.prectseg;

   prseg_ptr[0].elon = *elon;
   prseg_ptr[0].elat = *elat;
   prseg_ptr[0].nstk = nx;
   prseg_ptr[0].ndip = ny;
   prseg_ptr[0].flen = *fl;
   prseg_ptr[0].fwid = *fw;
   prseg_ptr[0].dlen = *dx;
   prseg_ptr[0].dwid = *dy;
   prseg_ptr[0].stk = *stk;
   prseg_ptr[0].dip = *dip;
   prseg_ptr[0].dtop = *dtop;
   prseg_ptr[0].shyp = *sh;
   prseg_ptr[0].dhyp = *dh;
   }

srf[0].srf_apnts.np = 0;
for(ig=0;ig<srf[0].srf_prect.nseg;ig++)
   srf[0].srf_apnts.np = srf[0].srf_apnts.np + (prseg_ptr[ig].nstk)*(prseg_ptr[ig].ndip);

srf[0].srf_apnts.apntvals = (struct srf_apointvalues *)check_malloc((srf[0].srf_apnts.np)*sizeof(struct srf_apointvalues));
}

void load_slip_srf(struct standrupformat *srf,struct stfpar *spar,struct pointsource *ps)
{
struct srf_allpoints *apnts_ptr;
struct srf_apointvalues *apval_ptr;
struct srf_prectsegments *prseg_ptr;
float area;
float *stf;
int i, j, ip, ip0, iseg, nseg, ioff, noff, ntot;

float dmin = 4.0;
float dmax = 6.0;
float rtfac, tzero;
float rtfac0, sabs;
float rtfac_min = 0.05;

dmin = spar->risetimedep - spar->risetimedep_range;
dmax = spar->risetimedep + spar->risetimedep_range;
rtfac0 = spar->risetimefac - 1.0;

apnts_ptr = &(srf->srf_apnts);
apval_ptr = apnts_ptr->apntvals;

nseg = srf[0].srf_prect.nseg;
prseg_ptr = srf[0].srf_prect.prectseg;

ntot = 0;
for(iseg=0;iseg<nseg;iseg++)
   ntot = ntot + prseg_ptr[iseg].nstk;

ioff = 0;
noff = 0;
for(iseg=0;iseg<nseg;iseg++)
   {
   for(j=0;j<prseg_ptr[iseg].ndip;j++)
      {
      for(i=0;i<prseg_ptr[iseg].nstk;i++)
         {
	 ip = noff + i + j*prseg_ptr[iseg].nstk;
	 ip0 = i + ioff + j*ntot;

         apval_ptr[ip].stf1 = (float *)check_malloc(spar->nt*sizeof(float));
         stf = apval_ptr[ip].stf1;

         apval_ptr[ip].dt = spar->dt;

	 sabs = sqrt(ps[ip0].slip*ps[ip0].slip);
         if(sabs > MINSLIP)
	    {
            if(ps[ip0].dep >= dmax)
               rtfac = 1.0;
            else if(ps[ip0].dep < dmax && ps[ip0].dep > dmin)
               rtfac = 1.0 + rtfac0*(dmax-(ps[ip0].dep))/(dmax-dmin);
            else
               rtfac = 1.0 + rtfac0;

            if(spar->rt_scalefac > 0.0)
               {
               rtfac = rtfac*exp(0.5*log(sabs))/spar->rt_scalefac;

               if(rtfac < rtfac_min)
                  rtfac = rtfac_min;
               }

            if(strcmp(spar->stype,"brune") == 0)
               {
               tzero = 0.1*exp(-1.0)*sqrt(ps[ip0].slip)/(1.2); /* assume slip in cm */     
               tzero = tzero*rtfac;

               apval_ptr[ip].nt1 = gen_brune_stf(&(ps[ip0].slip),&tzero,stf,spar->nt,&spar->dt);
               }
	    else if(strcmp(spar->stype,"urs") == 0)
               {
	       tzero = rtfac*spar->trise;

               apval_ptr[ip].nt1 = gen_2tri_stf(&(ps[ip0].slip),&tzero,stf,spar->nt,&spar->dt,&ps[ip0].dep);
               }
            else if(strcmp(spar->stype,"ucsb") == 0)
               {
               tzero = rtfac*spar->trise;

               apval_ptr[ip].nt1 = gen_ucsb_stf(&(ps[ip0].slip),&tzero,stf,spar->nt,&spar->dt);
               }
            else if(strcmp(spar->stype,"tri") == 0)
               {
               tzero = rtfac*spar->trise;

               apval_ptr[ip].nt1 = gen_tri_stf(&(ps[ip0].slip),&tzero,stf,spar->nt,&spar->dt);
               }
	    }
         else
            apval_ptr[ip].nt1 = 0;

         if(apval_ptr[ip].nt1)
            apval_ptr[ip].stf1 = (float *)check_realloc(apval_ptr[ip].stf1,(apval_ptr[ip].nt1)*sizeof(float));
         else
            {
            free(apval_ptr[ip].stf1);
            apval_ptr[ip].stf1 = NULL;
            }

         apval_ptr[ip].lon = ps[ip0].lon;
         apval_ptr[ip].lat = ps[ip0].lat;
         apval_ptr[ip].dep = ps[ip0].dep;
         apval_ptr[ip].stk = ps[ip0].stk;
         apval_ptr[ip].dip = ps[ip0].dip;
         apval_ptr[ip].area = ps[ip0].area;
         apval_ptr[ip].rake = ps[ip0].rak;
         apval_ptr[ip].slip1 = ps[ip0].slip;

         apval_ptr[ip].slip2 = 0.0;
         apval_ptr[ip].nt2 = 0;
         apval_ptr[ip].stf2 = NULL;
         apval_ptr[ip].slip3 = 0.0;
         apval_ptr[ip].nt3 = 0;
         apval_ptr[ip].stf3 = NULL;
         }
      }
   ioff = ioff + prseg_ptr[iseg].nstk;
   noff = noff + prseg_ptr[iseg].nstk*prseg_ptr[iseg].ndip;
   }
}

void load_rupt_srf(struct standrupformat *srf,struct pointsource *ps,float *sh,float *dh)
{
struct srf_allpoints *apnts_ptr;
struct srf_apointvalues *apval_ptr;
struct srf_prectsegments *prseg_ptr;
int i, j, ip, ip0, iseg, nseg, ioff, noff, ntot;

(srf->srf_prect).prectseg[0].shyp = *sh;
(srf->srf_prect).prectseg[0].dhyp = *dh;

apnts_ptr = &(srf->srf_apnts);
apval_ptr = apnts_ptr->apntvals;

nseg = srf[0].srf_prect.nseg;
prseg_ptr = srf[0].srf_prect.prectseg;

ntot = 0;
for(iseg=0;iseg<nseg;iseg++)
   ntot = ntot + prseg_ptr[iseg].nstk;

ioff = 0;
noff = 0;
for(iseg=0;iseg<nseg;iseg++)
   {
   for(j=0;j<prseg_ptr[iseg].ndip;j++)
      {
      for(i=0;i<prseg_ptr[iseg].nstk;i++)
         {
	 ip = noff + i + j*prseg_ptr[iseg].nstk;
	 ip0 = i + ioff + j*ntot;

         apval_ptr[ip].tinit = ps[ip0].rupt;
         }
      }
   ioff = ioff + prseg_ptr[iseg].nstk;
   noff = noff + prseg_ptr[iseg].nstk*prseg_ptr[iseg].ndip;
   }
}

void write2gsf(struct generic_slip *gslip,struct pointsource *ps,char *ifile,char *ofile)
{
FILE *fpr, *fpw, *fopfile();
int ip;
char str[1024];

if(strcmp(ofile,"stdout") == 0)
   fpw = stdout;
else
   fpw = fopfile(ofile,"w");

fpr = fopfile(ifile,"r");
fgets(str,1024,fpr);
while(strncmp(str,"#",1) == 0)
   {
   fprintf(fpw,"%s",str);
   fgets(str,1024,fpr);
   } 
fclose(fpr);
   
fprintf(fpw,"%d\n",gslip->np);
      
for(ip=0;ip<gslip->np;ip++)
   {
   fprintf(fpw,"%11.5f %11.5f %8.4f %8.4f %8.4f %6.1f %6.1f %6.1f %8.2f %8.3f %3d\n",
                                                               gslip->spar[ip].lon,
                                                               gslip->spar[ip].lat,
                                                               gslip->spar[ip].dep,
                                                               gslip->spar[ip].ds,
                                                               gslip->spar[ip].dw,
                                                               gslip->spar[ip].stk,
                                                               gslip->spar[ip].dip,
                                                               ps[ip].rak,
                                                               ps[ip].slip,
                                                               ps[ip].rupt,
                                                               gslip->spar[ip].segno);
   }
fclose(fpw);
}

void write2srf(struct standrupformat *srf,char *str,int outbin)
{
write_srf(srf,str,outbin);
}

void load_slip_srf_dd2(struct standrupformat *srf,struct stfpar2 *spar,struct pointsource *ps,long *seed,struct velmodel *vmod)
{
struct srf_allpoints *apnts_ptr;
struct srf_apointvalues *apval_ptr;
struct srf_prectsegments *prseg_ptr;
float area;
float *stf;
int i, j, ip, ip0, iseg, nseg, ioff, noff, ntot, k;

float dmin1, dmax1;
float dmin2, dmax2;
float rtfac, tzero;
float rtfac1, rtfac2, sabs, slip_max, slip_avg, slip_min;
float rtfac_min = 0.05;

float fzero = 0.0;

dmin1 = spar->risetimedep - spar->risetimedep_range;
dmax1 = spar->risetimedep + spar->risetimedep_range;
rtfac1 = spar->risetimefac - 1.0;

dmin2 = spar->deep_risetimedep - spar->deep_risetimedep_range;
dmax2 = spar->deep_risetimedep + spar->deep_risetimedep_range;
rtfac2 = spar->deep_risetimefac - 1.0;

apnts_ptr = &(srf->srf_apnts);
apval_ptr = apnts_ptr->apntvals;

nseg = srf[0].srf_prect.nseg;
prseg_ptr = srf[0].srf_prect.prectseg;

/* version 2.0 stuff, just in case */

srf[0].srf_hcmnt.nline = 0;
srf[0].nseg = 1;
srf[0].np_seg = (int *)check_malloc((srf[0].nseg)*sizeof(int));
srf[0].np_seg[0] = srf[0].srf_apnts.np;

/* end version 2.0 stuff, just in case */

ntot = 0;
for(iseg=0;iseg<nseg;iseg++)
   ntot = ntot + prseg_ptr[iseg].nstk;

slip_max = -1.0e+15;
slip_min = 1.0e+15;
slip_avg = 0.0;
for(ip0=0;ip0<ntot*prseg_ptr[0].ndip;ip0++)
   {
   sabs = sqrt(ps[ip0].slip*ps[ip0].slip);
   slip_avg = slip_avg + sabs;
   if(sabs > slip_max)
      slip_max = sabs;
   if(sabs < slip_min)
      slip_min = sabs;
   }
slip_min = 1.001*slip_min;
slip_avg = slip_avg/(ntot*prseg_ptr[0].ndip);
fprintf(stderr,"SLIP_AVG= %f\n",slip_avg);

ioff = 0;
noff = 0;
for(iseg=0;iseg<nseg;iseg++)
   {
   for(j=0;j<prseg_ptr[iseg].ndip;j++)
      {
      for(i=0;i<prseg_ptr[iseg].nstk;i++)
         {
	 ip = noff + i + j*prseg_ptr[iseg].nstk;
	 ip0 = i + ioff + j*ntot;

         apval_ptr[ip].stf1 = (float *)check_malloc(spar->nt*sizeof(float));
         stf = apval_ptr[ip].stf1;

         apval_ptr[ip].dt = spar->dt;

	 sabs = sqrt(ps[ip0].slip*ps[ip0].slip);
         if(sabs > MINSLIP)
	    {
            if(ps[ip0].dep <= dmin1)
               rtfac = 1.0 + rtfac1;
            else if(ps[ip0].dep < dmax1 && ps[ip0].dep > dmin1)
               rtfac = 1.0 + rtfac1*(dmax1-(ps[ip0].dep))/(dmax1-dmin1);
            else if(ps[ip0].dep >= dmax1 && ps[ip0].dep <= dmin2)
               rtfac = 1.0;
            else if(ps[ip0].dep < dmax2 && ps[ip0].dep > dmin2)
               rtfac = 1.0 + rtfac2*((ps[ip0].dep)-dmin2)/(dmax2-dmin2);
            else
               rtfac = 1.0 + rtfac2;

            if(spar->rt_scalefac > 0.0)
               {
	       if(sabs <= 1.0*slip_min)
                  rtfac = (1.0 + rtfac1)*exp(0.5*log(slip_max))/spar->rt_scalefac;
	       else
                  rtfac = rtfac*exp(0.5*log(sabs))/spar->rt_scalefac;

/* 20131119 RWG:

    Add random perturbations to rise time so that it is not 1:1 correlated with sqrt(slip).
    Perturbations are log normal with ln(sigma)=rt_rand*trise.  Default for rt_rand=0.5.

*/

               if(spar->rt_rand > 0.0)
                  rtfac = rtfac*exp(gaus_rand(&(spar->rt_rand),&fzero,seed));

               if(rtfac < rtfac_min)
                  rtfac = rtfac_min;
               }

            if(strcmp(spar->stype,"brune") == 0)
               {
               tzero = 0.1*exp(-1.0)*sqrt(ps[ip0].slip)/(1.2); /* assume slip in cm */     
               tzero = tzero*rtfac;

               apval_ptr[ip].nt1 = gen_brune_stf(&(ps[ip0].slip),&tzero,stf,spar->nt,&spar->dt);
               }
	    else if(strcmp(spar->stype,"urs") == 0)
               {
	       tzero = rtfac*spar->trise;

               apval_ptr[ip].nt1 = gen_2tri_stf(&(ps[ip0].slip),&tzero,stf,spar->nt,&spar->dt,&ps[ip0].dep);
               }
            else if(strcmp(spar->stype,"ucsb") == 0)
               {
               tzero = rtfac*spar->trise;

               apval_ptr[ip].nt1 = gen_ucsb_stf(&(ps[ip0].slip),&tzero,stf,spar->nt,&spar->dt);
               }
            else if(strcmp(spar->stype,"Mliu") == 0)
               {
               tzero = rtfac*spar->trise;

               apval_ptr[ip].nt1 = gen_Mliu_stf(&(ps[ip0].slip),&tzero,stf,spar->nt,&spar->dt);
               }
            else if(strcmp(spar->stype,"tri") == 0)
               {
               tzero = rtfac*spar->trise;

               apval_ptr[ip].nt1 = gen_tri_stf(&(ps[ip0].slip),&tzero,stf,spar->nt,&spar->dt);
               }
	    }
         else
            apval_ptr[ip].nt1 = 0;

         if(apval_ptr[ip].nt1)
            apval_ptr[ip].stf1 = (float *)check_realloc(apval_ptr[ip].stf1,(apval_ptr[ip].nt1)*sizeof(float));
         else
            {
            free(apval_ptr[ip].stf1);
            apval_ptr[ip].stf1 = NULL;
            }

         k = 0;
         while(ps[ip0].dep > vmod->dep[k] && k != vmod->nlay-1)
            k++;

         apval_ptr[ip].lon = ps[ip0].lon;
         apval_ptr[ip].lat = ps[ip0].lat;
         apval_ptr[ip].dep = ps[ip0].dep;
         apval_ptr[ip].stk = ps[ip0].stk;
         apval_ptr[ip].dip = ps[ip0].dip;
         apval_ptr[ip].area = ps[ip0].area;

         apval_ptr[ip].vs = 1.0e+05*vmod->vs[k];
         apval_ptr[ip].den = vmod->den[k];

         apval_ptr[ip].rake = ps[ip0].rak;
         apval_ptr[ip].slip1 = ps[ip0].slip;

         apval_ptr[ip].slip2 = 0.0;
         apval_ptr[ip].nt2 = 0;
         apval_ptr[ip].stf2 = NULL;
         apval_ptr[ip].slip3 = 0.0;
         apval_ptr[ip].nt3 = 0;
         apval_ptr[ip].stf3 = NULL;
         }
      }
   ioff = ioff + prseg_ptr[iseg].nstk;
   noff = noff + prseg_ptr[iseg].nstk*prseg_ptr[iseg].ndip;
   }
}

void load_vsden_srf(struct standrupformat *srf,struct velmodel *vmod)
{
struct srf_allpoints *apnts_ptr;
struct srf_apointvalues *apval_ptr;
struct srf_prectsegments *prseg_ptr;
int i, j, ip, k, iseg, nseg, ioff, noff, ntot;

apnts_ptr = &(srf->srf_apnts);
apval_ptr = apnts_ptr->apntvals;

nseg = srf[0].srf_prect.nseg;
prseg_ptr = srf[0].srf_prect.prectseg;

ntot = 0;
for(iseg=0;iseg<nseg;iseg++)
   ntot = ntot + prseg_ptr[iseg].nstk;

ioff = 0;
noff = 0;
for(iseg=0;iseg<nseg;iseg++)
   {
   for(j=0;j<prseg_ptr[iseg].ndip;j++)
      {
      for(i=0;i<prseg_ptr[iseg].nstk;i++)
         {
	 ip = noff + i + j*prseg_ptr[iseg].nstk;

         k = 0;
         while(apval_ptr[ip].dep > vmod->dep[k] && k != vmod->nlay-1)
            k++;

         apval_ptr[ip].vs = vmod->vs[k];
         apval_ptr[ip].den = vmod->den[k];
         }
      }
   ioff = ioff + prseg_ptr[iseg].nstk;
   noff = noff + prseg_ptr[iseg].nstk*prseg_ptr[iseg].ndip;
   }
}

void load_slip_srf_dd3(struct standrupformat *srf,struct stfpar2 *spar,struct pointsource *ps,float *rtime_r,struct velmodel *vmod)
{
struct srf_allpoints *apnts_ptr;
struct srf_apointvalues *apval_ptr;
struct srf_prectsegments *prseg_ptr;
float area;
float *stf;
int i, j, ip, ip0, iseg, nseg, ioff, noff, ntot, k;

float dmin1, dmax1;
float dmin2, dmax2;
float rtfac, tzero;
float rtfac1, rtfac2, sabs, slip_max, slip_avg, slip_min;
float rtfac_min = 0.05;

rtfac_min = spar->dt/spar->trise;

dmin1 = spar->risetimedep - spar->risetimedep_range;
dmax1 = spar->risetimedep + spar->risetimedep_range;
rtfac1 = spar->risetimefac - 1.0;

dmin2 = spar->deep_risetimedep - spar->deep_risetimedep_range;
dmax2 = spar->deep_risetimedep + spar->deep_risetimedep_range;
rtfac2 = spar->deep_risetimefac - 1.0;

apnts_ptr = &(srf->srf_apnts);
apval_ptr = apnts_ptr->apntvals;

nseg = srf[0].srf_prect.nseg;
prseg_ptr = srf[0].srf_prect.prectseg;

/* version 2.0 stuff, just in case */

srf[0].srf_hcmnt.nline = 0;
srf[0].nseg = 1;
srf[0].np_seg = (int *)check_malloc((srf[0].nseg)*sizeof(int));
srf[0].np_seg[0] = srf[0].srf_apnts.np;

/* end version 2.0 stuff, just in case */

ntot = 0;
for(iseg=0;iseg<nseg;iseg++)
   ntot = ntot + prseg_ptr[iseg].nstk;

slip_max = -1.0e+15;
slip_min = 1.0e+15;
slip_avg = 0.0;
for(ip0=0;ip0<ntot*prseg_ptr[0].ndip;ip0++)
   {
   sabs = sqrt(ps[ip0].slip*ps[ip0].slip);
   slip_avg = slip_avg + sabs;
   if(sabs > slip_max)
      slip_max = sabs;
   if(sabs < slip_min)
      slip_min = sabs;
   }
slip_min = 1.001*slip_min;
slip_avg = slip_avg/(ntot*prseg_ptr[0].ndip);
fprintf(stderr,"SLIP_AVG= %f\n",slip_avg);

ioff = 0;
noff = 0;
for(iseg=0;iseg<nseg;iseg++)
   {
   for(j=0;j<prseg_ptr[iseg].ndip;j++)
      {
      for(i=0;i<prseg_ptr[iseg].nstk;i++)
         {
	 ip = noff + i + j*prseg_ptr[iseg].nstk;
	 ip0 = i + ioff + j*ntot;

         apval_ptr[ip].stf1 = (float *)check_malloc(spar->nt*sizeof(float));
         stf = apval_ptr[ip].stf1;

         apval_ptr[ip].dt = spar->dt;

	 sabs = sqrt(ps[ip0].slip*ps[ip0].slip);
         if(sabs > MINSLIP)
	    {
            if(ps[ip0].dep <= dmin1)
               rtfac = 1.0 + rtfac1;
            else if(ps[ip0].dep < dmax1 && ps[ip0].dep > dmin1)
               rtfac = 1.0 + rtfac1*(dmax1-(ps[ip0].dep))/(dmax1-dmin1);
            else if(ps[ip0].dep >= dmax1 && ps[ip0].dep <= dmin2)
               rtfac = 1.0;
            else if(ps[ip0].dep < dmax2 && ps[ip0].dep > dmin2)
               rtfac = 1.0 + rtfac2*((ps[ip0].dep)-dmin2)/(dmax2-dmin2);
            else
               rtfac = 1.0 + rtfac2;

/* 20150324 RWG:
    New way using array of gaussian values roughly correlated with slip
    Note that rtime_r is already normalized to have average of 1
*/
            if(spar->rt_scalefac > 0.0)
               rtfac = rtfac*sqrt(rtime_r[ip0])/spar->rt_scalefac;

            if(rtfac < rtfac_min)
               rtfac = rtfac_min;

            if(strcmp(spar->stype,"brune") == 0)
               {
               tzero = 0.1*exp(-1.0)*sqrt(ps[ip0].slip)/(1.2); /* assume slip in cm */     
               tzero = tzero*rtfac;

               apval_ptr[ip].nt1 = gen_brune_stf(&(ps[ip0].slip),&tzero,stf,spar->nt,&spar->dt);
               }
	    else if(strcmp(spar->stype,"urs") == 0)
               {
	       tzero = rtfac*spar->trise;

               apval_ptr[ip].nt1 = gen_2tri_stf(&(ps[ip0].slip),&tzero,stf,spar->nt,&spar->dt,&ps[ip0].dep);
               }
            else if(strcmp(spar->stype,"ucsb") == 0)
               {
               tzero = rtfac*spar->trise;

               apval_ptr[ip].nt1 = gen_ucsb_stf(&(ps[ip0].slip),&tzero,stf,spar->nt,&spar->dt);
               }
            else if(strcmp(spar->stype,"Mliu") == 0)
               {
               tzero = rtfac*spar->trise;

               apval_ptr[ip].nt1 = gen_Mliu_stf(&(ps[ip0].slip),&tzero,stf,spar->nt,&spar->dt);
               }
            else if(strcmp(spar->stype,"tri") == 0)
               {
               tzero = rtfac*spar->trise;

               apval_ptr[ip].nt1 = gen_tri_stf(&(ps[ip0].slip),&tzero,stf,spar->nt,&spar->dt);
               }
	    }
         else
            apval_ptr[ip].nt1 = 0;

         if(apval_ptr[ip].nt1)
            apval_ptr[ip].stf1 = (float *)check_realloc(apval_ptr[ip].stf1,(apval_ptr[ip].nt1)*sizeof(float));
         else
            {
            free(apval_ptr[ip].stf1);
            apval_ptr[ip].stf1 = NULL;
            }

         k = 0;
         while(ps[ip0].dep > vmod->dep[k] && k != vmod->nlay-1)
            k++;

         apval_ptr[ip].lon = ps[ip0].lon;
         apval_ptr[ip].lat = ps[ip0].lat;
         apval_ptr[ip].dep = ps[ip0].dep;
         apval_ptr[ip].stk = ps[ip0].stk;
         apval_ptr[ip].dip = ps[ip0].dip;
         apval_ptr[ip].area = ps[ip0].area;

         apval_ptr[ip].vs = 1.0e+05*vmod->vs[k];
         apval_ptr[ip].den = vmod->den[k];

         apval_ptr[ip].rake = ps[ip0].rak;
         apval_ptr[ip].slip1 = ps[ip0].slip;

         apval_ptr[ip].slip2 = 0.0;
         apval_ptr[ip].nt2 = 0;
         apval_ptr[ip].stf2 = NULL;
         apval_ptr[ip].slip3 = 0.0;
         apval_ptr[ip].nt3 = 0;
         apval_ptr[ip].stf3 = NULL;
         }
      }
   ioff = ioff + prseg_ptr[iseg].nstk;
   noff = noff + prseg_ptr[iseg].nstk*prseg_ptr[iseg].ndip;
   }
}

void load_slip_srf_dd4(struct standrupformat *srf,struct stfpar2 *spar,struct pointsource *ps,float *rtime1_r,float *rtime2_r,struct velmodel *vmod)
{
struct srf_allpoints *apnts_ptr;
struct srf_apointvalues *apval_ptr;
struct srf_prectsegments *prseg_ptr;
float area;
float *stf;
int i, j, ip, ip0, iseg, nseg, ioff, noff, ntot, k;

float dmin1, dmax1;
float dmin2, dmax2;
float rtfac, tzero;
float rtfac1, rtfac2, sabs, slip_max, slip_avg, slip_min;
float rtfac_min = 0.05;

rtfac_min = spar->dt/spar->trise;

dmin1 = spar->risetimedep - spar->risetimedep_range;
dmax1 = spar->risetimedep + spar->risetimedep_range;
rtfac1 = spar->risetimefac - 1.0;

dmin2 = spar->deep_risetimedep - spar->deep_risetimedep_range;
dmax2 = spar->deep_risetimedep + spar->deep_risetimedep_range;
rtfac2 = spar->deep_risetimefac - 1.0;

apnts_ptr = &(srf->srf_apnts);
apval_ptr = apnts_ptr->apntvals;

nseg = srf[0].srf_prect.nseg;
prseg_ptr = srf[0].srf_prect.prectseg;

/* version 2.0 stuff, just in case */

srf[0].srf_hcmnt.nline = 0;
srf[0].nseg = 1;
srf[0].np_seg = (int *)check_malloc((srf[0].nseg)*sizeof(int));
srf[0].np_seg[0] = srf[0].srf_apnts.np;

/* end version 2.0 stuff, just in case */

ntot = 0;
for(iseg=0;iseg<nseg;iseg++)
   ntot = ntot + prseg_ptr[iseg].nstk;

slip_max = -1.0e+15;
slip_min = 1.0e+15;
slip_avg = 0.0;
for(ip0=0;ip0<ntot*prseg_ptr[0].ndip;ip0++)
   {
   sabs = sqrt(ps[ip0].slip*ps[ip0].slip);
   slip_avg = slip_avg + sabs;
   if(sabs > slip_max)
      slip_max = sabs;
   if(sabs < slip_min)
      slip_min = sabs;
   }
slip_min = 1.001*slip_min;
slip_avg = slip_avg/(ntot*prseg_ptr[0].ndip);
fprintf(stderr,"SLIP_AVG= %f\n",slip_avg);

ioff = 0;
noff = 0;
for(iseg=0;iseg<nseg;iseg++)
   {
   for(j=0;j<prseg_ptr[iseg].ndip;j++)
      {
      for(i=0;i<prseg_ptr[iseg].nstk;i++)
         {
	 ip = noff + i + j*prseg_ptr[iseg].nstk;
	 ip0 = i + ioff + j*ntot;

         apval_ptr[ip].stf1 = (float *)check_malloc(spar->nt*sizeof(float));
         stf = apval_ptr[ip].stf1;

         apval_ptr[ip].dt = spar->dt;

	 sabs = sqrt(ps[ip0].slip*ps[ip0].slip);
         if(sabs > MINSLIP)
	    {
            if(ps[ip0].dep <= dmin1)
               rtfac = 1.0 + rtfac1;
            else if(ps[ip0].dep < dmax1 && ps[ip0].dep > dmin1)
               rtfac = 1.0 + rtfac1*(dmax1-(ps[ip0].dep))/(dmax1-dmin1);
            else if(ps[ip0].dep >= dmax1 && ps[ip0].dep <= dmin2)
               rtfac = 1.0;
            else if(ps[ip0].dep < dmax2 && ps[ip0].dep > dmin2)
               rtfac = 1.0 + rtfac2*((ps[ip0].dep)-dmin2)/(dmax2-dmin2);
            else
               rtfac = 1.0 + rtfac2;

/* 20150324 RWG:
    New way using array of gaussian values roughly correlated with slip
    Note that rtime1_r is already normalized to have average of 1
*/
            if(spar->rt_scalefac > 0.0)
               rtfac = rtfac*sqrt(rtime1_r[ip0])/spar->rt_scalefac;

	    /*
            rtfac = rtfac*exp(rtime2_r[ip0]);
	    */
            rtfac = rtfac*(1.0+rtime2_r[ip0]);

            if(rtfac < rtfac_min)
               rtfac = rtfac_min;

            if(strcmp(spar->stype,"brune") == 0)
               {
               tzero = 0.1*exp(-1.0)*sqrt(ps[ip0].slip)/(1.2); /* assume slip in cm */     
               tzero = tzero*rtfac;

               apval_ptr[ip].nt1 = gen_brune_stf(&(ps[ip0].slip),&tzero,stf,spar->nt,&spar->dt);
               }
	    else if(strcmp(spar->stype,"urs") == 0)
               {
	       tzero = rtfac*spar->trise;

               apval_ptr[ip].nt1 = gen_2tri_stf(&(ps[ip0].slip),&tzero,stf,spar->nt,&spar->dt,&ps[ip0].dep);
               }
            else if(strcmp(spar->stype,"ucsb") == 0)
               {
               tzero = rtfac*spar->trise;

               apval_ptr[ip].nt1 = gen_ucsb_stf(&(ps[ip0].slip),&tzero,stf,spar->nt,&spar->dt);
               }
            else if(strcmp(spar->stype,"Mliu") == 0)
               {
               tzero = rtfac*spar->trise;

               apval_ptr[ip].nt1 = gen_Mliu_stf(&(ps[ip0].slip),&tzero,stf,spar->nt,&spar->dt);
               }
            else if(strcmp(spar->stype,"tri") == 0)
               {
               tzero = rtfac*spar->trise;

               apval_ptr[ip].nt1 = gen_tri_stf(&(ps[ip0].slip),&tzero,stf,spar->nt,&spar->dt);
               }
	    }
         else
            apval_ptr[ip].nt1 = 0;

         if(apval_ptr[ip].nt1)
            apval_ptr[ip].stf1 = (float *)check_realloc(apval_ptr[ip].stf1,(apval_ptr[ip].nt1)*sizeof(float));
         else
            {
            free(apval_ptr[ip].stf1);
            apval_ptr[ip].stf1 = NULL;
            }

         k = 0;
         while(ps[ip0].dep > vmod->dep[k] && k != vmod->nlay-1)
            k++;

         apval_ptr[ip].lon = ps[ip0].lon;
         apval_ptr[ip].lat = ps[ip0].lat;
         apval_ptr[ip].dep = ps[ip0].dep;
         apval_ptr[ip].stk = ps[ip0].stk;
         apval_ptr[ip].dip = ps[ip0].dip;
         apval_ptr[ip].area = ps[ip0].area;

         apval_ptr[ip].vs = 1.0e+05*vmod->vs[k];
         apval_ptr[ip].den = vmod->den[k];

         apval_ptr[ip].rake = ps[ip0].rak;
         apval_ptr[ip].slip1 = ps[ip0].slip;

         apval_ptr[ip].slip2 = 0.0;
         apval_ptr[ip].nt2 = 0;
         apval_ptr[ip].stf2 = NULL;
         apval_ptr[ip].slip3 = 0.0;
         apval_ptr[ip].nt3 = 0;
         apval_ptr[ip].stf3 = NULL;
         }
      }
   ioff = ioff + prseg_ptr[iseg].nstk;
   noff = noff + prseg_ptr[iseg].nstk*prseg_ptr[iseg].ndip;
   }
}
