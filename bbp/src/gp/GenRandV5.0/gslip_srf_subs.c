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

float flen_tot, flen_back, shyp_end;

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

   flen_tot = 0.0;
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
               dbldlen = 1.0*spar[ip].ds;
               dbldwid = 1.0*spar[ip].dw;
               dblstk = spar[ip].stk;
               dbldip = spar[ip].dip;
               }
            else
               {
               dbldlen = dbldlen + 1.0*spar[ip].ds;
               dbldwid = dbldwid + 1.0*spar[ip].dw;

               if(spar[ip].stk > dblstk/(float)(ip_seg-1.0) + 90.0)
                  {
                  dblstk = dblstk + 1.0*spar[ip].stk - 180.0;
                  dbldip = dbldip + 180.0 - 1.0*spar[ip].dip;
                  }
               else if(spar[ip].stk < dblstk/(float)(ip_seg-1.0) - 90.0)
                  {
                  dblstk = dblstk + 1.0*spar[ip].stk + 180.0;
                  dbldip = dbldip + 180.0 - 1.0*spar[ip].dip;
                  }
               else
                  {
                  dblstk = dblstk + 1.0*spar[ip].stk;
                  dbldip = dbldip + 1.0*spar[ip].dip;
                  }
               }
            }
         }

      prseg_ptr[i].dlen = dbldlen/(double)(1.0*ip_seg);
      prseg_ptr[i].dwid = dbldwid/(double)(1.0*ip_seg);

      prseg_ptr[i].stk = dblstk/(double)(1.0*ip_seg);
      prseg_ptr[i].dip = dbldip/(double)(1.0*ip_seg);

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

/* 2016-01-12: Discovered that se,sn computed immediately below are incorrect
   correct relations (replace sin with cos & multiply by 0.5) are given below these

      se = -prseg_ptr[i].dwid*sin(rperd*prseg_ptr[i].dip)*cos(rperd*prseg_ptr[i].stk);
      sn = prseg_ptr[i].dwid*sin(rperd*prseg_ptr[i].dip)*sin(rperd*prseg_ptr[i].stk);
*/
      se = -0.5*prseg_ptr[i].dwid*cos(rperd*prseg_ptr[i].dip)*cos(rperd*prseg_ptr[i].stk);
      sn = 0.5*prseg_ptr[i].dwid*cos(rperd*prseg_ptr[i].dip)*sin(rperd*prseg_ptr[i].stk);

      set_ll(&avglon,&avglat,&prseg_ptr[i].elon,&prseg_ptr[i].elat,&sn,&se);

      fprintf(stderr,"%3d: %11.5f %11.5f %2d %2d %10.4f %10.4f %.1f %.1f %10.4f\n",i,prseg_ptr[i].elon,prseg_ptr[i].elat,prseg_ptr[i].nstk,prseg_ptr[i].ndip,prseg_ptr[i].flen,prseg_ptr[i].fwid,prseg_ptr[i].stk,prseg_ptr[i].dip,prseg_ptr[i].dtop);

      flen_tot = flen_tot + prseg_ptr[i].flen;
      }

   shyp_end = *sh + 0.5*flen_tot;
   flen_back = 0.0;
fprintf(stderr,"%10.4f %10.4f %10.4f\n",*sh,flen_tot,shyp_end);
   for(i=0;i<nseg;i++)
      {
      prseg_ptr[i].shyp = shyp_end - flen_back - 0.5*prseg_ptr[i].flen;
      prseg_ptr[i].dhyp = *dh;

      flen_back = flen_back + prseg_ptr[i].flen;

      fprintf(stderr,"%3d: %10.4f %10.4f\n",i,prseg_ptr[i].shyp,flen_back);
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

/* 2025-02-20 Don't do this because init_plane_srf already does it */
/* load hypo here just to be sure it's in header
   (may have already been done in init_plane_srf if read_gsf=1)
(srf->srf_prect).prectseg[0].shyp = *sh;
(srf->srf_prect).prectseg[0].dhyp = *dh;
*/

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

void load_slip_srf_dd4_vsden(struct standrupformat *srf,struct stfpar2 *spar,struct pointsource *ps,float *rtime1_r,float *rtime2_r)
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

float beta, b_dmin, b_dmax, delta_b;

rtfac_min = 5.0*spar->dt/spar->trise;    /* so minimum rise time >= 5*dt */
rtfac_min = 10.0*spar->dt/spar->trise;    /* so minimum rise time >= 10*dt */
rtfac_min = spar->dt/spar->trise;    /* so minimum rise time >= dt */

dmin1 = spar->risetimedep - spar->risetimedep_range;
dmax1 = spar->risetimedep + spar->risetimedep_range;
rtfac1 = spar->risetimefac - 1.0;

dmin2 = spar->deep_risetimedep - spar->deep_risetimedep_range;
dmax2 = spar->deep_risetimedep + spar->deep_risetimedep_range;
rtfac2 = spar->deep_risetimefac - 1.0;

b_dmin = spar->beta_depth - spar->beta_depth_range;
b_dmax = spar->beta_depth + spar->beta_depth_range;
delta_b = spar->beta_deep - spar->beta_shal;

apnts_ptr = &(srf->srf_apnts);
apval_ptr = apnts_ptr->apntvals;

nseg = srf[0].srf_prect.nseg;
prseg_ptr = srf[0].srf_prect.prectseg;

/* version 2.0 stuff, just in case */

srf[0].srf_hcmnt.nline = 0;

/* 2018-04-02 updated for multiple segs */
srf[0].nseg = srf[0].srf_prect.nseg;
srf[0].np_seg = (int *)check_malloc((srf[0].nseg)*sizeof(int));

for(iseg=0;iseg<nseg;iseg++)
   srf[0].np_seg[iseg] = prseg_ptr[iseg].nstk*prseg_ptr[0].ndip;

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
            else if(strcmp(spar->stype,"OliuP") == 0) /* original form */
               {
               if(ps[ip0].dep <= b_dmin)
                  beta = spar->beta_shal;
               else if(ps[ip0].dep < b_dmax && ps[ip0].dep > b_dmin)
                  beta = spar->beta_shal + delta_b*(ps[ip0].dep - b_dmin)/(b_dmax-b_dmin);
               else
                  beta = spar->beta_deep;

               tzero = rtfac*spar->trise;

               apval_ptr[ip].nt1 = gen_OliuP_stf(&(ps[ip0].slip),&tzero,&beta,stf,spar->nt,&spar->dt);
               }
            else if(strcmp(spar->stype,"MliuP") == 0) /* modified form */
               {
               if(ps[ip0].dep <= b_dmin)
                  beta = spar->beta_shal;
               else if(ps[ip0].dep < b_dmax && ps[ip0].dep > b_dmin)
                  beta = spar->beta_shal + delta_b*(ps[ip0].dep - b_dmin)/(b_dmax-b_dmin);
               else
                  beta = spar->beta_deep;

               tzero = rtfac*spar->trise;

               apval_ptr[ip].nt1 = gen_MliuP_stf(&(ps[ip0].slip),&tzero,&beta,stf,spar->nt,&spar->dt);
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

         apval_ptr[ip].vs = 1.0e+05*ps[ip0].vs;
         apval_ptr[ip].den = ps[ip0].den;

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

float beta, b_dmin, b_dmax, delta_b;

rtfac_min = 5.0*spar->dt/spar->trise;    /* so minimum rise time >= 5*dt */
rtfac_min = 10.0*spar->dt/spar->trise;    /* so minimum rise time >= 10*dt */
rtfac_min = spar->dt/spar->trise;    /* so minimum rise time >= dt */

dmin1 = spar->risetimedep - spar->risetimedep_range;
dmax1 = spar->risetimedep + spar->risetimedep_range;
rtfac1 = spar->risetimefac - 1.0;

dmin2 = spar->deep_risetimedep - spar->deep_risetimedep_range;
dmax2 = spar->deep_risetimedep + spar->deep_risetimedep_range;
rtfac2 = spar->deep_risetimefac - 1.0;

b_dmin = spar->beta_depth - spar->beta_depth_range;
b_dmax = spar->beta_depth + spar->beta_depth_range;
delta_b = spar->beta_deep - spar->beta_shal;

apnts_ptr = &(srf->srf_apnts);
apval_ptr = apnts_ptr->apntvals;

nseg = srf[0].srf_prect.nseg;
prseg_ptr = srf[0].srf_prect.prectseg;

/* version 2.0 stuff, just in case */

srf[0].srf_hcmnt.nline = 0;

/* 2018-04-02 updated for multiple segs */
srf[0].nseg = srf[0].srf_prect.nseg;
srf[0].np_seg = (int *)check_malloc((srf[0].nseg)*sizeof(int));

for(iseg=0;iseg<nseg;iseg++)
   srf[0].np_seg[iseg] = prseg_ptr[iseg].nstk*prseg_ptr[0].ndip;

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
            else if(strcmp(spar->stype,"OliuP") == 0) /* original form */
               {
               if(ps[ip0].dep <= b_dmin)
                  beta = spar->beta_shal;
               else if(ps[ip0].dep < b_dmax && ps[ip0].dep > b_dmin)
                  beta = spar->beta_shal + delta_b*(ps[ip0].dep - b_dmin)/(b_dmax-b_dmin);
               else
                  beta = spar->beta_deep;

               tzero = rtfac*spar->trise;

               apval_ptr[ip].nt1 = gen_OliuP_stf(&(ps[ip0].slip),&tzero,&beta,stf,spar->nt,&spar->dt);
               }
            else if(strcmp(spar->stype,"MliuP") == 0) /* modified form */
               {
               if(ps[ip0].dep <= b_dmin)
                  beta = spar->beta_shal;
               else if(ps[ip0].dep < b_dmax && ps[ip0].dep > b_dmin)
                  beta = spar->beta_shal + delta_b*(ps[ip0].dep - b_dmin)/(b_dmax-b_dmin);
               else
                  beta = spar->beta_deep;

               tzero = rtfac*spar->trise;

               apval_ptr[ip].nt1 = gen_MliuP_stf(&(ps[ip0].slip),&tzero,&beta,stf,spar->nt,&spar->dt);
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

int get_seg_bounds(struct standrupformat *srf,float *gbnd)
{
struct srf_prectsegments *prseg_ptr;
int i, nb;

nb = srf[0].srf_prect.nseg - 1;
prseg_ptr = srf[0].srf_prect.prectseg;

gbnd[0] = prseg_ptr[0].flen;
for(i=1;i<nb;i++)
   gbnd[i] = gbnd[i-1] + prseg_ptr[i].flen;

return(nb);
}

void load_slip_srf_dd5_vsden(struct standrupformat *srf,struct stfpar2 *spar,struct pointsource *ps,float *rtime1_r,float *rtime2_r)
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

float beta, b_dmin, b_dmax, delta_b;

rtfac_min = 5.0*spar->dt/spar->trise;    /* so minimum rise time >= 5*dt */
rtfac_min = 50.0*spar->dt/spar->trise;    /* so minimum rise time >= 50*dt */
rtfac_min = spar->dt/spar->trise;    /* so minimum rise time >= dt */
rtfac_min = 3.0*spar->dt/spar->trise;    /* min. risetime = 3*dt; stf has 2 non-zero values */
rtfac_min = 5.0*spar->dt/spar->trise;    /* min. risetime = 5*dt; stf has 4 non-zero values */

fprintf(stderr,"rtfac_min= %5e\n",rtfac_min);

dmin1 = spar->risetimedep - spar->risetimedep_range;
dmax1 = spar->risetimedep + spar->risetimedep_range;
rtfac1 = spar->risetimefac - 1.0;

dmin2 = spar->deep_risetimedep - spar->deep_risetimedep_range;
dmax2 = spar->deep_risetimedep + spar->deep_risetimedep_range;
rtfac2 = spar->deep_risetimefac - 1.0;

b_dmin = spar->beta_depth - spar->beta_depth_range;
b_dmax = spar->beta_depth + spar->beta_depth_range;
delta_b = spar->beta_deep - spar->beta_shal;

apnts_ptr = &(srf->srf_apnts);
apval_ptr = apnts_ptr->apntvals;

nseg = srf[0].srf_prect.nseg;
prseg_ptr = srf[0].srf_prect.prectseg;

/* version 2.0 stuff, just in case */

srf[0].srf_hcmnt.nline = 0;

/* 2018-04-02 updated for multiple segs */
srf[0].nseg = srf[0].srf_prect.nseg;
srf[0].np_seg = (int *)check_malloc((srf[0].nseg)*sizeof(int));

for(iseg=0;iseg<nseg;iseg++)
   srf[0].np_seg[iseg] = prseg_ptr[iseg].nstk*prseg_ptr[0].ndip;

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

/* 20250520 RWG:
    Added scaling of risetime using ps array, see init_asperity_IO()
*/
            if(spar->rt_scalefac > 0.0 && rtime1_r[ip0] > 0.0)
               rtfac = ps[ip0].trise_fac*rtfac*rtime1_r[ip0]/spar->rt_scalefac;
            if(rtime1_r[ip0] <= 0.0)
               rtfac = 0.0;

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
            else if(strcmp(spar->stype,"OliuP2") == 0) /* original form, beta in ps array */
               {
               tzero = rtfac*spar->trise;

               apval_ptr[ip].nt1 = gen_OliuP_stf(&(ps[ip0].slip),&tzero,&(ps[ip0].beta),stf,spar->nt,&spar->dt);
               }
            else if(strcmp(spar->stype,"OliuP") == 0) /* original form */
               {
               if(ps[ip0].dep <= b_dmin)
                  beta = spar->beta_shal;
               else if(ps[ip0].dep < b_dmax && ps[ip0].dep > b_dmin)
                  beta = spar->beta_shal + delta_b*(ps[ip0].dep - b_dmin)/(b_dmax-b_dmin);
               else
                  beta = spar->beta_deep;

               tzero = rtfac*spar->trise;

               apval_ptr[ip].nt1 = gen_OliuP_stf(&(ps[ip0].slip),&tzero,&beta,stf,spar->nt,&spar->dt);
               }
            else if(strcmp(spar->stype,"MliuP") == 0) /* modified form */
               {
               if(ps[ip0].dep <= b_dmin)
                  beta = spar->beta_shal;
               else if(ps[ip0].dep < b_dmax && ps[ip0].dep > b_dmin)
                  beta = spar->beta_shal + delta_b*(ps[ip0].dep - b_dmin)/(b_dmax-b_dmin);
               else
                  beta = spar->beta_deep;

               tzero = rtfac*spar->trise;

               apval_ptr[ip].nt1 = gen_MliuP_stf(&(ps[ip0].slip),&tzero,&beta,stf,spar->nt,&spar->dt);
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

         apval_ptr[ip].vs = 1.0e+05*ps[ip0].vs;
         apval_ptr[ip].den = ps[ip0].den;

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

void srf_hfadj(struct standrupformat *srf,int nfine_stk,int nfine_dip,struct vruppars *vr,float *shyp,float *dhyp)
{
struct srf_prectsegments *prseg;
struct srf_apointvalues *apval;

float shyp_seg, dhyp_seg, flen_tot, flen_back;
float rfdep, rspd, dmin1, dmin2, dmax1, dmax2;
float dlen_f, dwid_f, tp0, *tpf, *rtf;
float tmin, tmax1, tmax2, tmax3, tend, sfac;
float xp, yp, xf, yf;

int it, ig, ip, ix, iy, ixs, iys, ips;
int it0, it_out, ntot, nts, nt_tmp;

float *stfp_i, *stfp_o;
float *y2vals, *utmp, *stfp_t;
float dt_faux;

long seed = 12345678;
float fzero = 0.0;
float fone = 1.0;

/* tsfac & rtime stuff */
float mom, mag;
struct complex *tsfac1_c, *rtime1_c;
float *tsfac1_r, *rtime1_r;

int nstk_tot, nstk_back, nstk_f, ndip_f;
int ipf, ipf_0;
float clen_s, clen_d, dstk_f, ddip_f, dks_f, dkd_f;
float wavelength_max, wavelength_min;
float tsf1_max, tsf1_avg, tsf1_min, tsf1_sig;

float tsfac_main = -1.0e+15;
float tsfac_bzero = -0.1;
float tsfac_slope = -0.5;
float tsfac1_sigma = 1.0;

float risetime_coef = 2.2;
float rtime1_sigma = 0.75;
float trise, rt1_max, rt1_avg, rt1_min, rt1_sig;

int kmodel = 2; /* MAI_FLAG */

dmin1 = vr[0].shal_vrup_dep - vr[0].shal_vrup_deprange;
dmax1 = vr[0].shal_vrup_dep + vr[0].shal_vrup_deprange;
dmin2 = vr[0].deep_vrup_dep - vr[0].deep_vrup_deprange;
dmax2 = vr[0].deep_vrup_dep + vr[0].deep_vrup_deprange;

tpf = (float *)check_malloc(nfine_stk*nfine_dip*sizeof(float));
rtf = (float *)check_malloc(nfine_stk*nfine_dip*sizeof(float));
stfp_i = NULL;
stfp_t = NULL;
y2vals = NULL;
utmp = NULL;

prseg = srf[0].srf_prect.prectseg;
apval = srf[0].srf_apnts.apntvals;

flen_back = 0.0;
flen_tot = 0.0;
nstk_tot = 0;
for(ig=0;ig<srf[0].nseg;ig++)
   {
   flen_tot = flen_tot + prseg[ig].flen;
   nstk_tot = nstk_tot + prseg[ig].nstk;
   }

/* get moment and correlation lengths */
mom = 0.0;
for(ip=0;ip<srf[0].srf_apnts.np;ip++)
   {
   sfac = sqrt(apval[ip].slip1*apval[ip].slip1 +
                     apval[ip].slip2*apval[ip].slip2 +
                     apval[ip].slip3*apval[ip].slip3);

   mom = mom + sfac*apval[ip].area*apval[ip].vs*apval[ip].vs*apval[ip].den;
   }
mag = 2.0*(log(mom)/log(10.0))/3.0 - 10.73;

fprintf(stderr,"mag= %.2f mom= %.5e\n",mag,mom);

clen_s = exp(log(10.0)*(0.5*mag - 2.50));
clen_d = exp(log(10.0)*(1.0*mag/3.0 - 1.50));

nstk_f = nfine_stk*nstk_tot;
if(nstk_f%2)
   nstk_f++;

ndip_f = nfine_dip*prseg[0].ndip;
if(ndip_f%2)
   ndip_f++;

dstk_f = flen_tot/(1.0*nstk_f);
ddip_f = prseg[0].fwid/(1.0*ndip_f);

dks_f = 1.0/(nstk_f*dstk_f);
dkd_f = 1.0/(ndip_f*ddip_f);

wavelength_max = 1.0e+15;
wavelength_min = 2.0*sqrt(dstk_f*ddip_f)/0.8;    /* 80% of nyquist */

/* start tsfac1 */

tsfac_main = tsfac_bzero + tsfac_slope*(1.0e-09*exp(log(mom)/3.0));

tsfac1_c = (struct complex *) check_malloc (nstk_f*ndip_f*sizeof(struct complex));
tsfac1_r = (float *) check_malloc (nstk_f*ndip_f*sizeof(float));

for(ip=0;ip<ndip_f*nstk_f;ip++)
   {
   tsfac1_c[ip].re = 1.0;
   tsfac1_c[ip].im = 0.0;
   }

fft2d_fftw(tsfac1_c,nstk_f,ndip_f,-1,&dstk_f,&ddip_f);

kfilt_gaus2(tsfac1_c,nstk_f,ndip_f,&dks_f,&dkd_f,&clen_s,&clen_d,&seed,kmodel,&wavelength_max,&wavelength_min);

fprintf(stderr,"\n     tsfac1_c[0].re= %14.6e tsfac1_c[0].im= %14.6e\n\n",tsfac1_c[0].re,tsfac1_c[0].im);

tsf1_avg = 0.0;
for(iy=0;iy<ndip_f;iy++)
   {
   for(ix=0;ix<nstk_f;ix++)
      {
      ip = ix + iy*nstk_f;
      tsfac1_r[ip] = tsfac1_c[ip].re;

      tsf1_avg = tsf1_avg + tsfac1_r[ip];
      }
   }
fprintf(stderr,"tsf1_avg= %13.5e, removed from field\n",tsf1_avg/(float)(nstk_f*ndip_f));

tsf1_avg = tsf1_avg/(float)(nstk_f*ndip_f);
for(ip=0;ip<ndip_f*nstk_f;ip++)
   tsfac1_r[ip] = tsfac1_r[ip] - tsf1_avg;

tsf1_avg = 0.0;

tsf1_sig = 0.0;
for(ip=0;ip<ndip_f*nstk_f;ip++)
   tsf1_sig = tsf1_sig + (tsfac1_r[ip])*(tsfac1_r[ip]);

tsf1_sig = sqrt(tsf1_sig/(ndip_f*nstk_f));

if(tsfac1_sigma < 0.0)
   tsfac1_sigma = 0.0;

sfac = tsfac1_sigma/tsf1_sig;
for(ip=0;ip<nstk_f*ndip_f;ip++)
   tsfac1_r[ip] = sfac*tsfac1_r[ip];

tsf1_min = 1.0e+15;
tsf1_max = -1.0e+15;
tsf1_sig = 0.0;
for(ip=0;ip<nstk_f*ndip_f;ip++)
   {
   if(tsfac1_r[ip] < tsf1_min)
      tsf1_min = tsfac1_r[ip];
   if(tsfac1_r[ip] > tsf1_max)
      tsf1_max = tsfac1_r[ip];

   tsf1_sig = tsf1_sig + (tsfac1_r[ip])*(tsfac1_r[ip]);
   }

tsf1_sig = sqrt(tsf1_sig/(ndip_f*nstk_f));

fprintf(stderr,"tsfac1: avg= %13.5e min= %13.5e max= %13.5e sig= %13.5e\n",tsf1_avg,tsf1_min,tsf1_max,tsf1_sig);

/* end tsfac1 */

/* start rtime1 */

trise = risetime_coef*(1.0e-09*exp(log(mom)/3.0));

rtime1_c = (struct complex *) check_malloc (nstk_f*ndip_f*sizeof(struct complex));
rtime1_r = (float *) check_malloc (nstk_f*ndip_f*sizeof(float));

for(ip=0;ip<ndip_f*nstk_f;ip++)
   {
   rtime1_c[ip].re = 1.0;
   rtime1_c[ip].im = 0.0;
   }

fft2d_fftw(rtime1_c,nstk_f,ndip_f,-1,&dstk_f,&ddip_f);

kfilt_gaus2(rtime1_c,nstk_f,ndip_f,&dks_f,&dkd_f,&clen_s,&clen_d,&seed,kmodel,&wavelength_max,&wavelength_min);

fprintf(stderr,"\n     rtime1_c[0].re= %14.6e rtime1_c[0].im= %14.6e\n\n",rtime1_c[0].re,rtime1_c[0].im);

rt1_avg = 0.0;
for(iy=0;iy<ndip_f;iy++)
   {
   for(ix=0;ix<nstk_f;ix++)
      {
      ip = ix + iy*nstk_f;
      rtime1_r[ip] = rtime1_c[ip].re;

      rt1_avg = rt1_avg + rtime1_r[ip];
      }
   }

fprintf(stderr,"rt1_avg= %13.5e, normalized to 1\n",rt1_avg/(float)(nstk_f*ndip_f));

rt1_avg = rt1_avg/(float)(nstk_f*ndip_f);
rt1_avg = sqrt(rt1_avg*rt1_avg);
for(ip=0;ip<ndip_f*nstk_f;ip++)
   rtime1_r[ip] = rtime1_r[ip]/rt1_avg;

rt1_avg = 1.0;

rt1_sig = 0.0;
for(ip=0;ip<ndip_f*nstk_f;ip++)
   rt1_sig = rt1_sig + (rtime1_r[ip]-rt1_avg)*(rtime1_r[ip]-rt1_avg);

rt1_sig = sqrt(rt1_sig/(ndip_f*nstk_f))/rt1_avg;

if(rtime1_sigma < 0.0)
   rtime1_sigma = 0.0;

fprintf(stderr,"rtime1_sigma(CoV)= %f\n",rtime1_sigma);

sfac = rtime1_sigma/rt1_sig;
for(ip=0;ip<nstk_f*ndip_f;ip++)
   rtime1_r[ip] = sfac*(rtime1_r[ip] - rt1_avg) + rt1_avg;

/* truncate any negative values, reset sig, min, max and avg */
sfac = 0.0;
rt1_avg = 0.0;
for(ip=0;ip<nstk_f*ndip_f;ip++)
   {
   if(rtime1_r[ip] <= 0.0)
      {
      sfac = sfac + 1.0;
      rtime1_r[ip] = 1.0e-15;
/*
*/
      }

   rt1_avg = rt1_avg + rtime1_r[ip];
   }
fprintf(stderr,"ratio (negative rt)/(positive rt)= %f\n",sfac/(nstk_f*ndip_f));

rt1_avg = rt1_avg/(float)(nstk_f*ndip_f);
for(ip=0;ip<ndip_f*nstk_f;ip++)
   rtime1_r[ip] = rtime1_r[ip]/rt1_avg;

rt1_avg = 1.0;

rt1_min = 1.0e+15;
rt1_max = -1.0e+15;
rt1_sig = 0.0;
for(ip=0;ip<nstk_f*ndip_f;ip++)
   {
   if(rtime1_r[ip] < rt1_min)
      rt1_min = rtime1_r[ip];
   if(rtime1_r[ip] > rt1_max)
      rt1_max = rtime1_r[ip];

   rt1_sig = rt1_sig + (rtime1_r[ip]-rt1_avg)*(rtime1_r[ip]-rt1_avg);
   }

rt1_sig = sqrt(rt1_sig/(ndip_f*nstk_f))/rt1_avg;

fprintf(stderr,"rtime1: avg= %13.5e min= %13.5e max= %13.5e sig= %13.5e\n",rt1_avg,rt1_min,rt1_max,rt1_sig);

/* end rtime1 */

ntot = 0;
nstk_back = 0;
for(ig=0;ig<srf[0].nseg;ig++)
   {
   if(*shyp > -1.0e+14)
      {
      shyp_seg = *shyp + 0.5*flen_tot - flen_back - 0.5*prseg[ig].flen;
      flen_back = flen_back + prseg[ig].flen;
      }
   else
      shyp_seg = prseg[ig].shyp;

   if(*dhyp > -1.0e+14)
      dhyp_seg = *dhyp;
   else
      dhyp_seg = prseg[ig].dhyp;

   prseg[ig].dlen = prseg[ig].flen/(1.0*prseg[ig].nstk);
   prseg[ig].dwid = prseg[ig].fwid/(1.0*prseg[ig].ndip);
   dlen_f = prseg[ig].dlen/(1.0*nfine_stk);
   dwid_f = prseg[ig].dwid/(1.0*nfine_dip);

/*
fprintf(stderr,"shyp_seg= %.5e dhyp_seg= %.5e dlen_f= %.5e dwid_f= %.5e\n",shyp_seg,dhyp_seg,dlen_f,dwid_f);
*/

   for(iy=0;iy<prseg[ig].ndip;iy++)
      {
      yp = (iy+0.5)*prseg[ig].dwid;

      for(ix=0;ix<prseg[ig].nstk;ix++)
         {
         xp = (ix+0.5)*prseg[ig].dlen - 0.5*prseg[ig].flen;
         ip = ix + iy*prseg[ig].nstk + ntot;

	 if(apval[ip].nt1 > 0 || apval[ip].nt2 > 0 || apval[ip].nt3 > 0)
	    {
            if(apval[ip].dep <= dmin1)
               rfdep = vr[0].shal_vrup;
            else if(apval[ip].dep < dmax1 && apval[ip].dep > dmin1)
               rfdep = (1.0 - (1.0 - vr[0].shal_vrup)*(dmax1-apval[ip].dep)/(dmax1-dmin1));
            else if(apval[ip].dep >= dmax1 && apval[ip].dep <= dmin2)
               rfdep = 1.0;
            else if(apval[ip].dep < dmax2 && apval[ip].dep > dmin2)
               rfdep = (1.0 - (1.0 - vr[0].deep_vrup)*(apval[ip].dep-dmin2)/(dmax2-dmin2));
            else
               rfdep = vr[0].deep_vrup;

            rspd = vr[0].rvfrac*rfdep*apval[ip].vs*1.0e-05;

/*
fprintf(stderr,"ip= %d xp= %.5e yp= %.5e rspd= %.5e\n",ip,xp,yp,rspd);
*/

	    tp0 = sqrt((xp-shyp_seg)*(xp-shyp_seg) + (yp-dhyp_seg)*(yp-dhyp_seg))/rspd;

	    ipf_0 = (ix+nstk_back)*nfine_stk + iy*nfine_dip*nstk_f;

	    rt1_avg = 0.0;
	    tmin = 1.0e+15;
	    tmax1 = -1.0e+15;
	    tmax2 = -1.0e+15;
	    tmax3 = -1.0e+15;
            for(iys=0;iys<nfine_dip;iys++)
               {
	       yf = yp - 0.5*prseg[ig].dwid + (iys+0.5)*dwid_f;

               for(ixs=0;ixs<nfine_stk;ixs++)
                  {
	          xf = xp - 0.5*prseg[ig].dlen + (ixs+0.5)*dlen_f;
                  ips = ixs + iys*nfine_stk;

	          tpf[ips] = sqrt((xf-shyp_seg)*(xf-shyp_seg) + (yf-dhyp_seg)*(yf-dhyp_seg))/rspd;

/* adjust timing to center point tinit value */

		  ipf = ipf_0 + ixs + iys*nstk_f;
	          tpf[ips] = tpf[ips] - tp0 + apval[ip].tinit + tsfac_main*tsfac1_r[ipf];

/* load random rtime values and compute mean over subfault */

		  ipf = ipf_0 + ixs + iys*nstk_f;
	          rtf[ips] = rtime1_r[ipf];
                  rt1_avg = rt1_avg + rtf[ips];

/*
fprintf(stderr,"%d) tpf= %.5e tp0= %.5e tinit= %.5e\n",ip,tpf[ips],tp0,apval[ip].tinit);
*/

	          if(tpf[ips] < tmin)
	             tmin = tpf[ips];
         
                  if(apval[ip].nt1 > 0)
	             {
	             tend = tpf[ips] + apval[ip].nt1*apval[ip].dt;
	             if(tend > tmax1)
	                tmax1 = tend;
		     }
         
                  if(apval[ip].nt2 > 0)
	             {
	             tend = tpf[ips] + apval[ip].nt2*apval[ip].dt;
	             if(tend > tmax2)
	                tmax2 = tend;
		     }
         
                  if(apval[ip].nt3 > 0)
	             {
	             tend = tpf[ips] + apval[ip].nt3*apval[ip].dt;
	             if(tend > tmax3)
	                tmax3 = tend;
		     }
if(ip == 100452)
   /*
   fprintf(stderr,"%10.4f",tpf[ips]);
   */
   fprintf(stderr,"%10.4f",rtf[ips]*apval[ip].nt1*apval[ip].dt);
	          }
if(ip == 100452)
   fprintf(stderr,"\n");
	       }
if(ip == 100452)
   fprintf(stderr,"rtf_avg= %.5e\n",rt1_avg/(nfine_stk*nfine_dip));

            apval[ip].tinit = tmin;

/* adjust mean of random rise time factors to unity and recompute tmax's */

	    rt1_avg = nfine_stk*nfine_dip/rt1_avg;
            for(ips=0;ips<nfine_stk*nfine_dip;ips++)
	       rtf[ips] = rtf[ips]*rt1_avg;

/*
for(ips=0;ips<nfine_stk*nfine_dip;ips++)
   rtf[ips] = 2.0*(rtf[ips]-1.0) + 1.0;
*/

	    tmax1 = -1.0e+15;
	    tmax2 = -1.0e+15;
	    tmax3 = -1.0e+15;
            for(ips=0;ips<nfine_stk*nfine_dip;ips++)
	       {
               if(apval[ip].nt1 > 0)
	          {
	          tend = tpf[ips] + rtf[ips]*apval[ip].nt1*apval[ip].dt;
	          if(tend > tmax1)
	             tmax1 = tend;
		  }
         
               if(apval[ip].nt2 > 0)
	          {
	          tend = tpf[ips] + rtf[ips]*apval[ip].nt2*apval[ip].dt;
	          if(tend > tmax2)
	             tmax2 = tend;
		  }
         
               if(apval[ip].nt3 > 0)
	          {
	          tend = tpf[ips] + rtf[ips]*apval[ip].nt3*apval[ip].dt;
	          if(tend > tmax3)
	             tmax3 = tend;
		  }
	       }

            if(apval[ip].nt1 > 0)
	       {
	       nts = apval[ip].nt1;
               stfp_i = (float *)check_realloc(stfp_i,nts*sizeof(float));

               stfp_o = apval[ip].stf1;
	       for(it=0;it<nts;it++)
	          stfp_i[it] = stfp_o[it];

               apval[ip].nt1 = (int)((tmax1-tmin)/apval[ip].dt + 0.5);
               apval[ip].stf1 = (float *)check_realloc(apval[ip].stf1,apval[ip].nt1*sizeof(float));
               stfp_o = apval[ip].stf1;
	       for(it=0;it<apval[ip].nt1;it++)
	          stfp_o[it] = 0.0;

               stfp_t = (float *)check_realloc(stfp_t,apval[ip].nt1*sizeof(float));
               y2vals = (float *)check_realloc(y2vals,nts*sizeof(float));
               utmp = (float *)check_realloc(utmp,nts*sizeof(float));

	       for(ips=0;ips<nfine_stk*nfine_dip;ips++)
	          {
	          it0 = (int)((tpf[ips] - apval[ip].tinit)/apval[ip].dt + 0.5);
	          if(it0 < 0)		/* shouldn't happen but ... */
	             it0 = 0;

		  dt_faux = rtf[ips]*apval[ip].dt;
                  spline_natural(nts,&(dt_faux),stfp_i,y2vals,utmp);

                  nt_tmp = (dt_faux*nts)/apval[ip].dt;
                  while(apval[ip].dt*nt_tmp > dt_faux*nts)
                     nt_tmp--;

	          for(it=0;it<nt_tmp;it++) /* just in case */
                     stfp_t[it] = 0.0;

                  splint_vector(nts,&(dt_faux),stfp_i,y2vals,nt_tmp,&(apval[ip].dt),stfp_t);
/*
XXXXX
*/

	          while((it0+nt_tmp) > apval[ip].nt1)	/* shouldn't happen but ... */
	             nt_tmp--;

/*
if(ip == 100452)
   fprintf(stderr,"%10.4f %d %d\n",dt_faux,nt_tmp,nts);
*/

/*
fprintf(stderr,"it0= %d nts= %d nt1= %d\n",it0,nts,apval[ip].nt1);
*/

		  sfac = apval[ip].dt/dt_faux;	/* adjust to match original slip */
		  sfac = 1.0;
	          for(it=0;it<nt_tmp;it++)
	             {
		     it_out = it + it0;
                     stfp_o[it_out] = stfp_o[it_out] + sfac*stfp_t[it];
                     }
		  }

/* ensure output slip-rate integrates to target slip */

	       sfac = 0.0;
               for(it=0;it<apval[ip].nt1;it++)
                  sfac = sfac + (apval[ip].dt)*stfp_o[it];

/*
fprintf(stderr,"%d) sfac= %.5e slip1= %.5e\n",ip,sfac,apval[ip].slip1);
exit(-1);
*/

	       sfac = apval[ip].slip1/sfac;
               for(it=0;it<apval[ip].nt1;it++)
                  stfp_o[it] = sfac*stfp_o[it];
	       }

            if(apval[ip].nt2 > 0)
	       {
	       nts = apval[ip].nt2;
               stfp_i = (float *)check_realloc(stfp_i,nts*sizeof(float));

               stfp_o = apval[ip].stf2;
	       for(it=0;it<nts;it++)
	          stfp_i[it] = stfp_o[it];

               apval[ip].nt2 = (int)((tmax2-tmin)/apval[ip].dt + 0.5);
               apval[ip].stf2 = (float *)check_realloc(apval[ip].stf2,apval[ip].nt2*sizeof(float));
               stfp_o = apval[ip].stf2;
	       for(it=0;it<apval[ip].nt2;it++)
	          stfp_o[it] = 0.0;

	       for(ips=0;ips<nfine_stk*nfine_dip;ips++)
	          {
	          it0 = (int)((tpf[ips] - apval[ip].tinit)/apval[ip].dt + 0.5);
	          if(it0 < 0)		/* shouldn't happen but ... */
	             it0 = 0;

		  nt_tmp = nts;
	          while((it0+nt_tmp) > apval[ip].nt2)	/* shouldn't happen but ... */
	             nt_tmp--;

	          for(it=0;it<nt_tmp;it++)
	             {
		     it_out = it + it0;
                     stfp_o[it_out] = stfp_o[it_out] + stfp_i[it];
                     }
		  }

/* ensure output slip-rate integrates to target slip */

	       sfac = 0.0;
               for(it=0;it<apval[ip].nt2;it++)
                  sfac = sfac + (apval[ip].dt)*stfp_o[it];

	       sfac = apval[ip].slip2/sfac;
               for(it=0;it<apval[ip].nt2;it++)
                  stfp_o[it] = sfac*stfp_o[it];
	       }

            if(apval[ip].nt3 > 0)
	       {
	       nts = apval[ip].nt3;
               stfp_i = (float *)check_realloc(stfp_i,nts*sizeof(float));

               stfp_o = apval[ip].stf3;
	       for(it=0;it<nts;it++)
	          stfp_i[it] = stfp_o[it];

               apval[ip].nt3 = (int)((tmax3-tmin)/apval[ip].dt + 0.5);
               apval[ip].stf3 = (float *)check_realloc(apval[ip].stf3,apval[ip].nt3*sizeof(float));
               stfp_o = apval[ip].stf3;
	       for(it=0;it<apval[ip].nt3;it++)
	          stfp_o[it] = 0.0;

	       for(ips=0;ips<nfine_stk*nfine_dip;ips++)
	          {
	          it0 = (int)((tpf[ips] - apval[ip].tinit)/apval[ip].dt + 0.5);
	          if(it0 < 0)		/* shouldn't happen but ... */
	             it0 = 0;

		  nt_tmp = nts;
	          while((it0+nt_tmp) > apval[ip].nt3)	/* shouldn't happen but ... */
	             nt_tmp--;

	          for(it=0;it<nt_tmp;it++)
	             {
		     it_out = it + it0;
                     stfp_o[it_out] = stfp_o[it_out] + stfp_i[it];
                     }
		  }

/* ensure output slip-rate integrates to target slip */

	       sfac = 0.0;
               for(it=0;it<apval[ip].nt3;it++)
                  sfac = sfac + (apval[ip].dt)*stfp_o[it];

	       sfac = apval[ip].slip3/sfac;
               for(it=0;it<apval[ip].nt3;it++)
                  stfp_o[it] = sfac*stfp_o[it];
	       }

if(ip == 100452)
   {
   xf = ((ix+0.5)/(1.0*prseg[ig].nstk) - 0.5)*prseg[ig].flen - prseg[ig].shyp;
   yf = (iy+0.5)*prseg[ig].fwid/(1.0*prseg[ig].ndip) - prseg[ig].dhyp;
   fprintf(stderr,"tmin= %10.4f nt1_new= %d nt1_orig= %d rhyp= %10.4f\n",tmin,apval[ip].nt1,nts,sqrt(xf*xf + yf*yf));
   fprintf(stderr,"%.6f %.6f %.5e\n",apval[ip].lon,apval[ip].lat,apval[ip].dep);
   }
            }
         }
      }

   ntot = ntot + srf[0].np_seg[ig];
   nstk_back = nstk_back + prseg[ig].nstk;
   }

free(stfp_i);
free(tpf);
free(tsfac1_c);
free(tsfac1_r);
free(rtf);
free(rtime1_c);
free(rtime1_r);
}

void load_slip_srf_dd4_v553(struct standrupformat *srf,struct stfpar2 *spar,struct pointsource *ps,float *rtime1_r,float *rtime2_r,struct velmodel *vmod)
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

float beta, b_dmin, b_dmax, delta_b;

rtfac_min = 5.0*spar->dt/spar->trise;    /* so minimum rise time >= 5*dt */
rtfac_min = 10.0*spar->dt/spar->trise;    /* so minimum rise time >= 10*dt */
rtfac_min = spar->dt/spar->trise;    /* so minimum rise time >= dt */

dmin1 = spar->risetimedep - spar->risetimedep_range;
dmax1 = spar->risetimedep + spar->risetimedep_range;
rtfac1 = spar->risetimefac - 1.0;

dmin2 = spar->deep_risetimedep - spar->deep_risetimedep_range;
dmax2 = spar->deep_risetimedep + spar->deep_risetimedep_range;
rtfac2 = spar->deep_risetimefac - 1.0;

b_dmin = spar->beta_depth - spar->beta_depth_range;
b_dmax = spar->beta_depth + spar->beta_depth_range;
delta_b = spar->beta_deep - spar->beta_shal;

apnts_ptr = &(srf->srf_apnts);
apval_ptr = apnts_ptr->apntvals;

nseg = srf[0].srf_prect.nseg;
prseg_ptr = srf[0].srf_prect.prectseg;

/* version 2.0 stuff, just in case */

srf[0].srf_hcmnt.nline = 0;

/* 2018-04-02 updated for multiple segs */
srf[0].nseg = srf[0].srf_prect.nseg;
srf[0].np_seg = (int *)check_malloc((srf[0].nseg)*sizeof(int));

for(iseg=0;iseg<nseg;iseg++)
   srf[0].np_seg[iseg] = prseg_ptr[iseg].nstk*prseg_ptr[0].ndip;

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
/* 20260504
               rtfac = rtfac*sqrt(rtime1_r[ip0])/spar->rt_scalefac;
*/
               rtfac = rtfac*rtime1_r[ip0]/spar->rt_scalefac;

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
            else if(strcmp(spar->stype,"OliuP") == 0) /* original form */
               {
               if(ps[ip0].dep <= b_dmin)
                  beta = spar->beta_shal;
               else if(ps[ip0].dep < b_dmax && ps[ip0].dep > b_dmin)
                  beta = spar->beta_shal + delta_b*(ps[ip0].dep - b_dmin)/(b_dmax-b_dmin);
               else
                  beta = spar->beta_deep;

               tzero = rtfac*spar->trise;

               apval_ptr[ip].nt1 = gen_OliuP_stf(&(ps[ip0].slip),&tzero,&beta,stf,spar->nt,&spar->dt);
               }
            else if(strcmp(spar->stype,"MliuP") == 0) /* modified form */
               {
               if(ps[ip0].dep <= b_dmin)
                  beta = spar->beta_shal;
               else if(ps[ip0].dep < b_dmax && ps[ip0].dep > b_dmin)
                  beta = spar->beta_shal + delta_b*(ps[ip0].dep - b_dmin)/(b_dmax-b_dmin);
               else
                  beta = spar->beta_deep;

               tzero = rtfac*spar->trise;

               apval_ptr[ip].nt1 = gen_MliuP_stf(&(ps[ip0].slip),&tzero,&beta,stf,spar->nt,&spar->dt);
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
