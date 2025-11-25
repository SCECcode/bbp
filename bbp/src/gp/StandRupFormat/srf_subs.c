#include "include.h"
#include "structure.h"
#include "function.h"
#include "defs.h"

void read_srf(struct standrupformat *srf,char *file,int bflag)
{
FILE *fpr, *fopfile();
struct srf_prectsegments *prseg_ptr;
struct srf_apointvalues *apval_ptr;
char *sptr, str[MAXLINE];

float *stf;
int i, j, k, it, ig;

char pword[32];
int fdr;

int ip, np_seg, npread, np_tot;
long foff;

int mrf1_flag, mrf6_flag;

if(bflag)
   {
   if(strcmp(file,"stdin") == 0)
      fdr = STDIN_FILENO;
   else
      fdr = opfile_ro(file);

   reed(fdr,srf->version,sizeof(srf->version));

   reed(fdr,pword,sizeof(pword));
   if(strcmp(pword,"PLANE") == 0)
      {
      sprintf(srf->type,"PLANE");

      reed(fdr,&(srf[0].srf_prect.nseg),sizeof(int));
      srf[0].srf_prect.prectseg = (struct srf_prectsegments *)check_malloc(srf[0].srf_prect.nseg*sizeof(struct srf_prectsegments));
      prseg_ptr = srf[0].srf_prect.prectseg;

      reed(fdr,prseg_ptr,(srf[0].srf_prect.nseg)*sizeof(struct srf_prectsegments));

      while(strncmp(pword,"POINTS",6) != 0)
         reed(fdr,pword,sizeof(pword));
      }

   if(strncmp(pword,"POINTS",6) == 0)
      {
      reed(fdr,&(srf[0].srf_apnts.np),sizeof(int));
      srf[0].srf_apnts.apntvals = (struct srf_apointvalues *)check_malloc((srf[0].srf_apnts.np)*sizeof(struct srf_apointvalues));

      apval_ptr = srf[0].srf_apnts.apntvals;

      for(i=0;i<srf[0].srf_apnts.np;i++)
         {
         reed(fdr,&(apval_ptr[i].lon),sizeof(float));
         reed(fdr,&(apval_ptr[i].lat),sizeof(float));
         reed(fdr,&(apval_ptr[i].dep),sizeof(float));
         reed(fdr,&(apval_ptr[i].stk),sizeof(float));
         reed(fdr,&(apval_ptr[i].dip),sizeof(float));
         reed(fdr,&(apval_ptr[i].area),sizeof(float));
         reed(fdr,&(apval_ptr[i].tinit),sizeof(float));
         reed(fdr,&(apval_ptr[i].dt),sizeof(float));
         reed(fdr,&(apval_ptr[i].rake),sizeof(float));
         reed(fdr,&(apval_ptr[i].slip1),sizeof(float));
         reed(fdr,&(apval_ptr[i].nt1),sizeof(int));
         reed(fdr,&(apval_ptr[i].slip2),sizeof(float));
         reed(fdr,&(apval_ptr[i].nt2),sizeof(int));
         reed(fdr,&(apval_ptr[i].slip3),sizeof(float));
         reed(fdr,&(apval_ptr[i].nt3),sizeof(int));

         apval_ptr[i].stf1 = (float *)check_malloc((apval_ptr[i].nt1)*sizeof(float));
         apval_ptr[i].stf2 = (float *)check_malloc((apval_ptr[i].nt2)*sizeof(float));
         apval_ptr[i].stf3 = (float *)check_malloc((apval_ptr[i].nt3)*sizeof(float));

         reed(fdr,apval_ptr[i].stf1,(apval_ptr[i].nt1)*sizeof(float));
         reed(fdr,apval_ptr[i].stf2,(apval_ptr[i].nt2)*sizeof(float));
         reed(fdr,apval_ptr[i].stf3,(apval_ptr[i].nt3)*sizeof(float));
         }
      }
   close(fdr);
   }
else
   {
   if(strcmp(file,"stdin") == 0)
      fpr = stdin;
   else
      fpr = fopfile(file,"r");

   fgets(str,MAXLINE,fpr);
   sscanf(str,"%s",srf[0].version);

   mrf1_flag = 0;
   mrf6_flag = 0;
   sprintf(srf[0].src_format,"SLIP");
   if(atof(srf[0].version) >= 3.0)
      {
      ip = sscanf(str,"%*s %s",srf[0].src_format);

      if(ip != 1 || (strcmp(srf[0].src_format,"SLIP") != 0 &&
                     strcmp(srf[0].src_format,"MOMENT") != 0 &&
                     strcmp(srf[0].src_format,"MOMENT-1MECH") != 0 &&
                     strcmp(srf[0].src_format,"MOMENT-6MECH") != 0))
         sprintf(srf[0].src_format,"SLIP");

      if(strcmp(srf[0].src_format,"MOMENT") == 0 || strcmp(srf[0].src_format,"MOMENT-1MECH") == 0)
         mrf1_flag = 1;

      if(strcmp(srf[0].src_format,"MOMENT-6MECH") == 0)
         mrf6_flag = 1;
      }

   /* reserve first 3 lines for writing out command that generated SRF */
   srf[0].srf_hcmnt.nline = 3;
   srf[0].srf_hcmnt.cbuf = (char *)check_malloc((srf[0].srf_hcmnt.nline)*MAXLINE*sizeof(char));

   /* initialize first 3 lines with NULL, if not reset, won't print out */
   srf[0].srf_hcmnt.cbuf[0] = '\0';
   srf[0].srf_hcmnt.cbuf[MAXLINE] = '\0';
   srf[0].srf_hcmnt.cbuf[2*MAXLINE] = '\0';

   fgets(str,MAXLINE,fpr);
   while(str[0] == '#')
      {
      srf[0].srf_hcmnt.nline++;
      srf[0].srf_hcmnt.cbuf = (char *)check_realloc(srf[0].srf_hcmnt.cbuf,(srf[0].srf_hcmnt.nline)*MAXLINE*sizeof(char));

      /* get rid of annoying newline */
      i = 0;
      while(str[i] != '\n' && str[i] != '\0')
         i++;
      str[i] = '\0';

      sptr = srf[0].srf_hcmnt.cbuf + (srf[0].srf_hcmnt.nline-1)*MAXLINE;
      sprintf(sptr,"%s",str);

      /* just to make sure string ends without issue */
      sptr[MAXLINE-1] = '\0';

      if(fgets(str,MAXLINE,fpr) == NULL) 
         break;            
      }                    

   sscanf(str,"%s",pword);

   if(strncmp(pword,"PLANE",5) == 0)
      {
      sscanf(str,"%s %d",srf[0].type,&(srf[0].srf_prect.nseg));

      srf[0].srf_prect.prectseg = (struct srf_prectsegments *)check_malloc(srf[0].srf_prect.nseg*sizeof(struct srf_prectsegments));
      prseg_ptr = srf[0].srf_prect.prectseg;

      for(ig=0;ig<srf[0].srf_prect.nseg;ig++)
         {
         fgets(str,MAXLINE,fpr);
         sscanf(str,"%f %f %d %d %f %f",&(prseg_ptr[ig].elon),
                                     &(prseg_ptr[ig].elat),
                                     &(prseg_ptr[ig].nstk),
                                     &(prseg_ptr[ig].ndip),
                                     &(prseg_ptr[ig].flen),
                                     &(prseg_ptr[ig].fwid));
         fgets(str,MAXLINE,fpr);
         sscanf(str,"%f %f %f %f %f",&(prseg_ptr[ig].stk),
                                  &(prseg_ptr[ig].dip),
                                  &(prseg_ptr[ig].dtop),
                                  &(prseg_ptr[ig].shyp),
                                  &(prseg_ptr[ig].dhyp));
         }
      }

   if(strncmp(pword,"POINTS",6) != 0)
      {
      fgets(str,MAXLINE,fpr);
      while(str[0] == '#')
         {
         srf[0].srf_hcmnt.nline++;
         srf[0].srf_hcmnt.cbuf = (char *)check_realloc(srf[0].srf_hcmnt.cbuf,(srf[0].srf_hcmnt.nline)*MAXLINE*sizeof(char));

         sptr = srf[0].srf_hcmnt.cbuf + (srf[0].srf_hcmnt.nline-1)*MAXLINE;
         sprintf(sptr,"%s",str);

         /* just to make sure string ends without issue */
         sptr[MAXLINE-2] = '\n';
         sptr[MAXLINE-1] = '\0';

         if(fgets(str,MAXLINE,fpr) == NULL)
            break;
         }

      sscanf(str,"%s",pword);
      }

   if(strncmp(pword,"POINTS",6) == 0)
      {
      srf[0].nseg = 0;
      srf[0].srf_apnts.np = 0;

      srf[0].np_seg = NULL;
      srf[0].srf_apnts.apntvals = NULL;

      while(strncmp(pword,"POINTS",6) == 0)
         {
         srf[0].nseg++;
         srf[0].np_seg = (int *)check_realloc(srf[0].np_seg,(srf[0].nseg)*sizeof(int));

         sscanf(str,"%*s %d",&np_seg);

	 np_tot = np_seg + srf[0].srf_apnts.np;

         srf[0].srf_apnts.apntvals = (struct srf_apointvalues *)check_realloc(srf[0].srf_apnts.apntvals,np_tot*sizeof(struct srf_apointvalues));

         npread = 0;
         for(i=0;i<np_seg;i++)
            {
	    ip = i + srf[0].srf_apnts.np;
            apval_ptr = &(srf[0].srf_apnts.apntvals[ip]);

	    foff = ftell(fpr);
            if(fgets(str,MAXLINE,fpr) == NULL)
	       break;
	    else
	       {
               sscanf(str,"%s",pword);
               if(strncmp(pword,"POINTS",6) == 0)
                  {
		  fseek(fpr,foff,SEEK_SET);
		  break;
		  }

	       npread++;

               if(atof(srf->version) < 2.0)
                  {
                  sscanf(str,"%f %f %f %f %f %f %f %f",&(apval_ptr->lon),
                                           &(apval_ptr->lat),
                                           &(apval_ptr->dep),
                                           &(apval_ptr->stk),
                                           &(apval_ptr->dip),
                                           &(apval_ptr->area),
                                           &(apval_ptr->tinit),
                                           &(apval_ptr->dt));

                  apval_ptr->vp = -1;
                  apval_ptr->vs = -1;
                  apval_ptr->den = -1;
	          }
               else if(atof(srf->version) < 3.0)
                  {
                  sscanf(str,"%f %f %f %f %f %f %f %f %f %f",&(apval_ptr->lon),
                                           &(apval_ptr->lat),
                                           &(apval_ptr->dep),
                                           &(apval_ptr->stk),
                                           &(apval_ptr->dip),
                                           &(apval_ptr->area),
                                           &(apval_ptr->tinit),
                                           &(apval_ptr->dt),
                                           &(apval_ptr->vs),
                                           &(apval_ptr->den));
                  apval_ptr->vp = -1;
	          }
               else if(atof(srf->version) >= 3.0)
                  {
                  sscanf(str,"%f %f %f %f %f %f %f %f %f %f %f",&(apval_ptr->lon),
                                           &(apval_ptr->lat),
                                           &(apval_ptr->dep),
                                           &(apval_ptr->stk),
                                           &(apval_ptr->dip),
                                           &(apval_ptr->area),
                                           &(apval_ptr->tinit),
                                           &(apval_ptr->dt),
                                           &(apval_ptr->vp),
                                           &(apval_ptr->vs),
                                           &(apval_ptr->den));
	          }

               if(mrf1_flag == 0 && mrf6_flag == 0)           /* expecting slip */
                  {
                  fgets(str,MAXLINE,fpr);
                  sscanf(str,"%f %f %d %f %d %f %d",&(apval_ptr->rake),
                                           &(apval_ptr->slip1),
                                           &(apval_ptr->nt1),
                                           &(apval_ptr->slip2),
                                           &(apval_ptr->nt2),
                                           &(apval_ptr->slip3),
                                           &(apval_ptr->nt3));

	          if(apval_ptr->nt1)
                     apval_ptr->stf1 = (float *)check_malloc((apval_ptr->nt1)*sizeof(float));
	          else
                     apval_ptr->stf1 = NULL;

                  stf = apval_ptr->stf1;

                  for(it=0;it<(apval_ptr->nt1);it++)
                     fscanf(fpr,"%f",&stf[it]);

	          if(apval_ptr->nt2)
                     apval_ptr->stf2 = (float *)check_malloc((apval_ptr->nt2)*sizeof(float));
	          else
                     apval_ptr->stf2 = NULL;

                  stf = apval_ptr->stf2;

                  for(it=0;it<(apval_ptr->nt2);it++)
                     fscanf(fpr,"%f",&stf[it]);

	          if(apval_ptr->nt3)
                     apval_ptr->stf3 = (float *)check_malloc((apval_ptr->nt3)*sizeof(float));
	          else
                     apval_ptr->stf3 = NULL;

                  stf = apval_ptr->stf3;

                  for(it=0;it<(apval_ptr->nt3);it++)
                     fscanf(fpr,"%f",&stf[it]);

                  /* get rouge newline character */
                  if((apval_ptr->nt1) || (apval_ptr->nt2) || (apval_ptr->nt3))
                     fgets(str,MAXLINE,fpr);
	          }

               else         /* expecting moment */
                  {
                  fgets(str,MAXLINE,fpr);
                  sscanf(str,"%lf %lf %lf %lf %lf %lf %d",
                                           &(apval_ptr->mnn),
                                           &(apval_ptr->mee),
                                           &(apval_ptr->mdd),
                                           &(apval_ptr->mne),
                                           &(apval_ptr->mnd),
                                           &(apval_ptr->med),
                                           &(apval_ptr->ntmr));

                  if(mrf1_flag == 1)         /* single time function */
		     {
	             if(apval_ptr->ntmr)
                        apval_ptr->mrf = (float *)check_malloc((apval_ptr->ntmr)*sizeof(float));
	             else
                        apval_ptr->mrf = NULL;

                     stf = apval_ptr->mrf;

                     for(it=0;it<(apval_ptr->ntmr);it++)
                        fscanf(fpr,"%f",&stf[it]);
		     }

                  else if(mrf6_flag == 1)         /* 6 time functions */
		     {
	             if(apval_ptr->ntmr)
                        apval_ptr->mr_nn = (float *)check_malloc((apval_ptr->ntmr)*sizeof(float));
	             else
                        apval_ptr->mr_nn = NULL;

                     stf = apval_ptr->mr_nn;

                     for(it=0;it<(apval_ptr->ntmr);it++)
                        fscanf(fpr,"%f",&stf[it]);

	             if(apval_ptr->ntmr)
                        apval_ptr->mr_ee = (float *)check_malloc((apval_ptr->ntmr)*sizeof(float));
	             else
                        apval_ptr->mr_ee = NULL;

                     stf = apval_ptr->mr_ee;

                     for(it=0;it<(apval_ptr->ntmr);it++)
                        fscanf(fpr,"%f",&stf[it]);

	             if(apval_ptr->ntmr)
                        apval_ptr->mr_dd = (float *)check_malloc((apval_ptr->ntmr)*sizeof(float));
	             else
                        apval_ptr->mr_dd = NULL;

                     stf = apval_ptr->mr_dd;

                     for(it=0;it<(apval_ptr->ntmr);it++)
                        fscanf(fpr,"%f",&stf[it]);

	             if(apval_ptr->ntmr)
                        apval_ptr->mr_ne = (float *)check_malloc((apval_ptr->ntmr)*sizeof(float));
	             else
                        apval_ptr->mr_ne = NULL;

                     stf = apval_ptr->mr_ne;

                     for(it=0;it<(apval_ptr->ntmr);it++)
                        fscanf(fpr,"%f",&stf[it]);

	             if(apval_ptr->ntmr)
                        apval_ptr->mr_nd = (float *)check_malloc((apval_ptr->ntmr)*sizeof(float));
	             else
                        apval_ptr->mr_nd = NULL;

                     stf = apval_ptr->mr_nd;

                     for(it=0;it<(apval_ptr->ntmr);it++)
                        fscanf(fpr,"%f",&stf[it]);

	             if(apval_ptr->ntmr)
                        apval_ptr->mr_ed = (float *)check_malloc((apval_ptr->ntmr)*sizeof(float));
	             else
                        apval_ptr->mr_ed = NULL;

                     stf = apval_ptr->mr_ed;

                     for(it=0;it<(apval_ptr->ntmr);it++)
                        fscanf(fpr,"%f",&stf[it]);

		     }

                  /* get rouge newline character */
                  if((apval_ptr->ntmr))
                     fgets(str,MAXLINE,fpr);
	          }

	       }
            }

	 srf[0].np_seg[(srf[0].nseg)-1] = npread;
         srf[0].srf_apnts.np = srf[0].srf_apnts.np + npread;

         if(fgets(str,MAXLINE,fpr) == NULL)
            break;
         else
	    {
            sscanf(str,"%s",pword);
            while(strncmp(pword,"POINTS",6) != 0)
               {
               if(fgets(str,MAXLINE,fpr) == NULL)
                  break;
               else
                  sscanf(str,"%s",pword);
               }
	    }
         }

      if(np_tot > srf[0].srf_apnts.np)
         srf[0].srf_apnts.apntvals = (struct srf_apointvalues *)check_realloc(srf[0].srf_apnts.apntvals,(srf[0].srf_apnts.np)*sizeof(struct srf_apointvalues));
      }
   else
      {
      fprintf(stderr,"*** No valid data found.  Last input line read is:\n");
      fprintf(stderr,"%s",str);
      exit(-99);
      }

   fclose(fpr);
   }
}

void write_srf(struct standrupformat *srf,char *file,int bflag)
{
if(atof(srf->version) < 2.0)
   write_srf1(srf,file,bflag);
else if(atof(srf->version) < 3.0)
   write_srf2(srf,file,bflag);
else if(atof(srf->version) >= 3.0)
   write_srf3(srf,file,bflag);
}

void write_srf1(struct standrupformat *srf,char *file,int bflag)
{
FILE *fpw, *fopfile();
struct srf_planerectangle *prect_ptr;
struct srf_prectsegments *prseg_ptr;
struct srf_allpoints *apnts_ptr;
struct srf_apointvalues *apval_ptr;

float area;
float *stf;
int i, j, k, nt6, it, ip, ig, ntot;

char pword[32];
int fdw;

prect_ptr = &(srf->srf_prect);
prseg_ptr = prect_ptr->prectseg;
apnts_ptr = &(srf->srf_apnts);
apval_ptr = apnts_ptr->apntvals;

if(bflag)
   {
   if(strcmp(file,"stdout") == 0)
      fdw = STDOUT_FILENO;
   else
      fdw = croptrfile(file);

   rite(fdw,srf->version,sizeof(srf->version));

   if(strcmp(srf->type,"PLANE") == 0)
      {
      rite(fdw,srf->type,sizeof(srf->type));
      rite(fdw,&(prect_ptr->nseg),sizeof(prect_ptr->nseg));
      rite(fdw,prseg_ptr,(prect_ptr->nseg)*sizeof(struct srf_prectsegments));
      }

   sprintf(pword,"POINTS");
   rite(fdw,pword,sizeof(pword));
   rite(fdw,&(apnts_ptr->np),sizeof(apnts_ptr->np));
   for(i=0;i<apnts_ptr->np;i++)
      {
      rite(fdw,&(apval_ptr[i].lon),sizeof(float));
      rite(fdw,&(apval_ptr[i].lat),sizeof(float));
      rite(fdw,&(apval_ptr[i].dep),sizeof(float));
      rite(fdw,&(apval_ptr[i].stk),sizeof(float));
      rite(fdw,&(apval_ptr[i].dip),sizeof(float));
      rite(fdw,&(apval_ptr[i].area),sizeof(float));
      rite(fdw,&(apval_ptr[i].tinit),sizeof(float));
      rite(fdw,&(apval_ptr[i].dt),sizeof(float));
      rite(fdw,&(apval_ptr[i].rake),sizeof(float));
      rite(fdw,&(apval_ptr[i].slip1),sizeof(float));
      rite(fdw,&(apval_ptr[i].nt1),sizeof(int));
      rite(fdw,&(apval_ptr[i].slip2),sizeof(float));
      rite(fdw,&(apval_ptr[i].nt2),sizeof(int));
      rite(fdw,&(apval_ptr[i].slip3),sizeof(float));
      rite(fdw,&(apval_ptr[i].nt3),sizeof(int));

      rite(fdw,apval_ptr[i].stf1,(apval_ptr[i].nt1)*sizeof(float));
      rite(fdw,apval_ptr[i].stf2,(apval_ptr[i].nt2)*sizeof(float));
      rite(fdw,apval_ptr[i].stf3,(apval_ptr[i].nt3)*sizeof(float));
      }
   close(fdw);
   }
else
   {
   if(strcmp(file,"stdout") == 0)
      fpw = stdout;
   else
      fpw = fopfile(file,"w");

   fprintf(fpw,"%s\n",srf->version);

   if(strcmp(srf->type,"PLANE") == 0)
      {
      fprintf(fpw,"%s %d\n",srf->type,prect_ptr->nseg);
      for(ig=0;ig<prect_ptr->nseg;ig++)
         {
         fprintf(fpw,"%12.6f %11.6f %5d %5d %10.4f %10.4f\n",prseg_ptr[ig].elon,
                                                        prseg_ptr[ig].elat,
                                                        prseg_ptr[ig].nstk,
                                                        prseg_ptr[ig].ndip,
                                                        prseg_ptr[ig].flen,
                                                        prseg_ptr[ig].fwid);
         fprintf(fpw,"%4.0f %4.0f %10.4f %10.4f %10.4f\n",prseg_ptr[ig].stk,
                                                    prseg_ptr[ig].dip,
                                                    prseg_ptr[ig].dtop,
                                                    prseg_ptr[ig].shyp,
                                                    prseg_ptr[ig].dhyp);
         }
      }

   fprintf(fpw,"POINTS %d\n",apnts_ptr->np);
   for(i=0;i<apnts_ptr->np;i++)
      {
      fprintf(fpw,"%12.6f %11.6f %12.5e %4.0f %4.0f %12.5e %10.4f %12.5e\n",
                                              apval_ptr[i].lon,
                                              apval_ptr[i].lat,
                                              apval_ptr[i].dep,
                                              apval_ptr[i].stk,
                                              apval_ptr[i].dip,
                                              apval_ptr[i].area,
                                              apval_ptr[i].tinit,
                                              apval_ptr[i].dt);
      fprintf(fpw,"%4.0f %10.4f %6d %10.4f %6d %10.4f %6d\n",
                                              apval_ptr[i].rake,
                                              apval_ptr[i].slip1,
                                              apval_ptr[i].nt1,
                                              apval_ptr[i].slip2,
                                              apval_ptr[i].nt2,
                                              apval_ptr[i].slip3,
                                              apval_ptr[i].nt3);

      stf = apval_ptr[i].stf1;
      nt6 = (apval_ptr[i].nt1)/6;
      for(k=0;k<nt6;k++)
         {
         for(j=0;j<6;j++)
            {
            it = 6*k + j;
            fprintf(fpw,"%13.5e",stf[it]);
            }
         fprintf(fpw,"\n");
         }

      if(6*nt6 != (apval_ptr[i].nt1))
         {
         for(j=6*nt6;j<(apval_ptr[i].nt1);j++)
            fprintf(fpw,"%13.5e",stf[j]);

         fprintf(fpw,"\n");
         }

      stf = apval_ptr[i].stf2;
      nt6 = (apval_ptr[i].nt2)/6;
      for(k=0;k<nt6;k++)
         {
         for(j=0;j<6;j++)
            {
            it = 6*k + j;
            fprintf(fpw,"%13.5e",stf[it]);
            }
         fprintf(fpw,"\n");
         }

      if(6*nt6 != (apval_ptr[i].nt2))
         {
         for(j=6*nt6;j<(apval_ptr[i].nt2);j++)
            fprintf(fpw,"%13.5e",stf[j]);

         fprintf(fpw,"\n");
         }

      stf = apval_ptr[i].stf3;
      nt6 = (apval_ptr[i].nt3)/6;
      for(k=0;k<nt6;k++)
         {
         for(j=0;j<6;j++)
            {
            it = 6*k + j;
            fprintf(fpw,"%13.5e",stf[it]);
            }
         fprintf(fpw,"\n");
         }

      if(6*nt6 != (apval_ptr[i].nt3))
         {
         for(j=6*nt6;j<(apval_ptr[i].nt3);j++)
            fprintf(fpw,"%13.5e",stf[j]);

         fprintf(fpw,"\n");
         }
      }
   fclose(fpw);
   }
}

void write_srf2(struct standrupformat *srf,char *file,int bflag)
{
FILE *fpw, *fopfile();
struct srf_planerectangle *prect_ptr;
struct srf_prectsegments *prseg_ptr;
struct srf_allpoints *apnts_ptr;
struct srf_apointvalues *apval_ptr;

float area;
float *stf;
int i, j, k, nt6, it, ig;

char *sptr, pword[32];
int fdw;

int ip, nprite;

prect_ptr = &(srf->srf_prect);
prseg_ptr = prect_ptr->prectseg;
apnts_ptr = &(srf->srf_apnts);
apval_ptr = apnts_ptr->apntvals;

if(bflag)
   {
   if(strcmp(file,"stdout") == 0)
      fdw = STDOUT_FILENO;
   else
      fdw = croptrfile(file);

   rite(fdw,srf->version,sizeof(srf->version));

   if(strcmp(srf->type,"PLANE") == 0)
      {
      rite(fdw,srf->type,sizeof(srf->type));
      rite(fdw,&(prect_ptr->nseg),sizeof(prect_ptr->nseg));
      rite(fdw,prseg_ptr,(prect_ptr->nseg)*sizeof(struct srf_prectsegments));
      }

   sprintf(pword,"POINTS");
   rite(fdw,pword,sizeof(pword));
   rite(fdw,&(apnts_ptr->np),sizeof(apnts_ptr->np));
   for(i=0;i<apnts_ptr->np;i++)
      {
      rite(fdw,&(apval_ptr[i].lon),sizeof(float));
      rite(fdw,&(apval_ptr[i].lat),sizeof(float));
      rite(fdw,&(apval_ptr[i].dep),sizeof(float));
      rite(fdw,&(apval_ptr[i].stk),sizeof(float));
      rite(fdw,&(apval_ptr[i].dip),sizeof(float));
      rite(fdw,&(apval_ptr[i].area),sizeof(float));
      rite(fdw,&(apval_ptr[i].tinit),sizeof(float));
      rite(fdw,&(apval_ptr[i].dt),sizeof(float));
      rite(fdw,&(apval_ptr[i].rake),sizeof(float));
      rite(fdw,&(apval_ptr[i].slip1),sizeof(float));
      rite(fdw,&(apval_ptr[i].nt1),sizeof(int));
      rite(fdw,&(apval_ptr[i].slip2),sizeof(float));
      rite(fdw,&(apval_ptr[i].nt2),sizeof(int));
      rite(fdw,&(apval_ptr[i].slip3),sizeof(float));
      rite(fdw,&(apval_ptr[i].nt3),sizeof(int));

      rite(fdw,apval_ptr[i].stf1,(apval_ptr[i].nt1)*sizeof(float));
      rite(fdw,apval_ptr[i].stf2,(apval_ptr[i].nt2)*sizeof(float));
      rite(fdw,apval_ptr[i].stf3,(apval_ptr[i].nt3)*sizeof(float));
      }
   close(fdw);
   }
else
   {
   if(strcmp(file,"stdout") == 0)
      fpw = stdout;
   else
      fpw = fopfile(file,"w");

   fprintf(fpw,"%s\n",srf->version);

   for(i=0;i<srf->srf_hcmnt.nline;i++)
      {
      sptr = (srf->srf_hcmnt.cbuf) + i*MAXLINE;
      if(sptr[0] != '\0')
         fprintf(fpw,"%s\n",sptr);
      }

   if(strcmp(srf->type,"PLANE") == 0)
      {
      fprintf(fpw,"%s %d\n",srf->type,prect_ptr->nseg);
      for(ig=0;ig<prect_ptr->nseg;ig++)
         {
         fprintf(fpw,"%12.6f %11.6f %5d %5d %10.4f %10.4f\n",prseg_ptr[ig].elon,
                                                        prseg_ptr[ig].elat,
                                                        prseg_ptr[ig].nstk,
                                                        prseg_ptr[ig].ndip,
                                                        prseg_ptr[ig].flen,
                                                        prseg_ptr[ig].fwid);
         fprintf(fpw,"%4.0f %4.0f %10.4f %10.4f %10.4f\n",prseg_ptr[ig].stk,
                                                    prseg_ptr[ig].dip,
                                                    prseg_ptr[ig].dtop,
                                                    prseg_ptr[ig].shyp,
                                                    prseg_ptr[ig].dhyp);
         }
      }

   nprite = 0;
   for(ig=0;ig<srf->nseg;ig++)
      {
      fprintf(fpw,"POINTS %d\n",srf->np_seg[ig]);
      for(ip=0;ip<srf->np_seg[ig];ip++)
         {
         i = ip + nprite;

         fprintf(fpw,"%12.6f %11.6f %12.5e %4.0f %4.0f %12.5e %13.6e %12.5e %13.5e %13.5e\n",
                                              apval_ptr[i].lon,
                                              apval_ptr[i].lat,
                                              apval_ptr[i].dep,
                                              apval_ptr[i].stk,
                                              apval_ptr[i].dip,
                                              apval_ptr[i].area,
                                              apval_ptr[i].tinit,
                                              apval_ptr[i].dt,
                                              apval_ptr[i].vs,
                                              apval_ptr[i].den);

         fprintf(fpw,"%4.0f %10.4f %6d %10.4f %6d %10.4f %6d\n",
                                              apval_ptr[i].rake,
                                              apval_ptr[i].slip1,
                                              apval_ptr[i].nt1,
                                              apval_ptr[i].slip2,
                                              apval_ptr[i].nt2,
                                              apval_ptr[i].slip3,
                                              apval_ptr[i].nt3);

         stf = apval_ptr[i].stf1;
         nt6 = (apval_ptr[i].nt1)/6;
         for(k=0;k<nt6;k++)
            {
            for(j=0;j<6;j++)
               {
               it = 6*k + j;
               fprintf(fpw,"%13.5e",stf[it]);
               }
            fprintf(fpw,"\n");
            }

         if(6*nt6 != (apval_ptr[i].nt1))
            {
            for(j=6*nt6;j<(apval_ptr[i].nt1);j++)
               fprintf(fpw,"%13.5e",stf[j]);

            fprintf(fpw,"\n");
            }

         stf = apval_ptr[i].stf2;
         nt6 = (apval_ptr[i].nt2)/6;
         for(k=0;k<nt6;k++)
            {
            for(j=0;j<6;j++)
               {
               it = 6*k + j;
               fprintf(fpw,"%13.5e",stf[it]);
               }
            fprintf(fpw,"\n");
            }

         if(6*nt6 != (apval_ptr[i].nt2))
            {
            for(j=6*nt6;j<(apval_ptr[i].nt2);j++)
               fprintf(fpw,"%13.5e",stf[j]);

            fprintf(fpw,"\n");
            }

         stf = apval_ptr[i].stf3;
         nt6 = (apval_ptr[i].nt3)/6;
         for(k=0;k<nt6;k++)
            {
            for(j=0;j<6;j++)
               {
               it = 6*k + j;
               fprintf(fpw,"%13.5e",stf[it]);
               }
            fprintf(fpw,"\n");
            }

         if(6*nt6 != (apval_ptr[i].nt3))
            {
            for(j=6*nt6;j<(apval_ptr[i].nt3);j++)
               fprintf(fpw,"%13.5e",stf[j]);

            fprintf(fpw,"\n");
            }
         }
      nprite = nprite + srf->np_seg[ig];
      }

   fclose(fpw);
   }
}

void write_srf3(struct standrupformat *srf,char *file,int bflag)
{
FILE *fpw, *fopfile();
struct srf_planerectangle *prect_ptr;
struct srf_prectsegments *prseg_ptr;
struct srf_allpoints *apnts_ptr;
struct srf_apointvalues *apval_ptr;

float area;
float *stf;
int i, j, k, nt6, it, ig;

char *sptr, pword[32];
int fdw;

int ip, nprite;

int mrf1_flag, mrf6_flag;

prect_ptr = &(srf->srf_prect);
prseg_ptr = prect_ptr->prectseg;
apnts_ptr = &(srf->srf_apnts);
apval_ptr = apnts_ptr->apntvals;

if(bflag)
   {
   if(strcmp(file,"stdout") == 0)
      fdw = STDOUT_FILENO;
   else
      fdw = croptrfile(file);

   rite(fdw,srf->version,sizeof(srf->version));

   if(strcmp(srf->type,"PLANE") == 0)
      {
      rite(fdw,srf->type,sizeof(srf->type));
      rite(fdw,&(prect_ptr->nseg),sizeof(prect_ptr->nseg));
      rite(fdw,prseg_ptr,(prect_ptr->nseg)*sizeof(struct srf_prectsegments));
      }

   sprintf(pword,"POINTS");
   rite(fdw,pword,sizeof(pword));
   rite(fdw,&(apnts_ptr->np),sizeof(apnts_ptr->np));
   for(i=0;i<apnts_ptr->np;i++)
      {
      rite(fdw,&(apval_ptr[i].lon),sizeof(float));
      rite(fdw,&(apval_ptr[i].lat),sizeof(float));
      rite(fdw,&(apval_ptr[i].dep),sizeof(float));
      rite(fdw,&(apval_ptr[i].stk),sizeof(float));
      rite(fdw,&(apval_ptr[i].dip),sizeof(float));
      rite(fdw,&(apval_ptr[i].area),sizeof(float));
      rite(fdw,&(apval_ptr[i].tinit),sizeof(float));
      rite(fdw,&(apval_ptr[i].dt),sizeof(float));
      rite(fdw,&(apval_ptr[i].rake),sizeof(float));
      rite(fdw,&(apval_ptr[i].slip1),sizeof(float));
      rite(fdw,&(apval_ptr[i].nt1),sizeof(int));
      rite(fdw,&(apval_ptr[i].slip2),sizeof(float));
      rite(fdw,&(apval_ptr[i].nt2),sizeof(int));
      rite(fdw,&(apval_ptr[i].slip3),sizeof(float));
      rite(fdw,&(apval_ptr[i].nt3),sizeof(int));

      rite(fdw,apval_ptr[i].stf1,(apval_ptr[i].nt1)*sizeof(float));
      rite(fdw,apval_ptr[i].stf2,(apval_ptr[i].nt2)*sizeof(float));
      rite(fdw,apval_ptr[i].stf3,(apval_ptr[i].nt3)*sizeof(float));
      }
   close(fdw);
   }
else
   {
   if(strcmp(file,"stdout") == 0)
      fpw = stdout;
   else
      fpw = fopfile(file,"w");

   fprintf(fpw,"%s %s\n",srf->version,srf->src_format);

   for(i=0;i<srf->srf_hcmnt.nline;i++)
      {
      sptr = (srf->srf_hcmnt.cbuf) + i*MAXLINE;
      if(sptr[0] != '\0')
         fprintf(fpw,"%s\n",sptr);
      }

   if(strcmp(srf->type,"PLANE") == 0)
      {
      fprintf(fpw,"%s %d\n",srf->type,prect_ptr->nseg);
      for(ig=0;ig<prect_ptr->nseg;ig++)
         {
         fprintf(fpw,"%12.6f %11.6f %5d %5d %10.4f %10.4f\n",prseg_ptr[ig].elon,
                                                        prseg_ptr[ig].elat,
                                                        prseg_ptr[ig].nstk,
                                                        prseg_ptr[ig].ndip,
                                                        prseg_ptr[ig].flen,
                                                        prseg_ptr[ig].fwid);
         fprintf(fpw,"%4.0f %4.0f %10.4f %10.4f %10.4f\n",prseg_ptr[ig].stk,
                                                    prseg_ptr[ig].dip,
                                                    prseg_ptr[ig].dtop,
                                                    prseg_ptr[ig].shyp,
                                                    prseg_ptr[ig].dhyp);
         }
      }

   mrf1_flag = 0;
   if(strcmp(srf->src_format,"MOMENT") == 0 || strcmp(srf->src_format,"MOMENT-1MECH") == 0)
      mrf1_flag = 1;

   mrf6_flag = 0;
   if(strcmp(srf->src_format,"MOMENT-6MECH") == 0)
      mrf6_flag = 1;

   nprite = 0;
   for(ig=0;ig<srf->nseg;ig++)
      {
      fprintf(fpw,"POINTS %d\n",srf->np_seg[ig]);
      for(ip=0;ip<srf->np_seg[ig];ip++)
         {
         i = ip + nprite;

         fprintf(fpw,"%12.6f %11.6f %12.5e %4.0f %4.0f %12.5e %13.6e %12.5e %13.5e %13.5e %13.5e\n",
                                              apval_ptr[i].lon,
                                              apval_ptr[i].lat,
                                              apval_ptr[i].dep,
                                              apval_ptr[i].stk,
                                              apval_ptr[i].dip,
                                              apval_ptr[i].area,
                                              apval_ptr[i].tinit,
                                              apval_ptr[i].dt,
                                              apval_ptr[i].vp,
                                              apval_ptr[i].vs,
                                              apval_ptr[i].den);

         if(mrf1_flag == 0 && mrf6_flag == 0)		/* expecting slip */
            {
            fprintf(fpw,"%4.0f %10.4f %6d %10.4f %6d %10.4f %6d\n",
                                                 apval_ptr[i].rake,
                                                 apval_ptr[i].slip1,
                                                 apval_ptr[i].nt1,
                                                 apval_ptr[i].slip2,
                                                 apval_ptr[i].nt2,
                                                 apval_ptr[i].slip3,
                                                 apval_ptr[i].nt3);

            stf = apval_ptr[i].stf1;
            nt6 = (apval_ptr[i].nt1)/6;
            for(k=0;k<nt6;k++)
               {
               for(j=0;j<6;j++)
                  {
                  it = 6*k + j;
                  fprintf(fpw,"%13.5e",stf[it]);
                  }
               fprintf(fpw,"\n");
               }

            if(6*nt6 != (apval_ptr[i].nt1))
               {
               for(j=6*nt6;j<(apval_ptr[i].nt1);j++)
                  fprintf(fpw,"%13.5e",stf[j]);

               fprintf(fpw,"\n");
               }

            stf = apval_ptr[i].stf2;
            nt6 = (apval_ptr[i].nt2)/6;
            for(k=0;k<nt6;k++)
               {
               for(j=0;j<6;j++)
                  {
                  it = 6*k + j;
                  fprintf(fpw,"%13.5e",stf[it]);
                  }
               fprintf(fpw,"\n");
               }

            if(6*nt6 != (apval_ptr[i].nt2))
               {
               for(j=6*nt6;j<(apval_ptr[i].nt2);j++)
                  fprintf(fpw,"%13.5e",stf[j]);

               fprintf(fpw,"\n");
               }

            stf = apval_ptr[i].stf3;
            nt6 = (apval_ptr[i].nt3)/6;
            for(k=0;k<nt6;k++)
               {
               for(j=0;j<6;j++)
                  {
                  it = 6*k + j;
                  fprintf(fpw,"%13.5e",stf[it]);
                  }
               fprintf(fpw,"\n");
               }

            if(6*nt6 != (apval_ptr[i].nt3))
               {
               for(j=6*nt6;j<(apval_ptr[i].nt3);j++)
                  fprintf(fpw,"%13.5e",stf[j]);

               fprintf(fpw,"\n");
               }
            }

         else		/* expecting moment */
            {
            fprintf(fpw,"%13.5e %13.5e %13.5e %13.5e %13.5e %13.5e %6d\n",
                                                 apval_ptr[i].mnn,
                                                 apval_ptr[i].mee,
                                                 apval_ptr[i].mdd,
                                                 apval_ptr[i].mne,
                                                 apval_ptr[i].mnd,
                                                 apval_ptr[i].med,
                                                 apval_ptr[i].ntmr);

            nt6 = (apval_ptr[i].ntmr)/6;

            if(mrf1_flag == 1)		/* single time function */
	       {
               stf = apval_ptr[i].mrf;
               for(k=0;k<nt6;k++)
                  {
                  for(j=0;j<6;j++)
                     {
                     it = 6*k + j;
                     fprintf(fpw,"%13.5e",stf[it]);
                     }
                  fprintf(fpw,"\n");
                  }

               if(6*nt6 != (apval_ptr[i].ntmr))
                  {
                  for(j=6*nt6;j<(apval_ptr[i].ntmr);j++)
                     fprintf(fpw,"%13.5e",stf[j]);

                  fprintf(fpw,"\n");
                  }
               }

            else if(mrf6_flag == 1)		/* 6 time functions */
               {
               stf = apval_ptr[i].mr_nn;
               for(k=0;k<nt6;k++)
                  {
                  for(j=0;j<6;j++)
                     {
                     it = 6*k + j;
                     fprintf(fpw,"%13.5e",stf[it]);
                     }
                  fprintf(fpw,"\n");
                  }

               if(6*nt6 != (apval_ptr[i].ntmr))
                  {
                  for(j=6*nt6;j<(apval_ptr[i].ntmr);j++)
                     fprintf(fpw,"%13.5e",stf[j]);

                  fprintf(fpw,"\n");
                  }

               stf = apval_ptr[i].mr_ee;
               for(k=0;k<nt6;k++)
                  {
                  for(j=0;j<6;j++)
                     {
                     it = 6*k + j;
                     fprintf(fpw,"%13.5e",stf[it]);
                     }
                  fprintf(fpw,"\n");
                  }

               if(6*nt6 != (apval_ptr[i].ntmr))
                  {
                  for(j=6*nt6;j<(apval_ptr[i].ntmr);j++)
                     fprintf(fpw,"%13.5e",stf[j]);

                  fprintf(fpw,"\n");
                  }

               stf = apval_ptr[i].mr_dd;
               for(k=0;k<nt6;k++)
                  {
                  for(j=0;j<6;j++)
                     {
                     it = 6*k + j;
                     fprintf(fpw,"%13.5e",stf[it]);
                     }
                  fprintf(fpw,"\n");
                  }

               if(6*nt6 != (apval_ptr[i].ntmr))
                  {
                  for(j=6*nt6;j<(apval_ptr[i].ntmr);j++)
                     fprintf(fpw,"%13.5e",stf[j]);

                  fprintf(fpw,"\n");
                  }

               stf = apval_ptr[i].mr_ne;
               for(k=0;k<nt6;k++)
                  {
                  for(j=0;j<6;j++)
                     {
                     it = 6*k + j;
                     fprintf(fpw,"%13.5e",stf[it]);
                     }
                  fprintf(fpw,"\n");
                  }

               if(6*nt6 != (apval_ptr[i].ntmr))
                  {
                  for(j=6*nt6;j<(apval_ptr[i].ntmr);j++)
                     fprintf(fpw,"%13.5e",stf[j]);

                  fprintf(fpw,"\n");
                  }

               stf = apval_ptr[i].mr_nd;
               for(k=0;k<nt6;k++)
                  {
                  for(j=0;j<6;j++)
                     {
                     it = 6*k + j;
                     fprintf(fpw,"%13.5e",stf[it]);
                     }
                  fprintf(fpw,"\n");
                  }

               if(6*nt6 != (apval_ptr[i].ntmr))
                  {
                  for(j=6*nt6;j<(apval_ptr[i].ntmr);j++)
                     fprintf(fpw,"%13.5e",stf[j]);

                  fprintf(fpw,"\n");
                  }

               stf = apval_ptr[i].mr_ed;
               for(k=0;k<nt6;k++)
                  {
                  for(j=0;j<6;j++)
                     {
                     it = 6*k + j;
                     fprintf(fpw,"%13.5e",stf[it]);
                     }
                  fprintf(fpw,"\n");
                  }

               if(6*nt6 != (apval_ptr[i].ntmr))
                  {
                  for(j=6*nt6;j<(apval_ptr[i].ntmr);j++)
                     fprintf(fpw,"%13.5e",stf[j]);

                  fprintf(fpw,"\n");
                  }
               }
            }

         }
      nprite = nprite + srf->np_seg[ig];
      }

   fclose(fpw);
   }
}

void free_srf_stf(struct standrupformat *srf)
{
struct srf_allpoints *apnts_ptr;
struct srf_apointvalues *apval_ptr;
int i;

apnts_ptr = &(srf->srf_apnts);
apval_ptr = apnts_ptr->apntvals;

for(i=0;i<apnts_ptr->np;i++)
   {
   free(apval_ptr[i].stf1);
   free(apval_ptr[i].stf2);
   free(apval_ptr[i].stf3);
   }
}

void sum_srf(struct standrupformat *srf0,struct standrupformat *srf1,struct standrupformat *srf2,float *new_rake)
{
struct srf_prectsegments *prseg_ptr0, *prseg_ptr2;
struct srf_apointvalues *apval_ptr[3];
float *stf, *stf1, *stf2, *stf3;
float tdel, rdif;
double cosR, sinR;
int i, j, k, it, ip, ig, newnt, itdel, ir;

double rperd = 0.017453293;

if(srf0[0].srf_apnts.np != srf1[0].srf_apnts.np)
   {
   fprintf(stderr,"*** number of points in srf1 (%d) not equal to number of points in srf2 (%d), exiting...\n",srf0[0].srf_apnts.np,srf1[0].srf_apnts.np);
   exit(-1);
   }

/* 1st, copy all header info from srf0 to srf2 */

strcpy(srf2[0].version,srf0[0].version);

srf2[0].type[0] = '\0';
if(strncmp(srf0[0].type,"PLANE",5) == 0)
   {
   strcpy(srf2[0].type,srf0[0].type);
   srf2[0].srf_prect.nseg = srf0[0].srf_prect.nseg;

   srf2[0].srf_prect.prectseg = (struct srf_prectsegments *)check_malloc(srf2[0].srf_prect.nseg*sizeof(struct srf_prectsegments));

   prseg_ptr0 = srf0[0].srf_prect.prectseg;
   prseg_ptr2 = srf2[0].srf_prect.prectseg;

   for(ig=0;ig<srf2[0].srf_prect.nseg;ig++)
      {
      prseg_ptr2[ig].elon = prseg_ptr0[ig].elon;
      prseg_ptr2[ig].elat = prseg_ptr0[ig].elat;
      prseg_ptr2[ig].nstk = prseg_ptr0[ig].nstk;
      prseg_ptr2[ig].ndip = prseg_ptr0[ig].ndip;
      prseg_ptr2[ig].flen = prseg_ptr0[ig].flen;
      prseg_ptr2[ig].fwid = prseg_ptr0[ig].fwid;
      prseg_ptr2[ig].stk = prseg_ptr0[ig].stk;
      prseg_ptr2[ig].dip = prseg_ptr0[ig].dip;
      prseg_ptr2[ig].dtop = prseg_ptr0[ig].dtop;
      prseg_ptr2[ig].shyp = prseg_ptr0[ig].shyp;
      prseg_ptr2[ig].dhyp = prseg_ptr0[ig].dhyp;
      }
   }

srf2[0].srf_apnts.np = srf0[0].srf_apnts.np;
srf2[0].srf_apnts.apntvals = (struct srf_apointvalues *)check_malloc((srf2[0].srf_apnts.np)*sizeof(struct srf_apointvalues));

for(i=0;i<srf2[0].srf_apnts.np;i++)
   {
   apval_ptr[0] = &(srf0[0].srf_apnts.apntvals[i]);
   apval_ptr[1] = &(srf1[0].srf_apnts.apntvals[i]);
   apval_ptr[2] = &(srf2[0].srf_apnts.apntvals[i]);

   /* these should be all the same for srf0 and srf1 */
   apval_ptr[2]->lon = apval_ptr[0]->lon;
   apval_ptr[2]->lat = apval_ptr[0]->lat;
   apval_ptr[2]->dep = apval_ptr[0]->dep;
   apval_ptr[2]->stk = apval_ptr[0]->stk;
   apval_ptr[2]->dip = apval_ptr[0]->dip;
   apval_ptr[2]->area = apval_ptr[0]->area;
   apval_ptr[2]->dt = apval_ptr[0]->dt;
   
   /* find earliest initiation time, reset pointers so ptr0 is earliest and ptr2 is latest */

   if(apval_ptr[1]->tinit < apval_ptr[0]->tinit) /* need to reset pointers */
      {
      apval_ptr[0] = &(srf1[0].srf_apnts.apntvals[i]);
      apval_ptr[1] = &(srf0[0].srf_apnts.apntvals[i]);
      }
   apval_ptr[2]->tinit = apval_ptr[0]->tinit;

   /* resolve all slip to rake_u1=new_rake, rake_u2=new_rake+90 */

   apval_ptr[2]->rake = *new_rake;

   apval_ptr[2]->slip1 = 0.0;
   apval_ptr[2]->nt1 = 0;
   apval_ptr[2]->stf1 = NULL;

   apval_ptr[2]->slip2 = 0.0;
   apval_ptr[2]->nt2 = 0;
   apval_ptr[2]->stf2 = NULL;

   apval_ptr[2]->slip3 = 0.0;
   apval_ptr[2]->nt3 = 0;
   apval_ptr[2]->stf3 = NULL;

   /* loop over u1, u2, u3 for srf0 and srf1 to build up stfs */

   for(ir=0;ir<2;ir++)
      {
      rdif = apval_ptr[ir]->rake - apval_ptr[2]->rake;
      cosR = cos(rperd*rdif);
      sinR = sin(rperd*rdif);

      apval_ptr[2]->slip1 = apval_ptr[2]->slip1 +
                            apval_ptr[ir]->slip1*cosR - apval_ptr[ir]->slip2*sinR;

      apval_ptr[2]->slip2 = apval_ptr[2]->slip2 +
                            apval_ptr[ir]->slip1*sinR + apval_ptr[ir]->slip2*cosR;

      apval_ptr[2]->slip3 = apval_ptr[2]->slip3 + apval_ptr[ir]->slip3;

      tdel = apval_ptr[ir]->tinit - apval_ptr[2]->tinit;
      itdel = (int)(tdel/(apval_ptr[2]->dt) + 0.5);

      stf1 = apval_ptr[2]->stf1;
      stf2 = apval_ptr[2]->stf2;
      stf3 = apval_ptr[2]->stf3;

      if(apval_ptr[ir]->nt1)
	 {
	 newnt = itdel + apval_ptr[ir]->nt1;

	 if(newnt > apval_ptr[2]->nt1)
	    {
            apval_ptr[2]->stf1 = (float *)check_realloc(apval_ptr[2]->stf1,newnt*sizeof(float));
            stf1 = apval_ptr[2]->stf1;

            for(it=(apval_ptr[2]->nt1);it<newnt;it++)
	       stf1[it] = 0.0;

	    apval_ptr[2]->nt1 = newnt;
	    }

	 if(newnt > apval_ptr[2]->nt2)
	    {
            apval_ptr[2]->stf2 = (float *)check_realloc(apval_ptr[2]->stf2,newnt*sizeof(float));
            stf2 = apval_ptr[2]->stf2;

            for(it=(apval_ptr[2]->nt2);it<newnt;it++)
	       stf2[it] = 0.0;

	    apval_ptr[2]->nt2 = newnt;
	    }

         stf = apval_ptr[ir]->stf1;
         for(it=itdel;it<newnt;it++)
	    {
            stf1[it] = stf1[it] + stf[it-itdel]*cosR;
            stf2[it] = stf2[it] + stf[it-itdel]*sinR;
	    }
	 }

      if(apval_ptr[ir]->nt2)
	 {
	 newnt = itdel + apval_ptr[ir]->nt2;

	 if(newnt > apval_ptr[2]->nt1)
	    {
            apval_ptr[2]->stf1 = (float *)check_realloc(apval_ptr[2]->stf1,newnt*sizeof(float));
            stf1 = apval_ptr[2]->stf1;

            for(it=(apval_ptr[2]->nt1);it<newnt;it++)
	       stf1[it] = 0.0;

	    apval_ptr[2]->nt1 = newnt;
	    }

	 if(newnt > apval_ptr[2]->nt2)
	    {
            apval_ptr[2]->stf2 = (float *)check_realloc(apval_ptr[2]->stf2,newnt*sizeof(float));
            stf2 = apval_ptr[2]->stf2;

            for(it=(apval_ptr[2]->nt2);it<newnt;it++)
	       stf2[it] = 0.0;

	    apval_ptr[2]->nt2 = newnt;
	    }

         stf = apval_ptr[ir]->stf2;
         for(it=itdel;it<newnt;it++)
	    {
            stf1[it] = stf1[it] - stf[it-itdel]*sinR;
            stf2[it] = stf2[it] + stf[it-itdel]*cosR;
	    }
	 }

      if(apval_ptr[ir]->nt3)
	 {
	 newnt = itdel + apval_ptr[ir]->nt3;

	 if(newnt > apval_ptr[2]->nt3)
	    {
            apval_ptr[2]->stf3 = (float *)check_realloc(apval_ptr[2]->stf3,newnt*sizeof(float));
            stf3 = apval_ptr[2]->stf3;

            for(it=(apval_ptr[2]->nt3);it<newnt;it++)
	       stf3[it] = 0.0;

	    apval_ptr[2]->nt3 = newnt;
	    }

         stf = apval_ptr[ir]->stf3;
         for(it=itdel;it<newnt;it++)
            stf3[it] = stf3[it] + stf[it-itdel];
	 }
      }

   if(apval_ptr[2]->nt1 == 0)
      apval_ptr[2]->slip1 = 0.0;
   if(apval_ptr[2]->nt2 == 0)
      apval_ptr[2]->slip2 = 0.0;
   if(apval_ptr[2]->nt3 == 0)
      apval_ptr[2]->slip3 = 0.0;
   }
}

void join_srf(struct standrupformat *srf0,struct standrupformat *srf1,struct standrupformat *srf2)
{
struct standrupformat *srfp_in;
struct srf_prectsegments *prseg_in, *prseg_out;
struct srf_apointvalues *apval_in, *apval_out;
char *sptr_in, *sptr_out;
float *stfin, *stfout;
int i, j, k, it, ip, ig;
int kg, npoff_in, npoff_out;

if(strcmp(srf0[0].version,srf1[0].version) != 0)
   {
   fprintf(stderr,"srf1.version= %s != srf2.version= %s, exiting ... \n",srf0[0].version,srf1[0].version);
   exit(-1);
   }

strcpy(srf2[0].version,srf0[0].version);

if(atof(srf2->version) < 2.0)
   {
   srf2[0].type[0] = '\0';
   if(strncmp(srf0[0].type,"PLANE",5) == 0 && strncmp(srf1[0].type,"PLANE",5) == 0)
      {
      strcpy(srf2[0].type,srf0[0].type);

      srf2[0].srf_prect.nseg = srf0[0].srf_prect.nseg + srf1[0].srf_prect.nseg;
      srf2[0].srf_prect.prectseg = (struct srf_prectsegments *)check_malloc(srf2[0].srf_prect.nseg*sizeof(struct srf_prectsegments));

      prseg_out = srf2[0].srf_prect.prectseg;
      for(ig=0;ig<srf2[0].srf_prect.nseg;ig++)
         {
         if(ig < srf0[0].srf_prect.nseg)
            {
            k = ig;
            prseg_in = srf0[0].srf_prect.prectseg;
	    }
         else
            {
            k = ig - srf0[0].srf_prect.nseg;
            prseg_in = srf1[0].srf_prect.prectseg;
	    }

         prseg_out[ig].elon = prseg_in[k].elon;
         prseg_out[ig].elat = prseg_in[k].elat;
         prseg_out[ig].nstk = prseg_in[k].nstk;
         prseg_out[ig].ndip = prseg_in[k].ndip;
         prseg_out[ig].flen = prseg_in[k].flen;
         prseg_out[ig].fwid = prseg_in[k].fwid;
         prseg_out[ig].stk = prseg_in[k].stk;
         prseg_out[ig].dip = prseg_in[k].dip;
         prseg_out[ig].dtop = prseg_in[k].dtop;
         prseg_out[ig].shyp = prseg_in[k].shyp;
         prseg_out[ig].dhyp = prseg_in[k].dhyp;
         }
      }

   srf2[0].srf_apnts.np = srf0[0].srf_apnts.np + srf1[0].srf_apnts.np;
   srf2[0].srf_apnts.apntvals = (struct srf_apointvalues *)check_malloc((srf2[0].srf_apnts.np)*sizeof(struct srf_apointvalues));

   for(i=0;i<srf2[0].srf_apnts.np;i++)
      {
      if(i < srf0[0].srf_apnts.np)
         {
         k = i;
         apval_in = &(srf0[0].srf_apnts.apntvals[k]);
         }
      else
         {
         k = i - srf0[0].srf_apnts.np;
         apval_in = &(srf1[0].srf_apnts.apntvals[k]);
         }

      apval_out = &(srf2[0].srf_apnts.apntvals[i]);

      apval_out->lon = apval_in->lon;
      apval_out->lat = apval_in->lat;
      apval_out->dep = apval_in->dep;
      apval_out->stk = apval_in->stk;
      apval_out->dip = apval_in->dip;
      apval_out->area = apval_in->area;
      apval_out->tinit = apval_in->tinit;
      apval_out->dt = apval_in->dt;
      apval_out->rake = apval_in->rake;

      apval_out->slip1 = apval_in->slip1;
      apval_out->nt1 = apval_in->nt1;
      apval_out->stf1 = NULL;

      if(apval_out->nt1)
         {
         apval_out->stf1 = (float *)check_realloc(apval_out->stf1,(apval_out->nt1)*sizeof(float));

         stfin = apval_in->stf1;
         stfout = apval_out->stf1;

         for(it=0;it<(apval_out->nt1);it++)
            stfout[it] = stfin[it];
         }

      apval_out->slip2 = apval_in->slip2;
      apval_out->nt2 = apval_in->nt2;
      apval_out->stf2 = NULL;

      if(apval_out->nt2)
         {
         apval_out->stf2 = (float *)check_realloc(apval_out->stf2,(apval_out->nt2)*sizeof(float));

         stfin = apval_in->stf2;
         stfout = apval_out->stf2;

         for(it=0;it<(apval_out->nt2);it++)
            stfout[it] = stfin[it];
         }

      apval_out->slip3 = apval_in->slip3;
      apval_out->nt3 = apval_in->nt3;
      apval_out->stf3 = NULL;

      if(apval_out->nt3)
         {
         apval_out->stf3 = (float *)check_realloc(apval_out->stf3,(apval_out->nt3)*sizeof(float));

         stfin = apval_in->stf3;
         stfout = apval_out->stf3;

         for(it=0;it<(apval_out->nt3);it++)
            stfout[it] = stfin[it];
         }
      }
   }

else if(atof(srf2->version) >= 2.0)
   {
   srf2[0].srf_hcmnt.nline = srf0[0].srf_hcmnt.nline + srf1[0].srf_hcmnt.nline;
   srf2[0].srf_hcmnt.cbuf = (char *)check_malloc((srf2[0].srf_hcmnt.nline)*MAXLINE*sizeof(char));

   for(i=0;i<srf2[0].srf_hcmnt.nline;i++)
      {
      if(i < srf0[0].srf_hcmnt.nline)
         {
         k = i;
         sptr_in = srf0[0].srf_hcmnt.cbuf + k*MAXLINE;
	 }
      else
         {
         k = i - srf0[0].srf_hcmnt.nline;
         sptr_in = srf1[0].srf_hcmnt.cbuf + k*MAXLINE;
	 }

      sptr_out = srf2[0].srf_hcmnt.cbuf + i*MAXLINE;
      strcpy(sptr_out,sptr_in);
      }

   srf2[0].type[0] = '\0';
   if(strncmp(srf0[0].type,"PLANE",5) == 0 && strncmp(srf1[0].type,"PLANE",5) == 0)
      {
      strcpy(srf2[0].type,srf0[0].type);

      srf2[0].srf_prect.nseg = srf0[0].srf_prect.nseg + srf1[0].srf_prect.nseg;
      srf2[0].srf_prect.prectseg = (struct srf_prectsegments *)check_malloc(srf2[0].srf_prect.nseg*sizeof(struct srf_prectsegments));

      prseg_out = srf2[0].srf_prect.prectseg;
      for(ig=0;ig<srf2[0].srf_prect.nseg;ig++)
         {
         if(ig < srf0[0].srf_prect.nseg)
            {
            k = ig;
            prseg_in = srf0[0].srf_prect.prectseg;
	    }
         else
            {
            k = ig - srf0[0].srf_prect.nseg;
            prseg_in = srf1[0].srf_prect.prectseg;
	    }

         prseg_out[ig].elon = prseg_in[k].elon;
         prseg_out[ig].elat = prseg_in[k].elat;
         prseg_out[ig].nstk = prseg_in[k].nstk;
         prseg_out[ig].ndip = prseg_in[k].ndip;
         prseg_out[ig].flen = prseg_in[k].flen;
         prseg_out[ig].fwid = prseg_in[k].fwid;
         prseg_out[ig].stk = prseg_in[k].stk;
         prseg_out[ig].dip = prseg_in[k].dip;
         prseg_out[ig].dtop = prseg_in[k].dtop;
         prseg_out[ig].shyp = prseg_in[k].shyp;
         prseg_out[ig].dhyp = prseg_in[k].dhyp;
         }
      }

   srf2[0].srf_apnts.np = srf0[0].srf_apnts.np + srf1[0].srf_apnts.np;
   srf2[0].srf_apnts.apntvals = (struct srf_apointvalues *)check_malloc((srf2[0].srf_apnts.np)*sizeof(struct srf_apointvalues));

   srf2[0].nseg = srf0[0].nseg + srf1[0].nseg;
   srf2[0].np_seg = (int *)check_malloc((srf2[0].nseg)*sizeof(int));

   npoff_out = 0;
   for(ig=0;ig<srf2[0].nseg;ig++)
      {
      if(ig < srf0[0].nseg)
         {
	 kg = ig;

	 if(ig == 0)
            npoff_in = 0;
	 else
	    npoff_in = npoff_in + srf0[0].np_seg[kg-1];

         srfp_in = srf0;
	 }
      else
         {
	 kg = ig - srf0[0].nseg;

	 if(kg == 0)
            npoff_in = 0;
	 else
	    npoff_in = npoff_in + srf1[0].np_seg[kg-1];

         srfp_in = srf1;
	 }

      if(ig == 0)
         npoff_in = 0;
      else
         npoff_out = npoff_out + srf2[0].np_seg[ig-1];

      srf2[0].np_seg[ig] = srfp_in[0].np_seg[kg];
      for(i=0;i<srf2[0].np_seg[ig];i++)
         {
         apval_in = &(srfp_in[0].srf_apnts.apntvals[i+npoff_in]);
         apval_out = &(srf2[0].srf_apnts.apntvals[i+npoff_out]);

         apval_out->lon = apval_in->lon;
         apval_out->lat = apval_in->lat;
         apval_out->dep = apval_in->dep;
         apval_out->stk = apval_in->stk;
         apval_out->dip = apval_in->dip;
         apval_out->area = apval_in->area;
         apval_out->tinit = apval_in->tinit;
         apval_out->dt = apval_in->dt;
         apval_out->vs = apval_in->vs;
         apval_out->den = apval_in->den;

         apval_out->rake = apval_in->rake;

         apval_out->slip1 = apval_in->slip1;
         apval_out->nt1 = apval_in->nt1;
         apval_out->stf1 = NULL;

         if(apval_out->nt1)
            {
            apval_out->stf1 = (float *)check_realloc(apval_out->stf1,(apval_out->nt1)*sizeof(float));

            stfin = apval_in->stf1;
            stfout = apval_out->stf1;

            for(it=0;it<(apval_out->nt1);it++)
               stfout[it] = stfin[it];
            }

         apval_out->slip2 = apval_in->slip2;
         apval_out->nt2 = apval_in->nt2;
         apval_out->stf2 = NULL;

         if(apval_out->nt2)
            {
            apval_out->stf2 = (float *)check_realloc(apval_out->stf2,(apval_out->nt2)*sizeof(float));

            stfin = apval_in->stf2;
            stfout = apval_out->stf2;

            for(it=0;it<(apval_out->nt2);it++)
               stfout[it] = stfin[it];
            }

         apval_out->slip3 = apval_in->slip3;
         apval_out->nt3 = apval_in->nt3;
         apval_out->stf3 = NULL;

         if(apval_out->nt3)
            {
            apval_out->stf3 = (float *)check_realloc(apval_out->stf3,(apval_out->nt3)*sizeof(float));

            stfin = apval_in->stf3;
            stfout = apval_out->stf3;

            for(it=0;it<(apval_out->nt3);it++)
               stfout[it] = stfin[it];
            }
         }
      }
   }
}

void select_depths_srf(struct standrupformat *srf0,struct standrupformat *srf1,float *dep0,float *dep1)
{
struct srf_prectsegments *prseg_in, *prseg_out;
struct srf_apointvalues *apval_in, *apval_out;
float *stfin, *stfout;
int i, j, k, it, ip, ig;

/* 1st, copy all header info from srf0 to srf1 */

strcpy(srf1[0].version,srf0[0].version);

srf1[0].type[0] = '\0';
if(strncmp(srf0[0].type,"PLANE",5) == 0)
   {
   strcpy(srf1[0].type,srf0[0].type);

   srf1[0].srf_prect.nseg = srf0[0].srf_prect.nseg;
   srf1[0].srf_prect.prectseg = (struct srf_prectsegments *)check_malloc(srf1[0].srf_prect.nseg*sizeof(struct srf_prectsegments));

   prseg_in = srf0[0].srf_prect.prectseg;
   prseg_out = srf1[0].srf_prect.prectseg;
   for(ig=0;ig<srf1[0].srf_prect.nseg;ig++)
      {
      prseg_out[ig].elon = prseg_in[ig].elon;
      prseg_out[ig].elat = prseg_in[ig].elat;
      prseg_out[ig].nstk = prseg_in[ig].nstk;
      prseg_out[ig].ndip = prseg_in[ig].ndip;
      prseg_out[ig].flen = prseg_in[ig].flen;
      prseg_out[ig].fwid = prseg_in[ig].fwid;
      prseg_out[ig].stk = prseg_in[ig].stk;
      prseg_out[ig].dip = prseg_in[ig].dip;
      prseg_out[ig].dtop = prseg_in[ig].dtop;
      prseg_out[ig].shyp = prseg_in[ig].shyp;
      prseg_out[ig].dhyp = prseg_in[ig].dhyp;
      }
   }

srf1[0].srf_apnts.apntvals = (struct srf_apointvalues *)check_malloc((srf0[0].srf_apnts.np)*sizeof(struct srf_apointvalues));

ip = 0;
for(i=0;i<srf0[0].srf_apnts.np;i++)
   {
   apval_in = &(srf0[0].srf_apnts.apntvals[i]);

   if(apval_in->dep >= *dep0 && apval_in->dep <= *dep1)
      {
      apval_out = &(srf1[0].srf_apnts.apntvals[ip]);

      apval_out->lon = apval_in->lon;
      apval_out->lat = apval_in->lat;
      apval_out->dep   = apval_in->dep;

      apval_out->stk   = apval_in->stk;
      apval_out->dip   = apval_in->dip;
      apval_out->area  = apval_in->area;

      apval_out->tinit = apval_in->tinit;
      apval_out->dt    = apval_in->dt;
      apval_out->rake  = apval_in->rake;

      apval_out->slip1 = apval_in->slip1;
      apval_out->nt1 = apval_in->nt1;
      apval_out->stf1 = NULL;

      if(apval_out->nt1)
         {
         apval_out->stf1 = (float *)check_realloc(apval_out->stf1,(apval_out->nt1)*sizeof(float));

         stfin = apval_in->stf1;
         stfout = apval_out->stf1;

         for(it=0;it<(apval_out->nt1);it++)
            stfout[it] = stfin[it];
         }

      apval_out->slip2 = apval_in->slip2;
      apval_out->nt2 = apval_in->nt2;
      apval_out->stf2 = NULL;

      if(apval_out->nt2)
         {
         apval_out->stf2 = (float *)check_realloc(apval_out->stf2,(apval_out->nt2)*sizeof(float));

         stfin = apval_in->stf2;
         stfout = apval_out->stf2;

         for(it=0;it<(apval_out->nt2);it++)
            stfout[it] = stfin[it];
         }

      apval_out->slip3 = apval_in->slip3;
      apval_out->nt3 = apval_in->nt3;
      apval_out->stf3 = NULL;

      if(apval_out->nt3)
         {
         apval_out->stf3 = (float *)check_realloc(apval_out->stf3,(apval_out->nt3)*sizeof(float));

         stfin = apval_in->stf3;
         stfout = apval_out->stf3;

         for(it=0;it<(apval_out->nt3);it++)
            stfout[it] = stfin[it];
         }

      ip++;
      }
   }
srf1[0].srf_apnts.np = ip;
}

void scale_srf(struct standrupformat *srf1,struct standrupformat *srf2,float *scale)
{
struct srf_prectsegments *prseg_ptr1, *prseg_ptr2;
struct srf_apointvalues *apval_ptr[2];
float *stfp1, *stfp2;
float tdel, rdif;
double cosR, sinR;
int i, j, k, it, ip, ig, newnt, itdel, ir;

double rperd = 0.017453293;

/* 1st, copy all header info from srf1 to srf2 */

strcpy(srf2[0].version,srf1[0].version);

srf2[0].type[0] = '\0';
if(strncmp(srf1[0].type,"PLANE",5) == 0)
   {
   strcpy(srf2[0].type,srf1[0].type);
   srf2[0].srf_prect.nseg = srf1[0].srf_prect.nseg;

   srf2[0].srf_prect.prectseg = (struct srf_prectsegments *)check_malloc(srf2[0].srf_prect.nseg*sizeof(struct srf_prectsegments));

   prseg_ptr1 = srf1[0].srf_prect.prectseg;
   prseg_ptr2 = srf2[0].srf_prect.prectseg;

   for(ig=0;ig<srf2[0].srf_prect.nseg;ig++)
      {
      prseg_ptr2[ig].elon = prseg_ptr1[ig].elon;
      prseg_ptr2[ig].elat = prseg_ptr1[ig].elat;
      prseg_ptr2[ig].nstk = prseg_ptr1[ig].nstk;
      prseg_ptr2[ig].ndip = prseg_ptr1[ig].ndip;
      prseg_ptr2[ig].flen = prseg_ptr1[ig].flen;
      prseg_ptr2[ig].fwid = prseg_ptr1[ig].fwid;
      prseg_ptr2[ig].stk = prseg_ptr1[ig].stk;
      prseg_ptr2[ig].dip = prseg_ptr1[ig].dip;
      prseg_ptr2[ig].dtop = prseg_ptr1[ig].dtop;
      prseg_ptr2[ig].shyp = prseg_ptr1[ig].shyp;
      prseg_ptr2[ig].dhyp = prseg_ptr1[ig].dhyp;
      }
   }

srf2[0].srf_apnts.np = srf1[0].srf_apnts.np;
srf2[0].srf_apnts.apntvals = (struct srf_apointvalues *)check_malloc((srf2[0].srf_apnts.np)*sizeof(struct srf_apointvalues));

for(i=0;i<srf2[0].srf_apnts.np;i++)
   {
   apval_ptr[0] = &(srf1[0].srf_apnts.apntvals[i]);
   apval_ptr[1] = &(srf2[0].srf_apnts.apntvals[i]);

   /* these should be all the same for srf0 and srf1 */
   apval_ptr[1]->lon = apval_ptr[0]->lon;
   apval_ptr[1]->lat = apval_ptr[0]->lat;
   apval_ptr[1]->dep = apval_ptr[0]->dep;
   apval_ptr[1]->stk = apval_ptr[0]->stk;
   apval_ptr[1]->dip = apval_ptr[0]->dip;
   apval_ptr[1]->area = apval_ptr[0]->area;
   apval_ptr[1]->dt = apval_ptr[0]->dt;
   
   apval_ptr[1]->tinit = apval_ptr[0]->tinit;
   apval_ptr[1]->rake = apval_ptr[0]->rake;

   apval_ptr[1]->slip1 = (*scale)*apval_ptr[0]->slip1;
   apval_ptr[1]->nt1 = apval_ptr[0]->nt1;
   apval_ptr[1]->stf1 = NULL;

   if(apval_ptr[1]->nt1)
      {
      apval_ptr[1]->stf1 = (float *)check_realloc(apval_ptr[1]->stf1,apval_ptr[1]->nt1*sizeof(float));
      stfp1 = apval_ptr[0]->stf1;
      stfp2 = apval_ptr[1]->stf1;
      for(it=0;it<apval_ptr[1]->nt1;it++)
         stfp2[it] = (*scale)*stfp1[it];
      }

   apval_ptr[1]->slip2 = (*scale)*apval_ptr[0]->slip2;
   apval_ptr[1]->nt2 = apval_ptr[0]->nt2;
   apval_ptr[1]->stf2 = NULL;

   if(apval_ptr[1]->nt2)
      {
      apval_ptr[1]->stf2 = (float *)check_realloc(apval_ptr[1]->stf2,apval_ptr[1]->nt2*sizeof(float));
      stfp1 = apval_ptr[0]->stf2;
      stfp2 = apval_ptr[1]->stf2;
      for(it=0;it<apval_ptr[1]->nt2;it++)
         stfp2[it] = (*scale)*stfp1[it];
      }

   apval_ptr[1]->slip3 = (*scale)*apval_ptr[0]->slip3;
   apval_ptr[1]->nt3 = apval_ptr[0]->nt3;
   apval_ptr[1]->stf3 = NULL;

   if(apval_ptr[1]->nt3)
      {
      apval_ptr[1]->stf3 = (float *)check_realloc(apval_ptr[1]->stf3,apval_ptr[1]->nt3*sizeof(float));
      stfp1 = apval_ptr[0]->stf3;
      stfp2 = apval_ptr[1]->stf3;
      for(it=0;it<apval_ptr[1]->nt3;it++)
         stfp2[it] = (*scale)*stfp1[it];
      }
   }
}

void copy_hcmnt(struct standrupformat *srf2,struct standrupformat *srf1)
{
int i;
char *sptr_in, *sptr_out;

if(atof(srf2->version) < 2.0)
   {
   srf2[0].srf_hcmnt.nline = 0;
   srf2[0].srf_hcmnt.cbuf = NULL;
   }
else if(atof(srf2->version) >= 2.0)
   {
   srf2[0].srf_hcmnt.nline = srf1[0].srf_hcmnt.nline;
   srf2[0].srf_hcmnt.cbuf = (char *)check_malloc((srf2[0].srf_hcmnt.nline)*MAXLINE*sizeof(char));

   for(i=0;i<srf2[0].srf_hcmnt.nline;i++)
      {
      sptr_in = srf1[0].srf_hcmnt.cbuf + i*MAXLINE;
      sptr_out = srf2[0].srf_hcmnt.cbuf + i*MAXLINE;
      strcpy(sptr_out,sptr_in);
      }
   }
}

void replace_sdr(struct standrupformat *srf0,struct standrupformat *srf1,struct standrupformat *srf2,int tflag)
{
struct standrupformat *srfp_in;
struct srf_prectsegments *prseg_in0, *prseg_in1, *prseg_out;
struct srf_apointvalues *apval_in0, *apval_in1, *apval_out;
char *sptr_in, *sptr_out;
float *stfin, *stfout;
int i, j, k, it, ip, ig;
int npoff;

strcpy(srf2[0].version,srf0[0].version);

   srf2[0].srf_hcmnt.nline = srf0[0].srf_hcmnt.nline;
   srf2[0].srf_hcmnt.cbuf = (char *)check_malloc((srf2[0].srf_hcmnt.nline)*MAXLINE*sizeof(char));

   for(i=0;i<srf2[0].srf_hcmnt.nline;i++)
      {
      sptr_in = srf0[0].srf_hcmnt.cbuf + i*MAXLINE;
      sptr_out = srf2[0].srf_hcmnt.cbuf + i*MAXLINE;
      strcpy(sptr_out,sptr_in);
      }

   srf2[0].type[0] = '\0';
   if(strncmp(srf0[0].type,"PLANE",5) == 0 && strncmp(srf1[0].type,"PLANE",5) == 0)
      {
      strcpy(srf2[0].type,srf0[0].type);

      srf2[0].srf_prect.nseg = srf0[0].srf_prect.nseg;
      srf2[0].srf_prect.prectseg = (struct srf_prectsegments *)check_malloc(srf2[0].srf_prect.nseg*sizeof(struct srf_prectsegments));

      prseg_in0 = srf0[0].srf_prect.prectseg;
      prseg_in1 = srf1[0].srf_prect.prectseg;
      prseg_out = srf2[0].srf_prect.prectseg;
      for(ig=0;ig<srf2[0].srf_prect.nseg;ig++)
         {
         prseg_out[ig].elon = prseg_in0[ig].elon;
         prseg_out[ig].elat = prseg_in0[ig].elat;
         prseg_out[ig].nstk = prseg_in0[ig].nstk;
         prseg_out[ig].ndip = prseg_in0[ig].ndip;
         prseg_out[ig].flen = prseg_in0[ig].flen;
         prseg_out[ig].fwid = prseg_in0[ig].fwid;
         prseg_out[ig].stk = prseg_in1[ig].stk;
         prseg_out[ig].dip = prseg_in1[ig].dip;
         prseg_out[ig].dtop = prseg_in0[ig].dtop;
         prseg_out[ig].shyp = prseg_in0[ig].shyp;
         prseg_out[ig].dhyp = prseg_in0[ig].dhyp;
         }
      }

   srf2[0].srf_apnts.np = srf0[0].srf_apnts.np;
   srf2[0].srf_apnts.apntvals = (struct srf_apointvalues *)check_malloc((srf2[0].srf_apnts.np)*sizeof(struct srf_apointvalues));

   srf2[0].nseg = srf0[0].nseg;
   srf2[0].np_seg = (int *)check_malloc((srf2[0].nseg)*sizeof(int));

   npoff = 0;
   for(ig=0;ig<srf2[0].nseg;ig++)
      {
      srf2[0].np_seg[ig] = srf0[0].np_seg[ig];
      for(i=0;i<srf2[0].np_seg[ig];i++)
         {
         apval_in0 = &(srf0[0].srf_apnts.apntvals[i+npoff]);
         apval_in1 = &(srf1[0].srf_apnts.apntvals[i+npoff]);
         apval_out = &(srf2[0].srf_apnts.apntvals[i+npoff]);

         apval_out->lon = apval_in0->lon;
         apval_out->lat = apval_in0->lat;
         apval_out->dep = apval_in0->dep;
         apval_out->stk = apval_in1->stk;
         apval_out->dip = apval_in1->dip;
         apval_out->area = apval_in0->area;

	 if(tflag)
            apval_out->tinit = apval_in1->tinit;
	 else
            apval_out->tinit = apval_in0->tinit;

         apval_out->dt = apval_in0->dt;
         apval_out->vs = apval_in0->vs;
         apval_out->den = apval_in0->den;

         apval_out->rake = apval_in1->rake;

         apval_out->slip1 = apval_in0->slip1;
         apval_out->nt1 = apval_in0->nt1;
         apval_out->stf1 = NULL;

         if(apval_out->nt1)
            {
            apval_out->stf1 = (float *)check_realloc(apval_out->stf1,(apval_out->nt1)*sizeof(float));

            stfin = apval_in0->stf1;
            stfout = apval_out->stf1;

            for(it=0;it<(apval_out->nt1);it++)
               stfout[it] = stfin[it];
            }

         apval_out->slip2 = apval_in0->slip2;
         apval_out->nt2 = apval_in0->nt2;
         apval_out->stf2 = NULL;

         if(apval_out->nt2)
            {
            apval_out->stf2 = (float *)check_realloc(apval_out->stf2,(apval_out->nt2)*sizeof(float));

            stfin = apval_in0->stf2;
            stfout = apval_out->stf2;

            for(it=0;it<(apval_out->nt2);it++)
               stfout[it] = stfin[it];
            }

         apval_out->slip3 = apval_in0->slip3;
         apval_out->nt3 = apval_in0->nt3;
         apval_out->stf3 = NULL;

         if(apval_out->nt3)
            {
            apval_out->stf3 = (float *)check_realloc(apval_out->stf3,(apval_out->nt3)*sizeof(float));

            stfin = apval_in0->stf3;
            stfout = apval_out->stf3;

            for(it=0;it<(apval_out->nt3);it++)
               stfout[it] = stfin[it];
            }
         }

      npoff = npoff + srf2[0].np_seg[ig];
      }
}

void load_command_srf(struct standrupformat *srf,int ac,char **av)
{
int i, nb;
char *sptr1, *sptr2;

if(srf[0].srf_hcmnt.nline == 0)
   {
   srf[0].srf_hcmnt.cbuf = (char *)check_malloc(2*MAXLINE*sizeof(char));
   srf[0].srf_hcmnt.nline = 2;
   }
else
   {
   srf[0].srf_hcmnt.cbuf = (char *)check_realloc(srf[0].srf_hcmnt.cbuf,(2+srf[0].srf_hcmnt.nline)*MAXLINE*sizeof(char));

   for(i=srf[0].srf_hcmnt.nline-1;i>=0;i--)
      {
      sptr1 = srf[0].srf_hcmnt.cbuf + (i)*MAXLINE;
      sptr2 = srf[0].srf_hcmnt.cbuf + (i+2)*MAXLINE;
      strcpy(sptr2,sptr1);
      }

   srf[0].srf_hcmnt.nline = srf[0].srf_hcmnt.nline + 2;
   }

sptr1 = srf[0].srf_hcmnt.cbuf;
nb = sprintf(sptr1,"# Command:");
sptr1 = sptr1 + nb;
for(i=0;i<ac;i++)
   {
   nb = sprintf(sptr1," %s",av[i]);
   sptr1 = sptr1 + nb;
   }

sptr1 = srf[0].srf_hcmnt.cbuf + MAXLINE;
sprintf(sptr1,"#");
}

void load_seed_srf(struct standrupformat *srf,int starting_seed,int ending_seed)
{
int i, nb;
char *sptr1;

if(srf[0].srf_hcmnt.nline == 0)
   {
   srf[0].srf_hcmnt.cbuf = (char *)check_malloc(2*MAXLINE*sizeof(char));
   srf[0].srf_hcmnt.nline = 3;
   }
else
   {
   srf[0].srf_hcmnt.nline = srf[0].srf_hcmnt.nline + 3;

   srf[0].srf_hcmnt.cbuf = (char *)check_realloc(srf[0].srf_hcmnt.cbuf,(srf[0].srf_hcmnt.nline)*MAXLINE*sizeof(char));
   }

sptr1 = srf[0].srf_hcmnt.cbuf + (srf[0].srf_hcmnt.nline-3)*MAXLINE;
nb = sprintf(sptr1,"# starting_seed= %lld",starting_seed);
sptr1 = sptr1 + MAXLINE;
nb = sprintf(sptr1,"# ending_seed= %lld",ending_seed);
sptr1 = sptr1 + MAXLINE;
sprintf(sptr1,"#");
}

void srf_to_mrf1(struct standrupformat *srf,struct standrupformat *mrf,struct velmodel *vm,int use_srf_lame,int pflag,int ac,char **av)
{
struct srf_prectsegments *prseg_in, *prseg_out;
struct srf_apointvalues *apval_in, *apval_out;
float *stfin, *stfout;
int i, j, it, ig;

double s_mom, m_mom;

float u1, u2, u3, sum;
double lam, l2m, mu;
double ux, uy, uz, vx, vy, vz;

double arg;
double cosS, sinS;
double cosD, sinD;
double cosL, sinL;

double rperd = 0.017453292519943;

if(atof(srf[0].version) < 2.0)
   {
   fprintf(stderr,"srf version= %s < 2.0, exiting ... \n",srf[0].version);
   exit(-1);
   }

s_mom = 0.0;
m_mom = 0.0;

sprintf(mrf[0].version,"3.0");

if(strncmp(mrf[0].src_format,"MOMENT-1MECH",12) != 0 && strncmp(mrf[0].src_format,"MOMENT-6MECH",12) != 0)
   sprintf(mrf[0].src_format,"MOMENT");

copy_hcmnt(mrf,srf);

if(pflag && atof(mrf[0].version) >= 2.0)
   load_command_srf(mrf,ac,av);

mrf[0].type[0] = '\0';
if(strncmp(srf[0].type,"PLANE",5) == 0)
   {
   strcpy(mrf[0].type,srf[0].type);

   mrf[0].srf_prect.nseg = srf[0].srf_prect.nseg;
   mrf[0].srf_prect.prectseg = (struct srf_prectsegments *)check_malloc(mrf[0].srf_prect.nseg*sizeof(struct srf_prectsegments));

   prseg_in = srf[0].srf_prect.prectseg;
   prseg_out = mrf[0].srf_prect.prectseg;
   for(ig=0;ig<mrf[0].srf_prect.nseg;ig++)
      {
      prseg_out[ig].elon = prseg_in[ig].elon;
      prseg_out[ig].elat = prseg_in[ig].elat;
      prseg_out[ig].nstk = prseg_in[ig].nstk;
      prseg_out[ig].ndip = prseg_in[ig].ndip;
      prseg_out[ig].flen = prseg_in[ig].flen;
      prseg_out[ig].fwid = prseg_in[ig].fwid;
      prseg_out[ig].stk = prseg_in[ig].stk;
      prseg_out[ig].dip = prseg_in[ig].dip;
      prseg_out[ig].dtop = prseg_in[ig].dtop;
      prseg_out[ig].shyp = prseg_in[ig].shyp;
      prseg_out[ig].dhyp = prseg_in[ig].dhyp;
      }
   }

mrf[0].srf_apnts.np = srf[0].srf_apnts.np;
mrf[0].srf_apnts.apntvals = (struct srf_apointvalues *)check_malloc((mrf[0].srf_apnts.np)*sizeof(struct srf_apointvalues));

mrf[0].nseg = srf[0].nseg;
mrf[0].np_seg = (int *)check_malloc((mrf[0].nseg)*sizeof(int));

for(i=0;i<mrf[0].nseg;i++)
   mrf[0].np_seg[i] = srf[0].np_seg[i];

apval_in = srf[0].srf_apnts.apntvals;
apval_out = mrf[0].srf_apnts.apntvals;

for(i=0;i<mrf[0].srf_apnts.np;i++)
   {
   apval_out[i].lon = apval_in[i].lon;
   apval_out[i].lat = apval_in[i].lat;
   apval_out[i].dep = apval_in[i].dep;
   apval_out[i].stk = apval_in[i].stk;
   apval_out[i].dip = apval_in[i].dip;
   apval_out[i].area = apval_in[i].area;
   apval_out[i].tinit = apval_in[i].tinit;
   apval_out[i].dt = apval_in[i].dt;
   apval_out[i].vp = apval_in[i].vp;
   apval_out[i].vs = apval_in[i].vs;
   apval_out[i].den = apval_in[i].den;

   apval_out[i].stk = apval_in[i].stk;
   apval_out[i].dip = apval_in[i].dip;
   apval_out[i].rake = apval_in[i].rake;

   apval_out[i].mrf = NULL;

   apval_out[i].ntmr = apval_in[i].nt1;
   if(apval_in[i].nt2 > apval_out[i].ntmr)
      apval_out[i].ntmr = apval_in[i].nt2;
   if(apval_in[i].nt3 > apval_out[i].ntmr)
      apval_out[i].ntmr = apval_in[i].nt3;

   if(apval_out[i].ntmr > 0)
      {
      apval_out[i].mrf = (float *)check_realloc(apval_out[i].mrf,(apval_out[i].ntmr)*sizeof(float));
      stfout = apval_out[i].mrf;

      for(it=0;it<(apval_out[i].ntmr);it++)
         stfout[it] = 0.0;

      if(apval_in[i].nt1)
         {
         stfin = apval_in[i].stf1;
         for(it=0;it<(apval_in[i].nt1);it++)
            stfout[it] = stfout[it] + stfin[it]*stfin[it];
         }

      if(apval_in[i].nt2)
         {
         stfin = apval_in[i].stf2;
         for(it=0;it<(apval_in[i].nt2);it++)
            stfout[it] = stfout[it] + stfin[it]*stfin[it];
         }

      if(apval_in[i].nt3)
         {
         stfin = apval_in[i].stf3;
         for(it=0;it<(apval_in[i].nt3);it++)
            stfout[it] = stfout[it] + stfin[it]*stfin[it];
         }

      sum = 0.0;
      for(it=0;it<(apval_out[i].ntmr);it++)
         {
         stfout[it] = sqrt(stfout[it]);
         sum = sum + (apval_out[i].dt)*stfout[it];
         }

      sum = 1.0/sum;
      for(it=0;it<(apval_out[i].ntmr);it++)
         stfout[it] = sum*stfout[it];
      }

   if(use_srf_lame != 0 && apval_out[i].vp > 0.0 && apval_out[i].vs > 0.0 && apval_out[i].den > 0.0)
      {
      l2m = apval_out[i].vp*apval_out[i].vp*apval_out[i].den;
      mu = apval_out[i].vs*apval_out[i].vs*apval_out[i].den;
      }
   else
      {
      j = 0;
      while(vm->dep[j] < apval_out[i].dep)
         j++;

      l2m = 1.0e+10*vm->vp[j]*vm->vp[j]*vm->den[j];
      mu = 1.0e+10*vm->vs[j]*vm->vs[j]*vm->den[j];

      if(apval_out[i].vp < 0.0)
         apval_out[i].vp = 1.0e+05*vm->vp[j];

      if(apval_out[i].vs < 0.0)
         apval_out[i].vs = 1.0e+05*vm->vs[j];

      if(apval_out[i].den < 0.0)
         apval_out[i].den = vm->den[j];
      }
   lam = l2m - 2.0*mu;

   u1 = apval_in[i].slip1;
   u2 = apval_in[i].slip2;
   u3 = apval_in[i].slip3;

   arg = apval_in[i].stk*rperd;
   cosS = cos(arg);
   sinS = sin(arg);

   arg = apval_in[i].dip*rperd;
   cosD = cos(arg);
   sinD = sin(arg);

   arg = apval_in[i].rake*rperd;
   cosL = cos(arg);
   sinL = sin(arg);

   vx = -sinD*sinS;
   vy =  sinD*cosS;
   vz = -cosD;

   ux = -(u3*sinD - cosD*(u1*sinL + u2*cosL))*sinS + (u1*cosL - u2*sinL)*cosS;
   uy =  (u3*sinD - cosD*(u1*sinL + u2*cosL))*cosS + (u1*cosL - u2*sinL)*sinS;
   uz = -u3*cosD - (u1*sinL + u2*cosL)*sinD;

   apval_out[i].mnn = (l2m*vx*ux + lam*(vy*uy + vz*uz))*apval_out[i].area;
   apval_out[i].mee = (l2m*vy*uy + lam*(vx*ux + vz*uz))*apval_out[i].area;
   apval_out[i].mdd = (l2m*vz*uz + lam*(vx*ux + vy*uy))*apval_out[i].area;

   apval_out[i].mne = mu*(vx*uy + vy*ux)*apval_out[i].area;
   apval_out[i].mnd = mu*(vx*uz + vz*ux)*apval_out[i].area;
   apval_out[i].med = mu*(vy*uz + vz*uy)*apval_out[i].area;

   s_mom = s_mom + sqrt(u1*u1 + u2*u2 + u3*u3)*mu*apval_in[i].area;
   m_mom = m_mom + sqrt(0.5*apval_out[i].mnn*apval_out[i].mnn
                      + 0.5*apval_out[i].mee*apval_out[i].mee
		      + 0.5*apval_out[i].mdd*apval_out[i].mdd
		      + apval_out[i].mne*apval_out[i].mne
		      + apval_out[i].mnd*apval_out[i].mnd
		      + apval_out[i].med*apval_out[i].med);

   if(strncmp(mrf[0].src_format,"MOMENT-6MECH",12) == 0)
      {
      apval_out[i].mr_nn = NULL;
      apval_out[i].mr_ee = NULL;
      apval_out[i].mr_dd = NULL;
      apval_out[i].mr_ne = NULL;
      apval_out[i].mr_nd = NULL;
      apval_out[i].mr_ed = NULL;

      if(apval_out[i].ntmr > 0)
         {
         stfin = apval_out[i].mrf;

         apval_out[i].mr_nn = (float *)check_realloc(apval_out[i].mr_nn,(apval_out[i].ntmr)*sizeof(float));
         stfout = apval_out[i].mr_nn;
         for(it=0;it<(apval_out[i].ntmr);it++)
            stfout[it] = apval_out[i].mnn*stfin[it];

         apval_out[i].mr_ee = (float *)check_realloc(apval_out[i].mr_ee,(apval_out[i].ntmr)*sizeof(float));
         stfout = apval_out[i].mr_ee;
         for(it=0;it<(apval_out[i].ntmr);it++)
            stfout[it] = apval_out[i].mee*stfin[it];

         apval_out[i].mr_dd = (float *)check_realloc(apval_out[i].mr_dd,(apval_out[i].ntmr)*sizeof(float));
         stfout = apval_out[i].mr_dd;
         for(it=0;it<(apval_out[i].ntmr);it++)
            stfout[it] = apval_out[i].mdd*stfin[it];

         apval_out[i].mr_ne = (float *)check_realloc(apval_out[i].mr_ne,(apval_out[i].ntmr)*sizeof(float));
         stfout = apval_out[i].mr_ne;
         for(it=0;it<(apval_out[i].ntmr);it++)
            stfout[it] = apval_out[i].mne*stfin[it];

         apval_out[i].mr_nd = (float *)check_realloc(apval_out[i].mr_nd,(apval_out[i].ntmr)*sizeof(float));
         stfout = apval_out[i].mr_nd;
         for(it=0;it<(apval_out[i].ntmr);it++)
            stfout[it] = apval_out[i].mnd*stfin[it];

         apval_out[i].mr_ed = (float *)check_realloc(apval_out[i].mr_ed,(apval_out[i].ntmr)*sizeof(float));
         stfout = apval_out[i].mr_ed;
         for(it=0;it<(apval_out[i].ntmr);it++)
            stfout[it] = apval_out[i].med*stfin[it];
         }
      }
   }
fprintf(stderr,"slip moment= %.5e\n",s_mom);
fprintf(stderr,"MT moment= %.5e\n",m_mom);
}

void srf_to_mrf6_dsamp(struct standrupformat *srf,struct standrupformat *mrf,struct velmodel *vm,int ncrs_stk,int ncrs_dip,int stk_off,int dip_off,int use_srf_lame,int pflag,int ac,char **av)
{
struct srf_prectsegments *prseg_in, *prseg_out;
struct srf_apointvalues *apval_in, *apval_out;
float *stfp, *stfout;

int it, ig, ip_in, ip_out, ix, iy, ixp, iyp;
int ntot_in, ntot_out, *nstk_out, *ndip_out, *nstk_in, *ndip_in;
int ix0, ixend, ixs, iy0, iyend, iys, nts, ips;
int id, it0, it_out, ntmr_back, nts_back;

float tmin, tmax, tend;

double Mnn, Mee, Mdd, Mne, Mnd, Med;
double *mr_nnD, *mr_eeD, *mr_ddD, *mr_neD, *mr_ndD, *mr_edD;
double *s_stfD;

double s_mom, m_mom_f, m_mom_c, sumD;

double u1, u2, u3;
double lam, l2m, mu;
double ux, uy, uz, vx, vy, vz;

double arg;
double cosS, sinS;
double cosD, sinD;
double cosL, sinL;

double rperd = 0.017453292519943;

if(atof(srf[0].version) < 2.0)
   {
   fprintf(stderr,"srf version= %s < 2.0, exiting ... \n",srf[0].version);
   exit(-1);
   }

s_mom = 0.0;
m_mom_f = 0.0;
m_mom_c = 0.0;

sprintf(mrf[0].version,"3.0");
sprintf(mrf[0].src_format,"MOMENT-6MECH");

copy_hcmnt(mrf,srf);

if(pflag && atof(mrf[0].version) >= 2.0)
   load_command_srf(mrf,ac,av);

mrf[0].type[0] = '\0';
if(strncmp(srf[0].type,"PLANE",5) == 0)
   {
   strcpy(mrf[0].type,srf[0].type);

   mrf[0].srf_prect.nseg = srf[0].srf_prect.nseg;
   mrf[0].srf_prect.prectseg = (struct srf_prectsegments *)check_malloc(mrf[0].srf_prect.nseg*sizeof(struct srf_prectsegments));

   prseg_in = srf[0].srf_prect.prectseg;
   prseg_out = mrf[0].srf_prect.prectseg;

   nstk_in = (int *)check_malloc((srf[0].srf_prect.nseg)*sizeof(int));
   ndip_in = (int *)check_malloc((srf[0].srf_prect.nseg)*sizeof(int));
   nstk_out = (int *)check_malloc((mrf[0].srf_prect.nseg)*sizeof(int));
   ndip_out = (int *)check_malloc((mrf[0].srf_prect.nseg)*sizeof(int));

   for(ig=0;ig<mrf[0].srf_prect.nseg;ig++)
      {
      prseg_out[ig].elon = prseg_in[ig].elon;
      prseg_out[ig].elat = prseg_in[ig].elat;

      prseg_out[ig].nstk = (int)(1.0*prseg_in[ig].nstk/ncrs_stk + 0.5);
      while(prseg_out[ig].nstk*ncrs_stk > prseg_in[ig].nstk)
         prseg_out[ig].nstk--;

      prseg_out[ig].ndip = (int)(1.0*prseg_in[ig].ndip/ncrs_dip + 0.5);
      while(prseg_out[ig].ndip*ncrs_dip > prseg_in[ig].ndip)
         prseg_out[ig].ndip--;

      prseg_out[ig].flen = prseg_out[ig].nstk*ncrs_stk*(prseg_in[ig].flen/prseg_in[ig].nstk);;
      prseg_out[ig].fwid = prseg_out[ig].ndip*ncrs_dip*(prseg_in[ig].fwid/prseg_in[ig].ndip);;

      prseg_out[ig].stk = prseg_in[ig].stk;
      prseg_out[ig].dip = prseg_in[ig].dip;
      prseg_out[ig].dtop = prseg_in[ig].dtop;
      prseg_out[ig].shyp = prseg_in[ig].shyp;
      prseg_out[ig].dhyp = prseg_in[ig].dhyp;

      nstk_in[ig] = prseg_in[ig].nstk;
      ndip_in[ig] = prseg_in[ig].ndip;
      nstk_out[ig] = prseg_out[ig].nstk;
      ndip_out[ig] = prseg_out[ig].ndip;
      }
   }

mrf[0].nseg = srf[0].nseg;
mrf[0].np_seg = (int *)check_malloc((mrf[0].nseg)*sizeof(int));

mrf[0].srf_apnts.np = 0;
for(ig=0;ig<mrf[0].nseg;ig++)
   {
   mrf[0].np_seg[ig] = nstk_out[ig]*ndip_out[ig];
   mrf[0].srf_apnts.np = mrf[0].srf_apnts.np + mrf[0].np_seg[ig];
   }

mrf[0].srf_apnts.apntvals = (struct srf_apointvalues *)check_malloc((mrf[0].srf_apnts.np)*sizeof(struct srf_apointvalues));

apval_in = srf[0].srf_apnts.apntvals;
apval_out = mrf[0].srf_apnts.apntvals;

mr_nnD = NULL;
mr_eeD = NULL;
mr_ddD = NULL;
mr_neD = NULL;
mr_ndD = NULL;
mr_edD = NULL;
s_stfD = NULL;

ntmr_back = 0;
nts_back = 0;

ntot_in = 0;
ntot_out = 0;
for(ig=0;ig<mrf[0].nseg;ig++)
   {
   for(iy=0;iy<ndip_out[ig];iy++)
      {
      iyp = iy*ncrs_dip + dip_off;
      for(ix=0;ix<nstk_out[ig];ix++)
         {
         ixp = ix*ncrs_stk + stk_off;

         ip_in = ixp + iyp*nstk_in[ig] + ntot_in;
         ip_out = ix + iy*nstk_out[ig] + ntot_out;

/* use center point values for following parameters */

         apval_out[ip_out].lon = apval_in[ip_in].lon;
         apval_out[ip_out].lat = apval_in[ip_in].lat;
         apval_out[ip_out].dep = apval_in[ip_in].dep;
         apval_out[ip_out].area = apval_in[ip_in].area*(ncrs_stk*ncrs_dip);
         apval_out[ip_out].dt = apval_in[ip_in].dt;
         apval_out[ip_out].vp = apval_in[ip_in].vp;
         apval_out[ip_out].vs = apval_in[ip_in].vs;
         apval_out[ip_out].den = apval_in[ip_in].den;

         apval_out[ip_out].stk = apval_in[ip_in].stk;
         apval_out[ip_out].dip = apval_in[ip_in].dip;
         apval_out[ip_out].rake = apval_in[ip_in].rake;

/* determine starting indicies of input points */

         ix0 = ixp - (ncrs_stk-1)/2;
	 ixend = ix0 + ncrs_stk;
	 if(ix0 < 0)
	    ix0 = 0;

         iy0 = iyp - (ncrs_dip-1)/2;
	 iyend = iy0 + ncrs_dip;
	 if(iy0 < 0)
	    iy0 = 0;

/* determine minimum tinit and total ntmr needed for this group of input points */

/*
fprintf(stderr,"ix0= %d ixend= %d iy0= %d iyend= %d\n",ix0,ixend,iy0,iyend);
*/

	 tmin = 1.0e+15;
	 tmax = -1.0e+15;
         for(iys=iy0;iys<iyend;iys++)
            {
            for(ixs=ix0;ixs<ixend;ixs++)
               {
               ips = ixs + iys*nstk_in[ig] + ntot_in;
         
	       if(apval_in[ips].tinit < tmin)
	          tmin = apval_in[ips].tinit;
         
               nts = apval_in[ips].nt1;
               if(apval_in[ips].nt2 > nts)
                  nts = apval_in[ips].nt2;
               if(apval_in[ips].nt3 > nts)
                  nts = apval_in[ips].nt3;
         
	       tend = apval_in[ips].tinit + nts*apval_in[ips].dt;
	       if(tend > tmax)
	          tmax = tend;

               if(use_srf_lame == 0 || apval_in[ips].vp < 0.0 || apval_in[ips].vs < 0.0 || apval_in[ips].den < 0.0)
                  {
                  id = 0;
                  while(vm->dep[id] < apval_in[ips].dep)
                     id++;

                  if(apval_in[ips].vp < 0.0)
                     apval_in[ips].vp = 1.0e+05*vm->vp[id];

                  if(apval_in[ips].vs < 0.0)
                     apval_in[ips].vs = 1.0e+05*vm->vs[id];

                  if(apval_in[ips].den < 0.0)
                     apval_in[ips].den = vm->den[id];

                  if(ips == ip_in) /* reset output center point values */
                     {
                     apval_out[ip_out].vp = apval_in[ip_in].vp;
                     apval_out[ip_out].vs = apval_in[ip_in].vs;
                     apval_out[ip_out].den = apval_in[ip_in].den;
	             }
                  }
	       }
            }

         apval_out[ip_out].tinit = tmin;

         apval_out[ip_out].mnn = 0.0;
         apval_out[ip_out].mee = 0.0;
         apval_out[ip_out].mdd = 0.0;
         apval_out[ip_out].mne = 0.0;
         apval_out[ip_out].mnd = 0.0;
         apval_out[ip_out].med = 0.0;

         apval_out[ip_out].mr_nn = NULL;
         apval_out[ip_out].mr_ee = NULL;
         apval_out[ip_out].mr_dd = NULL;
         apval_out[ip_out].mr_ne = NULL;
         apval_out[ip_out].mr_nd = NULL;
         apval_out[ip_out].mr_ed = NULL;

         apval_out[ip_out].ntmr = (int)((tmax-tmin)/apval_out[ip_out].dt + 0.5);

/*
fprintf(stderr,"ip_o= %d ip_i= %d ntmr= %d tmin= %.4f tmax= %.4f\n",ip_out,ip_in,apval_out[ip_out].ntmr,tmin,tmax);
*/

	 if(apval_out[ip_out].ntmr > 0)
	    {
	    if(apval_out[ip_out].ntmr > ntmr_back)
	       {
               mr_nnD = (double *)check_realloc(mr_nnD,(apval_out[ip_out].ntmr)*sizeof(double));
               mr_eeD = (double *)check_realloc(mr_eeD,(apval_out[ip_out].ntmr)*sizeof(double));
               mr_ddD = (double *)check_realloc(mr_ddD,(apval_out[ip_out].ntmr)*sizeof(double));
               mr_neD = (double *)check_realloc(mr_neD,(apval_out[ip_out].ntmr)*sizeof(double));
               mr_ndD = (double *)check_realloc(mr_ndD,(apval_out[ip_out].ntmr)*sizeof(double));
               mr_edD = (double *)check_realloc(mr_edD,(apval_out[ip_out].ntmr)*sizeof(double));

	       ntmr_back = apval_out[ip_out].ntmr;
	       }

	    for(it=0;it<apval_out[ip_out].ntmr;it++)
	       {
	       mr_nnD[it] = 0.0;
	       mr_eeD[it] = 0.0;
	       mr_ddD[it] = 0.0;
	       mr_neD[it] = 0.0;
	       mr_ndD[it] = 0.0;
	       mr_edD[it] = 0.0;
	       }

/* compute moment-rate functions for each input point and sum across group */

            for(iys=iy0;iys<iyend;iys++)
               {
               for(ixs=ix0;ixs<ixend;ixs++)
                  {
                  ips = ixs + iys*nstk_in[ig] + ntot_in;
         
                  nts = apval_in[ips].nt1;
                  if(apval_in[ips].nt2 > nts)
                     nts = apval_in[ips].nt2;
                  if(apval_in[ips].nt3 > nts)
                     nts = apval_in[ips].nt3;

	          if(nts > nts_back)
	             {
                     s_stfD = (double *)check_realloc(s_stfD,nts*sizeof(double));
	             nts_back = nts;
	             }

	          for(it=0;it<nts;it++)
	             s_stfD[it] = 0.0;

                  if(apval_in[ips].nt1)
                     {
                     stfp = apval_in[ips].stf1;
                     for(it=0;it<(apval_in[ips].nt1);it++)
                        s_stfD[it] = s_stfD[it] + stfp[it]*stfp[it];
                     }

                  if(apval_in[ips].nt2)
                     {
                     stfp = apval_in[ips].stf2;
                     for(it=0;it<(apval_in[ips].nt2);it++)
                        s_stfD[it] = s_stfD[it] + stfp[it]*stfp[it];
                     }

                  if(apval_in[ips].nt3)
                     {
                     stfp = apval_in[ips].stf3;
                     for(it=0;it<(apval_in[ips].nt3);it++)
                        s_stfD[it] = s_stfD[it] + stfp[it]*stfp[it];
                     }

                  sumD = 0.0;
                  for(it=0;it<nts;it++)
                     {
                     s_stfD[it] = sqrt(s_stfD[it]);
                     sumD = sumD + (apval_in[ips].dt)*s_stfD[it];
                     }

                  sumD = 1.0/sumD;
                  for(it=0;it<nts;it++)
                     s_stfD[it] = sumD*s_stfD[it];

                  l2m = apval_in[ips].vp*apval_in[ips].vp*apval_in[ips].den;
                  mu = apval_in[ips].vs*apval_in[ips].vs*apval_in[ips].den;
                  lam = l2m - 2.0*mu;

                  u1 = apval_in[ips].slip1;
                  u2 = apval_in[ips].slip2;
                  u3 = apval_in[ips].slip3;

                  arg = apval_in[ips].stk*rperd;
                  cosS = cos(arg);
                  sinS = sin(arg);

                  arg = apval_in[ips].dip*rperd;
                  cosD = cos(arg);
                  sinD = sin(arg);

                  arg = apval_in[ips].rake*rperd;
                  cosL = cos(arg);
                  sinL = sin(arg);

                  vx = -sinD*sinS;
                  vy =  sinD*cosS;
                  vz = -cosD;

                  ux = -(u3*sinD - cosD*(u1*sinL + u2*cosL))*sinS + (u1*cosL - u2*sinL)*cosS;
                  uy =  (u3*sinD - cosD*(u1*sinL + u2*cosL))*cosS + (u1*cosL - u2*sinL)*sinS;
                  uz = -u3*cosD - (u1*sinL + u2*cosL)*sinD;

                  Mnn = (l2m*vx*ux + lam*(vy*uy + vz*uz))*apval_in[ips].area;
                  Mee = (l2m*vy*uy + lam*(vx*ux + vz*uz))*apval_in[ips].area;
                  Mdd = (l2m*vz*uz + lam*(vx*ux + vy*uy))*apval_in[ips].area;

                  Mne = mu*(vx*uy + vy*ux)*apval_in[ips].area;
                  Mnd = mu*(vx*uz + vz*ux)*apval_in[ips].area;
                  Med = mu*(vy*uz + vz*uy)*apval_in[ips].area;

/* sum fine grid moment-rates into coarse grid arrays */
         
	          it0 = (int)((apval_in[ips].tinit - apval_out[ip_out].tinit)/apval_in[ips].dt + 0.5);
	          if(it0 < 0)
	             it0 = 0;

	          while((it0+nts) > apval_out[ip_out].ntmr)
	             nts--;

	          for(it=0;it<nts;it++)
	             {
		     it_out = it + it0;

	             mr_nnD[it_out] = mr_nnD[it_out] + Mnn*s_stfD[it];
	             mr_eeD[it_out] = mr_eeD[it_out] + Mee*s_stfD[it];
	             mr_ddD[it_out] = mr_ddD[it_out] + Mdd*s_stfD[it];
	             mr_neD[it_out] = mr_neD[it_out] + Mne*s_stfD[it];
	             mr_ndD[it_out] = mr_ndD[it_out] + Mnd*s_stfD[it];
	             mr_edD[it_out] = mr_edD[it_out] + Med*s_stfD[it];
	             }

                  s_mom = s_mom + sqrt(u1*u1 + u2*u2 + u3*u3)*mu*apval_in[ips].area;
                  m_mom_f = m_mom_f + sqrt(0.5*(Mnn*Mnn + Mee*Mee + Mdd*Mdd) + Mne*Mne + Mnd*Mnd + Med*Med);
	          }
               }

/* copy individual moment-rates to output structure */

            apval_out[ip_out].mr_nn = (float *)check_realloc(apval_out[ip_out].mr_nn,(apval_out[ip_out].ntmr)*sizeof(float));
            stfp = apval_out[ip_out].mr_nn;
            for(it=0;it<apval_out[ip_out].ntmr;it++)
	       {
	       stfp[it] = mr_nnD[it];
               apval_out[ip_out].mnn = apval_out[ip_out].mnn + (apval_out[ip_out].dt)*mr_nnD[it];
	       }

            apval_out[ip_out].mr_ee = (float *)check_realloc(apval_out[ip_out].mr_ee,(apval_out[ip_out].ntmr)*sizeof(float));
            stfp = apval_out[ip_out].mr_ee;
            for(it=0;it<apval_out[ip_out].ntmr;it++)
	       {
	       stfp[it] = mr_eeD[it];
               apval_out[ip_out].mee = apval_out[ip_out].mee + (apval_out[ip_out].dt)*mr_eeD[it];
	       }

            apval_out[ip_out].mr_dd = (float *)check_realloc(apval_out[ip_out].mr_dd,(apval_out[ip_out].ntmr)*sizeof(float));
            stfp = apval_out[ip_out].mr_dd;
            for(it=0;it<apval_out[ip_out].ntmr;it++)
	       {
	       stfp[it] = mr_ddD[it];
               apval_out[ip_out].mdd = apval_out[ip_out].mdd + (apval_out[ip_out].dt)*mr_ddD[it];
	       }

            apval_out[ip_out].mr_ne = (float *)check_realloc(apval_out[ip_out].mr_ne,(apval_out[ip_out].ntmr)*sizeof(float));
            stfp = apval_out[ip_out].mr_ne;
            for(it=0;it<apval_out[ip_out].ntmr;it++)
	       {
	       stfp[it] = mr_neD[it];
               apval_out[ip_out].mne = apval_out[ip_out].mne + (apval_out[ip_out].dt)*mr_neD[it];
	       }

            apval_out[ip_out].mr_nd = (float *)check_realloc(apval_out[ip_out].mr_nd,(apval_out[ip_out].ntmr)*sizeof(float));
            stfp = apval_out[ip_out].mr_nd;
            for(it=0;it<apval_out[ip_out].ntmr;it++)
	       {
	       stfp[it] = mr_ndD[it];
               apval_out[ip_out].mnd = apval_out[ip_out].mnd + (apval_out[ip_out].dt)*mr_ndD[it];
	       }

            apval_out[ip_out].mr_ed = (float *)check_realloc(apval_out[ip_out].mr_ed,(apval_out[ip_out].ntmr)*sizeof(float));
            stfp = apval_out[ip_out].mr_ed;
            for(it=0;it<apval_out[ip_out].ntmr;it++)
	       {
	       stfp[it] = mr_edD[it];
               apval_out[ip_out].med = apval_out[ip_out].med + (apval_out[ip_out].dt)*mr_edD[it];
	       }

            m_mom_c = m_mom_c + sqrt(0.5*apval_out[ip_out].mnn*apval_out[ip_out].mnn
                                   + 0.5*apval_out[ip_out].mee*apval_out[ip_out].mee
                                   + 0.5*apval_out[ip_out].mdd*apval_out[ip_out].mdd
                                   + apval_out[ip_out].mne*apval_out[ip_out].mne
                                   + apval_out[ip_out].mnd*apval_out[ip_out].mnd
                                   + apval_out[ip_out].med*apval_out[ip_out].med);
	    }
         }
      }

   ntot_in = ntot_in + srf[0].np_seg[ig];
   ntot_out = ntot_out + mrf[0].np_seg[ig];
   }

free(mr_nnD);
free(mr_eeD);
free(mr_ddD);
free(mr_neD);
free(mr_ndD);
free(mr_edD);
free(s_stfD);

fprintf(stderr,"slip moment= %.5e\n",s_mom);
fprintf(stderr,"MT moment (fine grid)= %.5e\n",m_mom_f);
fprintf(stderr,"MT moment (coarse grid)= %.5e\n",m_mom_c);
}

void srf_dwnsamp(struct standrupformat *srf_in,struct standrupformat *srf_out,int ncrs_stk,int ncrs_dip,int stk_off,int dip_off,int pflag,int ac,char **av)
{
struct srf_prectsegments *prseg_in, *prseg_out;
struct srf_apointvalues *apval_in, *apval_out;

int it, ig, ip_in, ip_out, ix, iy, ixp, iyp;
int ntot_in, ntot_out, *nstk_out, *ndip_out, *nstk_in, *ndip_in;
int ix0, ixend, ixs, iy0, iyend, iys, nts, ips, nf;
int id, it0, it_out;

float *stfp_i, *stfp_o;
float tmin, tmax1, tmax2, tmax3, tend, sfac;
float vp_avg, vs_avg, den_avg, stk_avg, dip_avg, rak_avg, slip1_avg, slip2_avg, slip3_avg;
float xx, yy;

double dperr, stk1, stk2, dip1, dip2, rak1, rak2;
double rperd = 0.017453292519943;

double u1, u2, u3, mu, s_mom_f, s_mom_c;

if(atof(srf_in[0].version) < 2.0)
   {
   fprintf(stderr,"srf version= %s < 2.0, exiting ... \n",srf_in[0].version);
   exit(-1);
   }

dperr = 1.0/rperd;

s_mom_f = 0.0;
s_mom_c = 0.0;

strcpy(srf_out[0].version,srf_in[0].version);
sprintf(srf_out[0].src_format,"SLIP");

copy_hcmnt(srf_out,srf_in);

if(pflag && atof(srf_out[0].version) >= 2.0)
   load_command_srf(srf_out,ac,av);

srf_out[0].type[0] = '\0';
if(strncmp(srf_in[0].type,"PLANE",5) == 0)
   {
   strcpy(srf_out[0].type,srf_in[0].type);

   srf_out[0].srf_prect.nseg = srf_in[0].srf_prect.nseg;
   srf_out[0].srf_prect.prectseg = (struct srf_prectsegments *)check_malloc(srf_out[0].srf_prect.nseg*sizeof(struct srf_prectsegments));

   prseg_in = srf_in[0].srf_prect.prectseg;
   prseg_out = srf_out[0].srf_prect.prectseg;

   nstk_in = (int *)check_malloc((srf_in[0].srf_prect.nseg)*sizeof(int));
   ndip_in = (int *)check_malloc((srf_in[0].srf_prect.nseg)*sizeof(int));
   nstk_out = (int *)check_malloc((srf_out[0].srf_prect.nseg)*sizeof(int));
   ndip_out = (int *)check_malloc((srf_out[0].srf_prect.nseg)*sizeof(int));

   for(ig=0;ig<srf_out[0].srf_prect.nseg;ig++)
      {
      prseg_out[ig].elon = prseg_in[ig].elon;
      prseg_out[ig].elat = prseg_in[ig].elat;

      prseg_out[ig].nstk = (int)(1.0*prseg_in[ig].nstk/ncrs_stk + 0.5);
      while(prseg_out[ig].nstk*ncrs_stk > prseg_in[ig].nstk)
         prseg_out[ig].nstk--;

      prseg_out[ig].ndip = (int)(1.0*prseg_in[ig].ndip/ncrs_dip + 0.5);
      while(prseg_out[ig].ndip*ncrs_dip > prseg_in[ig].ndip)
         prseg_out[ig].ndip--;

      prseg_out[ig].flen = prseg_out[ig].nstk*ncrs_stk*(prseg_in[ig].flen/prseg_in[ig].nstk);;
      prseg_out[ig].fwid = prseg_out[ig].ndip*ncrs_dip*(prseg_in[ig].fwid/prseg_in[ig].ndip);;

      prseg_out[ig].stk = prseg_in[ig].stk;
      prseg_out[ig].dip = prseg_in[ig].dip;
      prseg_out[ig].dtop = prseg_in[ig].dtop;
      prseg_out[ig].shyp = prseg_in[ig].shyp;
      prseg_out[ig].dhyp = prseg_in[ig].dhyp;

      nstk_in[ig] = prseg_in[ig].nstk;
      ndip_in[ig] = prseg_in[ig].ndip;
      nstk_out[ig] = prseg_out[ig].nstk;
      ndip_out[ig] = prseg_out[ig].ndip;
      }
   }

srf_out[0].nseg = srf_in[0].nseg;
srf_out[0].np_seg = (int *)check_malloc((srf_out[0].nseg)*sizeof(int));

srf_out[0].srf_apnts.np = 0;
for(ig=0;ig<srf_out[0].nseg;ig++)
   {
   srf_out[0].np_seg[ig] = nstk_out[ig]*ndip_out[ig];
   srf_out[0].srf_apnts.np = srf_out[0].srf_apnts.np + srf_out[0].np_seg[ig];
   }

srf_out[0].srf_apnts.apntvals = (struct srf_apointvalues *)check_malloc((srf_out[0].srf_apnts.np)*sizeof(struct srf_apointvalues));

apval_in = srf_in[0].srf_apnts.apntvals;
apval_out = srf_out[0].srf_apnts.apntvals;

ntot_in = 0;
ntot_out = 0;
for(ig=0;ig<srf_out[0].nseg;ig++)
   {
   for(iy=0;iy<ndip_out[ig];iy++)
      {
      iyp = iy*ncrs_dip + dip_off;
      for(ix=0;ix<nstk_out[ig];ix++)
         {
         ixp = ix*ncrs_stk + stk_off;

         ip_in = ixp + iyp*nstk_in[ig] + ntot_in;
         ip_out = ix + iy*nstk_out[ig] + ntot_out;

/* use center point values for following parameters */

         apval_out[ip_out].lon = apval_in[ip_in].lon;
         apval_out[ip_out].lat = apval_in[ip_in].lat;
         apval_out[ip_out].dep = apval_in[ip_in].dep;
         apval_out[ip_out].area = apval_in[ip_in].area*(ncrs_stk*ncrs_dip);
         apval_out[ip_out].dt = apval_in[ip_in].dt;

/* determine starting indicies of input points */

         ix0 = ixp - (ncrs_stk-1)/2;
	 ixend = ix0 + ncrs_stk;
	 if(ix0 < 0)
	    ix0 = 0;

         iy0 = iyp - (ncrs_dip-1)/2;
	 iyend = iy0 + ncrs_dip;
	 if(iy0 < 0)
	    iy0 = 0;

/* determine avg. values, minimum tinit and total nt needed for this group of input points */

/*
fprintf(stderr,"ix0= %d ixend= %d iy0= %d iyend= %d\n",ix0,ixend,iy0,iyend);
*/

	 nf = 0;
	 vp_avg = 0.0;
	 vs_avg = 0.0;
	 den_avg = 0.0;

	 stk_avg = 0.0;
	 dip_avg = 0.0;
	 rak_avg = 0.0;

	 stk1 = 0.0;
	 stk2 = 0.0;
	 dip1 = 0.0;
	 dip2 = 0.0;
	 rak1 = 0.0;
	 rak2 = 0.0;

	 slip1_avg = 0.0;
	 slip2_avg = 0.0;
	 slip3_avg = 0.0;

	 tmin = 1.0e+15;
	 tmax1 = -1.0e+15;
	 tmax2 = -1.0e+15;
	 tmax3 = -1.0e+15;
         for(iys=iy0;iys<iyend;iys++)
            {
            for(ixs=ix0;ixs<ixend;ixs++)
               {
               ips = ixs + iys*nstk_in[ig] + ntot_in;

	       nf++;
         
	       vp_avg = vp_avg + apval_in[ips].vp;
	       vs_avg = vs_avg + apval_in[ips].vs;
	       den_avg = den_avg + apval_in[ips].den;

	       stk_avg = stk_avg + apval_in[ips].stk;
	       dip_avg = dip_avg + apval_in[ips].dip;
	       rak_avg = rak_avg + apval_in[ips].rake;

	       stk1 = stk1 + cos(rperd*apval_in[ips].stk);
	       stk2 = stk2 + sin(rperd*apval_in[ips].stk);
	       dip1 = dip1 + cos(rperd*apval_in[ips].dip);
	       dip2 = dip2 + sin(rperd*apval_in[ips].dip);
	       rak1 = rak1 + cos(rperd*apval_in[ips].rake);
	       rak2 = rak2 + sin(rperd*apval_in[ips].rake);

	       slip1_avg = slip1_avg + apval_in[ips].slip1;
	       slip2_avg = slip2_avg + apval_in[ips].slip2;
	       slip3_avg = slip3_avg + apval_in[ips].slip3;

	       if(apval_in[ips].tinit < tmin)
	          tmin = apval_in[ips].tinit;
         
               if(apval_in[ips].nt1 > 0)
	          {
	          tend = apval_in[ips].tinit + apval_in[ips].nt1*apval_in[ips].dt;
	          if(tend > tmax1)
	             tmax1 = tend;
		  }
         
               if(apval_in[ips].nt2 > 0)
	          {
	          tend = apval_in[ips].tinit + apval_in[ips].nt2*apval_in[ips].dt;
	          if(tend > tmax2)
	             tmax2 = tend;
		  }
         
               if(apval_in[ips].nt3 > 0)
	          {
	          tend = apval_in[ips].tinit + apval_in[ips].nt3*apval_in[ips].dt;
	          if(tend > tmax3)
	             tmax3 = tend;
		  }
   /*
if(ip_out == 100452)
   fprintf(stderr,"%10.4f",apval_in[ips].tinit);
   fprintf(stderr,"%10.4f",apval_in[ips].dt*apval_in[ips].nt1);
   fprintf(stderr,"%10.4f",apval_in[ips].slip1);
   */
	       }
/*
if(ip_out == 100452)
   fprintf(stderr,"\n");
   */
            }

         sfac = 1.0/(1.0*nf);

         apval_out[ip_out].vp = vp_avg*sfac;
         apval_out[ip_out].vs = vs_avg*sfac;
         apval_out[ip_out].den = den_avg*sfac;

/* XXXX problem with unwrapping angles 
         apval_out[ip_out].stk = stk_avg*sfac;
         apval_out[ip_out].dip = dip_avg*sfac;
         apval_out[ip_out].rake = rak_avg*sfac;
*/
/* just use center point values
         apval_out[ip_out].stk = apval_in[ip_in].stk;
         apval_out[ip_out].dip = apval_in[ip_in].dip;
         apval_out[ip_out].rake = apval_in[ip_in].rake;
*/
/* use vector component values
*/
         apval_out[ip_out].stk = dperr*atan2(stk2,stk1);
	 while(apval_out[ip_out].stk > 360.0)
	    apval_out[ip_out].stk = apval_out[ip_out].stk - 360.0;
	 while(apval_out[ip_out].stk < 0.0)
	    apval_out[ip_out].stk = apval_out[ip_out].stk + 360.0;

         apval_out[ip_out].dip = dperr*atan2(dip2,dip1);
	 while(apval_out[ip_out].dip < 0.0)
	    apval_out[ip_out].dip = apval_out[ip_out].dip + 180.0;

         apval_out[ip_out].rake = dperr*atan2(rak2,rak1);
	 while(apval_out[ip_out].rake > 180.0)
	    apval_out[ip_out].rake = apval_out[ip_out].rake - 360.0;
	 while(apval_out[ip_out].rake < -180.0)
	    apval_out[ip_out].rake = apval_out[ip_out].rake + 360.0;

         apval_out[ip_out].tinit = tmin;

         apval_out[ip_out].stf1 = NULL;
         if(tmax1 > 0.0)
	    {
            apval_out[ip_out].slip1 = slip1_avg*sfac;
            apval_out[ip_out].nt1 = (int)((tmax1-tmin)/apval_out[ip_out].dt + 0.5);
	    }
	 else
	    {
            apval_out[ip_out].slip1 = 0.0;
            apval_out[ip_out].nt1 = 0;
	    }

         apval_out[ip_out].stf2 = NULL;
         if(tmax2 > 0.0)
	    {
            apval_out[ip_out].slip2 = slip2_avg*sfac;
            apval_out[ip_out].nt2 = (int)((tmax2-tmin)/apval_out[ip_out].dt + 0.5);
	    }
	 else
	    {
            apval_out[ip_out].slip2 = 0.0;
            apval_out[ip_out].nt2 = 0;
	    }

         apval_out[ip_out].stf3 = NULL;
         if(tmax3 > 0.0)
	    {
            apval_out[ip_out].slip3 = slip3_avg*sfac;
            apval_out[ip_out].nt3 = (int)((tmax3-tmin)/apval_out[ip_out].dt + 0.5);
	    }
	 else
	    {
            apval_out[ip_out].slip3 = 0.0;
            apval_out[ip_out].nt3 = 0;
	    }

/*
if(ip_out == 100452)
   {
   xx = ((ix+0.5)/(1.0*prseg_out[ig].nstk) - 0.5)*prseg_out[ig].flen - prseg_out[ig].shyp;
   yy = (iy+0.5)*prseg_out[ig].fwid/(1.0*prseg_out[ig].ndip) - prseg_out[ig].dhyp;
   fprintf(stderr,"tmin= %10.4f nt1_new= %d nt1_orig= %d rhyp= %10.4f\n",tmin,apval_out[ip_out].nt1,apval_in[ip_in].nt1,sqrt(xx*xx + yy*yy));
   }
*/

/*
fprintf(stderr,"ip_o= %d ip_i= %d nt1= %d tmin= %.4f tmax1= %.4f\n",ip_out,ip_in,apval_out[ip_out].nt1,tmin,tmax1);
*/

	 if(apval_out[ip_out].nt1 > 0 || apval_out[ip_out].nt2 > 0 || apval_out[ip_out].nt3 > 0)
	    {
	    if(apval_out[ip_out].nt1 > 0)
	       {
               apval_out[ip_out].stf1 = (float *)check_realloc(apval_out[ip_out].stf1,(apval_out[ip_out].nt1)*sizeof(float));

               stfp_o = apval_out[ip_out].stf1;
	       for(it=0;it<apval_out[ip_out].nt1;it++)
	          stfp_o[it] = 0.0;
	       }

	    if(apval_out[ip_out].nt2 > 0)
	       {
               apval_out[ip_out].stf2 = (float *)check_realloc(apval_out[ip_out].stf2,(apval_out[ip_out].nt2)*sizeof(float));

               stfp_o = apval_out[ip_out].stf2;
	       for(it=0;it<apval_out[ip_out].nt2;it++)
	          stfp_o[it] = 0.0;
	       }

	    if(apval_out[ip_out].nt3 > 0)
	       {
               apval_out[ip_out].stf3 = (float *)check_realloc(apval_out[ip_out].stf3,(apval_out[ip_out].nt3)*sizeof(float));

               stfp_o = apval_out[ip_out].stf3;
	       for(it=0;it<apval_out[ip_out].nt1;it++)
	          stfp_o[it] = 0.0;
	       }

/* compute slip-rate functions for each input point and sum across group */

            for(iys=iy0;iys<iyend;iys++)
               {
               for(ixs=ix0;ixs<ixend;ixs++)
                  {
                  ips = ixs + iys*nstk_in[ig] + ntot_in;
         
	          it0 = (int)((apval_in[ips].tinit - apval_out[ip_out].tinit)/apval_in[ips].dt + 0.5);
	          if(it0 < 0)		/* shouldn't happen but ... */
	             it0 = 0;

                  if(apval_in[ips].nt1)
                     {
		     nts = apval_in[ips].nt1;
	             while((it0+nts) > apval_out[ip_out].nt1)	/* shouldn't happen but ... */
	                nts--;

                     stfp_i = apval_in[ips].stf1;
                     stfp_o = apval_out[ip_out].stf1;
	             for(it=0;it<nts;it++)
	                {
		        it_out = it + it0;
                        stfp_o[it_out] = stfp_o[it_out] + stfp_i[it];
                        }
                     }

                  if(apval_in[ips].nt2)
                     {
		     nts = apval_in[ips].nt2;
	             while((it0+nts) > apval_out[ip_out].nt2)	/* shouldn't happen but ... */
	                nts--;

                     stfp_i = apval_in[ips].stf2;
                     stfp_o = apval_out[ip_out].stf2;
	             for(it=0;it<nts;it++)
	                {
		        it_out = it + it0;
                        stfp_o[it_out] = stfp_o[it_out] + stfp_i[it];
                        }
                     }

                  if(apval_in[ips].nt3)
                     {
		     nts = apval_in[ips].nt3;
	             while((it0+nts) > apval_out[ip_out].nt3)	/* shouldn't happen but ... */
	                nts--;

                     stfp_i = apval_in[ips].stf3;
                     stfp_o = apval_out[ip_out].stf3;
	             for(it=0;it<nts;it++)
	                {
		        it_out = it + it0;
                        stfp_o[it_out] = stfp_o[it_out] + stfp_i[it];
                        }
                     }

                  u1 = apval_in[ips].slip1;
                  u2 = apval_in[ips].slip2;
                  u3 = apval_in[ips].slip3;
                  mu = apval_in[ips].vs*apval_in[ips].vs*apval_in[ips].den;

                  s_mom_f = s_mom_f + sqrt(u1*u1 + u2*u2 + u3*u3)*mu*apval_in[ips].area;
	          }
               }

/* ensure coarse slip-rates integrate to target slip */

	    if(apval_out[ip_out].nt1)
	       {
               stfp_o = apval_out[ip_out].stf1;

	       sfac = 0.0;
               for(it=0;it<apval_out[ip_out].nt1;it++)
                  sfac = sfac + (apval_out[ip_out].dt)*stfp_o[it];

	       sfac = apval_out[ip_out].slip1/sfac;
               for(it=0;it<apval_out[ip_out].nt1;it++)
                  stfp_o[it] = sfac*stfp_o[it];
	       }

	    if(apval_out[ip_out].nt2)
	       {
               stfp_o = apval_out[ip_out].stf2;

	       sfac = 0.0;
               for(it=0;it<apval_out[ip_out].nt2;it++)
                  sfac = sfac + (apval_out[ip_out].dt)*stfp_o[it];

	       sfac = apval_out[ip_out].slip2/sfac;
               for(it=0;it<apval_out[ip_out].nt2;it++)
                  stfp_o[it] = sfac*stfp_o[it];
	       }

	    if(apval_out[ip_out].nt3)
	       {
               stfp_o = apval_out[ip_out].stf3;

	       sfac = 0.0;
               for(it=0;it<apval_out[ip_out].nt3;it++)
                  sfac = sfac + (apval_out[ip_out].dt)*stfp_o[it];

	       sfac = apval_out[ip_out].slip3/sfac;
               for(it=0;it<apval_out[ip_out].nt3;it++)
                  stfp_o[it] = sfac*stfp_o[it];
	       }


            u1 = apval_out[ip_out].slip1;
            u2 = apval_out[ip_out].slip2;
            u3 = apval_out[ip_out].slip3;
            mu = apval_out[ip_out].vs*apval_out[ip_out].vs*apval_out[ip_out].den;

            s_mom_c = s_mom_c + sqrt(u1*u1 + u2*u2 + u3*u3)*mu*apval_out[ip_out].area;
	    }
         }
      }

   ntot_in = ntot_in + srf_in[0].np_seg[ig];
   ntot_out = ntot_out + srf_out[0].np_seg[ig];
   }

fprintf(stderr,"slip moment (fine grid)= %.5e\n",s_mom_f);
fprintf(stderr,"slip moment (coarse grid)= %.5e\n",s_mom_c);
}

void srf_to_mrf6_dsampXXX(struct standrupformat *srf,struct standrupformat *mrf,struct velmodel *vm,int ncrs_stk,int ncrs_dip,int use_srf_lame,int pflag,int ac,char **av)
{
struct srf_prectsegments *prseg_in, *prseg_out;
struct srf_apointvalues *apval_in, *apval_out;
float *stfp, *stfout;

int it, ig, ip_in, ip_out, ix, iy, ixp, iyp;
int ntot_in, ntot_out, *nstk_out, *ndip_out, *nstk_in, *ndip_in, *stk_off, *dip_off;
int ix0, ixend, ixs, iy0, iyend, iys, nts, ips;
int id, it0, it_out, ntmr_back, nts_back;

float tmin, tmax, tend;

double Mnn, Mee, Mdd, Mne, Mnd, Med;
double *mr_nnD, *mr_eeD, *mr_ddD, *mr_neD, *mr_ndD, *mr_edD;
double *s_stfD;

double s_mom, m_mom_f, m_mom_c, sumD;

double u1, u2, u3;
double lam, l2m, mu;
double ux, uy, uz, vx, vy, vz;

double arg;
double cosS, sinS;
double cosD, sinD;
double cosL, sinL;

double rperd = 0.017453292519943;

if(atof(srf[0].version) < 2.0)
   {
   fprintf(stderr,"srf version= %s < 2.0, exiting ... \n",srf[0].version);
   exit(-1);
   }

s_mom = 0.0;
m_mom_f = 0.0;
m_mom_c = 0.0;

sprintf(mrf[0].version,"3.0");
sprintf(mrf[0].src_format,"MOMENT-6MECH");

copy_hcmnt(mrf,srf);

if(pflag && atof(mrf[0].version) >= 2.0)
   load_command_srf(mrf,ac,av);

mrf[0].type[0] = '\0';
if(strncmp(srf[0].type,"PLANE",5) == 0)
   {
   strcpy(mrf[0].type,srf[0].type);

   mrf[0].srf_prect.nseg = srf[0].srf_prect.nseg;
   mrf[0].srf_prect.prectseg = (struct srf_prectsegments *)check_malloc(mrf[0].srf_prect.nseg*sizeof(struct srf_prectsegments));

   prseg_in = srf[0].srf_prect.prectseg;
   prseg_out = mrf[0].srf_prect.prectseg;

   nstk_in = (int *)check_malloc((srf[0].srf_prect.nseg)*sizeof(int));
   ndip_in = (int *)check_malloc((srf[0].srf_prect.nseg)*sizeof(int));
   nstk_out = (int *)check_malloc((mrf[0].srf_prect.nseg)*sizeof(int));
   ndip_out = (int *)check_malloc((mrf[0].srf_prect.nseg)*sizeof(int));
   stk_off = (int *)check_malloc((mrf[0].srf_prect.nseg)*sizeof(int));
   dip_off = (int *)check_malloc((mrf[0].srf_prect.nseg)*sizeof(int));

   for(ig=0;ig<mrf[0].srf_prect.nseg;ig++)
      {
      prseg_out[ig].elon = prseg_in[ig].elon;
      prseg_out[ig].elat = prseg_in[ig].elat;

      stk_off[ig] = (int)(0.5*(ncrs_stk-1.0) + 0.5);
      prseg_out[ig].nstk = 1;
      ixend = ncrs_stk;
      while(ixend < prseg_in[ig].nstk)
	 {
         prseg_out[ig].nstk++;
         ixend = ncrs_stk*prseg_out[ig].nstk;
	 }

      dip_off[ig] = (int)(0.5*(ncrs_dip-1.0) + 0.5);
      prseg_out[ig].ndip = 1;
      iyend = ncrs_dip;
      while(iyend < prseg_in[ig].ndip)
	 {
         prseg_out[ig].ndip++;
         iyend = ncrs_dip*prseg_out[ig].ndip;
	 }

      /*
      prseg_out[ig].nstk = (int)(1.0*prseg_in[ig].nstk/ncrs_stk + 0.5);
      while(prseg_out[ig].nstk*ncrs_stk > prseg_in[ig].nstk)
         prseg_out[ig].nstk--;

      prseg_out[ig].ndip = (int)(1.0*prseg_in[ig].ndip/ncrs_dip + 0.5);
      while(prseg_out[ig].ndip*ncrs_dip > prseg_in[ig].ndip)
         prseg_out[ig].ndip--;
      */


      prseg_out[ig].flen = prseg_out[ig].nstk*ncrs_stk*(prseg_in[ig].flen/prseg_in[ig].nstk);
      prseg_out[ig].fwid = prseg_out[ig].ndip*ncrs_dip*(prseg_in[ig].fwid/prseg_in[ig].ndip);

      prseg_out[ig].stk = prseg_in[ig].stk;
      prseg_out[ig].dip = prseg_in[ig].dip;
      prseg_out[ig].dtop = prseg_in[ig].dtop;
      prseg_out[ig].shyp = prseg_in[ig].shyp;
      prseg_out[ig].dhyp = prseg_in[ig].dhyp;

      nstk_in[ig] = prseg_in[ig].nstk;
      ndip_in[ig] = prseg_in[ig].ndip;
      nstk_out[ig] = prseg_out[ig].nstk;
      ndip_out[ig] = prseg_out[ig].ndip;
      }
   }

mrf[0].nseg = srf[0].nseg;
mrf[0].np_seg = (int *)check_malloc((mrf[0].nseg)*sizeof(int));

mrf[0].srf_apnts.np = 0;
for(ig=0;ig<mrf[0].nseg;ig++)
   {
   mrf[0].np_seg[ig] = nstk_out[ig]*ndip_out[ig];
   mrf[0].srf_apnts.np = mrf[0].srf_apnts.np + mrf[0].np_seg[ig];
   }

mrf[0].srf_apnts.apntvals = (struct srf_apointvalues *)check_malloc((mrf[0].srf_apnts.np)*sizeof(struct srf_apointvalues));

apval_in = srf[0].srf_apnts.apntvals;
apval_out = mrf[0].srf_apnts.apntvals;

mr_nnD = NULL;
mr_eeD = NULL;
mr_ddD = NULL;
mr_neD = NULL;
mr_ndD = NULL;
mr_edD = NULL;
s_stfD = NULL;

ntmr_back = 0;
nts_back = 0;

ntot_in = 0;
ntot_out = 0;
for(ig=0;ig<mrf[0].nseg;ig++)
   {
   for(iy=0;iy<ndip_out[ig];iy++)
      {
      iyp = iy*ncrs_dip + dip_off[ig];
      for(ix=0;ix<nstk_out[ig];ix++)
         {
         ixp = ix*ncrs_stk + stk_off[ig];

         ip_in = ixp + iyp*nstk_in[ig] + ntot_in;
         ip_out = ix + iy*nstk_out[ig] + ntot_out;

/* use center point values for following parameters */

         apval_out[ip_out].lon = apval_in[ip_in].lon;
         apval_out[ip_out].lat = apval_in[ip_in].lat;
         apval_out[ip_out].dep = apval_in[ip_in].dep;
         apval_out[ip_out].area = apval_in[ip_in].area*(ncrs_stk*ncrs_dip);
         apval_out[ip_out].dt = apval_in[ip_in].dt;
         apval_out[ip_out].vp = apval_in[ip_in].vp;
         apval_out[ip_out].vs = apval_in[ip_in].vs;
         apval_out[ip_out].den = apval_in[ip_in].den;

         apval_out[ip_out].stk = apval_in[ip_in].stk;
         apval_out[ip_out].dip = apval_in[ip_in].dip;
         apval_out[ip_out].rake = apval_in[ip_in].rake;

/* determine starting indicies of input points */

         ix0 = ixp - (ncrs_stk-1)/2;
	 if(ix0 < 0)
	    ix0 = 0;

	 ixend = ix0 + ncrs_stk;
/* not sure if this 'if' works */
	 if(ixend > nstk_in[ig] || (ixend + stk_off[ig] + 1) > nstk_in[ig])
	    ixend = nstk_in[ig];

         iy0 = iyp - (ncrs_dip-1)/2;
	 if(iy0 < 0)
	    iy0 = 0;

	 iyend = iy0 + ncrs_dip;
/* not sure if this 'if' works */
	 if(iyend > ndip_in[ig] || (iyend + dip_off[ig] + 1) > ndip_in[ig])
	    iyend = ndip_in[ig];


/* determine minimum tinit and total ntmr needed for this group of input points */

/*
fprintf(stderr,"ix0= %d ixend= %d iy0= %d iyend= %d\n",ix0,ixend,iy0,iyend);
*/

	 tmin = 1.0e+15;
	 tmax = -1.0e+15;
         for(iys=iy0;iys<iyend;iys++)
            {
            for(ixs=ix0;ixs<ixend;ixs++)
               {
               ips = ixs + iys*nstk_in[ig] + ntot_in;
         
	       if(apval_in[ips].tinit < tmin)
	          tmin = apval_in[ips].tinit;
         
               nts = apval_in[ips].nt1;
               if(apval_in[ips].nt2 > nts)
                  nts = apval_in[ips].nt2;
               if(apval_in[ips].nt3 > nts)
                  nts = apval_in[ips].nt3;
         
	       tend = apval_in[ips].tinit + nts*apval_in[ips].dt;
	       if(tend > tmax)
	          tmax = tend;

               if(use_srf_lame == 0 || apval_in[ips].vp < 0.0 || apval_in[ips].vs < 0.0 || apval_in[ips].den < 0.0)
                  {
                  id = 0;
                  while(vm->dep[id] < apval_in[ips].dep)
                     id++;

                  if(apval_in[ips].vp < 0.0)
                     apval_in[ips].vp = 1.0e+05*vm->vp[id];

                  if(apval_in[ips].vs < 0.0)
                     apval_in[ips].vs = 1.0e+05*vm->vs[id];

                  if(apval_in[ips].den < 0.0)
                     apval_in[ips].den = vm->den[id];

                  if(ips == ip_in) /* reset output center point values */
                     {
                     apval_out[ip_out].vp = apval_in[ip_in].vp;
                     apval_out[ip_out].vs = apval_in[ip_in].vs;
                     apval_out[ip_out].den = apval_in[ip_in].den;
	             }
                  }
	       }
            }

         apval_out[ip_out].tinit = tmin;

         apval_out[ip_out].mnn = 0.0;
         apval_out[ip_out].mee = 0.0;
         apval_out[ip_out].mdd = 0.0;
         apval_out[ip_out].mne = 0.0;
         apval_out[ip_out].mnd = 0.0;
         apval_out[ip_out].med = 0.0;

         apval_out[ip_out].mr_nn = NULL;
         apval_out[ip_out].mr_ee = NULL;
         apval_out[ip_out].mr_dd = NULL;
         apval_out[ip_out].mr_ne = NULL;
         apval_out[ip_out].mr_nd = NULL;
         apval_out[ip_out].mr_ed = NULL;

         apval_out[ip_out].ntmr = (int)((tmax-tmin)/apval_out[ip_out].dt + 0.5);

/*
*/
fprintf(stderr,"ip_o= %d ip_i= %d ntmr= %d tmin= %.4f tmax= %.4f\n",ip_out,ip_in,apval_out[ip_out].ntmr,tmin,tmax);

	 if(apval_out[ip_out].ntmr > 0)
	    {
	    if(apval_out[ip_out].ntmr > ntmr_back)
	       {
               mr_nnD = (double *)check_realloc(mr_nnD,(apval_out[ip_out].ntmr)*sizeof(double));
               mr_eeD = (double *)check_realloc(mr_eeD,(apval_out[ip_out].ntmr)*sizeof(double));
               mr_ddD = (double *)check_realloc(mr_ddD,(apval_out[ip_out].ntmr)*sizeof(double));
               mr_neD = (double *)check_realloc(mr_neD,(apval_out[ip_out].ntmr)*sizeof(double));
               mr_ndD = (double *)check_realloc(mr_ndD,(apval_out[ip_out].ntmr)*sizeof(double));
               mr_edD = (double *)check_realloc(mr_edD,(apval_out[ip_out].ntmr)*sizeof(double));

	       ntmr_back = apval_out[ip_out].ntmr;
	       }

	    for(it=0;it<apval_out[ip_out].ntmr;it++)
	       {
	       mr_nnD[it] = 0.0;
	       mr_eeD[it] = 0.0;
	       mr_ddD[it] = 0.0;
	       mr_neD[it] = 0.0;
	       mr_ndD[it] = 0.0;
	       mr_edD[it] = 0.0;
	       }

/* compute moment-rate functions for each input point and sum across group */

            for(iys=iy0;iys<iyend;iys++)
               {
               for(ixs=ix0;ixs<ixend;ixs++)
                  {
                  ips = ixs + iys*nstk_in[ig] + ntot_in;
         
                  nts = apval_in[ips].nt1;
                  if(apval_in[ips].nt2 > nts)
                     nts = apval_in[ips].nt2;
                  if(apval_in[ips].nt3 > nts)
                     nts = apval_in[ips].nt3;

	          if(nts > nts_back)
	             {
                     s_stfD = (double *)check_realloc(s_stfD,nts*sizeof(double));
	             nts_back = nts;
	             }

	          for(it=0;it<nts;it++)
	             s_stfD[it] = 0.0;

                  if(apval_in[ips].nt1)
                     {
                     stfp = apval_in[ips].stf1;
                     for(it=0;it<(apval_in[ips].nt1);it++)
                        s_stfD[it] = s_stfD[it] + stfp[it]*stfp[it];
                     }

                  if(apval_in[ips].nt2)
                     {
                     stfp = apval_in[ips].stf2;
                     for(it=0;it<(apval_in[ips].nt2);it++)
                        s_stfD[it] = s_stfD[it] + stfp[it]*stfp[it];
                     }

                  if(apval_in[ips].nt3)
                     {
                     stfp = apval_in[ips].stf3;
                     for(it=0;it<(apval_in[ips].nt3);it++)
                        s_stfD[it] = s_stfD[it] + stfp[it]*stfp[it];
                     }

                  sumD = 0.0;
                  for(it=0;it<nts;it++)
                     {
                     s_stfD[it] = sqrt(s_stfD[it]);
                     sumD = sumD + (apval_in[ips].dt)*s_stfD[it];
                     }

                  sumD = 1.0/sumD;
                  for(it=0;it<nts;it++)
                     s_stfD[it] = sumD*s_stfD[it];

                  l2m = apval_in[ips].vp*apval_in[ips].vp*apval_in[ips].den;
                  mu = apval_in[ips].vs*apval_in[ips].vs*apval_in[ips].den;
                  lam = l2m - 2.0*mu;

                  u1 = apval_in[ips].slip1;
                  u2 = apval_in[ips].slip2;
                  u3 = apval_in[ips].slip3;

                  arg = apval_in[ips].stk*rperd;
                  cosS = cos(arg);
                  sinS = sin(arg);

                  arg = apval_in[ips].dip*rperd;
                  cosD = cos(arg);
                  sinD = sin(arg);

                  arg = apval_in[ips].rake*rperd;
                  cosL = cos(arg);
                  sinL = sin(arg);

                  vx = -sinD*sinS;
                  vy =  sinD*cosS;
                  vz = -cosD;

                  ux = -(u3*sinD - cosD*(u1*sinL + u2*cosL))*sinS + (u1*cosL - u2*sinL)*cosS;
                  uy =  (u3*sinD - cosD*(u1*sinL + u2*cosL))*cosS + (u1*cosL - u2*sinL)*sinS;
                  uz = -u3*cosD - (u1*sinL + u2*cosL)*sinD;

                  Mnn = (l2m*vx*ux + lam*(vy*uy + vz*uz))*apval_in[ips].area;
                  Mee = (l2m*vy*uy + lam*(vx*ux + vz*uz))*apval_in[ips].area;
                  Mdd = (l2m*vz*uz + lam*(vx*ux + vy*uy))*apval_in[ips].area;

                  Mne = mu*(vx*uy + vy*ux)*apval_in[ips].area;
                  Mnd = mu*(vx*uz + vz*ux)*apval_in[ips].area;
                  Med = mu*(vy*uz + vz*uy)*apval_in[ips].area;

/* sum fine grid moment-rates into coarse grid arrays */
         
	          it0 = (int)((apval_in[ips].tinit - apval_out[ip_out].tinit)/apval_in[ips].dt + 0.5);
	          if(it0 < 0)
	             it0 = 0;

	          while((it0+nts) > apval_out[ip_out].ntmr)
	             nts--;

	          for(it=0;it<nts;it++)
	             {
		     it_out = it + it0;

	             mr_nnD[it_out] = mr_nnD[it_out] + Mnn*s_stfD[it];
	             mr_eeD[it_out] = mr_eeD[it_out] + Mee*s_stfD[it];
	             mr_ddD[it_out] = mr_ddD[it_out] + Mdd*s_stfD[it];
	             mr_neD[it_out] = mr_neD[it_out] + Mne*s_stfD[it];
	             mr_ndD[it_out] = mr_ndD[it_out] + Mnd*s_stfD[it];
	             mr_edD[it_out] = mr_edD[it_out] + Med*s_stfD[it];
	             }

                  s_mom = s_mom + sqrt(u1*u1 + u2*u2 + u3*u3)*mu*apval_in[ips].area;
                  m_mom_f = m_mom_f + sqrt(0.5*(Mnn*Mnn + Mee*Mee + Mdd*Mdd) + Mne*Mne + Mnd*Mnd + Med*Med);
	          }
               }

/* copy individual moment-rates to output structure */

            apval_out[ip_out].mr_nn = (float *)check_realloc(apval_out[ip_out].mr_nn,(apval_out[ip_out].ntmr)*sizeof(float));
            stfp = apval_out[ip_out].mr_nn;
            for(it=0;it<apval_out[ip_out].ntmr;it++)
	       {
	       stfp[it] = mr_nnD[it];
               apval_out[ip_out].mnn = apval_out[ip_out].mnn + (apval_out[ip_out].dt)*mr_nnD[it];
	       }

            apval_out[ip_out].mr_ee = (float *)check_realloc(apval_out[ip_out].mr_ee,(apval_out[ip_out].ntmr)*sizeof(float));
            stfp = apval_out[ip_out].mr_ee;
            for(it=0;it<apval_out[ip_out].ntmr;it++)
	       {
	       stfp[it] = mr_eeD[it];
               apval_out[ip_out].mee = apval_out[ip_out].mee + (apval_out[ip_out].dt)*mr_eeD[it];
	       }

            apval_out[ip_out].mr_dd = (float *)check_realloc(apval_out[ip_out].mr_dd,(apval_out[ip_out].ntmr)*sizeof(float));
            stfp = apval_out[ip_out].mr_dd;
            for(it=0;it<apval_out[ip_out].ntmr;it++)
	       {
	       stfp[it] = mr_ddD[it];
               apval_out[ip_out].mdd = apval_out[ip_out].mdd + (apval_out[ip_out].dt)*mr_ddD[it];
	       }

            apval_out[ip_out].mr_ne = (float *)check_realloc(apval_out[ip_out].mr_ne,(apval_out[ip_out].ntmr)*sizeof(float));
            stfp = apval_out[ip_out].mr_ne;
            for(it=0;it<apval_out[ip_out].ntmr;it++)
	       {
	       stfp[it] = mr_neD[it];
               apval_out[ip_out].mne = apval_out[ip_out].mne + (apval_out[ip_out].dt)*mr_neD[it];
	       }

            apval_out[ip_out].mr_nd = (float *)check_realloc(apval_out[ip_out].mr_nd,(apval_out[ip_out].ntmr)*sizeof(float));
            stfp = apval_out[ip_out].mr_nd;
            for(it=0;it<apval_out[ip_out].ntmr;it++)
	       {
	       stfp[it] = mr_ndD[it];
               apval_out[ip_out].mnd = apval_out[ip_out].mnd + (apval_out[ip_out].dt)*mr_ndD[it];
	       }

            apval_out[ip_out].mr_ed = (float *)check_realloc(apval_out[ip_out].mr_ed,(apval_out[ip_out].ntmr)*sizeof(float));
            stfp = apval_out[ip_out].mr_ed;
            for(it=0;it<apval_out[ip_out].ntmr;it++)
	       {
	       stfp[it] = mr_edD[it];
               apval_out[ip_out].med = apval_out[ip_out].med + (apval_out[ip_out].dt)*mr_edD[it];
	       }

            m_mom_c = m_mom_c + sqrt(0.5*apval_out[ip_out].mnn*apval_out[ip_out].mnn
                                   + 0.5*apval_out[ip_out].mee*apval_out[ip_out].mee
                                   + 0.5*apval_out[ip_out].mdd*apval_out[ip_out].mdd
                                   + apval_out[ip_out].mne*apval_out[ip_out].mne
                                   + apval_out[ip_out].mnd*apval_out[ip_out].mnd
                                   + apval_out[ip_out].med*apval_out[ip_out].med);
	    }
         }
      }

   ntot_in = ntot_in + srf[0].np_seg[ig];
   ntot_out = ntot_out + mrf[0].np_seg[ig];
   }

free(mr_nnD);
free(mr_eeD);
free(mr_ddD);
free(mr_neD);
free(mr_ndD);
free(mr_edD);
free(s_stfD);

fprintf(stderr,"slip moment= %.5e\n",s_mom);
fprintf(stderr,"MT moment (fine grid)= %.5e\n",m_mom_f);
fprintf(stderr,"MT moment (coarse grid)= %.5e\n",m_mom_c);
}
