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

                  apval_ptr->vs = -1;
                  apval_ptr->den = -1;
	          }
               else if(atof(srf->version) >= 2.0)
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
	          }

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
else if(atof(srf->version) >= 2.0)
   write_srf2(srf,file,bflag);
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
