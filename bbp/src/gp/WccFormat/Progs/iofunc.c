#include        "include.h"
#include        "structure.h"
#include        "function.h"

float *read_wccseis(char *ifile,struct statdata *shead,float *s,int bflag)
{
FILE *fpr;
int fdr, nt6, i, j;
char header1[1024];

shead->hr = 0;
shead->min = 0;
shead->sec = 0;
shead->edist = 0;
shead->az = 0;
shead->baz = 0;

if(bflag)
   {
   if(strcmp(ifile,"stdin") == 0)
      fdr = STDIN_FILENO;
   else
      fdr = opfile_ro(ifile);

   reed(fdr,shead,sizeof(struct statdata));
   s = (float *) check_realloc(s,shead->nt*sizeof(float));
   reed(fdr,s,shead->nt*sizeof(float));
   close(fdr);
   }
else
   {
   if(strcmp(ifile,"stdin") == 0)
      fpr = stdin;
   else
      fpr = fopfile(ifile,"r");

   fgets(header1,1024,fpr);
   getheader(header1,shead);

   fgets(header1,1024,fpr);
   sscanf(header1,"%d %f %d %d %f %f %f %f",&shead->nt,
                                            &shead->dt,
                                            &shead->hr,
                                            &shead->min,
                                            &shead->sec,
                                            &shead->edist,
                                            &shead->az,
                                            &shead->baz);

   s = (float *) check_realloc(s,shead->nt*sizeof(float));

   for(i=0;i<shead->nt;i++)
      fscanf(fpr,"%f",&s[i]);

   fclose(fpr);
   }

return(s);
}

void getheader(char *str,struct statdata *hd)
{
int i;

i = 0;
while(str[i] != ' ' && str[i] != '\t' && str[i] != '\0')
   i++;

if(i < STATCHAR-1)
   {
   strncpy(hd->stat,str,i);
   hd->stat[i] = '\0';
   }
else
   {
   strncpy(hd->stat,str,STATCHAR-1);
   hd->stat[STATCHAR-1] = '\0';
   }

while(str[i] == ' ' || str[i] == '\t')
   i++;
str = str + i;

i = 0;
while(str[i] != ' ' && str[i] != '\t' && str[i] != '\n' && str[i] != '\0')
   i++;

if(i < COMPCHAR-1)
   {
   strncpy(hd->comp,str,i);
   hd->comp[i] = '\0';
   }
else
   {
   strncpy(hd->comp,str,COMPCHAR-1);
   hd->comp[COMPCHAR-1] = '\0';
   }

while(str[i] == ' ' || str[i] == '\t')
   i++;
str = str + i;
 
i = 0;
while(str[i] != '\n' && str[i] != '\0')
   i++;

if(i < TITLCHAR-1)
   {
   if(str[i] != '\0')
      str[i] = '\0';
   }
else
   str[TITLCHAR-1] = '\0';
 
strcpy(hd->stitle,str);
}

void write_wccseis(char *ofile,struct statdata *shead,float *s,int bflag)
{
FILE *fpw;
int fdw, nt6, i, j;

if(bflag)
   {
   if(strcmp(ofile,"stdout") == 0)
      fdw = STDOUT_FILENO;
   else
      fdw = croptrfile(ofile);

   rite(fdw,shead,sizeof(struct statdata));
   rite(fdw,s,shead->nt*sizeof(float));
   close(fdw);
   }
else
   {
   if(strcmp(ofile,"stdout") == 0)
      fpw = stdout;
   else
      fpw = fopfile(ofile,"w");

   fprintf(fpw,"%-10s %3s %s\n",shead->stat,shead->comp,shead->stitle);
   fprintf(fpw,"%d %12.5e %d %d %12.5e %12.5e %12.5e %12.5e\n",shead->nt,
                                           shead->dt,
                                           shead->hr,
                                           shead->min,
                                           shead->sec,
                                           shead->edist,
                                           shead->az,
                                           shead->baz);

   nt6 = shead->nt/6;                                   
   for(i=0;i<nt6;i++)
      {
      for(j=0;j<6;j++)
         fprintf(fpw,"%13.5e",s[6*i + j]);

      fprintf(fpw,"\n");                   
      }

   if(6*nt6 != shead->nt)
      {
      for(i=6*nt6;i<shead->nt;i++)
         fprintf(fpw,"%13.5e",s[i]);

      fprintf(fpw,"\n");             
      }
   fclose(fpw);
   }
}

FILE *fopfile(char *name,char *mode)
{
FILE *fp;

if((fp = fopen(name,mode)) == NULL)
   {
   fprintf(stderr,"CAN'T FOPEN FILE = %s, MODE = %s\n", name, mode);
   exit(-1);
   }
return(fp);
}

int opfile_ro(char *name)
{
int fd;
if ((fd = open (name, RDONLY_FLAGS, 0444)) == -1)
   fprintf (stderr, "CAN'T OPEN FILE %s\n", name);
return (fd);
}

int opfile(char *name)
{
int fd;
if ((fd = open (name, RDWR_FLAGS, 0664)) == -1)
   fprintf (stderr, "CAN'T OPEN FILE %s\n", name);
return (fd);
}

int croptrfile(char *name)
{
int fd;
if ((fd = open (name, CROPTR_FLAGS, 0664)) == -1)
   fprintf (stderr, "CAN'T OPEN FILE %s\n", name);
return (fd);
}

int reed(int fd, void *pntr, int length)
{
int temp;
if ((temp = read(fd, pntr, length)) < length)
   {
   fprintf (stderr, "READ ERROR\n");
   fprintf (stderr, "%d attempted  %d read\n", length, temp);
   exit(-1);
   }
return(temp);
}

int rite(int fd, void *pntr, int length)
{
int temp;
if ((temp = write(fd, pntr, length)) < length)
   {
   fprintf (stderr, "WRITE ERROR\n");
   fprintf (stderr, "%d attempted  %d written\n", length, temp);
   exit(-1);
   }
return(temp);
}

void *check_realloc(void *ptr,size_t len)
{
ptr = (char *) realloc (ptr,len);

if(ptr == NULL)
   {
   fprintf(stderr,"*****  memory reallocation error\n");
   exit(-1);
   }

return(ptr);
}

void *check_malloc(size_t len)
{
char *ptr;

ptr = (char *) malloc (len);
 
if(ptr == NULL)
   {
   fprintf(stderr,"*****  memory allocation error\n");
   exit(-1);
   }
 
return(ptr);
}
