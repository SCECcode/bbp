#include        "include.h"
#include        "structure.h"
#include        "function.h"

FILE *fopfile(char *name,char *mode)
{
FILE *fp;

if((fp = fopen(name,mode)) == NULL)
   {
   fprintf(stderr,"CAN'T FOPEN FILE = %s, MODE = %s\n", name, mode);

   fprintf(stderr,"errno= %d\n",errno);

   if(errno == EACCES)
      fprintf(stderr,"EACCES: errno= %d\n",errno);
   if(errno == EEXIST)
      fprintf(stderr,"EEXIST: errno= %d\n",errno);
   if(errno == EFAULT)
      fprintf(stderr,"EFAULT: errno= %d\n",errno);
   if(errno == EISDIR)
      fprintf(stderr,"EISDIR: errno= %d\n",errno);
   if(errno == ELOOP)
      fprintf(stderr,"ELOOP: errno= %d\n",errno);
   if(errno == EMFILE)
      fprintf(stderr,"EMFILE: errno= %d\n",errno);
   if(errno == ENAMETOOLONG)
      fprintf(stderr,"ENAMETOOLONG: errno= %d\n",errno);
   if(errno == ENFILE)
      fprintf(stderr,"ENFILE: errno= %d\n",errno);
   if(errno == ENODEV)
      fprintf(stderr,"ENODEV: errno= %d\n",errno);
   if(errno == ENOENT)
      fprintf(stderr,"ENOENT: errno= %d\n",errno);
   if(errno == ENOMEM)
      fprintf(stderr,"ENOMEM: errno= %d\n",errno);
   if(errno == ENOSPC)
      fprintf(stderr,"ENOSPC: errno= %d\n",errno);
   if(errno == ENOTDIR)
      fprintf(stderr,"ENOTDIR: errno= %d\n",errno);
   if(errno == ENXIO)
      fprintf(stderr,"ENXIO: errno= %d\n",errno);
   if(errno == EOVERFLOW)
      fprintf(stderr,"EOVERFLOW: errno= %d\n",errno);
   if(errno == EPERM)
      fprintf(stderr,"EPERM: errno= %d\n",errno);
   if(errno == EROFS)
      fprintf(stderr,"EROFS: errno= %d\n",errno);
   if(errno == ETXTBSY)
      fprintf(stderr,"ETXTBSY: errno= %d\n",errno);
   if(errno == EWOULDBLOCK)
      fprintf(stderr,"EWOULDBLOCK: errno= %d\n",errno);

   exit(-1);
   }
return(fp);
}

int opfile_ro(char *name)
{
int fd;
if ((fd = open (name, O_RDONLY, 0444)) == -1)
   {
   fprintf (stderr, "CAN'T OPEN FILE %s\n", name);
   exit(-1);
   }
return (fd);
}

int opfile(char *name)
{
int fd;
if ((fd = open (name, O_RDWR, 0664)) == -1)
   {
   fprintf (stderr, "CAN'T OPEN FILE %s\n", name);
   exit(-1);
   }
return (fd);
}

int croptrfile(char *name)
{
int fd;
if ((fd = open (name, O_CREAT | O_TRUNC | O_RDWR, 0664)) == -1)
   {
   fprintf (stderr, "CAN'T OPEN FILE %s\n", name);
   exit(-1);
   }
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

void fortran_rite(int fd,int nargs, ...)
{
va_list ap;
void *ptr[MAX_VAR_LIST];
int len[MAX_VAR_LIST];
int totlen = 0;
int i;

va_start(ap,nargs);
for(i=0;i<nargs;i++)
   {
   ptr[i] = va_arg(ap,void *);
   len[i] = va_arg(ap,int);
   totlen = totlen + len[i];
   }
va_end(ap);

rite(fd,&totlen,sizeof(int));

for(i=0;i<nargs;i++)
   rite(fd,ptr[i],len[i]);

rite(fd,&totlen,sizeof(int));
}

void write_seis(char *dir,char *stat,char *sname,char *comp,float *st,float *dt,int nt,float *ts)
{
FILE *fopfile(), *fpw;
int i, j, nt6;
char outfile[2048], header[128], stitle[128];

sprintf(outfile,"%s/%s.%s",dir,stat,comp);

fpw = fopfile(outfile,"w");

sprintf(stitle,"%-10s%3s %s\n",sname,comp,"TITLE");
fprintf(fpw,"%s",stitle);

sprintf(header,"%d %12.5e 0 0 %12.5e 0.0 0.0 0.0\n",nt,*dt,*ts);
fprintf(fpw,"%s",header);

nt6 = nt/6;
for(i=0;i<nt6;i++)
   {
   for(j=0;j<6;j++)
      fprintf(fpw,"%13.5e",st[6*i + j]);

   fprintf(fpw,"\n");
   }
if(6*nt6 != nt)
   {
   for(i=6*nt6;i<nt;i++)
      fprintf(fpw,"%13.5e",st[i]);

   fprintf(fpw,"\n");
   }

fclose(fpw);
}

char *skipval(int j,char *str)
{
while(str[0] == ' ' || str[0] == '\t' || str[0] == '\b' || str[0] == '\n')
   str++;

while(j--)
   {
   while(str[0] != ' ' && str[0] != '\t' && str[0] != '\b' && str[0] != '\n')
      str++;

   while(str[0] == ' ' || str[0] == '\t' || str[0] == '\b' || str[0] == '\n')
      str++;
   }

return(str);
}

void makedir(char *ipath)
{
struct stat sbuf;
char stmp[2048], str[2048], path[2048];
int rtn, j;
mode_t mode = 00777;

strcpy(path,ipath);

j = 0;
while(path[j] != '\0')
   j++;

j--;
while(path[j] == '/')
   j--;
path[j+1] = '\0';
path[j+2] = '\0';

j = 0;
while(path[j] != '\0')
   {
   while(path[j] != '/' && path[j] != '\0')
      j++;

   if(j != 0)
      {
      strncpy(stmp,path,j);
      stmp[j] = '\0';

      rtn = stat(stmp,&sbuf); /* stat directory path to see if it already exists */

      if(rtn == -1 && errno == ENOENT) /* try to make the directory path */
         {
         rtn = mkdir(stmp,mode);

         if(rtn == -1)
            {
            if(errno != EEXIST)
               {
               sprintf(str,"makedir() cannot make directory %s, exiting",stmp);
               perror(str);
               exit(-1);
               }
            }
         }

      else if(rtn == -1 && errno != ENOENT) /* some other problem */
         {
         sprintf(str,"problem with stat() on %s, exiting",stmp);
         perror(str);
         exit(-1);
         }
      }
   j++;
   }

/*
   Double check to make sure directory exists.  This is a brute-force
   method, but I ran inot problems with automounted directories using the
   error-checking above.  RWG 9/20/99
*/

rtn = mkdir(stmp,mode);
if(rtn == -1 && errno != EEXIST)
   {
   sprintf(str,"makedir() cannot make directory %s, exiting",stmp);
   perror(str);
   exit(-1);
   }
}
