#include        "include.h"
#include        "structure.h"
#include        "function.h"

int size_float = sizeof(float);
float float_swap(char *);

main(int ac, char **av)
{
FILE *fpr, *fopfile();
struct seisheader seishead;
int fd, i, j, k, n, ngf[10000], nfile;

int swap_bytes = 0;

char infile[128];
char filelist[128];
char *sptr1, *sptr2, *statlist[10000];

setpar(ac, av);
mstpar("filelist","s",filelist);
getpar("swap_bytes","d",&swap_bytes);
endpar();

fpr = fopfile(filelist,"r");

j = 0;
while(fscanf(fpr,"%s",infile) != EOF)
   {
   fd = opfile_ro(infile);

   reed(fd,&ngf[j],sizeof(int));
   if(swap_bytes)
      swap_in_place(1,(char *)(&ngf[j]));

   statlist[j] = check_malloc (ngf[j]*STATCHAR*sizeof(char));
   sptr1 = statlist[j];

   for(i=0;i<ngf[j];i++)
      {
      reed(fd,&seishead,sizeof(struct seisheader));
      strcpy(sptr1+i*STATCHAR,seishead.name);
      }

   close(fd);
   j++;
   }

nfile = j;

for(i=0;i<nfile;i++)
   {
   sptr1 = statlist[i];
   fprintf(stderr,"file %d of %d\n",i+1,nfile);
   for(j=i+1;j<nfile;j++)
      {
      sptr2 = statlist[j];

      for(k=0;k<ngf[i];k++)
         {
         for(n=0;n<ngf[j];n++)
            {
	    if(strcmp(sptr1+k*STATCHAR,sptr2+n*STATCHAR) == 0)
	       fprintf(stderr,"STAT= %s, FILES %d and %d\n",sptr1+k*STATCHAR,i,j);
	    }
	 }

      }
   }
}

long long_swap(char *cbuf)
{
union
   {
   char cval[4];
   long lval;
   } l_union;

l_union.cval[3] = cbuf[0];
l_union.cval[2] = cbuf[1];
l_union.cval[1] = cbuf[2];
l_union.cval[0] = cbuf[3];

return(l_union.lval);
}

float float_swap(char *cbuf)
{
union
   {
   char cval[4];
   float fval;
   } f_union;

f_union.cval[3] = cbuf[0];
f_union.cval[2] = cbuf[1];
f_union.cval[1] = cbuf[2];
f_union.cval[0] = cbuf[3];

return(f_union.fval);
}

void swap_in_place(int n,char *cbuf)
{
char cv;

while(n--)
   {
   cv = cbuf[0];
   cbuf[0] = cbuf[3];
   cbuf[3] = cv;

   cv = cbuf[1];
   cbuf[1] = cbuf[2];
   cbuf[2] = cv;

   cbuf = cbuf + 4;
   }
}
