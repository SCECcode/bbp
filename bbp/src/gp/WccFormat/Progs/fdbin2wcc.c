#include        "include.h"
#include        "structure.h"
#include        "function.h"

int size_float = sizeof(float);
float float_swap(char *);

main(int ac, char **av)
{
FILE *fpr, *fopfile();
struct statdata shead;
struct tsheader tshead;
struct seisheader seishead;
struct mtheader mthead;
float *s;
int i, ix;

int nx = -1;

int swap_bytes = 0;
char cbuf[512];

float scale = 1.0;
int flip = 1;

int read_tsheader = 0;

int outbin = 0;
int momten = 0;
int mtindx;
off_t off, cur_off;
int ngf;

int all_in_one = 0;
int seisindx = -1;
int notfound = 1;
int search_by_stat = 1;
int dcomp;

int fd;
char infile[128];
char filelist[128];
char outfile[128];

int itshft, itend;
int tst2zero = 0;
float tpadfront = 0.0;

shead.stitle[0] = '\0';
shead.stat[0] = '\0';
shead.comp[0] = '\0';
shead.nt = -1;
shead.dt = -1;
shead.hr = 0;
shead.min = 0;
shead.sec = 0.0;
shead.edist = 0.0;
shead.az = 0.0;
shead.baz = 0.0;

setpar(ac, av);

getpar("all_in_one","d",&all_in_one);
getpar("read_tsheader","d",&read_tsheader);

if(all_in_one)
   {
   mstpar("filelist","s",filelist);
   getpar("seisindx","d",&seisindx);

   search_by_stat = 1;
   if(seisindx >= 0)
      search_by_stat = 0;

   if(search_by_stat)
      mstpar("stat","s",shead.stat);
   else
      getpar("stat","s",shead.stat);

   getpar("nt","d",&shead.nt);
   getpar("comp","s",shead.comp);
   }
else if(read_tsheader)
   {
   mstpar("infile","s",infile);
   mstpar("stat","s",shead.stat);
   mstpar("comp","s",shead.comp);

   getpar("nx","d",&nx);
   getpar("nt","d",&shead.nt);
   getpar("dt","f",&shead.dt);
   }
else
   {
   mstpar("infile","s",infile);
   mstpar("nx","d",&nx);
   mstpar("nt","d",&shead.nt);
   mstpar("dt","f",&shead.dt);
   mstpar("stat","s",shead.stat);
   mstpar("comp","s",shead.comp);
   }

mstpar("ix","d",&ix);
mstpar("outfile","s",outfile);

getpar("title","s",shead.stitle);
getpar("epi","f",&shead.edist);
getpar("tst","f",&shead.sec);
getpar("tpadfront","f",&tpadfront);
getpar("tst2zero","d",&tst2zero);
getpar("scale","f",&scale);
getpar("flip","d",&flip);
getpar("outbin","d",&outbin);
getpar("swap_bytes","d",&swap_bytes);
getpar("momten","d",&momten);
if(momten)
   mstpar("mtindx","d",&mtindx);

endpar();

/*  read in input data filenames */

if(momten == 0 && all_in_one == 0)
   {
   fd = opfile_ro(infile);

   off = (off_t)(ix*size_float);
   if(read_tsheader == 1)
      {
      reed(fd,&tshead,sizeof(struct tsheader));

      nx = 3*tshead.nx*tshead.ny*tshead.nz;
      if(shead.nt < 0 || shead.nt > tshead.nt)
         shead.nt = tshead.nt;
      if(shead.dt < 0)
         shead.dt = tshead.dt;

      lseek(fd,off,SEEK_CUR);
      }
   else
      lseek(fd,off,SEEK_SET);

   itshft = 0;
   itend = shead.nt;
   if(tst2zero == 1)
      {
      itshft = (int)(shead.sec/shead.dt);
      shead.nt = itshft + shead.nt;
      shead.sec = 0.0;
      }
 
   s = (float *) check_malloc (shead.nt*size_float);

   if(itshft < 0)
      {
      off = (off_t)(-itshft*nx*size_float);
      lseek(fd,off,SEEK_CUR);
      itshft = 0;
      itend = shead.nt;
      }
   else if(itshft > 0)
      {
      for(i=0;i<itshft;i++)
         s[i] = 0.0;

      itend = shead.nt - itshft;
      }

   for(i=0;i<itend;i++)
      {
      reed(fd,&s[i+itshft],size_float);
      if(swap_bytes)
         swap_in_place(1,(char *)(&s[i+itshft]));

      off = (off_t)((nx-1)*size_float);
      lseek(fd,off,SEEK_CUR);
      s[i+itshft] = scale*flip*s[i+itshft];
      }
   }
else if(momten == 1)
   {
   fd = opfile_ro(infile);

   reed(fd,&ngf,sizeof(int));
   reed(fd,&mthead,sizeof(struct mtheader));
   shead.nt = mthead.nt;
   shead.dt = mthead.dt;

   itshft = 0;
   itend = shead.nt;
   if(tst2zero == 1)
      {
      itshft = (int)(shead.sec/shead.dt);
      shead.nt = itshft + shead.nt;
      shead.sec = 0.0;
      }

   s = (float *) check_malloc (shead.nt*size_float);

   off = sizeof(int) + (off_t)(mtindx)*sizeof(struct mtheader);
   lseek(fd,off,0);
   reed(fd,&mthead,sizeof(struct mtheader));

   fprintf(stderr,"Moment tensor header information-\n");
   fprintf(stderr,"   Title: %s\n",mthead.title);
   fprintf(stderr,"   FD model rotation=%.1f deg. clockwise from North\n",mthead.modelrot);
   fprintf(stderr,"   source type: xmom=%13.5e dyne-cm\n",mthead.xmom);
   fprintf(stderr,"                ymom=%13.5e dyne-cm\n",mthead.ymom);
   fprintf(stderr,"                zmom=%13.5e dyne-cm\n",mthead.zmom);
   fprintf(stderr,"   source location:  (x,y,z)=(%d,%d,%d)\n",mthead.xsrc,mthead.ysrc,mthead.zsrc);
   fprintf(stderr,"   station index=%5d\n",mthead.indx);
   fprintf(stderr,"   station location: (x,y,z)=(%d,%d,%d)\n",mthead.xsta,mthead.ysta,mthead.zsta);

/* now read GF for 'mtindx' */

   off = sizeof(int) + (off_t)(ngf*sizeof(struct mtheader) + 6*mtindx*size_float);
   lseek(fd,off,SEEK_SET);

   off = (off_t)(ix*size_float);
   lseek(fd,ix*size_float,SEEK_CUR);
   for(i=0;i<shead.nt;i++)
      {
      reed(fd,&s[i],size_float);
      lseek(fd,5*size_float,1);
      s[i] = scale*flip*s[i];
      }
   }
else if(all_in_one == 1)
   {
   fpr = fopfile(filelist,"r");

   notfound = 1;
   while(fscanf(fpr,"%s",infile) != EOF && notfound)
      {
      fd = opfile_ro(infile);

      reed(fd,&ngf,sizeof(int));
      if(swap_bytes)
	 swap_in_place(1,(char *)(&ngf));

      i = 0;
      while(i<ngf && notfound)
	 {
         reed(fd,&seishead,sizeof(struct seisheader));

         if(swap_bytes)
            {
	    swap_in_place(1,(char *)(&seishead.indx));
	    swap_in_place(1,(char *)(&seishead.ix));
	    swap_in_place(1,(char *)(&seishead.iy));
	    swap_in_place(1,(char *)(&seishead.iz));
	    swap_in_place(1,(char *)(&seishead.nt));
	    swap_in_place(1,(char *)(&seishead.dt));
	    swap_in_place(1,(char *)(&seishead.h));
	    swap_in_place(1,(char *)(&seishead.modelrot));
	    swap_in_place(1,(char *)(&seishead.modellat));
	    swap_in_place(1,(char *)(&seishead.modellon));
	    }

/*
fprintf(stderr,"name= %s\n",seishead.name);
*/

	 if(search_by_stat)
	    {
	    if(strcmp(shead.stat,seishead.name) == 0)
	       {
	       seisindx = i;
	       notfound = 0;
	       }
	    }
	 else
	    {
	    if(seishead.indx == seisindx)
	       {
	       seisindx = i;
	       notfound = 0;
	       }
	    }
	 i++;
	 }

      if(notfound)
	 close(fd);
      }

   if(notfound)
      {
      fprintf(stderr,"Unable to find time history, exiting...\n");
      exit(-1);
      }
   else
      fprintf(stderr,"Found stat= %s in file %s\n\n",seishead.name,infile);

   if(shead.stat[0] == '\0')
      strcpy(shead.stat,seishead.name);
   if(shead.comp[0] == '\0')
      {
      if(ix == 2)
	 {
	 if(flip == 1)
	    sprintf(shead.comp,"dwn");
	 else
	    sprintf(shead.comp,"ver");
	 }
      else
	 {
         if(ix == 0)
	    dcomp = 90 + seishead.modelrot;
         if(ix == 1)
	    dcomp = 180 + seishead.modelrot;

	 if(flip == -1)
	    dcomp = dcomp + 180;

	 while(dcomp < 0)
	    dcomp = dcomp + 360;
	 while(dcomp >= 360)
	    dcomp = dcomp - 360;

	 sprintf(shead.comp,"%d",dcomp);
	 }
      }
   if(shead.nt < 0 || shead.nt > seishead.nt)
      shead.nt = seishead.nt;

   shead.dt = seishead.dt;

   itshft = 0;
   itend = shead.nt;
   if(tst2zero == 1)
      {
      itshft = (int)(shead.sec/shead.dt);
      shead.nt = itshft + shead.nt;
      shead.sec = 0.0;
      }

   s = (float *) check_malloc (shead.nt*size_float);

   fprintf(stderr,"Seismogram File header information-\n");
   fprintf(stderr,"   station name= %s\n",seishead.name);
   fprintf(stderr,"   station index=%5d\n",seishead.indx);
   fprintf(stderr,"   station location: (x,y,z)=(%d,%d,%d)\n",seishead.ix,seishead.iy,seishead.iz);

/* now read seismo for 'seisindx' */

   off = sizeof(int) + (off_t)(ngf*sizeof(struct seisheader) + 3*seisindx*size_float);
   lseek(fd,off,SEEK_SET);

   off = (off_t)(ix*size_float);
   lseek(fd,off,SEEK_CUR);

   if(itshft < 0)
      {
      off = (off_t)(-itshft*3*ngf*size_float);
      lseek(fd,off,SEEK_CUR);
      itshft = 0;
      itend = shead.nt;
      }
   else if(itshft > 0)
      {
      for(i=0;i<itshft;i++)
         s[i] = 0.0;

      itend = shead.nt - itshft;
      }

   off = (off_t)((3*ngf - 1)*size_float);
   for(i=0;i<itend;i++)
      {
      reed(fd,&s[i+itshft],size_float);
      if(swap_bytes)
         swap_in_place(1,(char *)(&s[i+itshft]));

      lseek(fd,off,SEEK_CUR);
      s[i+itshft] = scale*flip*s[i+itshft];
      }
   }

close(fd);

if(tpadfront > 0.0)
   {
   itshft = (int)(tpadfront/shead.dt + 0.5);
   shead.nt = shead.nt + itshft;
   shead.sec = shead.sec - itshft*shead.dt;

   s = (float *) check_realloc(s,shead.nt*sizeof(float));

   for(i=shead.nt-1;i>=itshft;i--)
      s[i] = s[i-itshft];
   for(i=0;i<itshft;i++)
      s[i] = 0.0;
   }

write_wccseis(outfile,&shead,s,outbin);
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
