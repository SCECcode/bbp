#include "include.h"
#include "structure.h"
#include "function.h"
#include "defs.h"
#include "getpar.h"

int main(int ac,char **av)
{
int nseg;
char infile[4096], type[4096], outfile[4096];

struct standrupformat srf;

int inbin = 0;
int calc_xy = 0;

float xoff = 0.0;
float yoff = 0.0;

float svmin = 0.0;
float slipmin = -1.0;

int tsstr = 0;
int tsend = -1;
int tsinc = 1;

int keepsign = 0;
int dump_slip = 0;

int lonlatdep = 0;
float depmin = -1.0e+15;
float depmax = 1.0e+15;

sprintf(infile,"stdin");
sprintf(outfile,"stdout");
sprintf(type,"slip");
nseg = 0;

setpar(ac,av);
getpar("infile","s",infile);
getpar("outfile","s",outfile);
getpar("type","s",type);
getpar("nseg","d",&nseg);
getpar("inbin","d",&inbin);
getpar("calc_xy","d",&calc_xy);
getpar("keepsign","d",&keepsign);
getpar("xoff","f",&xoff);
getpar("yoff","f",&yoff);
getpar("tsstr","d",&tsstr);
getpar("tsend","d",&tsend);
getpar("tsinc","d",&tsinc);
getpar("svmin","f",&svmin);
getpar("slipmin","f",&slipmin);

getpar("lonlatdep","d",&lonlatdep);
getpar("depmin","f",&depmin);
getpar("depmax","f",&depmax);

getpar("dump_slip","d",&dump_slip);
endpar();

read_srf(&srf,infile,inbin);

if(lonlatdep)
   write_lld(outfile,&srf,nseg,&depmin,&depmax,type);
else
   write_xyz(outfile,&srf,type,nseg,calc_xy,&xoff,&yoff,tsstr,tsend,tsinc,&svmin,&slipmin,keepsign,dump_slip);
}
