#include "include.h"
#include "structure.h"
#include "function.h"
#include "getpar.h"

void copy_vpvs(struct velmodel *vsm,struct velmodel *vpm);

int main(int ac,char **av)
{
struct velmodel vpmod, vsmod;
FILE *fpr, *fopfile();
int nr, nd, i, j, ip, ishall;
float *rng, *dep, *ptime, *stime;
double rayp, rad;
char locfile[2048], velfile[2048];

double eps = 1.0e-12;
float recd = 0.0;
float shallowest_depth = 0.0;
float htol = 0.1;
float rloc = 0.0;

sprintf(locfile,"locations.dat");

setpar(ac,av);
mstpar("velfile","s",velfile);
getpar("locfile","s",locfile);
getpar("shallowest_depth","f",&shallowest_depth);
endpar();

fpr = fopfile(locfile,"r");

fscanf(fpr,"%d",&nr);
rng = (float *) check_malloc (nr*sizeof(float));
for(i=0;i<nr;i++)
   fscanf(fpr,"%f",&rng[i]);

fscanf(fpr,"%d",&nd);
dep = (float *) check_malloc (nd*sizeof(float));
for(i=0;i<nd;i++)
   fscanf(fpr,"%f",&dep[i]);

fclose(fpr);

ptime = (float *) check_malloc (nr*nd*sizeof(float));
stime = (float *) check_malloc (nr*nd*sizeof(float));

read_velmodel(velfile,&vsmod);
copy_vpvs(&vsmod,&vpmod);

for(i=0;i<nd;i++)
   {
   for(j=0;j<nr;j++)
      {
      get_rupt(&vpmod,&htol,&dep[i],&recd,&rloc,&rng[j],&rayp,&rad,&ptime[j+i*nr]);
      get_rupt(&vsmod,&htol,&dep[i],&recd,&rloc,&rng[j],&rayp,&rad,&stime[j+i*nr]);
      }
   }

for(i=0;i<nd;i++)
   {
   ishall = i;
   while(dep[ishall] < shallowest_depth)
      ishall++;

   for(j=0;j<nr;j++)
      {
      ip = j + ishall*nr;

      printf("%13.7f %13.7f %13.7f\n",dep[i],rng[j],ptime[ip]);
      printf("%13.7f %13.7f %13.7f\n",dep[i],rng[j],stime[ip]);
      }
   }
}

void copy_vpvs(struct velmodel *vsm,struct velmodel *vpm)
{
int i;

vpm->nlay = vsm->nlay;

vpm->vp = (float *)check_malloc(vpm->nlay*sizeof(float));
vpm->vs = (double *)check_malloc(vpm->nlay*sizeof(double));
vpm->den = (float *)check_malloc(vpm->nlay*sizeof(float));
vpm->th = (float *)check_malloc(vpm->nlay*sizeof(float));
vpm->dep = (float *)check_malloc(vpm->nlay*sizeof(float));
vpm->mu = (float *)check_malloc(vpm->nlay*sizeof(float));
vpm->invb2 = (double *)check_malloc(vpm->nlay*sizeof(double));

for(i=0;i<vpm->nlay;i++)
   {
   vpm->th[i] = vsm->th[i];
   vpm->vp[i] = vsm->vp[i];
   vpm->den[i] = vsm->den[i];
   vpm->dep[i] = vsm->dep[i];
   vpm->mu[i] = vsm->mu[i];

   vpm->vs[i] = vsm->vp[i];
   vpm->invb2[i] = 1.0/(vsm->vp[i]*vsm->vp[i]);
   }
}
