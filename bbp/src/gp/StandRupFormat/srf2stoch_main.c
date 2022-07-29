#include "include.h"
#include "structure.h"
#include "function.h"
#include "defs.h"
#include "getpar.h"

#define         RPERD           0.017453293
#define         DPERR           57.29577951

int main(int ac,char **av) {
	struct standrupformat srf;
	char infile[256], outfile[256];
	int inbin = 0;
	int i, ix, iy;
	int debug = 0;
	FILE* fpw;

	sprintf(infile,"stdin");
	sprintf(outfile,"stdout");

	setpar(ac,av);
	getpar("infile","s",infile);
	getpar("outfile","s",outfile);
	getpar("inbin","d",&inbin);
	endpar();

	read_srf(&srf,infile,inbin);

	if(strcmp(outfile,"stdout") == 0)
	   fpw = stdout;
	else
	   fpw = fopfile(outfile,"w");


	fprintf(fpw,"%d\n",srf.srf_prect.nseg);
	fflush(fpw);

	struct slipfile sfile;

	srf2stoch(ac, av, &srf, &sfile, debug);

	for(i=0;i<srf.srf_prect.nseg;i++) {
		fprintf(fpw,"%10.4f %10.4f %5d %5d %8.2f %8.2f\n",sfile.elon[i],sfile.elat[i],sfile.nx[i],sfile.ny[i],sfile.dx[i],sfile.dy[i]);
		fprintf(fpw,"%4.0f %4.0f %4.0f %8.2f %8.2f %8.2f\n",sfile.strike[i],sfile.dip[i],sfile.ravg[i],sfile.dtop[i],sfile.shypo[i],sfile.dhypo[i]);
		for(iy=0; iy<sfile.ny[i]; iy++) {
			for(ix=0; ix<sfile.nx[i]; ix++) {
				fprintf(fpw,"%13.5e",sfile.sp[i*NQ*NP + iy*sfile.nx[i] + ix]);
			}
			fprintf(fpw,"\n");
		}
		for(iy=0; iy<sfile.ny[i]; iy++) {
                        for(ix=0; ix<sfile.nx[i]; ix++) {
//                                fprintf(fpw,"%13.5e",sfile.tr[i*NQ*NP + iy*NP + ix]);
								fprintf(fpw,"%13.5e",sfile.tr[i*NQ*NP + iy*sfile.nx[i] + ix]);
                        }
						fprintf(fpw,"\n");
                }
		for(iy=0; iy<sfile.ny[i]; iy++) {
                        for(ix=0; ix<sfile.nx[i]; ix++) {
                                //fprintf(fpw,"%13.5e",sfile.ti[i*NQ*NP + iy*NP + ix]);
								fprintf(fpw,"%13.5e",sfile.ti[i*NQ*NP + iy*sfile.nx[i] + ix]);
                        }
						fprintf(fpw,"\n");
                }
	}
	fflush(fpw);
	fclose(fpw);
}
