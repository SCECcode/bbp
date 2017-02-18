/**************************
*routines for 3 comp seism*
*by Jan Schmedes          *
***************************/

#include"utlFault.h"

typedef struct{
		char name[32];
		int num;
		double lon,lat;
		double dt;
		double *ns;
		double *ew;
		double *ud;
		}SEISMOGRAM;

typedef struct{
		int numStat;
		SEISMOGRAM *seis;
		}ALLSEIS;

void readStatList(ALLSEIS *se);
void read3comps(ALLSEIS *se);
void read3compsBB(ALLSEIS *se);
void write3Comp(ALLSEIS *se);
