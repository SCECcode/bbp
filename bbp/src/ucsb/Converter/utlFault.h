//version 2, November 6, 2007
#include"include.h"

//use top center lon lat as center, that is as origin of local COSY

typedef struct{
		double topCenLon,topCenLat;
		double topCenX,topCenY,topCenZ;
		double ruplen,ddwidth;
		double strike,dip,rake;
		double hypoStrike,hypoDip;
		double M0;
		double Mw; //new v2
		int nx,ny;
		double dx,dy,dx2; //new v2
		int randomSeed;
		double dt,fc;
		char velmod[32];
		char name[132];
		}globalIn; //global input

typedef struct{
		double hypoX,hypoY,hypoZ;
		double ruplen,ddwidth;
		double strike,dip,rake;
		double hypoStrike,hypoDip;
		double M0;
		int nx,ny;
		int randomSeed;
		double dt,fc,dx;
		char velmod[32];
		char name[32];
		}pcOut; //output for output_LAH

void readIn(globalIn *in,char filename[]);
void convert2PC(globalIn *in,pcOut *out);
void writePC(pcOut *out);
