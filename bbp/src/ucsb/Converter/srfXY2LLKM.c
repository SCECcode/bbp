/*******************************
*convert x y str to lonlat str *
*written by J. Schmedes        *
********************************/
//v2, November 6, 2007
//use Rob's subroutines below now
//exchanged v2
#include"utlFault.h"

#define         ERAD            6378.139
extern void gcproj(float *xf,float *yf,float *rlon,float *rlat,float *ref_rad,double *g0,double *b0,double *amat,double *ainv,int gflag);
extern void gen_matrices(double *amat,double *ainv,float *alpha,float *ref_lon,float *ref_lat);
//end exchanged v2


int main(void)
{
	FILE *fin;
	FILE *fout;
	char in[132],out[132],in2[132],line[132];
	float lon,lat,lonc,latc,x,y,xc,yc;
	double junk1,junk2,junk3,junk4,junk5,junk6;
	int nt,npoi,ncol,i,j,ndo;
	int j1,j2;
	int nFile,n;
	globalIn inp;
//new v2
//parameter for Rob's code
	double g0, b0;
	double amat[9], ainv[9];
	int xy2ll = 0;
	int ll2xy = 1;
	float erad = ERAD;
	float rot;
//end new v2

	//read in global parameters to get center
	readIn(&inp,"faultGlobal.in");
	
	//set origin
	lonc=inp.topCenLon;
	latc=inp.topCenLat;
//set origin
//exchanged v2
	rot=0-90; //no rotation of COSY
   	gen_matrices(amat,ainv,&rot,&lonc,&latc);
	g0=0.0;
	b0=0.0;
//end exchanged v2

	ncol=6; //should have 6 cols for source time function

	sprintf(in2,"%s",inp.name);
	sprintf(out,"%s.srf",inp.name);
	fout=fopen(out,"w");
	fin=fopen(in2,"r");
	if(fin==NULL)
	{
		printf("can't find %s\n",in2);
		exit(0);
	}
	if(fin==NULL)
	{
		printf("can't open %s\n",out);
		exit(0);
	}

	fgets(line,132,fin);
	fputs(line,fout);
	fgets(line,132,fin);
	fputs(line,fout);
	fgets(line,132,fin);
	sscanf(line,"%f %f %d %d %lf %lf",&xc,&yc,&j1,&j2,&junk1,&junk2);
	xc*=0.001;
	yc*=0.001;
//exchanged v2
      		gcproj(&xc,&yc,&lon,&lat,&erad,&g0,&b0,amat,ainv,xy2ll);
//end exchanged v2
	fprintf(fout,"%f %f %d %d %lf %lf\n",lon,lat,j1,j2,junk1,junk2);
	fgets(line,132,fin);
	sscanf(line,"%lf %lf %lf %lf %lf",&junk1,&junk2,&junk3,&junk4,&junk5);
	j1=(int)junk1;	
	j2=(int)junk2;	
	//fputs(line,fout);
	fprintf(fout,"%d %d %lf %lf %lf\n",j1,j2,junk3,junk4,junk5);
	fgets(line,132,fin);
	fputs(line,fout);


	sscanf(line,"%*s %d",&npoi);
	
	fgets(line,132,fin);
	for(i=0;i<npoi;i++)
	{
		sscanf(line,"%f %f %lf %lf %lf %lf %lf %lf",&x,&y,&junk1,&junk2,&junk3,&junk4,&junk5,&junk6);
		j1=(int)junk2;	
		j2=(int)junk3;	
//exchanged v2
		x*=0.001;
		y*=0.001;
      		gcproj(&x,&y,&lon,&lat,&erad,&g0,&b0,amat,ainv,xy2ll);
//end exchanged v2
		fprintf(fout,"%7.4f %7.4f %7.4lf %d %d %6.5e %7.4lf %6.5e\n",lon,lat,junk1/1000.,j1,j2,junk4,junk5,junk6);
		fgets(line,132,fin);
		sscanf(line,"%lf %lf %d",&junk2,&junk1,&nt);
		//fputs(line,fout);
		j1=(int)junk2;
		fprintf(fout,"%d %lf %d 0.0 0 0.0 0\n",j1,junk1,nt); 
		//printf("%d\n",nt);
		ndo=floor((double)nt/((double)ncol))+1;
		if((ndo-1)*ncol==nt) ndo--;
		for(j=0;j<ndo;j++)
		{
			fgets(line,132,fin);
			fputs(line,fout);
		}
		fgets(line,132,fin);
	}

	fclose(fin);
	fclose(fout);
	
	return 0;
}
