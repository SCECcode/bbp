/*******************************
*convert lon lat str to xy str *
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
	FILE *fout,*sourceout;
	char in[132],out[132],in2[132],line[132];
	float lon,lat,lonc,latc,x,y;
	double junk1,junk2,junk3,junk4,junk5,junk6;
	int nt,npoi,ncol,i,j,ndo;
	int j1,j2;
	int nFile,n;
	globalIn inp;
	float rot;
//new v2
//parameter for Rob's code
	double g0, b0;
	double amat[9], ainv[9];
	int xy2ll = 0;
	int ll2xy = 1;
	float erad = ERAD;
//end new v2

	readIn(&inp,"faultGlobal.in");
	ncol=6;
	
	sourceout=fopen("source_model.list","w");
	if(sourceout==NULL)
	{
		printf("can't open source_model.list\n");
		exit(0);
	}

	fprintf(sourceout,"1 0.0\n");

	sprintf(in,"%s.srf",inp.name);
//.if for internal format
	sprintf(out,"%s.001.srf",inp.name);
	fin=fopen(in,"r");
	fout=fopen(out,"w");
	fprintf(sourceout,"%s\n",out);
	if(fin==NULL)
	{
		sprintf(in,"%s",inp.name);
		fin=fopen(in,"r");
		if(fin==NULL)
		{
			printf("can't find %s.srf or %s\n",inp.name,in);
			exit(0);
		}
	}
	if(fout==NULL)
	{
		printf("can't open %s\n",out);
		exit(0);
	}

	fgets(line,132,fin);
	fputs(line,fout);
	fgets(line,132,fin);
	fputs(line,fout);
	fgets(line,132,fin);
	sscanf(line,"%f %f %d %d %lf %lf",&lonc,&latc,&j1,&j2,&junk1,&junk2);
	//as top center is origin of local cosy
	fprintf(fout,"0. 0.  %d %d %lf %lf\n",j1,j2,junk1,junk2);
	fgets(line,132,fin);
	fputs(line,fout);
	fgets(line,132,fin);
	fputs(line,fout);

//set origin
//exchanged v2
	rot=0-90; //no rotation of COSY
   	gen_matrices(amat,ainv,&rot,&lonc,&latc);
	g0=0.0;
	b0=0.0;
//end exchanged v2

	sscanf(line,"%*s %d",&npoi);

	fgets(line,132,fin);
	for(i=0;i<npoi;i++)
	{
		sscanf(line,"%f %f %lf %lf %lf %lf %lf %lf",&lon,&lat,&junk1,&junk2,&junk3,&junk4,&junk5,&junk6);
		//printf("%f %f %f %f\n",lon,lat,lonc,latc);
//exchanged v2
      		gcproj(&x,&y,&lon,&lat,&erad,&g0,&b0,amat,ainv,ll2xy);
//end exchanged v2
		fprintf(fout,"%lf %lf %lf %lf %lf %lf %lf %lf\n",x*1000.,y*1000.,junk1*1000.,junk2,junk3,junk4,junk5,junk6);
		fgets(line,132,fin);
		sscanf(line,"%*s %*s %d",&nt);
		fputs(line,fout);
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
	fclose(sourceout);

	return 0;
}
