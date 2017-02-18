/*******************************
*convert lon lat to xy station *
*format used by PC code        *
*written by J. Schmedes        *
*using subroutine provided by  *
*R. Graves                     *
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
	float lon,lat,lonc,latc,x,y;
	FILE *fin;
	FILE *fout;
	char in[132],out[132],line[132],name[32];
	int i,nstat;
	double hx,hy,strike,comp1,comp2;
//parameter for Rob's code
	double g0, b0;
	double amat[9], ainv[9];
	int xy2ll = 0;
	int ll2xy = 1;
	float erad = ERAD;
	float rot;


	printf("lon and lat of center\n");
	scanf("%f %f",&lonc,&latc);
	hx=0.;
	hy=0.;
	strike=0; //no rotation of COSY
	comp1=0.; //N-S
	comp2=90.; //E-W

	sprintf(in,"stations.ll");
	sprintf(out,"stations.xy");	
	fin=fopen(in,"r");
        fout=fopen(out,"w");
        if(fin==NULL)
        {
                printf("can't find %s\n",in);
                exit(0);
        }
        if(fin==NULL)
        {
                printf("can't open %s\n",out);
                exit(0);
        }

//set origin
	rot=0-90; //no rotation of COSY
   	gen_matrices(amat,ainv,&rot,&lonc,&latc);
	g0=0.0;
	b0=0.0;

	fgets(line,132,fin);
	sscanf(line,"%d",&nstat);
	fprintf(fout,"%d %lf %lf %lf\n",nstat,hx,hy,strike);
	for(i=0;i<nstat;i++)	
	{
		fgets(line,132,fin);
		sscanf(line,"%f %f %s",&lon,&lat,name);
      		gcproj(&x,&y,&lon,&lat,&erad,&g0,&b0,amat,ainv,ll2xy);
		printf("%f %f\n",x,y);
		fprintf(fout,"%s\n",name);
		fprintf(fout,"%lf %lf 0.000 %lf %lf\n",x*1000,y*1000,comp1,comp2);
	}
	fclose(fout);
	fclose(fin);
	
	return 0;
}
