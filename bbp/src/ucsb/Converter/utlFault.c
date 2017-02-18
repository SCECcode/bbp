#include"utlFault.h"
void readIn(globalIn *in,char filename[])
{
	FILE *fin;
	int i;
	char line[132];
	int nx,ny;

	fin=fopen(filename,"r");
	if(fin==NULL)
	{
		printf("file %s does not exist!\n", filename);
		exit(0);
	}

	fgets(line,132,fin);
	sscanf(line,"%lf %lf %lf",&in->topCenLon,&in->topCenLat,&in->topCenZ);
	//top center as cosy origin
	in->topCenX=0.0;
	in->topCenY=0.0;
	fgets(line,132,fin);
	sscanf(line,"%lf %lf",&in->ruplen,&in->ddwidth);
	fgets(line,132,fin);
	sscanf(line,"%lf %lf %lf",&in->strike,&in->dip,&in->rake);
	fgets(line,132,fin);
	sscanf(line,"%lf %lf",&in->hypoStrike,&in->hypoDip);
	fgets(line,132,fin);
	sscanf(line,"%lf",&in->Mw); //changed v2
	fgets(line,132,fin);
	sscanf(line,"%lf %lf",&in->dx,&in->dy); //changed v2
	fgets(line,132,fin);
	sscanf(line,"%d",&in->randomSeed);
	fgets(line,132,fin);
	sscanf(line,"%lf",&in->dt);
	fgets(line,132,fin);
	sscanf(line,"%lf",&in->fc);
	fgets(line,132,fin);
	sscanf(line,"%s",in->velmod);
	fgets(line,132,fin);
	sscanf(line,"%s",in->name);
	in->dx2=in->dx*1000.;

//new v2
	if(in->randomSeed>0)
		in->randomSeed*=-1;
	in->M0=pow(10,1.5*in->Mw+9.05);
	nx=floor(in->ruplen/in->dx);
	ny=floor(in->ddwidth/in->dy);
	printf("%d %d\n",nx,ny);
	in->nx=2;
	while(in->nx<nx)
		in->nx*=2.;
	in->ny=2;
	while(in->ny<ny)
		in->ny*=2.;
	printf("%d %d\n",in->nx,in->ny);
	in->dx=in->ruplen/((double)in->nx);
	in->dy=in->ddwidth/((double)in->ny);
//end new v2
	
	fclose(fin);
}

void convert2PC(globalIn *in,pcOut *out)
{
	double rad=M_PI/180.;
	out->hypoStrike=(in->hypoStrike+in->ruplen/2.);
	out->hypoDip=in->hypoDip;
	out->hypoX=(in->topCenX+cos(in->strike*rad)*in->hypoStrike);	
	out->hypoX-=cos(in->dip*rad)*sin(in->strike*rad)*in->hypoDip;
	out->hypoY=in->topCenY+sin(in->strike*rad)*in->hypoStrike;	
	out->hypoY+=cos(in->dip*rad)*cos(in->strike*rad)*in->hypoDip;
	out->hypoZ=(in->topCenZ+sin(in->dip*rad)*in->hypoDip);
	out->ruplen=in->ruplen;
	out->ddwidth=in->ddwidth;
	out->strike=in->strike;
	out->dip=in->dip;
	out->rake=in->rake;
	out->randomSeed=in->randomSeed;
	out->M0=in->M0;
	out->fc=in->fc;
	out->nx=in->nx;
	out->ny=in->ny;
	out->dt=in->dt;
	out->dx=in->dx2;
	strcpy(out->velmod,in->velmod);
	strcpy(out->name,in->name);
}

void writePC(pcOut *out)
{
	FILE *fout;

	fout=fopen("output_LAH.inp","w");

	if(fout==NULL)
	{
		printf("cannot open file for writing!\n");
		exit(0);
	}
	
	fprintf(fout,"4\n"); //source time function #4
	fprintf(fout,"%lf %lf\n",out->ruplen,out->ddwidth);
	fprintf(fout,"%lf %lf %lf\n",out->hypoStrike,out->hypoDip,out->hypoZ);
	fprintf(fout,"%lf %lf\n",out->hypoX,out->hypoY);
	fprintf(fout,"%g %lf\n",out->M0,out->fc);
	fprintf(fout,"%lf %lf %lf\n",out->strike,out->dip,out->rake);
	fprintf(fout,"%d %d\n",out->nx,out->ny);
	fprintf(fout,"%d %d %d\n",out->randomSeed,out->randomSeed-6,out->randomSeed-4);
	fprintf(fout,"1,1\n");
	fprintf(fout,"%s\n",out->velmod);
	fprintf(fout,"0.0\n");
	fprintf(fout,"1\n");
	fprintf(fout,"%s\n",out->name);
	fprintf(fout,"%lf\n",out->dt);
}
	
void writeKM(pcOut *out)
{
	FILE *fout;

	fout=fopen("KinModel.inp","w");

	if(fout==NULL)
	{
		printf("cannot open file for writing!\n");
		exit(0);
	}
	
	if(out->randomSeed<0)
		out->randomSeed*=-1;
	fprintf(fout,"%lf %lf\n",out->ruplen*1000.,out->ddwidth*1000.);
	fprintf(fout,"%lf %lf\n",out->hypoStrike*1000.,out->hypoDip*1000.);
	fprintf(fout,"%lf %lf %lf\n",out->hypoX*1000.,out->hypoY*1000.,out->hypoZ*1000.);
	fprintf(fout,"%g %lf\n",out->M0,out->fc);						//JORGE Here I make the program write the corner frequency into KinModel.inp
	fprintf(fout,"%lf %lf %lf\n",out->strike,out->dip,out->rake);
	fprintf(fout,"%lf %lf\n",out->dx,out->dt);
	fprintf(fout,"1\n");
	fprintf(fout,"%d %d %d\n",out->randomSeed+23123,out->randomSeed+123213,out->randomSeed+454532);		//JORGE I modified this such that there are no longer negative random seeds
	fprintf(fout,"%s\n",out->velmod);
}
	
	
