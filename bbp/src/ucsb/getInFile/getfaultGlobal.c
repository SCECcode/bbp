/******************************
*creates faultGlobal.in from
*a srf file
******************************/

#include<stdio.h>
#include<stdlib.h>
#include<string.h>

int main(void)
{
	FILE *fin;
	int i;
	FILE *fout;
	double lat,lon,depthCenter;
	double length,width;
	double strike,dip,rake;
	double posS,posD;
	double mag;
	double dx,dy;
	long rand=1234567890;
	double dt;
	double corner=0;
	char velName[132];
	char filename[132];
	char outname[132]="faultGlobal.in";
	char line[132];

	printf("name of velocity model?\n");
	scanf("%s",velName);
	printf("name of srf file (has to end with .srf)?\n");
	scanf("%s",filename);
	fin=fopen(filename,"r");
	fgets(line,132,fin);
	fgets(line,132,fin);
	fgets(line,132,fin);
	sscanf(line,"%lf %lf %*d %*d %lf %lf",&lon,&lat,&length,&width);
	fgets(line,132,fin);
	sscanf(line,"%lf %lf %lf %lf %lf",&strike,&dip,&depthCenter,&posS,&posD);

	fclose(fin);
	fout=fopen(outname,"w");
	fprintf(fin,"%lf %lf %lf\n",lon,lat,depthCenter);
	fprintf(fin,"%lf %lf\n",length,width);
	fprintf(fin,"%lf %lf -1\n",strike,dip);
	fprintf(fin,"%lf %lf\n",posS,posD);
	fprintf(fin,"-1\n");
	fprintf(fin,"-1 -1\n");
	fprintf(fin,"1234567890\n");
	fprintf(fin,"-1\n");
	fprintf(fin,"-1\n");
	fprintf(fin,"%s\n",velName);	
	for(i=0;i<strlen(filename)-4;i++)
		fprintf(fin,"%c",filename[i]);
	fprintf(fin,"\n");

	fclose(fout);

	return 0;
}
	
