/**************************
*routines for 3 comp seism*
*by Jan Schmedes          *
***************************/
#include"utl3Comp.h"

void readStatList(ALLSEIS *se)
{
	FILE *fin;
	char line[132];
	int i;

	fin=fopen("stations.ll","r");
	if(fin==NULL)
	{
		printf("cannot find stations.ll\n");
		exit(0);
	}
	fgets(line,132,fin);
	sscanf(line,"%d",&se->numStat);
	(*se).seis=(SEISMOGRAM *)calloc(se->numStat,sizeof(SEISMOGRAM));
	for(i=0;i<se->numStat;i++)
	{
		fgets(line,132,fin);
		sscanf(line,"%lf %lf %s",&(*se).seis[i].lon,&(*se).seis[i].lat,(*se).seis[i].name);
	}
	fclose(fin);
}
		
void readStatListxy(ALLSEIS *se)
{
	FILE *fin;
	char line[132],name[132];
	int i;
	
	printf("name of stations file?\n");
	scanf("%s",name);
	fin=fopen(name,"r");
	if(fin==NULL)
	{
		printf("cannot find %s\n",name);
		exit(0);
	}
	fgets(line,132,fin);
	sscanf(line,"%d",&se->numStat);
	(*se).seis=(SEISMOGRAM *)calloc(se->numStat,sizeof(SEISMOGRAM));
	for(i=0;i<se->numStat;i++)
	{
		fgets(line,132,fin);
		sscanf(line," %s",(*se).seis[i].name);
//here lon and lat refer to x and y in cartesian
		fgets(line,132,fin);
		sscanf(line,"%lf %lf",&(*se).seis[i].lon,&(*se).seis[i].lat);
	}
	fclose(fin);
}
void read3comps(ALLSEIS *se)
{
	FILE *ns,*ew,*ud;
	int i,j;
	char nsN[32],ewN[32],udN[32];
	char line[132];
	
	for(i=0;i<se->numStat;i++)
	{
		sprintf(nsN,"%s.000.gm1D.001",(*se).seis[i].name);
		sprintf(ewN,"%s.090.gm1D.001",(*se).seis[i].name);
		sprintf(udN,"%s.ver.gm1D.001",(*se).seis[i].name);
		ns=fopen(nsN,"r");
		ew=fopen(ewN,"r");
		ud=fopen(udN,"r");
		if(ns==NULL || ew==NULL || ud==NULL)
		{
			printf("component(s) of %s missing\n",(*se).seis[i].name);
			printf("%s %s %s\n", nsN, ewN, udN);
			exit(0);
		}

		fgets(line,132,ns);
		fgets(line,132,ew);
		fgets(line,132,ud);
		sscanf(line,"%d %lf",&(*se).seis[i].num,&(*se).seis[i].dt);
		(*se).seis[i].ns=(double*)calloc((*se).seis[i].num,sizeof(double));
		(*se).seis[i].ew=(double*)calloc((*se).seis[i].num,sizeof(double));
		(*se).seis[i].ud=(double*)calloc((*se).seis[i].num,sizeof(double));
		for(j=0;j<(*se).seis[i].num;j++)
			fscanf(ns,"%lf",&(*se).seis[i].ns[j]);
		fclose(ns);
		for(j=0;j<(*se).seis[i].num;j++)
			fscanf(ew,"%lf",&(*se).seis[i].ew[j]);
		fclose(ew);
		for(j=0;j<(*se).seis[i].num;j++)
			fscanf(ud,"%lf",&(*se).seis[i].ud[j]);
		fclose(ud);
	}
}

void read3compsBB(ALLSEIS *se)
{
	FILE *ns,*ew,*ud;
	int i,j;
	char nsN[32],ewN[32],udN[32];
	char line[132];
	
	for(i=0;i<se->numStat;i++)
	{
		sprintf(nsN,"%s.000.gmBB.001",(*se).seis[i].name);
		sprintf(ewN,"%s.090.gmBB.001",(*se).seis[i].name);
		sprintf(udN,"%s.ver.gmBB.001",(*se).seis[i].name);
		ns=fopen(nsN,"r");
		ew=fopen(ewN,"r");
		ud=fopen(udN,"r");
		if(ns==NULL || ew==NULL || ud==NULL)
		{
			printf("component(s) of %s missing\n",(*se).seis[i].name);
			printf("%s %s %s\n", nsN, ewN, udN);
			exit(0);
		}

		fgets(line,132,ns);
		fgets(line,132,ew);
		fgets(line,132,ud);
		sscanf(line,"%d %lf",&(*se).seis[i].num,&(*se).seis[i].dt);
		(*se).seis[i].ns=(double*)calloc((*se).seis[i].num,sizeof(double));
		(*se).seis[i].ew=(double*)calloc((*se).seis[i].num,sizeof(double));
		(*se).seis[i].ud=(double*)calloc((*se).seis[i].num,sizeof(double));
		for(j=0;j<(*se).seis[i].num;j++)
			fscanf(ns,"%lf",&(*se).seis[i].ns[j]);
		fclose(ns);
		for(j=0;j<(*se).seis[i].num;j++)
			fscanf(ew,"%lf",&(*se).seis[i].ew[j]);
		fclose(ew);
		for(j=0;j<(*se).seis[i].num;j++)
			fscanf(ud,"%lf",&(*se).seis[i].ud[j]);
		fclose(ud);
	}
}
void write3Comp(ALLSEIS *se)
{
	char name[32];
	FILE *fout;
	int i,j;
	globalIn in;

	readIn(&in,"faultGlobal.in");

	for(i=0;i<se->numStat;i++)
	{
		sprintf(name,"%s.3comp",(*se).seis[i].name);
		fout=fopen(name,"w");
		if(fout==NULL)
		{
			printf("cannot open file for writing!\n");
			exit(0);
		}	
		fprintf(fout,"%% lon and lat: %lf %lf\n",(*se).seis[i].lon,(*se).seis[i].lat);
		fprintf(fout,"%% number of points: %d\n",(*se).seis[i].num);
		fprintf(fout,"%% name of source model used for computation: %s\n",in.name);
		fprintf(fout,"%% name of velocity model used: %s\n",in.velmod);
		fprintf(fout,"%% seismic moment of event: %g Nm\n",in.M0);
		fprintf(fout,"%% dt(s)   n-s(cm/s)   e-w(cm/s)  u-d(cm/s)\n");
		for(j=0;j<(*se).seis[i].num;j++)
		        /* FS/JC: Feb 2015: Multiplying by 100 to convert from m/s to cm/s */
			fprintf(fout,"%lf %.9lf %.9lf %.9lf\n",(double)j*(*se).seis[i].dt,(*se).seis[i].ns[j]*100,(*se).seis[i].ew[j]*100,(*se).seis[i].ud[j]*100);
		fclose(fout);
	}
}
void write3Compxy(ALLSEIS *se)
{
	char name[32];
	FILE *fout;
	int i,j;


	for(i=0;i<se->numStat;i++)
	{
		sprintf(name,"%s.3comp",(*se).seis[i].name);
		fout=fopen(name,"w");
		if(fout==NULL)
		{
			printf("cannot open file for writing!\n");
			exit(0);
		}	
		fprintf(fout,"%% x and y: %lf %lf\n",(*se).seis[i].lon,(*se).seis[i].lat);
		fprintf(fout,"%% number of points: %d\n",(*se).seis[i].num);
		fprintf(fout,"%% dummy line\n");
		fprintf(fout,"%% dummy line\n");
		fprintf(fout,"%% dummy line\n");
		fprintf(fout,"%% dt(s)   n-s(cm/s)   e-w(cm/s)  u-d(cm/s)\n");
		for(j=0;j<(*se).seis[i].num;j++)
			fprintf(fout,"%lf %lf %lf %lf\n",(double)j*(*se).seis[i].dt,(*se).seis[i].ns[j],(*se).seis[i].ew[j],(*se).seis[i].ud[j]);
		fclose(fout);
	}
}
