/*******************************
*program to find right duration*
*for greens functions          *
*author: Jan Schmede           *
*date: Oct. 1, 2007            *
*version: 1                    *
*******************************/
#include"utlDuration.h"

void readGreen(GREEN **gr)
{
	GREEN *green;
	FILE *fin;
	char line[512];
	int i,j;
	green=(GREEN*)calloc(8,sizeof(GREEN));

	fin=fopen("test","r");
	if(fin==NULL)
	{
		printf("cannot open file test containing Greens functions!\n");
		exit(0);
	}
	//first 2 junk
	fgets(line,512,fin);	
	fgets(line,512,fin);	
	for(i=0;i<8;i++)
	{
		green[i].max=0.;
		fgets(line,512,fin);
		fgets(line,512,fin);
		sscanf(line,"%d %lf",&green[i].num,&green[i].dt);
		green[i].absval=(double*)calloc(green[i].num,sizeof(double));
		for(j=0;j<green[i].num;j++)
		{
			fscanf(fin,"%lf",&green[i].absval[j]);
			green[i].absval[j]=fabs(green[i].absval[j]);
			if(green[i].absval[j]>green[i].max)
			{
				green[i].max=green[i].absval[j];
			}
		}
		fgets(line,512,fin);
	}
	
	fclose(fin);
	*gr=green;
}	

void getThreshold(GREEN **gr)
{
	FILE *fin;
	double denom;
	double tlength;
	int i;

	fin=fopen("durationParam.in","r");
	if(fin==NULL)
	{
		printf("cannot find file durationParam.in, use 1/40 of maximum for threshold (denominator=40) and 2 sec for windowlength\n");
		denom=40.;
		tlength=2.;
	}
	else
	{
		fscanf(fin,"%lf %lf",&denom,&tlength);
		fclose(fin);
	}
	for(i=0;i<8;i++)
	{
		(*gr+i)->thresh=(*gr+i)->max/denom;
		(*gr+i)->window=floor(tlength/((*gr+i)->dt));
	}
}

void getMinLength(GREEN **gr)
{
	int i,pos;
	int j;	
	double mean;
	double factor;

	for(i=0;i<8;i++)
	{
		factor=1./((double)(*gr+i)->window);
		pos=(*gr+i)->num;
		mean=0.;
		for(j=0;j<(*gr+i)->window;j++)
			mean+=(*gr+i)->absval[pos-j]*factor;
		//printf("%lf %lf %d %d\n",mean,(*gr+i)->thresh,pos,i);
		while(mean<(*gr+i)->thresh)
		{
			mean+=-(*gr+i)->absval[pos]*factor;
			mean+=(*gr+i)->absval[pos-(*gr+i)->window];
			pos--;
			printf("%lf %lf %d %d\n",mean,(*gr+i)->thresh,pos,i);
		}		
		(*gr+i)->minNum=pos;
	}
}
