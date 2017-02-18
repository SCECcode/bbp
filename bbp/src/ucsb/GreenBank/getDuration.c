/*******************************
*program to find right duration*
*for greens functions          *
*author: Jan Schmede           *
*date: Oct. 1, 2007            *
*version: 1                    *
*******************************/
#include"utlDuration.h"

int main(void)
{
	GREEN *g;
	int i;
	int longer=0;

	readGreen(&g);
	getThreshold(&g);
	getMinLength(&g);
	for(i=0;i<8;i++)
	{
		printf("%d\n",g[i].minNum);
		if(g[i].minNum==g[i].num)
			longer=1;
	}
	if(longer==1)
	{
		printf("Need longer Greensfunctions!!!\n");
	}
	
	return 0;
}
