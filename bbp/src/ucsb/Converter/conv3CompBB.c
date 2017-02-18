/*******************************************
*convert three one component time histories*
*into one file with all three components   *
********************************************/

#include"utl3Comp.h"

int main(void)
{
	ALLSEIS se;

	readStatList(&se);
	printf("stations read\n");
	read3compsBB(&se);
	printf("input read\n");
	write3Comp(&se);
}
