/*******************************************
*convert three one component time histories*
*into one file with all three components   *
********************************************/

#include"utl3Comp.h"

int main(void)
{
	ALLSEIS se;

	readStatListxy(&se);
	read3comps(&se);
	write3Compxy(&se);
}
