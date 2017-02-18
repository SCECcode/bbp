/**********************************
*convert global input to PC input *
***********************************/

#include"utlFault.h"

int main(void)
{
	globalIn in;
	pcOut out;

	readIn(&in,"faultGlobal.in");
	convert2PC(&in,&out);
	writeKM(&out);

	return 0;
}

