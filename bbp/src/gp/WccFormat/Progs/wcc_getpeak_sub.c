#include "include.h"
#include "structure.h"
#include "function.h"
#include "getpar.h"

float wcc_getpeak(int param_string_len, char** param_string, float* s1, struct statdata* head1) {
	float amax;
	int i;

	float max = -1.0e+20;
	float min = 1.0e+20;
	int keepsign = 0;
	float scale = 1.0;

	setpar(param_string_len, param_string);
	getpar("keepsign","d",&keepsign);
	getpar("scale","f",&scale);
	endpar();
 
	for(i=0;i<head1->nt;i++)
	   {
	   if(s1[i] > max)
	      max = s1[i];
	   if(s1[i] < min)
	      min = s1[i];
	   }
	 
	if(max >= 0.0 && min < 0.0)
	   {
	   if(max > -min)
	      amax = max;
	   else
	      {
	      amax = -min;
	      if(keepsign)
	         amax = -amax;
	      }
	   }
	
	else if(max >= 0.0 && min >= 0.0)
	   amax = max;
	
	else if(max < 0.0 && min < 0.0)
	   {
	   amax = -min;
	   if(keepsign)
	      amax = -amax;
	   }
	return scale*amax;
}
