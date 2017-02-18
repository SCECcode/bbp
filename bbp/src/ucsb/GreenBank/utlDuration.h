/*******************************
*program to find right duration*
*for greens functions          *
*author: Jan Schmede           *
*date: Oct. 1, 2007            *
*version: 1                    *
*******************************/

#include<math.h>
#include<string.h>
#include<stdio.h>
#include<stdlib.h>

typedef struct{
		int num; //number points
		double dt; //sampling 
		double *absval; //absolute values
		double max; //maximum
		double thresh; //threshold
		int window; //window length
		int minNum; //minimum number necessary to not have energy at end
		}GREEN;

void readGreen(GREEN **gr);
void getThreshold(GREEN **gr);
void getMinLength(GREEN **gr);
