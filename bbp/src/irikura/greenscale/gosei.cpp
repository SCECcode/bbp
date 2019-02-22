#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <complex>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

using namespace std;

void GetGousei( std::vector<double> *sum, std::vector<double> *wv, int midshift , int vrshift)
{
	int i;
	for( i = 0 ; i < midshift ; i++ ){
	  (*sum)[i+vrshift] += (*wv)[i];
	}
	return;
}
