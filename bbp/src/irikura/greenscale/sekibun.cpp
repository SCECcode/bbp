#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <complex>
#include <algorithm>
#include <string>

void Ssekibun(std::vector<double> *indata, std::vector<double> *outdata,
              int n, double s_time, double term, double gain)
{
	int		i;
	double	PI, OMEGA, OMSET;
	double	tmp1, tmp2, tmp3, tmp4;

	PI = 3.14159265358979;
	OMEGA = 2.0 * PI / term;
	OMSET = OMEGA * s_time;

	tmp4 = 12.0 * (1.0 + (gain * OMSET)) + pow(OMSET,2.0);

	tmp1 = 10.0 * pow(OMSET, 2.0) - 24.0;
	tmp1 = tmp1 / tmp4;
	tmp2 = 12.0 * (1.0 - (gain * OMSET)) + pow(OMSET,2.0);
	tmp2 = tmp2 / tmp4;

	tmp3 = 6.0 * s_time / tmp4;

	(*outdata)[0] = tmp3 * (*indata)[0];
    (*outdata)[1] = -tmp1 * (*outdata)[0] + tmp3 * (*indata)[1];
	
	for (i = 2; i < n; i++) {
		(*outdata)[i] = - tmp1 * (*outdata)[i-1]
						- tmp2 * (*outdata)[i-2]
						+ tmp3 * ( (*indata)[i] - (*indata)[i-2] );
	}
}
