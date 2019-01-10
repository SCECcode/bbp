#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define JMAX 40

void wvfast(int n, double* ar, double* ai, int ind);
double a0_func(double crit, double *x, int nt, double dt, double offset);
double rtbis(double (*func)(double, double *, int, double, double), double x1, double x2, double xacc, double *xx, int nt, double dt, double offset);

int seismic_intensity(int nt, double dt, double xr[], double& eji)
{
	double *fx1r, *fx2r;
	double *fx1i, *fx2i;
	double *xx;
	double y, y2, y3, y4, y6, y8, y10, y12, aa;
	double dur, dur2, eom, dom, freq, a0, amax, xtmp;
	double f1, f2, f3;
	int na, i, nf;
	int nom;

	eji = 0.0;
	nf = 2;
	do {
		nf *= 2;
	} while (nf < nt);

// allocate memory and initialize

	if ((fx1r = (double*)calloc(nf, sizeof(double))) == NULL) return(1);
	if ((fx2r = (double*)calloc(nf, sizeof(double))) == NULL) return(1);
	if ((fx1i = (double*)calloc(nf, sizeof(double))) == NULL) return(1);
	if ((fx2i = (double*)calloc(nf, sizeof(double))) == NULL) return(1);
	if ((xx   = (double*)calloc(nf, sizeof(double))) == NULL) return(1);

    for( int j = 0 ; j < nt ; j++)
	{
	    fx1r[j] = xr[j];
	    fx2r[j] = xr[j];
	}

// FOURIER TRANSFORM

	wvfast(nf, fx1r, fx1i, -1);
	wvfast(nf, fx2r, fx2i, -1);

	dur = nf * dt;
	if (dur == 0.0) return(1);

	for (i = 0; i < nf; i++) {
		fx1r[i] *= dt;
		fx2r[i] *= dt;
		fx1i[i] *= dt;
		fx2i[i] *= dt;
	}

	dom = 6.2831853070e0 / dur;
	nom = nf / 2 + 1;
	eom = dom * nom;

// BAND PASS filter

	for (i = 1; i <= nom; i++) {
		freq = i / dur;
		f1 = sqrt(1.0 / freq);
		y = freq / 10.0;
		y2 = y * y;
		y4 = y2 * y2;
		y6 = y2 * y4;
		y8 =y4 * y4;
		y10 = y4 * y6;
		y12 = y4 * y8;
		f2 = 1.0 / sqrt(1.0 + 0.694 * y2 
				+ 0.241 * y4 + 0.0557 * y6 
				+ 0.009664 * y8 + 0.00134 * y10 
				+ 0.000155 * y12);
		y = 2.0 * freq;
		y3 = y * y * y;
		f3 = sqrt(1.0 - exp(-y3));
		aa = f1 * f2 * f3;
		fx1r[i - 1] *= aa;
		fx1i[i - 1] *= aa;
		fx2r[i - 1] *= aa;
		fx2i[i - 1] *= aa;
	}

// INVERSE FOURIER TRANSFORM

	for (i = 2; i <= nf / 2; i++) {
		fx1r[nf + 2 - i - 1] = fx1r[i - 1];
		fx1i[nf + 2 - i - 1] = -fx1i[i - 1];
		fx2r[nf + 2 - i - 1] = fx2r[i - 1];
		fx2i[nf + 2 - i - 1] = -fx2i[i - 1];	
	}

	wvfast(nf, fx1r, fx1i, 1);
	wvfast(nf, fx2r, fx2i, 1);

// filtered WAVE DATA

	amax = 0.0;
	dur2 = dur * dur;
	for (i = 0; i < nt; i++) {
		xtmp = fx1r[i] * fx1r[i];
		xtmp += fx2r[i] * fx2r[i];
		xtmp /= dur2;
		xx[i] = sqrt(xtmp);
		if (xx[i] > amax) amax = xx[i];
	}

// Calculation of the A0-Tau relation

	na = 1000;
	a0 = rtbis(&a0_func, 0.0, amax, amax / na, xx, nt, dt, 0.299);	
	eji = 2.0 * log10(a0) + 0.943;

// Memory Free

	free(fx1r);
	free(fx2r);
	free(fx1i);
	free(fx2i);
    free(xx);

	return(0);
} 

void wvfast(int n, double* ar, double* ai,  int ind) 
{
	double tempr, tempi, thetai;
	static double pai = 3.1415926535897932384626434;
	int i, j, m, k;
	int kmax, istep;

//calucuration of fft

	j = 1;
	for (i = 1; i <= n; i++) {
		if (i < j) {
			tempr = ar[j - 1];
			tempi = ai[j - 1];
			ar[j - 1] = ar[i - 1];
			ai[j - 1] = ai[i - 1];
			ar[i - 1] = tempr;
			ai[i - 1] = tempi;
		} 
		m = n / 2;
		do {
			if (j <= m) break;
			j -= m;
			m /= 2;
		} while (m >= 2);
		j += m;
	}
	kmax = 1;
	while (kmax < n) {
		istep = kmax * 2;
		for (k = 1; k <= kmax; k++) {
			thetai = pai * (double)(ind * (k - 1)) 
					/ (double) (kmax);
			for (i = k; i <= n; i += istep) {
				j = i + kmax;
				tempr = ar[j - 1] * cos(thetai) 
					- ai[j - 1] * sin(thetai);
				tempi = ar[j - 1] * sin (thetai) 
					+ ai[j - 1] * cos(thetai);
				ar[j - 1] = ar[i - 1] - tempr;
				ai[j - 1] = ai[i - 1] - tempi;
				ar[i - 1] = ar[i - 1] + tempr;
				ai[i - 1] = ai[i - 1] + tempi;
			}
		}
		kmax = istep;
	}
}

double
rtbis(
	double (*func)(double, double *, int, double, double),
	double x1,
	double x2,
	double xacc,
	double *xx,
	int nt,
	double dt,
	double offset
)
{
	int j;
	double dx, f, fmid, xmid, rtb;

	f = (*func)(x1, xx, nt, dt, offset);
	fmid = (*func)(x2, xx, nt, dt, offset);

	if (f * fmid >= 0.0) return 0.0;
	rtb = f < 0.0 ? (dx = x2 - x1, x1) : (dx = x1 - x2, x2);
	for (j = 1; j <= JMAX; j++) {
		fmid = (*func)(xmid = rtb + (dx *= 0.5), xx, nt, dt, offset);
		if (fmid <= 0.0 ) rtb = xmid;
		if (fabs(dx) < xacc || fmid == 0.0) return rtb;
	}
	return 0.0;
}
	
double a0_func(double crit, double* x, int nt, double dt, double offset)
{
	int i;
	int cnt;
	double ans;

	for (i = 0, cnt = 0; i < nt; i++) {
		if (x[i] > crit) cnt++;
	}
	ans = cnt * dt - offset;
	return ans;
}
