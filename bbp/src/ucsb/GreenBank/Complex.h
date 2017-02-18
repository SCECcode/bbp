/*********************************************************
*			Complex.h
*  header file for complex data type, source codes are in
*  Complex.c and fft.c
*********************************************************/

#ifndef __MY_COMPLEX__
  #define __MY_COMPLEX__


/* data type */
typedef struct { float x; float y;} complex;

/* constants */
#define PI	3.1415926
#define Imag cmplx(0., 1.)
#define One  cmplx(1., 0.)
#define Zero cmplx(0., 0.)

/* basic operations */
complex cplus(complex a, complex b);
complex cmltp(complex a, complex b);
complex cngtv(complex a);
complex cinvs(complex a);
complex conjg(complex a);
complex dmltp(float a, complex b);
complex csqrt(complex a);
complex cmplx(float x, float y);
complex cphase(complex w);
double  ccabs(complex a);

/* fft */
void    fft(complex *a, int n, float dt);	/* dt>0: forw.; dt<0: inv */
void    fftr(complex *x, int n, float dt);

/* convolution and correlation */
void	cor(complex *a, complex *b, float dt, int nft);
void	conv(float *, int, float *, int);
float	*crscrl(int,float *,float *,int);
float	maxCor(float *, float *, int, int *, float *);

/* integration */
float amp(float t1, float t2, float *data, int n);
void cumsum(float *a, int n, float dt);

/* windowing */
float *coswndw(int, float);

/* high-pass filtering */
void	filter(complex *, int, float, float, float, int);

/* find max. absolute values in an array */
int findMaxAbs(float *cor, int n, float *amp);

/* remove rerend a+b*x */
void rtrend(float *, int);

/* some operation on spectrum */
void fltGauss(complex *u, int n, float gauss);
void shiftSpec(complex *u, int n, float shift);
void specAdd(complex *a, complex *b, int n);
void specMul(complex *a, complex *b, int n);
void specScale(complex *a, float c, int n);
float specPwr(complex *u, int n);

#endif
