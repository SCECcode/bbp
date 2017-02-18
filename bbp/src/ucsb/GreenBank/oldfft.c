/*****************************************************************
-			fft.c
-	subroutines doing fast fourier transform, correlation, etc
-
- Includes:
-	fft()	- for complex sequence
-	fftr()	- for real sequence
-	fftr_() - wrap of fftr for fortran calling
-	cor()	- cross-correlation of two seqs. using their spectrums
-	crscrl()- cross-correlation of two seq., returns a portion
-	maxCor()- max cross-correlation
-       conv()  - convolve two time seq. in the time domain
-	amp()	- integrate a time seq. between two points
-	cumsum()- cummulative sume of a time seq.
-	coswndw() - generate a cosine taper window
-	filter()	- high-pass filtering data in frequency domain
-	findMax()	- find maxium in an array
-	rtrend()	- remove trend (a+b*t)
-	fltGauss()	- apply the Gaussian filter
-	shiftSpec()	- do time shift in freq. domain
-	specAdd()	- spectrum addition
-	specMul()	- spectrum multiplication
-	specScale()	- multiply a constant to spectrum
-	specPwr()	- compute auto-corr. using spectrum
-
- Revision History
-	Lupei Zhu	06/20/94	Initial revision
-	Lupei Zhu	12/02/99	conv() now can handle ns>n
-	Lupei Zhu	01/01/00	add more subroutines
******************************************************************/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include "Complex.h"

/*---------------------------------------------------------------
*   fft()
*   discrete fourier transform of complex sequence x[i], i=0,1,...,n-1.
*
*	fft{x}[i] = dt*SUM x[k]*exp(-j*i*k*pi/n) over k=0 to n-1
*   or
*       inv_fft{x}[i] = (1/(n*dt)) SUM x[k]*exp(j*i*k*pi/n) over k=0 to n-1
*   This should agree with analog Fourier transform in amplitude
*   Input arguments:
*   x (complex *)	- array for FFT (IN/OUT)
*   n (int)	- dimension of x[], n=2^N, N>0.
*   dt (float)	- time sampling interval, forward (>0) or inverse (<0) FFT;
*--------------------------------------------------------------*/
void	fft(complex *a, int n, float dt) {
  int i, j, k, step, m;
  complex u, w, t;
  double pi;
  pi = -PI;
  if (dt<0.) pi = PI;
  for (m=n/2,j=0,i=1; i<n-1; i++) {
    for (k=m; k<=j; k/=2) j -= k;
    j += k;
    if(i<j) {
      t = a[i];
      a[i] = a[j];
      a[j] = t;
    }
  }
  for (m=1,step=2; m<n; m=step, step*=2) {
    for (u=One,w=cmplx(cos(pi/m),sin(pi/m)),j=0; j<m; j++) {
      for (i=j; i<n; i+=step) {
	k = i+m;
	t = cmltp(a[k], u);
	a[k] = cplus(a[i], cngtv(t));
	a[i] = cplus(a[i], t);
      }
      u = cmltp(u, w);
    }
  }
  if (dt<0.) dt=-1./(n*dt);
  for (i=0; i<n; i++) a[i] = dmltp(dt, a[i]);
}


/*---------------------------------------------------------------
*   fftr()
*   fast fourier transform of real sequence
*   Input arguments:
*   x (complex *)	- array for FFT (IN/OUT). For the forward
*	transf., x is the real sequence stored as complex array;
*	for the inverse transf., x is the half of spectrum f and
*	f(n) is in Im(x[0]).
*   n (int)		- dimension of x[], n=2^N, N>0.
*   dt (float)		- forward (>0) or inverse (<0) FFT;
*--------------------------------------------------------------*/
void	fftr(complex *x, int n, float dt) {
  int	i, j, n2;
  float	delw, w;
  complex t, g, h, isg;
  n2 = n/2;
  delw = PI/n;
  isg = Imag;
  if (dt>0.) {
     delw = -delw;
     isg = cngtv(isg);
     fft(x, n, dt);
  }
  x[0] = cmplx(x[0].x+x[0].y, x[0].x-x[0].y);
  for (i=1, w=delw; i<n2; i++, w+=delw) {
    j = n-i;
    t = conjg(x[j]);
    g = cplus(x[i], t);
    h = cplus(x[i], cngtv(t));
    h = cmltp(cmplx(cos(w), sin(w)), h);
    x[i] = dmltp(0.5, cplus(g, cmltp(isg,h)));
    x[j] = dmltp(0.5, cplus(conjg(g), cmltp(isg, conjg(h))));
  }
  x[n2] = conjg(x[n2]);
  if (dt<0.) {
    x[0] = dmltp(0.5, x[0]);
    fft(x, n, dt);
  }
}


/*---------------------------------------------------------------
*   fftr_()
*   a wrap of fftr() for fortran calling
*   Input arguments:
*   x (complex *)	- array for FFT (IN/OUT)
*   n (int)	- dimension of the x[], n=2^N, N>0.
*   dt (float)	- forward (>0.) or inverse (<0.) FFT;
*--------------------------------------------------------------*/
void	fftr_(float *x, int n, float dt) {
  fftr((complex *) x, n/2, dt);
}


/*
correlation, IFFT{data[w]*conjugate(src[w])} = int(data(tau)*src(t-tau),tau)
the zero-lag is at data[nft/2].
*/
void	cor(
	    complex	*src,		/* In: source function */
	    complex     *data,		/* In: data */
	                                /* Out: cross-correlation */
	    float	dt,		/* In: dt */
	    int 	nft		/* In: number of pts */
	    )
{
  int	j;
  float	aa;

  aa = -1.;
  data[0]=cmplx(data[0].x*src[0].x, aa*data[0].y*src[0].y);
  for (j=1; j<nft; j++) {
    data[j]=cmltp(data[j], conjg(src[j]));
    data[j]=dmltp(aa, data[j]);
    aa = -aa;
  }
  fftr(data, nft, -dt);
}


/*
   Convolving s[] with f[] in time domain, the result is
   brought back in f[]. so the result will be good for the
   case that s[] is shorter than f[]
*/
void    conv(float *s, int ns, float *f, int n) {
  int   i,j,k,m;
  float *g,*pt;
  m = n+ns;
  g=(float *) malloc(m*sizeof(float));
  for(i=0;i<ns;i++) g[i]=0.;
  for(pt=f,k=0;k<n;k++,i++,pt++){
     g[i]=*pt;
     for(*pt=0.,j=0;j<ns;j++) *pt += g[i-j]*s[j];
  }
  free(g);
}


/*
 cross-correlate rec with syn: sum(rec[i]*syn[j-i],j). Note no dt.
 only return a m+1 points of cross-correlation around the zero-lag
 (m=2*k; the dedays are -k, -k+1, ..., 0, ..., k-1, k)
*/
float	*crscrl(int npt,float *rec,float *syn,int m) {
  int	i,nft,nft2;
  float	*ss1,*ss2;
  nft=2;while(nft<npt)nft*=2;nft2=nft;nft*=2;
  ss1=(float *)calloc(nft, sizeof(float));
  ss2=(float *)calloc(nft, sizeof(float));
  memcpy((char *)ss1, (char *) rec, npt*sizeof(float));
  memcpy((char *)ss2, (char *) syn, npt*sizeof(float));
  for(i=npt;i<nft;i++) {ss1[i]=0.;ss2[i]=0.;}
  fftr((complex *) ss1,nft2,1.);
  fftr((complex *) ss2,nft2,1.);

  cor((complex *) ss2, (complex *) ss1, 1., nft2);

  nft2 -= m/2;
  i = (m+1)*sizeof(float);
  ss2=(float *)realloc(ss2, i);
  memcpy(ss2, ss1+nft2, i);
  free(ss1);
  return(ss2);
}


/* data(t) = amp*syn(t+delay), return max. cross-correlation */
float maxCor(float *data, float *syn, int n, int *delay, float *amp) {
  int	i,m;
  float	*crs,c,dataAuto,synAuto;
  for(dataAuto=0.,synAuto=0.,i=0;i<n;i++) {
     dataAuto += data[i]*data[i];
     synAuto += syn[i]*syn[i];
  }
  m = 2*n;
  crs = crscrl(n,data,syn,m);
  for(c=-FLT_MAX,i=0;i<=m;i++) {
     if (c<crs[i]) {
        c=crs[i];
        *delay=i;
     }
  }
  free(crs);
  *delay -= n;
  *amp = c/synAuto;
  return(c/sqrt(dataAuto*synAuto));
}


/* integrate data[n] between t1 and t2 (normalized time by dt)*/
float amp(float t1, float t2, float *data, int n) {
  int i, it1, it2;
  float dd, am;
  if ( t1 < 0 ) t1=0;
  if ( t2 > n-1 ) t2=n-1;
  if ( t1 > n-1 || t2 < t1) return 0.;
  it1 = floor(t1);
  i = it1 + 1;
  dd = i-t1;
  am = dd*(dd*data[it1]+(2.-dd)*data[i]);
/* return data[it1]*dd + data[i]*(1.-dd);*/
  it2 = ceil(t2);
  while (i<it2) {
    am += data[i]+data[i+1];
    i++;
  }
  dd = i-t2;
  am -= dd*(dd*data[i-1]+(2.-dd)*data[i]);
  return 0.5*am;
}


/*  cummulative sum of a(t) */
void cumsum(float *a, int n, float dt) {
  int i;
  float u;
  for(u=0,i=0;i<n;i++) {
     u+=a[i]*dt;
     a[i] = u;
  }
}


/* return a symetric taper function of length n
	f(i) = 0.5*(1-cos(i*pi/n1)),	i=0, n1
	f(i) = 1, 			i=n1, n/2
  where n1=w*n, 0<w<0.5
*/
float	*coswndw(int n, float w) {
     int   j, n1;
     float t, dt, *wndw;
     if ( (wndw=(float *)malloc(n*sizeof(float))) == NULL ) return NULL;
     for(j=0;j<n;j++) wndw[j]=1.;
     if (w>0.5) w=0.5;
     n1 = rint(w*n); if (n1<1) n1=1;
     t = 0.;
     dt = PI/n1;
     for(j=0;j<n1;j++,t+=dt) {
        wndw[j]=wndw[n-j-1]=0.5*(1-cos(t));
     }
     return wndw;
}


/* windowing spectrum d[i], sgn=1 -> high-pass; sgn=-1 -> low-pass */
void	filter(complex *d, int n, float f1, float f2, float dt, int sgn) {
     int i, if1, if2;
     float a;
     dt = 0.5/dt/n;
     if1 = rint(f1/dt);
     if2 = rint(f2/dt); if (if2>=n) if2=n-1;
     if (if2<=if1) {
	fprintf(stderr, "filter freq. wrong f2<f1\n");
	return;
     }
     dt = PI/(if2-if1);
     for(a=0.,i=if1;i<if2;i++,a+=dt) d[i] = dmltp(0.5*(1-sgn*cos(a)), d[i]);
     if (sgn<0) {
	for(i=if2;i<n;i++) d[i] = Zero;
	d[0].y = 0.;
     } else  {
        d[0].x = 0.;
        for(i=1;i<if1;i++) d[i] = Zero;
     }
}


/* find max. absolute value and location in an array */
int findMaxAbs(float *cor, int n, float *amp) { 
    int i,shift=0;
    *amp = 0.;
    for(i=0;i<n;i++) {
       if (fabs(cor[i])>fabs(*amp)) {
          *amp = cor[i]; 
          shift = i;
       }
    }
    return shift;
}


/* remove trend a*i + b */
void	rtrend(float *y, int n) {
     int i;
     double a, b, a11, a12, a22, y1, y2;
     y1 = y2 = 0.;
     for(i=0;i<n;i++) {
       y1 += i*y[i];
       y2 += y[i];
     }
     a12 = 0.5*n*(n-1);
     a11 = a12*(2*n-1)/3.;
     a22 = n;
     b = a11*a22-a12*a12;
     a = (a22*y1-a12*y2)/b; 
     b = (a11*y2-a12*y1)/b;
     for(i=0;i<n;i++) {
       y[i] = y[i] - a*i - b;
     }
}


/* multiply exp(-w^2/(4*gauss^2)) to spectrum u, which is equivalent to
   convolve f(t)=(gauss/sqrt(pi))*exp(-(gauss*t)^2) to u(t).
   f(t) has unit area. The input gauss is actually gauss*dt */
void fltGauss(complex *u, int n, float gauss) {
    int j;
    float w, delw, agg;
    delw=PI/n;
    for(w=delw,j=1;j<n;j++,w+=delw) {
	agg = 0.5*w/gauss;
	u[j] = dmltp(exp(-agg*agg),u[j]);
    }
    agg = 0.5*w/gauss;
    u[0].y = (exp(-agg*agg))*u[0].y;
}


/* compute f(t)=u(t-shift*dt), which is equ. to multiply
exp(-j*w*shift) to its spectrum u*/
void shiftSpec(complex *u, int n, float shift) {
    int j;
    float w, delw;
    delw=PI/n;
    for(w=delw,j=1;j<n;j++,w+=delw) {
	u[j] = cmltp(cmplx(cos(w*shift),-sin(w*shift)),u[j]);
    }
    u[0].y = cos(w*shift)*u[0].y;
}


/* compute spectrum power which equals auto_correlation*dt */
float specPwr(complex *u, int n) {
  int j;
  float a, temp;
  for(a=0.,j=0;j<n;j++) {
     temp = ccabs(u[j]);
     a += temp*temp;
  }
  return (a/n);
}


/* add spectrums a=a+b */
void specAdd(complex *a, complex *b, int n) {
  int j;
  for(j=0;j<n;j++) a[j] = cplus(a[j], b[j]);
}


/* specMul a=a*b */
void specMul(complex *a, complex *b, int n) {
  int	j;
  a[0]=cmplx(a[0].x*b[0].x, a[0].y*b[0].y);
  for (j=1; j<n; j++) a[j]=cmltp(a[j], b[j]);
}


/* multiply spectrum a by a constant c */
void specScale(complex *a, float c, int n) {
  int j;
  for(j=0;j<n;j++) a[j] = dmltp(c, a[j]);
}
