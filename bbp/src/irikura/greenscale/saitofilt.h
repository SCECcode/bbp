#include <stdio.h>
#include <math.h>
#include <complex>
#include <vector>
#include <stdlib.h>
#include <string.h>

//void butlop(h,m,gn,n,fp,fs,ap,as);
//void buthip(h,m,gn,n,fp,fs,ap,as);
//void butpas(h,m,gn,n,fl,fh,fs,ap,as);
//void tandem(x,y,n,h,m,nml);
//void recfil(x,y,n,h,nml);

int butlop(double *h,int *m,double *gn,int *n,double fp,double fs,double ap,double as);
int recfil(double *x,double *y,int n,double *h,int nml);
int tandem(double *x,double *y,int n,double *h,int m,int nml);
int butpas(double *h,int *m,double *gn,int *n,double fl,double fh,double fs,double ap,double as);
int buthip(double *h,int *m,double *gn,int *n,double fp,double fs,double ap,double as);
