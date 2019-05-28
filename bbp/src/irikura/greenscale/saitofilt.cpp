#define		PI		3.141593
#define		HP		1.570796
#include	<math.h>

#include "saitofilt.h"
/*
+      BUTTERWORTH LOW PASS FILTER COEFFICIENT
+
+      ARGUMENTS
+        H      : FILTER COEFFICIENTS
+        M      : ORDER OF FILTER  (M=(N+1)/2)
+        GN     : GAIN FACTOR
+        N      : ORDER OF BUTTERWORTH FUNCTION
+        FP     : PASS BAND FREQUENCY  (NON-DIMENSIONAL)
+        FS     : STOP BAND FREQUENCY
+        AP     : MAX. ATTENUATION IN PASS BAND
+        AS     : MIN. ATTENUATION IN STOP BAND
+
+      M. SAITO  (17/XII/75)
*/
int butlop(double *h,int *m,double *gn,int *n,double fp,double fs,double ap,double as)
{
	//double *h,fp,fs,ap,as,*gn;
	//int *m,*n;
	{
	double wp,ws,tp,ts,pa,sa,cc,c,dp,g,fj,c2,sj,tj,a;
	int k,j;
	if(fabs(fp)<fabs(fs)) wp=fabs(fp)*PI;
	else wp=fabs(fs)*PI;
	if(fabs(fp)>fabs(fs)) ws=fabs(fp)*PI;
	else ws=fabs(fs)*PI;
	if(wp==0.0 || wp==ws || ws>=HP)
		{
		printf("? (butlop) invalid input : fp=%14.6e fs=%14.6e ?\n",fp,fs);
		return(1);
		}
/****  DETERMINE N & C */
	tp=tan(wp);
	ts=tan(ws);
	if(fabs(ap)<fabs(as)) pa=fabs(ap);
	else pa=fabs(as);
	if(fabs(ap)>fabs(as)) sa=fabs(ap);
	else sa=fabs(as);
	if(pa==0.0) pa=0.5;
	if(sa==0.0) sa=5.0;
    int nn;
    double ddd;
    ddd = fabs(log(pa/sa)/log(tp/ts))+0.5;
    nn = (int)(fabs(log(pa/sa)/log(tp/ts))+0.5);
	if((*n=(int)(fabs(log(pa/sa)/log(tp/ts))+0.5))<2){
        *n=2;
    }
	cc=exp(log(pa*sa)/(double)(*n))/(tp*ts);
	c=sqrt(cc);

	dp=HP/(double)(*n);
	*m=(*n)/2;
	k=(*m)*4;
	g=fj=1.0;
	c2=2.0*(1.0-c)*(1.0+c);

	for(j=0;j<k;j+=4){
		sj=pow(cos(dp*fj),2.0);
		tj=sin(dp*fj);
		fj=fj+2.0;
		a=1.0/(pow(c+tj,2.0)+sj);
		g=g*a;
		h[j  ]=2.0;
		h[j+1]=1.0;
		h[j+2]=c2*a;
		h[j+3]=(pow(c-tj,2.0)+sj)*a;
    }
/****  EXIT */
	*gn=g;
	if(*n%2==0) return(0);
/****  FOR ODD N */
	*m=(*m)+1;
	*gn=g/(1.0+c);
	h[k  ]=1.0;
	h[k+1]=0.0;
	h[k+2]=(1.0-c)/(1.0+c);
	h[k+3]=0.0;
	return(0);
	}
}
/*
+      BUTTERWORTH HIGH PASS FILTER COEFFICIENT
+
+      ARGUMENTS
+        H      : FILTER COEFFICIENTS
+        M      : ORDER OF FILTER  (M=(N+1)/2)
+        GN     : GAIN FACTOR
+        N      : ORDER OF BUTTERWORTH FUNCTION
+        FP     : PASS BAND FREQUENCY  (NON-DIMENSIONAL)
+        FS     : STOP BAND FREQUENCY
+        AP     : MAX. ATTENUATION IN PASS BAND
+        AS     : MIN. ATTENUATION IN STOP BAND
+
+      M. SAITO  (7/I/76)
*/
int buthip(double *h,int *m,double *gn,int *n,double fp,double fs,double ap,double as)
{
	double wp,ws,tp,ts,pa,sa,cc,c,dp,g,fj,c2,sj,tj,a;
	int k,j;
	if(fabs(fp)>fabs(fs)) wp=fabs(fp)*PI;
	else wp=fabs(fs)*PI;
	if(fabs(fp)<fabs(fs)) ws=fabs(fp)*PI;
	else ws=fabs(fs)*PI;
	if(wp==0.0 || wp==ws || wp>=HP)
		{
		printf("? (buthip) invalid input : fp=%14.6e fs=%14.6e ?\n",fp,fs);
		return(1);
		}
/****  DETERMINE N & C */
	tp=tan(wp);
	ts=tan(ws);
	if(fabs(ap)<fabs(as)) pa=fabs(ap);
	else pa=fabs(as);
	if(fabs(ap)>fabs(as)) sa=fabs(ap);
	else sa=fabs(as);
	if(pa==0.0) pa=0.5;
	if(sa==0.0) sa=5.0;
	if((*n=(int)(fabs(log(sa/pa)/log(tp/ts))+0.5))<2) *n=2;
	cc=exp(log(pa*sa)/(double)(*n))*(tp*ts);
	c=sqrt(cc);

	dp=HP/(double)(*n);
	*m=(*n)/2;
	k=(*m)*4;
	g=fj=1.0;
	c2=(-2.0)*(1.0-c)*(1.0+c);

	for(j=0;j<k;j+=4){
		sj=pow(cos(dp*fj),2.0);
		tj=sin(dp*fj);
		fj=fj+2.0;
		a=1.0/(pow(c+tj,2.0)+sj);
		g=g*a;
		h[j  ]=(-2.0);
		h[j+1]=1.0;
		h[j+2]=c2*a;
		h[j+3]=(pow(c-tj,2.0)+sj)*a;
		}
/****  EXIT */
	*gn=g;
	if(*n%2==0) return(0);
/****  FOR ODD N */
	*m=(*m)+1;
	*gn=g/(c+1.0);
	h[k  ]=(-1.0);
	h[k+1]=0.0;
	h[k+2]=(c-1.0)/(c+1.0);
	h[k+3]=0.0;
	return(0);
}
/*
+      BUTTERWORTH BAND PASS FILTER COEFFICIENT
+
+      ARGUMENTS
+        H      : FILTER COEFFICIENTS
+        M      : ORDER OF FILTER
+        GN     : GAIN FACTOR
+        N      : ORDER OF BUTTERWORTH FUNCTION
+        FL     : LOW  FREQUENCY CUT-OFF  (NON-DIMENSIONAL)
+        FH     : HIGH FREQUENCY CUT-OFF
+        FS     : STOP BAND FREQUENCY
+        AP     : MAX. ATTENUATION IN PASS BAND
+        AS     : MIN. ATTENUATION IN STOP BAND
+
+      M. SAITO  (7/I/76)
*/
int butpas(double *h,int *m,double *gn,int *n,double fl,double fh,double fs,double ap,double as)
{
	//double *h,fl,fh,fs,ap,as,*gn;
	//int *m,*n;
	double wl,wh,ws,clh,op,ww,ts,os,pa,sa,cc,c,dp,g,fj,rr,tt,
		re,ri,a,wpc,wmc;
	int k,l,j,i;
	struct {
		double r;
		double c;
    } oj,aa,cq,r[2];
	if(fabs(fl)<fabs(fh)) wl=fabs(fl)*PI;
	else wl=fabs(fh)*PI;
	if(fabs(fl)>fabs(fh)) wh=fabs(fl)*PI;
	else wh=fabs(fh)*PI;
	ws=fabs(fs)*PI;
	if(wl==0.0 || wl==wh || wh>=HP || ws==0.0 || ws>=HP ||
			(ws-wl)*(ws-wh)<=0.0){
		printf("? (butpas) invalid input : fl=%14.6e fh=%14.6e fs=%14.6e ?\n",
			fl,fh,fs);
		*m=0;
		*gn=1.0;
		return(1);
    }
/****  DETERMINE N & C */
	clh=1.0/(cos(wl)*cos(wh));
	op=sin(wh-wl)*clh;
	ww=tan(wl)*tan(wh);
	ts=tan(ws);
	os=fabs(ts-ww/ts);
	if(fabs(ap)<fabs(as)) pa=fabs(ap);
	else pa=fabs(as);
	if(fabs(ap)>fabs(as)) sa=fabs(ap);
	else sa=fabs(as);
	if(pa==0.0) pa=0.5;
	if(sa==0.0) sa=5.0;
	if((*n=(int)(fabs(log(pa/sa)/log(op/os))+0.5))<2) *n=2;
	cc=exp(log(pa*sa)/(double)(*n))/(op*os);
	c=sqrt(cc);
	ww=ww*cc;

	dp=HP/(double)(*n);
	k=(*n)/2;
	*m=k*2;
	l=0;
	g=fj=1.0;

	for(j=0;j<k;j++){
		oj.r=cos(dp*fj)*0.5;
		oj.c=sin(dp*fj)*0.5;
		fj=fj+2.0;
		aa.r=oj.r*oj.r-oj.c*oj.c+ww;
		aa.c=2.0*oj.r*oj.c;
		rr=sqrt(aa.r*aa.r+aa.c*aa.c);
		tt=atan(aa.c/aa.r);
		cq.r=sqrt(rr)*cos(tt/2.0);
		cq.c=sqrt(rr)*sin(tt/2.0);
		r[0].r=oj.r+cq.r;
		r[0].c=oj.c+cq.c;
		r[1].r=oj.r-cq.r;
		r[1].c=oj.c-cq.c;
		g=g*cc;

		for(i=0;i<2;i++){
			re=r[i].r*r[i].r;
			ri=r[i].c;
			a=1.0/((c+ri)*(c+ri)+re);
			g=g*a;
			h[l  ]=0.0;
			h[l+1]=(-1.0);
			h[l+2]=2.0*((ri-c)*(ri+c)+re)*a;
			h[l+3]=((ri-c)*(ri-c)+re)*a;
			l=l+4;
        }
    }
/****  EXIT */
	*gn=g;
	if(*n==(*m)) return(0);
/****  FOR ODD N */
	*m=(*m)+1;
	wpc=  cc *cos(wh-wl)*clh;
	wmc=(-cc)*cos(wh+wl)*clh;
	a=1.0/(wpc+c);
	*gn=g*c*a;
	h[l  ]=0.0;
	h[l+1]=(-1.0);
	h[l+2]=2.0*wmc*a;
	h[l+3]=(wpc-c)*a;
	return(0);
}
/*
+      RECURSIVE FILTERING : F(Z) = (1+A*Z+AA*Z**2)/(1+B*Z+BB*Z**2)
+
+      ARGUMENTS
+        X      : INPUT TIME SERIES
+        Y      : OUTPUT TIME SERIES  (MAY BE EQUIVALENT TO X)
+        N      : LENGTH OF X & Y
+        H      : FILTER COEFFICIENTS ; H(1)=A, H(2)=AA, H(3)=B, H(4)=BB
+        NML    : >0 ; FOR NORMAL  DIRECTION FILTERING
+                 <0 ; FOR REVERSE DIRECTION FILTERING
+		 uv		: past data and results saved
+
+      M. SAITO  (6/XII/75)
*/
int recfil(double *x,double *y,int n,double *h,int nml)
{
	//int n,nml;
	//double *x,*y,*h;
	int i,j,jd;
	double a,aa,b,bb,u1,u2,u3,v1,v2,v3;
	if(n<=0){
		printf("? (recfil) invalid input : n=%d ?\n",n);
		return(1);
    }
	if(nml>=0){
		j=0;
		jd=1;
    }else{
		j=n-1;
		jd=(-1);
    }
	a =h[0];
	aa=h[1];
	b =h[2];
	bb=h[3];
	u1 = u2 = v1 = v2 = 0.0;
/****  FILTERING */
	for(i=0;i<n;i++){
		u3=u2;
		u2=u1;
		u1=x[j];
		v3=v2;
		v2=v1;
		v1=u1+a*u2+aa*u3-b*v2-bb*v3;
		y[j]=v1;
		j+=jd;
    }
	return(0);
}
/*
+      RECURSIVE FILTERING IN SERIES
+
+      ARGUMENTS
+        X      : INPUT TIME SERIES
+        Y      : OUTPUT TIME SERIES  (MAY BE EQUIVALENT TO X)
+        N      : LENGTH OF X & Y
+        H      : COEFFICIENTS OF FILTER
+        M      : ORDER OF FILTER
+        NML    : >0 ; FOR NORMAL  DIRECTION FILTERING
+                 <0 ;     REVERSE DIRECTION FILTERING
+		 uv		: past data and results saved
+
+      SUBROUTINE REQUIRED : RECFIL
+
+      M. SAITO  (6/XII/75)
*/
int tandem(double *x,double *y,int n,double *h,int m,int nml)
{
	//double *x,*y,*h;
	//int n,m,nml;
	int i;
	if(n<=0 || m<=0){
		printf("? (tandem) invalid input : n=%d m=%d ?\n",n,m);
		return(1);
    }
/****  1-ST CALL */
	recfil(x,y,n,h,nml);
/****  2-ND AND AFTER */
	if(m>1){
        for(i=1;i<m;i++){
            recfil(y,y,n,&h[i*4],nml);
        }
    }
	return(0);
}
