// FastFourierTransform.cpp: CFastFourierTransform クラスのインプリメンテーション
//
//////////////////////////////////////////////////////////////////////

#include "FastFourierTransform.h"

//////////////////////////////////////////////////////////////////////
// 構築/消滅
//////////////////////////////////////////////////////////////////////

CFastFourierTransform::CFastFourierTransform()
{

}

CFastFourierTransform::~CFastFourierTransform()
{

}

void CFastFourierTransform::Trans( std::vector< std::complex<double> >& cx, int ind )
{
	int i, j, k, m;
	std::complex<double> ctemp, ctheta;
	static const double PI = 3.141592653589793238;

	int n = cx.size();

	for( i=1, j=1; i<=n; i++, j+=m ){
		if (i<j) std::swap( cx[j-1], cx[i-1] );
		m = n/2;
  		do{
			if( j <= m ) break;
			j -= m;
			m /= 2;
		}while( m >= 2 );
	}

	int istep;
	for( int kmax = 1; kmax < n; ){
		istep = kmax*2;
		for( k=1; k<=kmax; k++ ){
			ctheta = std::complex<double> ( 0.0, PI*ind*(k-1)/kmax );
			for( i=k; i<=n; i+=istep ){
				j=i+kmax;
				if(( j>n )|| (i>n)){
				}
				ctemp=cx[j-1]*exp(ctheta);
				cx[j-1] = cx[i-1]-ctemp;
				cx[i-1] = cx[i-1]+ctemp;
			}
		}
		kmax=istep;
	}
}

void CFastFourierTransform::Inverse( std::vector< std::complex<double> >& x )
{
	int nt = x.size();
	std::vector< std::complex<double> > cwt(x.size());

	cwt[0]= x[0];
	for( int istep=1; istep<nt/2; istep++ ){
		cwt[istep]=x[istep];
		cwt[nt-istep]=std::conj( x[istep] );
	}
	cwt[nt/2] = x[nt/2];

	CFastFourierTransform::Trans(cwt,-1);

	x = cwt;
}
