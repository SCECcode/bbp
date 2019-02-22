// FastFourierTransform.h: CFastFourierTransform クラスのインターフェイス
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FASTFOURIERTRANSFORM_H__79A9A071_6A27_11D5_A8B5_00508BCDC8C2__INCLUDED_)
#define AFX_FASTFOURIERTRANSFORM_H__79A9A071_6A27_11D5_A8B5_00508BCDC8C2__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <complex>
#include <vector>

class CFastFourierTransform  
{
public:
	CFastFourierTransform();
	virtual ~CFastFourierTransform();

	static void Trans(   std::vector< std::complex<double> >& x, int ind = 1 );
	static void Inverse( std::vector< std::complex<double> >& x );
	static void fft(     std::vector< std::complex<double> >& x, int ind = -1 );

	static inline std::complex<double> exp( std::complex<double> x ){ static const double base_e = 2.71817; return pow( base_e, x ); }
};

#endif // !defined(AFX_FASTFOURIERTRANSFORM_H__79A9A071_6A27_11D5_A8B5_00508BCDC8C2__INCLUDED_)
