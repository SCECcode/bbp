// VerticalSHResponse.cpp: CVerticalSHResponse クラスのインプリメンテーション
//
//////////////////////////////////////////////////////////////////////

#include "VerticalSHResponse.h"
#include "BasenSpectrum.h"
#include "complex.h"
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <numeric>
#include <complex>

//////////////////////////////////////////////////////////////////////
// 構築/消滅
//////////////////////////////////////////////////////////////////////

CVerticalSHResponse::CVerticalSHResponse( 
	const std::vector<double>& velocity,
	const std::vector<double>& density,
	const std::vector<double>& thickness) : m_velocity( velocity ), m_density( density ), m_thickness( thickness )    
{
	m_a.resize( velocity.size() );
	m_b.resize( velocity.size() );
}

CVerticalSHResponse::~CVerticalSHResponse()
{

}

void CVerticalSHResponse::CalcResponse( unsigned int numofSampling, double delta_time )
{
	static const double PI = 3.141592653589793238; //円周率
	static const std::complex<double> IMAG = std::complex<double>(0.0, 1.0); //虚数単位

	int                  i, k;          // ループカウンタ
	double               w, x;          // テンポラリ
	std::complex<double> c, d;          // テンポラリ
	double               theta,impd1,impd2; // インピーダンス
    
	//係数ベクタをリサイズ
	for( i = 0; i < m_a.size(); i++ ) m_a[i].resize( numofSampling, std::complex<double>(0.0,0.0) );
	for( i = 0; i < m_b.size(); i++ ) m_b[i].resize( numofSampling, std::complex<double>(0.0,0.0) );

	//計算開始
	int ln = m_a.size();
	for( i = 0;  i <= numofSampling/2;  i++ )
	{
		w = (double)( 2.0 * PI * i ) / ( numofSampling * delta_time );
		double pi = 4.0 * atan(1.0);

			//地表面
		m_a[0][i] = std::complex<double>( 1.0, 0.0 );
		m_b[0][i] = std::complex<double>( 0.0, 0.0 );

		for( k = 1; k < ln ; k++ )
		{
			theta = w / (m_velocity[k-1]) * m_thickness[k-1];
			impd1 = m_density[k-1] * m_velocity[k-1] * 1000 ;
			impd2 = m_density[k  ] * m_velocity[k  ] * 1000 ;
			m_a[k][i] =          m_a[k-1][i] * cos( theta )  + IMAG * m_b[k-1][i] * sin( theta );
			m_b[k][i] = ( IMAG * m_a[k-1][i] * sin( theta )  +        m_b[k-1][i] * cos( theta ) ) * impd1 / impd2;		
		}

	//cos, sinをexpに変換
		for( k = 0; k < ln ; k++ )
		{
			c = m_a[k][i];
			d = m_b[k][i];
			m_a[k][i] = (c+d)/2.0;
			m_b[k][i] = (c-d)/2.0;
		}

		//入射波で規格化
		for( k = 0; k < ln ; k++ )
		{
			m_a[k][i] /= m_b[ln-1][i];
			m_b[k][i] /= m_b[ln-1][i];
		}

	for( k = 0, x = 0.0; k < ln ; k++ )
	{
			//負の周波数に共役な値をセット
			if( i != 0 )
			{
				m_a[k][numofSampling-i] = std::conj( m_a[k][i] );
				m_b[k][numofSampling-i] = std::conj( m_b[k][i] );
			}
	}
}
}

void CVerticalSHResponse::GetSpectrum( unsigned int layerNo, std::vector< std::complex<double> >& spectrum, double delta_time )
{
	int size =  m_a[layerNo].size();
    int i;
	spectrum.resize( size );
	for( i = 0 ; i < size ; i++ ){
		spectrum[i] = ( m_a[layerNo][i] + m_b[layerNo][i] )*( 8.0 + log10( delta_time * (double)i+1.0 ))/10.0;
	}
    //printf("!!!");
}
