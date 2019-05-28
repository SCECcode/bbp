// VerticalSHResponse.h: CVerticalSHResponse クラスのインターフェイス
//
// 多層構造に垂直入射する波に対する周波数応答を計算するクラスです。
// 層の境界に置いて、変位と応力の連続境界条件を解くことによって求めています。
//
// コンストラクタで媒質構造を与えてオブジェクトを構築し、
// CalcResponse()でΔt、サンプル数を与えて周波数応答を計算します。
// その後、任意の層番号を指定して、GetSpectrum()を通してstd::vector< std::complex<double> >に
// 応答を取得してください。
// 一度構築した周波数応答は破棄するまで有効です。
//
// 媒質構造は波の位相速度、密度、層厚で規定されます。物性値の単位系は揃えてください。
//
// 参考文献：
// 佐藤泰夫、1978、弾性波動論、岩波書店、118--120
//
// Ver. 1.02  2002 S.Senna(MSS)
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_VERTICALSHRESPONSE_H__79A9A070_6A27_11D5_A8B5_00508BCDC8C2__INCLUDED_)
#define AFX_VERTICALSHRESPONSE_H__79A9A070_6A27_11D5_A8B5_00508BCDC8C2__INCLUDED_

#if _MSC_VER > 1000
//#pragma once
#endif // _MSC_VER > 1000

//#pragma warning(disable:4786)

#include <vector>
#include <complex>


class CVerticalSHResponse
{
public:
	CVerticalSHResponse( const std::vector<double>& velocity, const std::vector<double>& density, const std::vector<double>& thickness);
	virtual ~CVerticalSHResponse();

	void CalcResponse( unsigned int numofSampling, double delta_time );
	void GetSpectrum( unsigned int layerNo, std::vector< std::complex<double> >& spectrum, double delta_time );
	void GetWave( unsigned int layerNo, const std::vector<double>& inWave, std::vector<double>& outWave );
    
	static inline std::complex<double> exp( std::complex<double> x ){ static const double base_e = 2.718281828459045; return pow( base_e, x ); }
	//static inline std::complex<double> exp( const std::complex<double> *x ){ static const double base_e = 2.718281828459045; return std::pow( base_e, x ); }

	double DeltaT(){ return m_dt; }
	unsigned int NumofSample(){ return m_a[0].size(); }

private:
	const std::vector<double>& m_velocity;
	const std::vector<double>& m_density;
	const std::vector<double>& m_thickness;

	std::vector< std::vector< std::complex<double> > > m_a; //上向き成分フーリエ係数　m_a[層番号][時刻or周波数]
	std::vector< std::vector< std::complex<double> > > m_b; //下向き成分フーリエ係数　m_b[層番号][時刻or周波数]

	double m_dt;
};

#endif // !defined(AFX_VERTICALSHRESPONSE_H__79A9A070_6A27_11D5_A8B5_00508BCDC8C2__INCLUDED_)
