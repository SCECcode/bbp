// BasenSpectrum.h: CBasenSpectrum クラスのインターフェイス
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_BASENSPECTRUM_H__E25839DF_6A11_11D5_B427_FFFFFF000000__INCLUDED_)
#define AFX_BASENSPECTRUM_H__E25839DF_6A11_11D5_B427_FFFFFF000000__INCLUDED_

#if _MSC_VER > 1000
//#pragma once
#endif // _MSC_VER > 1000

#include <complex>
#include <vector>

class CQvalue
{
private:
	double	m_dDepend_1, m_dDepend_2, m_dDepend_3;
	double	m_dF;
public:
	CQvalue();
	~CQvalue();
	void Set(double b1, double b2, double f, double b3);
	double GetQvalue( double f );
	void GetQvalueVector(std::vector<double> fvec, std::vector<double>& qval );
};

class CBasenSpectrum  
{
private:
	double	m_dRo_pq;	// 要素断層における地殻の密度(g/cm^3)
	double	m_dRo_sb;	// 要素断層における地震基盤の密度(g/cm^3)
	double	m_dBeta_pq;	// 要素断層における剪断波速度(km/s)
	double	m_dBeta_sb;	// 要素断層における地震基盤の剪断波速度(km/s)
	double	m_dM0_pq;	// 地震モーメント(dyne-cm)
    double  m_dRs;      // 震源アスペリティ―からの観測距離(km)
	double	m_dSig_pq;	// 実効応力(bar)
	double	m_dFrad;	// ラディエーションパターン係数
	double	m_dFcor_pq;	// ブルーンモデルのコーナ周波数
	double	m_dFmax_pq;	// fmax
	double	m_dM;		// 係数
	double  m_dMw_pq;   // モーメントマグニチュード
    double  m_dAvs;     // 地表から30mの平均S波
    double  m_dDslip;   // 最終すべり量(m)
	double  m_dVr;      // 破壊伝播速度(km/s)
    double  m_dBu;      // 媒質の剛性率(GPa)
	double  m_dVpq;     // 最大すべり速度(m/s)

private:
	void init();

public:
	CQvalue	m_cQval;

	CBasenSpectrum();
   ~CBasenSpectrum();

	void GetEnvelope(double dR_pq, int nSampleFrequency, int nSamplePoint, std::vector<double>& env ,double& td ,double& envv );
	void GetBasenSpectrum(double dR_pq, int nSampleFrequency, int nSamplePoint, std::vector<double>& frq, std::vector<double>& spc );
	void GetBasenWave( double dR_pq, int nSampleFrequency, int nSamplePoint, std::vector<double>& wave, std::vector<double>& wave2, std::vector<double> SAngle, int Ssize );
	void GetBasenCalc(	std::vector< std::complex<double> > wwave, std::vector< std::complex<double> >& wwave2, std::vector<double>& wave2, std::vector<double> acc, double td, int nSampleFrequency, int nSamplePoint );

public:
	// setter
	void SetRoPQ(double dRo_pq){m_dRo_pq=dRo_pq*1000.0;}
	void SetRoSB(double dRo_sb){m_dRo_sb=dRo_sb*1000.0;}
	void SetBetaPQ(double dBeta_pq){m_dBeta_pq=dBeta_pq*1000.0;}
	void SetBetaSB(double dBeta_sb){m_dBeta_sb=dBeta_sb*1000.0;}
	void SetM0PQ(double dM0_pq){m_dM0_pq=dM0_pq;}
	void SetRs(double dRs){m_dRs=dRs*1000.0;}
	void SetStress(double dSig){m_dSig_pq=dSig;}
	void SetRadiation(double dFrad){m_dFrad=dFrad;}
	void SetFMax(double dFMax){m_dFmax_pq=dFMax;}
	void SetFCor(double dFCor){m_dFcor_pq=dFCor;}
	void SetMw(double dMw){m_dMw_pq = dMw;}
	void SetMavs(double dAvs){m_dAvs = dAvs;}
	void SetM(double m){m_dM=m;}
    void SetDslip(double dDslip){m_dDslip = dDslip;}
	void SetVr(double dVr){m_dVr = dVr;}
    void SetBu(double dBu){m_dBu = dBu;}
	void SetVpq(double dVpq){m_dVpq = dVpq;}

	// getter
	double GetRoPQ(){return m_dRo_pq/1000.0;}
	double GetRoSB(){return m_dRo_sb/1000.0;}
	double GetBetaPQ(){return m_dBeta_pq/1000.0;}
	double GetBetaSB(){return m_dBeta_sb/1000.0;}
	double GetM0PQ(){return m_dM0_pq;}
	double GetRs(){return m_dRs/1000.0;}
	double GetStress(){return m_dSig_pq;}
	double GetRadiation(){return m_dFrad;}
	double GetFCor(){return m_dFcor_pq;}
	double GetFMax(){return m_dFmax_pq;}
    double GetMw(){return m_dMw_pq;}
	double Getavs(){return m_dAvs;}
	double GetM(){return m_dM;}
    double GetDslip(){return m_dDslip;}
	double GetVr(){return m_dVr;}
    double GetBu(){return m_dBu;}
	double GetVpq(){return m_dVpq;}
};

#endif // !defined(AFX_BASENSPECTRUM_H__E25839DF_6A11_11D5_B427_FFFFFF000000__INCLUDED_)
