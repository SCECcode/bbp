// BasenSpectrum.cpp: CBasenSpectrum クラスのインプリメンテーション
//
//////////////////////////////////////////////////////////////////////
#include "knetascii.h"
#include "BasenSpectrum.h"
//#include "four1.h"
#include "FastFourierTransform.h"
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <numeric>
#include <complex>

#include <vector>
#include <algorithm>

using namespace std;

double td = 0;
double envv = 0;

//////////////////////////////////////////////////////////////////////
// 構築/消滅
//////////////////////////////////////////////////////////////////////

CQvalue::CQvalue()
{
}

CQvalue::~CQvalue()
{
}

void CQvalue::Set(double b1, double b2, double f, double b3)
{
	m_dDepend_1 = b1;
	m_dDepend_2 = b2;
	m_dDepend_3 = b3;
	m_dF = f;
}

double CQvalue::GetQvalue( double f )
{

	if( m_dF < f )
		return m_dDepend_1*pow(f, m_dDepend_2);
	else
		return m_dDepend_3;
}

void CQvalue::GetQvalueVector(std::vector<double> fvec, std::vector<double>& qval )
{
	int	i;
	qval.clear();

	for(i=0;i<qval.size();i++){
		if( m_dF < fvec[i] ) qval.push_back( m_dDepend_1 * pow(fvec[i], m_dDepend_2));
		else qval.push_back( m_dDepend_3 );
		}

}

CBasenSpectrum::CBasenSpectrum()
{
	init();
}

CBasenSpectrum::~CBasenSpectrum()
{

}

void CBasenSpectrum::init()
{
  // 全てをSI単位系に直す
	m_dRo_pq *= 1000.0;
	m_dRo_sb *= 1000.0;
	m_dBeta_pq *= 1000.0;
	m_dBeta_sb *= 1000.0;
    m_dVr *= 1000.0;
	m_dBu *= 1e9;

}

void CBasenSpectrum::GetBasenSpectrum(
									  double dR_pq,				// km
									  int nSampleFrequency,
									  int nSamplePoint,
									  std::vector<double>& frq,
									  std::vector<double>& spc
									  )
{
	frq.clear();
	spc.clear();

	// 距離をSI単位に直す
	double	dR = dR_pq * 1000.0;
	double	pi = 4.0 * atan(1.0);
	double	pi_2= 2.0 * pi;
	double	sqrpi = sqrt(pi);

	double	dNyquist = (double)nSampleFrequency/2.0;
	double	df = dNyquist/(double)nSamplePoint;
	int	i;
	// コーナ周波数は剪断応力とモーメントの経験式で求める
	double	c1 = pow(7.0/16.0, 1.0/6.0);
	double	fcpq = c1/sqrpi*m_dBeta_pq*pow(m_dSig_pq/m_dM0_pq, 1.0/3.0);
	//double	fcpq2 = c1/sqrpi*m_dBeta_pq*pow(m_dSig_pq*10/m_dM0_pq, 1.0/3.0);

	for(i=0;i<nSamplePoint;i++){
		// 周波数
		double	f = df * (double)i;
		// Q値
		double dQ;
		dQ = m_cQval.GetQvalue(f);

		// ラディエーションパターンの影響
		double	a1 = m_dFrad/(4*pi*m_dRo_pq*pow(m_dBeta_pq, 3.0));
		// ブルーンの震源加速度スペクトル
		double	a2 = m_dM0_pq*pow(pi_2*f, 2.0)/(1.0+pow(f/fcpq , 2.0));
		// fmaxの影響
		double	a3 = 1.0 / sqrt(1.0 + pow(f/m_dFmax_pq,m_dM));
		// 距離減衰
		double	a4 = 1.0 / dR;

		double	a5 = exp( -pi * f * dR /( dQ * m_dBeta_pq ));

		double	a6 = 2.0 * sqrt((m_dRo_pq*m_dBeta_pq)/(m_dRo_sb*m_dBeta_sb));

                double a7 = 0.0;
		if( f > 50.0 ) a7 = 0.57;
		else if( f <= 50.0 && f > 11.1 ) a7 = 0.57;
		else if( f <= 11.1 && f > 7.69 ) a7 = 0.57;
		else if( f <= 7.69 && f > 3.33 ) a7 = 0.57;
		else if( f <= 3.33 && f > 1.66 ) a7 = 0.57;
		else if( f <= 1.66 && f > 1.0 ) a7 = 0.57;
		else if( f <= 1.0 && f > 0.5 ) a7 = 0.57;
		else if( f <= 0.5 && f > 0.2 ) a7 = 0.57;
		else if( f <= 0.2 ) a7 = 0.57;

		/* TODO : 100? 加速度スペクトルだからSI m -> cm　にした。*/
		double	s  = 100 * a1 * a2 * a3 * a4 * a5 * a6 * a7;

        //////////////
        //printf("freq= %f  :  spec= %f \n",f,s );
        //////////////
        
		frq.push_back(f);
		spc.push_back(s);

	}

}

void CBasenSpectrum::GetBasenWave(
								  double dR_pq,
								  int nSampleFrequency,
								  int nSamplePoint,
								  std::vector<double>& wave,
								  std::vector<double>& wave2,
								  std::vector<double> SAngle,
								  int Ssize
								  )
{
	wave.clear();

	std::vector<double> frq, acc, env, wavv1;
//	double	pi = 4.0 * atan(1.0);
	//double* wav = (double*)malloc(sizeof(double)*nSamplePoint*2);
	//double* wav2 = (double*)malloc(sizeof(double)*nSamplePoint*2);
	std::vector<std::complex<double> > wwave;
	std::vector<std::complex<double> > wwave2;
    //<-add
	std::vector<std::complex<double> > wwav(nSamplePoint);
    int k;
    //2006/01/17 ->
	int i,j;
	int nSample = nSamplePoint/2 + 1;
	/*----------------------------------------------------------------*/

	// 理論的な地震基盤でのスペクトルを計算する
	GetBasenSpectrum(dR_pq, nSampleFrequency, nSample, frq, acc);

	// エンベロープ関数を求める
	GetEnvelope(dR_pq, nSampleFrequency, nSamplePoint, env , td , envv );

	// 加速度フーリエスペクトルに割り当てる
	for(i=0;i<nSample;i++)
	{
		//wav[2*i] = acc[i]*cos(SAngle[i]);
		//wav[2*i+1]= acc[i]*sin(SAngle[i]);
        wwav[i] = complex<double>(acc[i]*cos(SAngle[i]),acc[i]*sin(SAngle[i]));
	}
    //wwav[nSample-1]=complex<double>(1.0,0.0);
	// 負の周波数成分には共役複素数を割り当てる
    /*
	for(i=nSample-2, j=nSample*2 ;i>=1;i--, j+=2)
	{
		wav[j  ] = wav[2*i  ];
		wav[j+1]= -wav[2*i+1];
    }
	for(i=nSample-2, k=0 ; i>=1 ; i--, k++)
    {
        wwav[nSample + k] = conj(wwav[i]);
	}*/

	// 時間領域に戻す
	//four1(wav-1, nSamplePoint, -1);
    /* (wav-1) ? 引数をあわせるため　Fortran->C */
    CFastFourierTransform::Inverse( wwav );
	// 実部を波形変数に格納する
	for(i=0;i<nSamplePoint;i++)
	{
		//wave.push_back((env[i]*wav[2*i])/(double)nSamplePoint );
		wave.push_back((env[i]*real(wwav[i]))/(double)nSamplePoint );
		wwave.push_back(wave[i]);
        //printf(" %lf\n ",wave[i]);
	}
        //printf(" end\n ");
/////////////
    //出力
    //double maxval1 = *std::max_element( wave.begin(), wave.end() );
    //double minval1 = *std::min_element( wave.begin(), wave.end() );
    //double absmax1 = 0;

    //if(fabs( maxval1 ) > fabs( minval1 )) absmax1 = fabs( maxval1 );
    //else absmax1 = fabs( minval1 );

    //printf("SB-WAVE MAX1(gal) = %f\n ",absmax1 );
/////////////


	// 位相割り当ての繰り返し
	for(int iii = 0 ; iii < 1 ; iii++ ) //100回繰り返し
	{
		wave2.clear();
		GetBasenCalc( wwave, wwave2, wave2, acc, td, nSampleFrequency, nSamplePoint );
/////////////
        // 地震基盤の波形ファイル出力
		//double maxval = *std::max_element( wave2.begin(), wave2.end() );
		//double minval = *std::min_element( wave2.begin(), wave2.end() );
		//double absmax = 0;

		//if(fabs( maxval ) > fabs( minval )) absmax = fabs( maxval );
		//else absmax = fabs( minval );

		//printf("SB-WAVE MAX2(gal) = %f\n ",absmax );
/////////////
		wwave = wwave2;
		wwave2.clear();
	}

//	free(wav);

}

void CBasenSpectrum::GetBasenCalc(
									std::vector< std::complex<double> > wwave,
									std::vector< std::complex<double> >& wwave2,
									std::vector<double>& wave2,
									std::vector<double> acc,
									double td,
									int nSampleFrequency,
									int nSamplePoint
									)
{
		int nSample = nSamplePoint/2 + 1;
		std::vector<double> wavv1;
//		double* wav2 = (double*)malloc(sizeof(double)*nSamplePoint*2);
    	std::vector<std::complex<double> > wwav2(nSamplePoint);

		CFastFourierTransform::Trans( wwave );

//////////
        //std::vector<std::complex<double> > inputWave( nSamplePoint );

//////////
		// 位相を取り出し再度加速度フーリエスペクトルに割り当てる。
		for(int jj=0 ; jj < nSample ; jj++ )
		{
			wavv1.push_back(arg( wwave[jj] ));//位相取得（振幅情報は捨てる）
			//wav2[2*jj]   = acc[jj]*cos(wavv1[jj]);//* (double)nSampleFrequency;
			//wav2[2*jj+1] = acc[jj]*sin(wavv1[jj]);//* (double)nSampleFrequency;
            wwav2[jj] = complex<double>(acc[jj]*cos(wavv1[jj]),acc[jj]*sin(wavv1[jj]));
		}
        //<-2006/01/17
        //Inverse 関数で共役複素数を割り当てるのでこの部分はコメントアウト
		// 再度負の周波数成分には共役複素数を割り当てる
        /*
		int ii,ll,k;
		for(ii=nSample-2,ll=nSample*2 ;ii >= 1 ; ii--, ll+=2 )
		{
			wav2[ll  ] = wav2[2*ii  ];
			wav2[ll+1]= -wav2[2*ii+1];
		}
    	for(ii=nSample-2, k=0 ; ii>=1 ; ii--, k++)
        {
            wwav2[nSample + k] = conj(wwav2[ii]);
    	}*/
		// 時間領域に戻す
		//four1(wav2-1, nSamplePoint, -1);
        CFastFourierTransform::Inverse( wwav2 );
        //<-

        /*
        for(int kk=0 ; kk < nSamplePoint ; kk++ ){
            inputWave[kk] = inputWave[kk]/ (double)nSamplePoint * (double)nSampleFrequency;
            printf("%d, %lf\n",kk,real(inputWave[kk]));
        }
        */
		// 再度実部を波形変数に格納する
		for(int kk = 0 ; kk < nSamplePoint ; kk++ )
		{
			if( kk < td * nSampleFrequency )// (td * nSampleFrequency) = tdまでのデータ数
			{
				//wave2.push_back(wav2[2*kk] / (double)nSamplePoint * (double)nSampleFrequency);  //  <=?
				wave2.push_back(real(wwav2[kk]) / (double)nSamplePoint * (double)nSampleFrequency);  //  <=?
				wwave2.push_back(wave2[kk]);
			}
			else
			{
				wave2.push_back(0.0);
				wwave2.push_back(wave2[kk]);
			}
		}
		//for(int kk = 0 ; kk < nSamplePoint ; kk++ )
		//{
		//    printf("%lf\n",wave2[kk]);
        //}
		// mallocで確保したメモリを開放する
//		free(wav2);
}

void CBasenSpectrum::GetEnvelope(double dR_pq,
								 int nSampleFrequency,
								 int nSamplePoint,
								 vector<double>& env ,
								 double& td,
								 double& envv
								 )
{

	double t;
    double dMjma_pq;
	// 経験式
	dMjma_pq = (log10(m_dM0_pq*1e7)-16.2)/1.5;

	double ta = 0;
	double tb = pow(10,(0.229*dMjma_pq-1.112));
	double tc = (tb + pow(10,(0.433*dMjma_pq-1.936)));
	       td = (tc + pow(10,(0.778*log10(dR_pq)-0.34)));
	double envsum = 0;
	envv = 0;

	tb += ta;
	tc += ta;
	td += ta;

	env.clear();
	env.resize(nSamplePoint,0.0);

	int i;
	for(i=0;i<nSamplePoint;i++){

		t = (double)i/(double)nSampleFrequency;

		if     (t <= 0  )	env[i] = 0.0;
		else if(t <= ta )	env[i] = 0.0;
		else if(t <= tb )	env[i] = pow((t-ta)/(tb-ta), 2.0);
		else if(t <= tc )   env[i] = 1.0;
		else if(t <= td )	env[i] = exp(-(log(10.0))*((t-tc)/(td-tc)));
		else				env[i] = 0.0;

		if( t > tc  && env[i] < 0.1 )  env[i] = 0.0;
	}   
		//面積比
		envsum = accumulate(env.begin(),env.end(),0.0f);
		envv = (td*nSampleFrequency)/(envsum);

}
