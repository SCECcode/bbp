// Copyright (c) 2011-2018 National Research Institute for Earth
// Science and Disaster Resilience (NIED)
//
// All rights reserved.
//
// Programs statgreen, statgreen2, and greenscale computes ground
// motion time series by the stochastic Green's function (SGF) method.
//
// Disclaimer:
/* Users take full responsibility for the use of the software. Neither
   NIED nor any of the contributors shall be responsible for any loss,
   damages, or costs which may be incurred by the users or third
   parties arising in the use of the software, whatever the reason may
   be, even if advised of the possibility of such damages. */
//
// References:
//
// (1) Reference for the SGF software
/*
- Senna, S. and H. Fujiwara (2011): Development of Estimation Tools
  for Earthquake Ground Motion, Technical Note of the National
  Research Institute for Earth Science and Disaster Prevention,
  No. 354 (in Japanese)
*/
//
// (2) Stochastic Green's function (SGF) method
/*
- Dan, K. and T. Sato (1998): Strong Motion Prediction by
  Semi-empirical Method Based on Variable-slip Rupture Model of
  Earthquake Fault, Journal of Structural and Construction Engineering
  (Transactions of AIJ), No.509, 49-60. (in Japanese)

- Dan, K., M. Watanabe, T. Sato, J. Miyakoshi, and T. Sato (2000):
  Isoseismal Map of Strong Motions for the 1923 Kanto Earthquake
  (MJMA7.9) by Stochastic Green's Function Method, Journal of
  Structural and Construction Engineering (Transactions of AIJ),
  No.530, 53-62. (in Japanese)

- Satoh, T., H. Kawase, and T. Sato (1994): Engineering bedrock waves
  obtained through the identification analysis based on borehole
  records and their statistical envelope characteristics, Journal of
  Structural and Construction Engineering (Transactions of AIJ),
  No.461, 19-28. (in Japanese)
*/
//
// (3) NIED broadband ground motion computation by recipe
/*
- Fujiwara, H. et al. (2009). Technical Reports on National Seismic
  Hazard Maps for Japan, Technical Note of the National Res. Inst. for
  Earth Science and Disaster Prevention 336, 528pp (in English)

- Fujiwara, H. et al. (2012): Some Improvements of Seismic Hazard
  Assessment based on the 2011 Tohoku Earthquake, Technical Note of
  the National Research Institute for Earth Science and Disaster
  Prevention, No.379. (in Japanese)

- They can be obtained from:

  http://www.j-shis.bosai.go.jp/map/JSHIS2/download.html?lang=en
*/

#include "greenscale.h"
#include "FastFourierTransform.h"
#include "knetascii.h"
#include "knetout.h"
#include "filter.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <complex>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

using namespace std;

int main(int argc, char* argv[])
{

	FLT mnfault,subf1;

	double omega,fhz,Rs,Rpq,minRpq;
	double Qvalue,tmp;
	double EEnumber,ELatitude,ELongitude, EDepth, EDpq, EVpq, ESigpq ,ELagTime;
	double OBBserno ,OBLatitude, OBLongitude, OBDepth;
	int Samplefrq;
    double pi = atan( 1.0 ) * 4.0;
    double intensity = 0;
	double xxi,ddi,XXDi,XXEQ;
	std::complex<double> sekibunc;
    int nShiftFlg; //シフト時間を決めるフラグ

	minRpq = 30000000;//最大観測距離(m)

	//--------------------- 読み込みファイル数エラー処理------------------------
    if( argc == 11 || argc == 9 ){
		;
	}else{
		std::cerr << "Usage: greenscale dir_name station elem_param fault_param "<<
                     "log filt_flg(HP LP NO) freq shiftflg start_num end_num" << std::endl;
		std::cerr << "       or " << std::endl;
		std::cerr << "Usage: greenscale dir_name station elem_param fault_param "<<
                     "log filt_flg(HP LP NO) freq shiftflg" << std::endl;
		return 1;
    }
	//------------------------- 要素地震波形の読み込み--------------------------

	std::vector<double> data;

	//char* filename = argv[1];

	//------------------ 要素断層データをファイルから読み込む。 ----------------

	char* filename1 = argv[3];
	//要素値(elementファイル)を読み込む

	//fieldにvectorを登録
	std::vector<double> Enumber;
	std::vector<double> latitude;
	std::vector<double> longitude;
	std::vector<double> depth;
	std::vector<double> ElementDpq;
	std::vector<double> ElementSigpq;
	std::vector<double> dElementLagTime;

	std::vector< std::vector<double>* > field;
    field.push_back( &Enumber );
	field.push_back( &latitude );
	field.push_back( &longitude  );
	field.push_back( &depth );
	field.push_back( &ElementDpq );
	field.push_back( &ElementSigpq );
	field.push_back( &dElementLagTime );
	
	unsigned int numofElement = read_csv( filename1, field );
	if( numofElement == 0 ){
		std::cerr << filename1 << " not found." << std::endl;
		return 1;
	}

	//------------------ 観測点データをファイルから読み込む。 ----------------

	char* filename2 = argv[2];
    string folder = argv[1];

	// 観測点(observファイル)を読み込む
	// observにvectorを登録

	std::vector<string> OBserno;
	std::vector<string> OLatitude;
	std::vector<string> OLongitude;
	std::vector<string> ODepth;
	std::vector<string> OFileName;

	//std::vector< std::vector<double>* > observ;
	std::vector< std::vector<string>* > observ;
    observ.push_back( &OBserno );
	observ.push_back( &OLatitude );
	observ.push_back( &OLongitude  );
	observ.push_back( &ODepth );
	observ.push_back( &OFileName );

	unsigned int numofObpoint = read_string( filename2, observ );
	if( numofObpoint == 0 ){
		std::cerr << filename2 << " not found." << std::endl;
		return 1;
	}

	//------------ 断層パラメータ，破壊開始等パラメータ読み込み ----------------

    char* efaultparam = argv[4];

	std::vector<string> RLatitude;
	std::vector<string> RLongitude;
	std::vector<string> RDepth;
	std::vector<string> DDLatitude;
	std::vector<string> DDLongitude;
	std::vector<string> DDDepth;
	std::vector<string> ESamplefrq;
	std::vector<string> ESampletime;
	std::vector<string> Em_dRo_pq;
	std::vector<string> Em_dRo_sb;
	std::vector<string> EBetas;
	std::vector<string> Em_dBeta_sb;
	std::vector<string> EM0s;
	std::vector<string> Esigmas;
	std::vector<string> Em_dFrad;
	std::vector<string> Efmax;
	std::vector<string> Em_dM;
	std::vector<string> Qdepend1;
	std::vector<string> Qdepend2;
	std::vector<string> Qdepend3;
	std::vector<string> Qdepend4;
    //add green scale param
    std::vector<string> ELs;
    std::vector<string> EWs;
    std::vector<string> EDs;
    std::vector<string> Emyu;
    std::vector<string> Evr;
    std::vector<string> Eoffset;

	//std::vector< std::vector<double>* > efaultparamset;
	std::vector< std::vector<string>* > efaultparamset;

	efaultparamset.push_back( &RLatitude );
	efaultparamset.push_back( &RLongitude );
	efaultparamset.push_back( &RDepth );
	efaultparamset.push_back( &DDLatitude );
	efaultparamset.push_back( &DDLongitude );
	efaultparamset.push_back( &DDDepth );
	efaultparamset.push_back( &ESamplefrq );
	efaultparamset.push_back( &ESampletime );
	efaultparamset.push_back( &Em_dRo_pq );
	efaultparamset.push_back( &Em_dRo_sb );
	efaultparamset.push_back( &EBetas );
	efaultparamset.push_back( &Em_dBeta_sb );
	efaultparamset.push_back( &EM0s );
	efaultparamset.push_back( &Esigmas );
	efaultparamset.push_back( &Em_dFrad );
	efaultparamset.push_back( &Efmax );
	efaultparamset.push_back( &Em_dM );
	efaultparamset.push_back( &Qdepend1 );
	efaultparamset.push_back( &Qdepend2 );
	efaultparamset.push_back( &Qdepend3 );
	efaultparamset.push_back( &Qdepend4 );
	efaultparamset.push_back( &ELs );
	efaultparamset.push_back( &EWs );
	efaultparamset.push_back( &EDs );
	efaultparamset.push_back( &Emyu );
	efaultparamset.push_back( &Evr );
	efaultparamset.push_back( &Eoffset );

	unsigned int numofObefault = read_fault( efaultparam, efaultparamset );
	if( numofObefault == 0 ){
		std::cerr << efaultparam << " not found." << std::endl;
		return 1;
	}

	//------------------------------ 出力ログファイル名入力 ----------------------------

	char* Outputlogfile = argv[5];
    //ハイパスかローパスかのフラグ
    string sHL_Flg = argv[6];
    if(sHL_Flg == "HP" || sHL_Flg == "LP" || sHL_Flg == "NO"){
        ;
    }else{
        std::cerr << "Filter Flg error" << std::endl;
        return 1;
    }
    //接続周波数の指定
    string sC_Freq = argv[7];
    //ループのはじめ
    unsigned int lSta,lEnd;
    if(argc == 9){
        lSta = 0;
        lEnd = numofObpoint;
    }else{
        lSta = atoi(argv[9]);
        if(lSta > 0 ){
            //指定した番号とファイルの番号の整合を取る。直感的にわかりやすくしておく。
            lSta = lSta -1;
        }
        lEnd = atoi(argv[10])  ;
        if(lEnd > numofObpoint){
            //最大値を超える場合は，最大値に戻す。
            lEnd = numofObpoint;
        }
    }

    //震源からのトラベルタイムを引くかどうかのフラグ
    //nShiftFlg==1 のとき従来どおり引く。0のとき引かない
    nShiftFlg = atoi(argv[8]);
    if(nShiftFlg == 0 || nShiftFlg == 1){
        ;//OK
    }else{
		std::cerr << "shift flg error nShiftFlg = (0 or 1) : input="<< nShiftFlg << std::endl;
        return 1;
    }
	//-----------------------------------波形合成計算--------------------------------------

	double sss2[65535];

	//ifstream *wavefile = new ifstream();
	//wavefile->open(filename,ios::in);

	unsigned int l;
    bool bFlg=false;
	//while(!wavefile->eof())
	//{
        //bFlg=false;
        //引数の指定により計算する観測点のループを振り分ける。
        //
    for( l = lSta ; l < lEnd ; l++ ) {// 観測点ループ
        //ｌは観測点のカウンタ
        bFlg=false;
        if(numofObpoint != 1){
    		if( l == numofObpoint)
                break;
        }else{
            if( l == 1 )
                break;
        }
		//char line1[65535];
		//wavefile->getline(line1,65535);
		//printf(" ELEMENT-WAVE-FILENAME = %s\n",line1);
        char line1[65535];
        //wavefile->getline(line1,65535);
        //OFileName
        //if( strcmp(line1 , "\0") == 0){
        if( strcmp(OFileName[l].c_str() , "\0") == 0){
            bFlg=true;
            break;
        }
        //printf(" ELEMENT-WAVE-FILENAME = %s\n",line1);
        printf(" ELEMENT-WAVE-FILENAME = %s\n",("./"+folder+"/"+OFileName[l]).c_str());

		knetascii knet1;
		ifstream ifs1;
		//ifs1.open(line1);
		ifs1.open(("./"+folder+"/"+OFileName[l]).c_str());
		ifs1 >> knet1;
		ifs1.close();
		knet1.get_data(data);
		Samplefrq = knet1.get_sampling_frequency();
		
		//フーリエ用データリサイズ
		int oridatasize = (int)data.size();
		int datasize;
		for(int nn = 0 ;  nn < 16 ; nn++ )
		{
			if( oridatasize > (int)std::pow(2.0,nn) )
                datasize = (int)std::pow(2.0,nn+1);
		}
		// 計算用
		std::vector< std::complex<double> > wdat(datasize);
		std::vector< std::complex<double> > wtmp(datasize);
		std::vector< std::complex<double> > wtmp2(datasize);
		std::vector< std::complex<double> > wtmps(datasize);
		std::vector< std::complex<double> > wtmpp(datasize);
		std::vector< std::complex<double> > wdats2(datasize);
		std::vector< std::complex<double> > wdatp2(datasize);

		// 計算結果用
		std::vector<double> sss( datasize,0 ); // S波合成波形
		std::vector<double> sss1( datasize,0 ); // S波要素波形
		std::vector<double> ppp( datasize,0 );  // P波合成波形
		std::vector<double> ppp1( datasize,0 ); // P波要素波形
		std::vector<double> sekibuns( datasize,0 ); // 積分速度波形S波
		std::vector<double> sekibunp( datasize,0 ); // 積分速度波形P波
		std::complex<double> Sscalefactor; //S波増幅率
		std::complex<double> Pscalefactor; //P波増幅率
		std::vector<double> Xi( numofObpoint,0 );
		std::vector<double> Di( numofObpoint,0 );
		std::vector<double> XDi( numofObpoint,0 );
		std::vector<double> XEQ( numofObpoint,0 );

		//データサイズの変更
		data.resize(datasize);
		wdat.resize(datasize);
		wtmp.resize(datasize);
		wtmp2.resize(datasize);
		wtmps.resize(datasize);
		wtmpp.resize(datasize);
		wdats2.resize(datasize);
		wdatp2.resize(datasize);
        sss2[datasize];//配列のCHECK?

		subf1.LengthElement = atof(ELs[0].c_str())*1000.0;
		subf1.WidthElement  = atof(EWs[0].c_str())*1000.0;
		subf1.Ds            = atof(EDs[0].c_str());
		subf1.myu           = atof(Emyu[0].c_str());
		subf1.M0s           = atof(EM0s[0].c_str());
		//subf1.Beta_s        = atof(EBetas[0].c_str());
		subf1.Beta_s        = atof(EBetas[0].c_str()) * 1000.0; //(km/s) -> (m/s)
		subf1.Sigma_s       = atof(Esigmas[0].c_str());
		subf1.vr            = atof(Evr[0].c_str()) * 1000.0;    //(km/s) -> (m/s)
		subf1.fmax          = atof(Efmax[0].c_str());
		subf1.Qdepend1      = atof(Qdepend1[0].c_str());
		subf1.Qdepend2      = atof(Qdepend2[0].c_str());
		subf1.Qdepend3      = atof(Qdepend3[0].c_str());
		subf1.Qdepend4      = atof(Qdepend4[0].c_str());
		subf1.offsettime    = atof(Eoffset[0].c_str());

		init( &mnfault, &subf1 );
		MakeFaultPara( &subf1 );

		double dt = 1.0/(double)Samplefrq;

		// データ内容表示
		printf(" STATION-POINTS = %d , ELEMENT-POINTS = %d\n", numofObpoint , numofElement );
		printf(" ORIGINAL-DATASIZE = %d , SET-DATASIZE = %d , SAMPLING-FREQUENCY(Hz) = %d\n\n", oridatasize , datasize , Samplefrq );

		for( int i = 0 ; i < datasize ; i++ ) wtmp[i] = data[i];

		CFastFourierTransform::Trans( wtmp,1 );

		OBBserno     =  atof(OBserno[l].c_str()) ;
		OBLatitude   =  atof(OLatitude[l].c_str());
		OBLongitude  =  atof(OLongitude[l].c_str());
		OBDepth      =  atof(ODepth[l].c_str());

		Rs = 0;
		minRpq = 30000000;

		double RsLatitude   =  atof(RLatitude[0].c_str());       // 断層中央点緯度
		double RsLongitude  =  atof(RLongitude[0].c_str());      // 断層中央点経度
		double RsDepth		=  atof(RDepth[0].c_str());		   // 断層中央深度
		double DLatitude	=  atof(DDLatitude[0].c_str());      // 破壊開始点緯度
		double DLongitude	=  atof(DDLongitude[0].c_str());     // 破壊開始点経度
		double DDepth		=  atof(DDDepth[0].c_str());         // 破壊開始点深度

		Distance( OBLongitude, OBLatitude, OBDepth ,RsLongitude, RsLatitude, RsDepth , Rs );

	//------------------------- 重ね合わせ計算部(壇・佐藤(1998)）----------------------
        //TDistanceは、要素断層にかかわらず一定であるので、ループのそとに
        //出した。破壊開始点と観測点によって決まる数値である。
        double TDistance = 0; // 破壊開始点から観測点までの距離(m)
		// 破壊開始点-観測点距離計算
		Distance( DLongitude, DLatitude, DDepth ,OBLongitude, OBLatitude, OBDepth , TDistance );
		double  tvrt     =  TDistance  / subf1.Beta_s; // 破壊開始点から観測点までの時間(sec)

		for( unsigned int k = 0 ; k < numofElement ; k++ ){
            //ｋは要素番号のカウンタ
			EEnumber   = Enumber[k];
			ELatitude  = latitude[k];
			ELongitude = longitude[k];
			EDepth     = depth[k];
			EDpq       = ElementDpq[k];
			ESigpq     = ElementSigpq[k];
            ELagTime   = dElementLagTime[k];

			wtmp2 = wtmp;

			Rpq = 0;

			Distance( OBLongitude, OBLatitude, OBDepth ,ELongitude, ELatitude, EDepth , Rpq );

			EVpq = (2.0*ESigpq*1000000*subf1.Beta_s)/subf1.myu;
			double EomegaDpq = ( EVpq / EDpq );
			double Sigmapq   = ESigpq*1000000;                 // MPa -> Pa
			double dNyquist = (double)Samplefrq/2.0;
			double df = dNyquist/((double)datasize/2.0);

			Xi[l] += Rpq/1000.0;
			//Di[l] += EDpq;                                      //org
			Di[l] += std::pow( EDpq,2);
			XDi[l] += std::pow( EDpq/(Rpq/1000.0) , 2 );        //org
			//XDi[l] += (EDpq / std::pow( (Rpq/1000.0) , 2 ));
			XEQ[l] = sqrt(Di[l]/XDi[l]);

			wtmp2[0] = 0.0;
	//----------------------------------- 周波数ループ内 ----------------------------------------
			for( int iii = 1 ; iii < (int)(datasize/2) ; iii++ ){

				fhz = df * (double)iii;
				omega = 2 * pi * fhz;                         //ωの計算

				if( subf1.Qdepend3 <= fhz )
                    Qvalue = (subf1.Qdepend1 * pow(fhz, subf1.Qdepend2));
				else
                    Qvalue = subf1.Qdepend4;

				Sscalefactor = ScaleFactorS( omega, Qvalue, EDpq, EomegaDpq, Rpq, Rs, Sigmapq, mnfault, subf1 ); // S波スケールファクターFpqを求める
				Pscalefactor = ScaleFactorP( omega, Qvalue, EDpq, EomegaDpq, Rpq, Rs, Sigmapq, mnfault, subf1 ); // P波スケールファクターFpqを求める
				wtmps[iii] = wtmp2[iii] * Sscalefactor;
			}
			wdats2 = wtmps;
            ///////
			CFastFourierTransform::Inverse( wdats2 ); //波形を周波数領域から時間領域に変換する。（S波）
            ////////
			tmp = 1.0 / (double) datasize;

			for( int m = 0 ; m < datasize ; m++ ){
				wdats2[m] = wdats2[m] * tmp;
				sss1[m] = wdats2[m].real();
			}

			double DDistance = 0; // 破壊開始点からの距離(m)
			//double TDistance = 0; // 破壊開始点から観測点までの距離(m)   -> ループの外へ
			// 破壊開始点-要素断層距離計算
			//Distance( DLongitude, DLatitude, DDepth ,ELongitude, ELatitude, EDepth , DDistance );
			// 破壊開始点-観測点距離計算  -> ループの外へ
			//Distance( DLongitude, DLatitude, DDepth ,OBLongitude, OBLatitude, OBDepth , TDistance );
			double  vrt;//      =  DDistance  / subf1.vr; // 破壊開始点からの各要素の発震時間(sec)
            vrt = ELagTime;
			//double  tvrt     =  TDistance  / subf1.Beta_s; // 破壊開始点から観測点までの時間(sec)　-> ループの外へ
			double  travelts =  Rpq / subf1.Beta_s;   // 要素地震毎のS波の観測点への到達時間(sec)
			double  traveltp =  Rpq / subf1.Alpha_p;  // 要素地震毎のP波の観測点への到達時間(sec)

			if( Rpq < minRpq ) minRpq = Rpq;// 断層最短距離(m)

			int vrshifts ;
            if( nShiftFlg == 1){
                vrshifts = (int)( ( subf1.offsettime + travelts + vrt - tvrt ) / (1.0/(double)Samplefrq) );
            }else{
                vrshifts = (int)( ( subf1.offsettime + travelts + vrt ) / (1.0/(double)Samplefrq) );
            }
			int midsift = datasize - vrshifts;//S波伝播時間シフト時間長

			GetGousei( &sss , &sss1 , midsift , vrshifts );// S波合成(NS,EW)

			//for( unsigned int pp = 0 ; pp < 	wdats2.size() ; pp++ ){
			//	sss2[pp] = sss[pp];
			//}   -> ループの外へ
		}
        for( unsigned int pp = 0 ; pp <	wdats2.size() ; pp++ ){
			sss2[pp] = sss[pp];
		}

		// 等価震源距離計算用関数
        //xxi = Xi[l]/numofElement;
		//xxi = Xi[l];
		//ddi = Di[l];
		//XXDi = XDi[l];
		XXEQ = XEQ[l];

		// 積分速度波形計算
		Ssekibun( &sss, &sekibuns, datasize, dt, 10, 0.6321);

		// 計測震度計算
		seismic_intensity( datasize, dt, sss2, intensity );

		// S波加速度最大値
		double amax = *max_element(sss.begin(),sss.end());
		double amin = *min_element(sss.begin(),sss.end());
		double absamax = 0;
		if( fabs(amax) > fabs(amin) ) absamax = fabs(amax);
		else absamax = fabs(amin);

		// S波速度最大値
		double vmax = *max_element(sekibuns.begin(),sekibuns.end());
		double vmin = *min_element(sekibuns.begin(),sekibuns.end());
		double absvmax = 0;
		if( fabs(vmax) > fabs(vmin) ) absvmax = fabs(vmax);
		else absvmax = fabs(vmin);

		// 加速度波形出力
		std::ostringstream filename_strm;
		filename_strm << "ac" << std::setw(6) << std::setfill('0') << OBBserno << ".dat";
		std::ofstream strm( filename_strm.str().c_str() );
		if( strm.is_open() ){
			output_in_kneta( strm, sss, dt, oridatasize, Samplefrq );
			strm.close();
		}else{
			std::cerr << "file out error." << std::endl;
		}

		// 速度波形出力
		std::ostringstream filename_strm1;
		filename_strm1 << "ve" << std::setw(6) << std::setfill('0') << OBBserno << ".dat";
		std::ofstream strm1( filename_strm1.str().c_str() );
		if( strm1.is_open() ){
			output_in_knetv( strm1, sekibuns, dt, oridatasize, Samplefrq );
			strm1.close();
		}else{
			std::cerr << "file out error." << std::endl;
		}
        //速度波形ファイル名をフィルタールーチンに渡す。
        double dFilteredMax=0.0;
        if(sHL_Flg != "NO"){
            filter_main(filename_strm1.str(),sHL_Flg,sC_Freq,&dFilteredMax);
        }
		// 解析結果ログファイル出力
        //フィルター後の結果を出力するためここに移動
		FILE *fpo1;
		if( (fpo1 = fopen( Outputlogfile , "a" )) == NULL )exit(1);
		fprintf( fpo1 , "%8.0lf\t%8.4lf\t%8.4lf\t%8.4lf\t%9.5lf\t%8.4lf\t%8.4lf\t%8.4lf\t%s\t%9.5lf\n"
			, OBBserno , OBLatitude , OBLongitude , absamax , absvmax ,
              minRpq/1000.0 , XXEQ , intensity , sHL_Flg.c_str(), dFilteredMax);
		fclose( fpo1 );
		// 解析結果画面出力
		printf(" CALC No. = %d  MAXACC(gal) = %f \n MAXVEL(kine) = %f FilteredMAXVEL(kine) = %f\n",
                             l+1 , absamax , absvmax, dFilteredMax );
		printf(" MINR(km) = %f Xeq(km) = %f INTENSITY = %f\n\n", minRpq/1000.0 , XXEQ , intensity );

		// データクリア
		std::fill( data.begin(), data.end(), 0);
		std::fill( wtmp.begin(), wtmp.end(), 0);
		std::fill( wtmp2.begin(), wtmp2.end(), 0);
		std::fill( wtmps.begin(), wtmps.end(), 0);
		std::fill( wtmpp.begin(), wtmpp.end(), 0);
		std::fill( wdats2.begin(), wdats2.end(), 0);
		std::fill( wdatp2.begin(), wdatp2.end(), 0);
		std::fill( wdat.begin(), wdat.end(), 0);
		std::fill( sss.begin(), sss.end(), 0 );
		std::fill( sss1.begin(), sss1.end(), 0 );
		std::fill( sekibuns.begin(), sekibuns.end(), 0 );
		std::fill( Xi.begin(), Xi.end(), 0 );
		std::fill( Di.begin(), Di.end(), 0 );
		std::fill( XDi.begin(), XDi.end(), 0 );
		std::fill( XEQ.begin(), XEQ.end(), 0 );
        //if( l == numofObpoint -1) break;
	}
    //if(bFlg) break;
	//	if(wavefile->eof()) break;
	//}
	//	wavefile->close();
	//	delete wavefile;
	return 0;
}

