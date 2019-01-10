#include "filter.h"
#include "knetout.h"
#include "knetascii.h"
#include "saitofilt.h"

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
#include <math.h>


using namespace std;

void filter_main(string sInFileName, string sHL_Flg, string sC_Freq,double* dMax)
{
    //sHL_Flg = "HP" or "LP"
    //HP: high pass
    //LP: low pass
    std::vector<double> data1;
    double samplingf;

	//WAVEDATA
    //char* filename1 = sInFileName.c_str();

    knetascii knet1;
    knet1.kind_of_data_ =  KNET_VEL;
    ifstream ifs1;
    ifs1.open(sInFileName.c_str());
    ifs1 >> knet1;
    ifs1.close();
    knet1.get_data(data1);
	samplingf = knet1.get_sampling_frequency();
	int wavenum = data1.size();
    std::vector<double> data2( data1.size() );
	double dt = 1.0/samplingf;
	//引数から接続周波数をセットする
	double fmatch = atof(sC_Freq.c_str());
    double flow  = fmatch / 1.50;
	double fhigh = fmatch * 1.50;
	double ap = 0.10;
	double as = 10.0;
	double gn;
    double h[200];
    for(int i=0 ; i<200 ; i++){
        h[i]=0.0;
    }
    int m,n;
    int nml1 = -1;
    int nml2 = 1;
	double *xx,*yy,*filop,*fihip;
    xx = new double[wavenum];
    yy = new double[wavenum];
    filop = new double[wavenum];
    fihip = new double[wavenum];

    string stradd;
    for( int i = 0 ; i < wavenum ; i++ ){
        xx[i] = data1[i];
        //cout << data1[i] << "\n";

    }
    //フィルター種別の振り分けと文字列の作成
    if(sHL_Flg == "HP"){
	    buthip(h,&m,&gn,&n,fhigh*dt,flow*dt,ap,as);
        stradd = "_HP";
    }else if(sHL_Flg == "LP" ){
    	butlop(h,&m,&gn,&n,flow*dt,fhigh*dt,ap,as);
        stradd = "_LP";
    }else{
        cout << "Please Set HL_Flg \"HP\":high pass \"LP\":low pass\n";
        return ;
    }
    //
	tandem(xx,yy,wavenum,h,m,nml1);

    if(sHL_Flg == "HP" ){
	    tandem(yy,fihip,wavenum,h,m,nml2);
    }else if(sHL_Flg == "LP" ){
    	tandem(yy,filop,wavenum,h,m,nml2);
    }
	for( int j = 0 ; j < wavenum ; j++ ){
        if(sHL_Flg == "HP" ){
            data2[j] = fihip[j]*gn*gn;
        }else if(sHL_Flg == "LP"){
            data2[j] = filop[j]*gn*gn;
        }
    }
	//FILTERWAVE
    string strall,strFreq,strFileName,strFileExt;

    strall = sInFileName;
    strFileName = fnGetFileName(strall);
    strFileExt  = fnGetExtName(strall);
    //接続周波数指定文字列に'.'があれば，'_'に置き換える。
    string sFreq = sC_Freq;
    fnReplace(sFreq, string("."), string("_") );
    //出力用ファイル名を作成する
    strall = strFileName + stradd + sFreq + strFileExt;

    std::ofstream strm1( strall.c_str() );

    if( strm1.is_open() ){
        output_in_knet( strm1, data2, dt,knet1 );
        strm1.close();
    }else{
        std::cerr << "file out error." << std::endl;
    }
	double maxval = *std::max_element( data2.begin(), data2.end() );
	double minval = *std::min_element( data2.begin(), data2.end() );
	double absmax = 0;
	if(fabs( maxval ) > fabs( minval )) absmax = fabs( maxval );
	else absmax = fabs( minval );
    *dMax = absmax;
    delete xx;
    delete yy;
    delete fihip;
    delete filop;
	return ;

}
