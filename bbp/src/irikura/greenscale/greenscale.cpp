// Copyright (c) 2011-2018 National Research Institute for Earth Science and Disaster Resilience (NIED)
// All rights reserved.
//
// Programs statgreen, statgreen2, and greenscale computes ground motion time series by the stochastic Green's function (SGF) method.
//
// Disclaimer:
/* Users take full responsibility for the use of the software. Neither NIED nor any of the contributors shall be responsible for any loss, damages, or costs which may be incurred by the users or third parties arising in the use of the software, whatever the reason may be, even if advised of the possibility of such damages. */
//
// References:
// (1) Reference for the SGF software
/*
- Senna, S. and H. Fujiwara (2011): Development of Estimation Tools for Earthquake Ground Motion, Technical Note of the National Research Institute for Earth Science and Disaster Prevention, No. 354 (in Japanese)
*/
// (2) Stochastic Green's function (SGF) method
/*
- Dan, K. and T. Sato (1998): Strong Motion Prediction by Semi-empirical Method Based on Variable-slip Rupture Model of Earthquake Fault, Journal of Structural and Construction Engineering (Transactions of AIJ), No.509, 49-60. (in Japanese)
- Dan, K., M. Watanabe, T. Sato, J. Miyakoshi, and T. Sato (2000): Isoseismal Map of Strong Motions for the 1923 Kanto Earthquake (MJMA7.9) by Stochastic Green's Function Method, Journal of Structural and Construction Engineering (Transactions of AIJ), No.530, 53-62. (in Japanese)
- Satoh, T., H. Kawase, and T. Sato (1994): Engineering bedrock waves obtained through the identification analysis based on borehole records and their statistical envelope characteristics, Journal of Structural and Construction Engineering (Transactions of AIJ), No.461, 19-28. (in Japanese)
*/
// (3) NIED broadband ground motion computation by recipe
/*
- Fujiwara, H. et al. (2009). Technical Reports on National Seismic Hazard Maps for Japan, Technical Note of the National Res. Inst. for Earth Science and Disaster Prevention 336, 528pp (in English)
- Fujiwara, H. et al. (2012): Some Improvements of Seismic Hazard Assessment based on the 2011 Tohoku Earthquake, Technical Note of the National Research Institute for Earth Science and Disaster Prevention, No.379. (in Japanese)
- They can be obtained from http://www.j-shis.bosai.go.jp/map/JSHIS2/download.html?lang=en .
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
    int nShiftFlg; //�V�t�g���Ԃ����߂��t���O

	minRpq = 30000000;//�ő��ϑ�����(m)

	//--------------------- �ǂݍ��݃t�@�C�����G���[����------------------------
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
	//------------------------- �v�f�n�k�g�`�̓ǂݍ���--------------------------

	std::vector<double> data;

	//char* filename = argv[1];

	//------------------ �v�f�f�w�f�[�^���t�@�C�������ǂݍ��ށB ----------------

	char* filename1 = argv[3];
	//�v�f�l(element�t�@�C��)���ǂݍ���

	//field��vector���o�^
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

	//------------------ �ϑ��_�f�[�^���t�@�C�������ǂݍ��ށB ----------------

	char* filename2 = argv[2];
    string folder = argv[1];

	// �ϑ��_(observ�t�@�C��)���ǂݍ���
	// observ��vector���o�^

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

	//------------ �f�w�p�����[�^�C�j���J�n���p�����[�^�ǂݍ��� ----------------

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

	//------------------------------ �o�̓��O�t�@�C�������� ----------------------------

	char* Outputlogfile = argv[5];
    //�n�C�p�X�����[�p�X���̃t���O
    string sHL_Flg = argv[6];
    if(sHL_Flg == "HP" || sHL_Flg == "LP" || sHL_Flg == "NO"){
        ;
    }else{
        std::cerr << "Filter Flg error" << std::endl;
        return 1;
    }
    //�ڑ����g���̎w��
    string sC_Freq = argv[7];
    //���[�v�̂͂���
    unsigned int lSta,lEnd;
    if(argc == 9){
        lSta = 0;
        lEnd = numofObpoint;
    }else{
        lSta = atoi(argv[9]);
        if(lSta > 0 ){
            //�w�肵���ԍ��ƃt�@�C���̔ԍ��̐����������B�����I�ɂ킩���₷�����Ă����B
            lSta = lSta -1;
        }
        lEnd = atoi(argv[10])  ;
        if(lEnd > numofObpoint){
            //�ő��l�𒴂����ꍇ�́C�ő��l�ɖ߂��B
            lEnd = numofObpoint;
        }
    }

    //�k�������̃g���x���^�C�����������ǂ����̃t���O
    //nShiftFlg==1 �̂Ƃ��]���ǂ��������B0�̂Ƃ������Ȃ�
    nShiftFlg = atoi(argv[8]);
    if(nShiftFlg == 0 || nShiftFlg == 1){
        ;//OK
    }else{
		std::cerr << "shift flg error nShiftFlg = (0 or 1) : input="<< nShiftFlg << std::endl;
        return 1;
    }
	//-----------------------------------�g�`�����v�Z--------------------------------------

	double sss2[65535];

	//ifstream *wavefile = new ifstream();
	//wavefile->open(filename,ios::in);

	unsigned int l;
    bool bFlg=false;
	//while(!wavefile->eof())
	//{
        //bFlg=false;
        //�����̎w���ɂ����v�Z�����ϑ��_�̃��[�v���U�蕪�����B
        //
    for( l = lSta ; l < lEnd ; l++ ) {// �ϑ��_���[�v
        //���͊ϑ��_�̃J�E���^
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

		//�t�[���G�p�f�[�^���T�C�Y
		int oridatasize = (int)data.size();
		int datasize;
		for(int nn = 0 ;  nn < 16 ; nn++ )
		{
			if( oridatasize > (int)std::pow(2.0,nn) )
                datasize = (int)std::pow(2.0,nn+1);
		}
		// �v�Z�p
		std::vector< std::complex<double> > wdat(datasize);
		std::vector< std::complex<double> > wtmp(datasize);
		std::vector< std::complex<double> > wtmp2(datasize);
		std::vector< std::complex<double> > wtmps(datasize);
		std::vector< std::complex<double> > wtmpp(datasize);
		std::vector< std::complex<double> > wdats2(datasize);
		std::vector< std::complex<double> > wdatp2(datasize);

		// �v�Z���ʗp
		std::vector<double> sss( datasize,0 ); // S�g�����g�`
		std::vector<double> sss1( datasize,0 ); // S�g�v�f�g�`
		std::vector<double> ppp( datasize,0 );  // P�g�����g�`
		std::vector<double> ppp1( datasize,0 ); // P�g�v�f�g�`
		std::vector<double> sekibuns( datasize,0 ); // �ϕ����x�g�`S�g
		std::vector<double> sekibunp( datasize,0 ); // �ϕ����x�g�`P�g
		std::complex<double> Sscalefactor; //S�g������
		std::complex<double> Pscalefactor; //P�g������
		std::vector<double> Xi( numofObpoint,0 );
		std::vector<double> Di( numofObpoint,0 );
		std::vector<double> XDi( numofObpoint,0 );
		std::vector<double> XEQ( numofObpoint,0 );

		//�f�[�^�T�C�Y�̕ύX
		data.resize(datasize);
		wdat.resize(datasize);
		wtmp.resize(datasize);
		wtmp2.resize(datasize);
		wtmps.resize(datasize);
		wtmpp.resize(datasize);
		wdats2.resize(datasize);
		wdatp2.resize(datasize);
        sss2[datasize];//�z����CHECK?

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

		// �f�[�^���e�\��
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

		double RsLatitude   =  atof(RLatitude[0].c_str());       // �f�w�����_�ܓx
		double RsLongitude  =  atof(RLongitude[0].c_str());      // �f�w�����_�o�x
		double RsDepth		=  atof(RDepth[0].c_str());		   // �f�w�����[�x
		double DLatitude	=  atof(DDLatitude[0].c_str());      // �j���J�n�_�ܓx
		double DLongitude	=  atof(DDLongitude[0].c_str());     // �j���J�n�_�o�x
		double DDepth		=  atof(DDDepth[0].c_str());         // �j���J�n�_�[�x

		Distance( OBLongitude, OBLatitude, OBDepth ,RsLongitude, RsLatitude, RsDepth , Rs );

	//------------------------- �d�ˍ��킹�v�Z��(�d�E����(1998)�j----------------------
        //TDistance�́A�v�f�f�w�ɂ������炸�����ł����̂ŁA���[�v�̂��Ƃ�
        //�o�����B�j���J�n�_�Ɗϑ��_�ɂ����Č��܂鐔�l�ł����B
        double TDistance = 0; // �j���J�n�_�����ϑ��_�܂ł̋���(m)
		// �j���J�n�_-�ϑ��_�����v�Z
		Distance( DLongitude, DLatitude, DDepth ,OBLongitude, OBLatitude, OBDepth , TDistance );
		double  tvrt     =  TDistance  / subf1.Beta_s; // �j���J�n�_�����ϑ��_�܂ł̎���(sec)

		for( unsigned int k = 0 ; k < numofElement ; k++ ){
            //���͗v�f�ԍ��̃J�E���^
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
	//----------------------------------- ���g�����[�v�� ----------------------------------------
			for( int iii = 1 ; iii < (int)(datasize/2) ; iii++ ){

				fhz = df * (double)iii;
				omega = 2 * pi * fhz;                         //�ւ̌v�Z

				if( subf1.Qdepend3 <= fhz )
                    Qvalue = (subf1.Qdepend1 * pow(fhz, subf1.Qdepend2));
				else
                    Qvalue = subf1.Qdepend4;

				Sscalefactor = ScaleFactorS( omega, Qvalue, EDpq, EomegaDpq, Rpq, Rs, Sigmapq, mnfault, subf1 ); // S�g�X�P�[���t�@�N�^�[Fpq�����߂�
				Pscalefactor = ScaleFactorP( omega, Qvalue, EDpq, EomegaDpq, Rpq, Rs, Sigmapq, mnfault, subf1 ); // P�g�X�P�[���t�@�N�^�[Fpq�����߂�
				wtmps[iii] = wtmp2[iii] * Sscalefactor;
			}
			wdats2 = wtmps;
            ///////
			CFastFourierTransform::Inverse( wdats2 ); //�g�`�����g���̈悩�玞�ԗ̈��ɕϊ������B�iS�g�j
            ////////
			tmp = 1.0 / (double) datasize;

			for( int m = 0 ; m < datasize ; m++ ){
				wdats2[m] = wdats2[m] * tmp;
				sss1[m] = wdats2[m].real();
			}

			double DDistance = 0; // �j���J�n�_�����̋���(m)
			//double TDistance = 0; // �j���J�n�_�����ϑ��_�܂ł̋���(m)   -> ���[�v�̊O��
			// �j���J�n�_-�v�f�f�w�����v�Z
			//Distance( DLongitude, DLatitude, DDepth ,ELongitude, ELatitude, EDepth , DDistance );
			// �j���J�n�_-�ϑ��_�����v�Z  -> ���[�v�̊O��
			//Distance( DLongitude, DLatitude, DDepth ,OBLongitude, OBLatitude, OBDepth , TDistance );
			double  vrt;//      =  DDistance  / subf1.vr; // �j���J�n�_�����̊e�v�f�̔��k����(sec)
            vrt = ELagTime;
			//double  tvrt     =  TDistance  / subf1.Beta_s; // �j���J�n�_�����ϑ��_�܂ł̎���(sec)�@-> ���[�v�̊O��
			double  travelts =  Rpq / subf1.Beta_s;   // �v�f�n�k����S�g�̊ϑ��_�ւ̓��B����(sec)
			double  traveltp =  Rpq / subf1.Alpha_p;  // �v�f�n�k����P�g�̊ϑ��_�ւ̓��B����(sec)

			if( Rpq < minRpq ) minRpq = Rpq;// �f�w�ŒZ����(m)

			int vrshifts ;
            if( nShiftFlg == 1){
                vrshifts = (int)( ( subf1.offsettime + travelts + vrt - tvrt ) / (1.0/(double)Samplefrq) );
            }else{
                vrshifts = (int)( ( subf1.offsettime + travelts + vrt ) / (1.0/(double)Samplefrq) );
            }
			int midsift = datasize - vrshifts;//S�g�`�d���ԃV�t�g���Ԓ�

			GetGousei( &sss , &sss1 , midsift , vrshifts );// S�g����(NS,EW)

			//for( unsigned int pp = 0 ; pp < 	wdats2.size() ; pp++ ){
			//	sss2[pp] = sss[pp];
			//}   -> ���[�v�̊O��
		}
        for( unsigned int pp = 0 ; pp <	wdats2.size() ; pp++ ){
			sss2[pp] = sss[pp];
		}

		// �����k�������v�Z�p�֐�
        //xxi = Xi[l]/numofElement;
		//xxi = Xi[l];
		//ddi = Di[l];
		//XXDi = XDi[l];
		XXEQ = XEQ[l];

		// �ϕ����x�g�`�v�Z
		Ssekibun( &sss, &sekibuns, datasize, dt, 10, 0.6321);

		// �v���k�x�v�Z
		seismic_intensity( datasize, dt, sss2, intensity );

		// S�g�����x�ő��l
		double amax = *max_element(sss.begin(),sss.end());
		double amin = *min_element(sss.begin(),sss.end());
		double absamax = 0;
		if( fabs(amax) > fabs(amin) ) absamax = fabs(amax);
		else absamax = fabs(amin);

		// S�g���x�ő��l
		double vmax = *max_element(sekibuns.begin(),sekibuns.end());
		double vmin = *min_element(sekibuns.begin(),sekibuns.end());
		double absvmax = 0;
		if( fabs(vmax) > fabs(vmin) ) absvmax = fabs(vmax);
		else absvmax = fabs(vmin);

		// �����x�g�`�o��
		std::ostringstream filename_strm;
		filename_strm << "ac" << std::setw(6) << std::setfill('0') << OBBserno << ".dat";
		std::ofstream strm( filename_strm.str().c_str() );
		if( strm.is_open() ){
			output_in_kneta( strm, sss, dt, oridatasize, Samplefrq );
			strm.close();
		}else{
			std::cerr << "file out error." << std::endl;
		}

		// ���x�g�`�o��
		std::ostringstream filename_strm1;
		filename_strm1 << "ve" << std::setw(6) << std::setfill('0') << OBBserno << ".dat";
		std::ofstream strm1( filename_strm1.str().c_str() );
		if( strm1.is_open() ){
			output_in_knetv( strm1, sekibuns, dt, oridatasize, Samplefrq );
			strm1.close();
		}else{
			std::cerr << "file out error." << std::endl;
		}
        //���x�g�`�t�@�C�������t�B���^�[���[�`���ɓn���B
        double dFilteredMax=0.0;
        if(sHL_Flg != "NO"){
            filter_main(filename_strm1.str(),sHL_Flg,sC_Freq,&dFilteredMax);
        }
		// ���͌��ʃ��O�t�@�C���o��
        //�t�B���^�[���̌��ʂ��o�͂��邽�߂����Ɉړ�
		FILE *fpo1;
		if( (fpo1 = fopen( Outputlogfile , "a" )) == NULL )exit(1);
		fprintf( fpo1 , "%8.0lf\t%8.4lf\t%8.4lf\t%8.4lf\t%9.5lf\t%8.4lf\t%8.4lf\t%8.4lf\t%s\t%9.5lf\n"
			, OBBserno , OBLatitude , OBLongitude , absamax , absvmax ,
              minRpq/1000.0 , XXEQ , intensity , sHL_Flg.c_str(), dFilteredMax);
		fclose( fpo1 );
		// ���͌��ʉ��ʏo��
		printf(" CALC No. = %d  MAXACC(gal) = %f \n MAXVEL(kine) = %f FilteredMAXVEL(kine) = %f\n",
                             l+1 , absamax , absvmax, dFilteredMax );
		printf(" MINR(km) = %f Xeq(km) = %f INTENSITY = %f\n\n", minRpq/1000.0 , XXEQ , intensity );

		// �f�[�^�N���A
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
