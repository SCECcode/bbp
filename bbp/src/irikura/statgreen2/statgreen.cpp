//#pragma warning(disable:4786)

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

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <vector>
#include <complex>
#include <algorithm>

#include "VerticalSHResponse.h"
#include "FastFourierTransform.h"
#include "knetascii.h"
#include "BasenSpectrum.h"
#include "calcdistance.h"
#include "statgreen.h"
#include "knetout.h"
#include "csvread.h"

using namespace std;

// CVerticalSHResponse�N���X�̃T���v���v���O����
int main( int argc, char* argv[] )
{
    vector<double> frq, acc, wave, wave2;
	CBasenSpectrum spec;
	CQvalue qspec;
	knetascii knet;
	ofstream ost1;
    unsigned int lStNum,lEndNum;

	if( argc == 5 ){
        lStNum = lEndNum = 0;
	}else if(argc == 7){
        lStNum = lEndNum = 0;
    }else{
    	std::cerr << "Usage: statgreen station faultparam soil phase start_num end_num " << std::endl;
    	std::cerr << "       or " << std::endl;
    	std::cerr << "Usage: statgreen station faultparam soil phase " << std::endl;
		return 1;
	}

    //------------------- �X�e�[�V�����̍s�ǂ� ---------------------
	//char* stnumber  = argv[1];
	//int stationnum  = atoi(stnumber);
    //------------------- �X�e�[�V�����t�@�C���̓ǂݍ��� -----------
	char* station_filename  = argv[1];
	//std::vector<double> stationNo;
	//std::vector<double> slatitude;
	//std::vector<double> slongitude;
	//std::vector<double> sdepth;
    //<- 060130 add
	vector<string> stationNo;
	vector<string> slatitude;
	vector<string> slongitude;
	vector<string> sdepth;
    vector<string> sSASFileName; //
    //->

	double sslatitude;
	double sslongitude;
	double ssdepth;
	double dR_s;

	//�����l���ǂݍ���
	//field��vector���o�^
	//std::vector< std::vector<double>* > station;
	std::vector< std::vector<string>* > station;
	station.push_back( &stationNo );
	station.push_back( &slatitude );
	station.push_back( &slongitude  );
	station.push_back( &sdepth );
	station.push_back( &sSASFileName );

	//�ǂݍ���
	unsigned int numofStation = read_sta( station_filename, station );
	if( argc == 5 ){
        lStNum  = 0;
        lEndNum = numofStation;
	}else if(argc == 7){
        lStNum  = atoi(argv[5]);
        if(lStNum > 0){
            //�ϑ��_�ԍ��ɂ��킹��
            lStNum = lStNum - 1;
        }
        lEndNum = atoi(argv[6]) ;
        if(lEndNum > numofStation){
            lEndNum = numofStation;
        }
    }


	if( numofStation == 0 ){
		std::cerr << station_filename << " not found." << std::endl;
		return 1;
	}
	std::cerr << "station " << numofStation << std::endl;

	//-------------------- �f�w�p�����[�^�̓ǂݍ��� --------------------
	char* fault_filename  = argv[2];

    /*
	std::vector<double> Elatitude;
	std::vector<double> Elongitude;
	std::vector<double> Edepth;
	std::vector<double> ESamplefrq;
	std::vector<double> ESampletime;
	std::vector<double> Em_dRo_pq;
	std::vector<double> Em_dRo_sb;
	std::vector<double> Em_dBeta_pq;
	std::vector<double> Em_dBeta_sb;
	std::vector<double> Em_dM0_pq;
	std::vector<double> Em_dSig_pq;
	std::vector<double> Em_dFrad;
	std::vector<double> Em_dFmax_pq;
	std::vector<double> Em_dM;
	std::vector<double> Em_dDepend_1;
	std::vector<double> Em_dDepend_2;
	std::vector<double> Em_dF;
	std::vector<double> Em_dDepend_3;
    */
	std::vector<string> Elatitude;
	std::vector<string> Elongitude;
	std::vector<string> Edepth;
	std::vector<string> Deslatitude;
	std::vector<string> Deslongitude;
	std::vector<string> Desdepth;
	std::vector<string> ESamplefrq;
	std::vector<string> ESampletime;
	std::vector<string> Em_dRo_pq;
	std::vector<string> Em_dRo_sb;
	std::vector<string> Em_dBeta_pq;
	std::vector<string> Em_dBeta_sb;
	std::vector<string> Em_dM0_pq;
	std::vector<string> Em_dSig_pq;
	std::vector<string> Em_dFrad;
	std::vector<string> Em_dFmax_pq;
	std::vector<string> Em_dM;
	std::vector<string> Em_dDepend_1;
	std::vector<string> Em_dDepend_2;
	std::vector<string> Em_dF;
	std::vector<string> Em_dDepend_3;
    //add green scale param
    std::vector<string> Em_sLength;
    std::vector<string> Em_sWidth;
    std::vector<string> Em_sSlip;
    std::vector<string> Em_sRigidity;
    std::vector<string> Em_sVr;
    std::vector<string> Em_sOffset;

	//�����l���ǂݍ���
	//field��vector���o�^
	//std::vector< std::vector<double>* > faultparam;
	std::vector< std::vector<string>* > faultparam;
	faultparam.push_back( &Elatitude );
	faultparam.push_back( &Elongitude );
	faultparam.push_back( &Edepth  );
	faultparam.push_back( &Deslatitude );
	faultparam.push_back( &Deslongitude );
	faultparam.push_back( &Desdepth  );
	faultparam.push_back( &ESamplefrq );
	faultparam.push_back( &ESampletime );
	faultparam.push_back( &Em_dRo_pq );
	faultparam.push_back( &Em_dRo_sb );
	faultparam.push_back( &Em_dBeta_pq );
	faultparam.push_back( &Em_dBeta_sb );
	faultparam.push_back( &Em_dM0_pq );
	faultparam.push_back( &Em_dSig_pq );
	faultparam.push_back( &Em_dFrad );
	faultparam.push_back( &Em_dFmax_pq );
	faultparam.push_back( &Em_dM );
	faultparam.push_back( &Em_dDepend_1 );
	faultparam.push_back( &Em_dDepend_2 );
	faultparam.push_back( &Em_dF );
	faultparam.push_back( &Em_dDepend_3 );
    // add green scale param
    faultparam.push_back( &Em_sLength );
    faultparam.push_back( &Em_sWidth );
    faultparam.push_back( &Em_sSlip );
    faultparam.push_back( &Em_sRigidity );
    faultparam.push_back( &Em_sVr );
    faultparam.push_back( &Em_sOffset );

	// �ǂݍ���
	//unsigned int numofefault = read_csv( fault_filename, faultparam );
	unsigned int numofefault = read_fault( fault_filename, faultparam );

	if( numofefault == 0 ){
		std::cerr << fault_filename << " not found." << std::endl;
		return 1;
	}
	std::cerr << "fault_eparam " << numofefault << std::endl;
    //<- 20060131 �t�H�[�}�b�g�ύX�ɂ����ǂݍ��݃��[�`���̕ύX
    //�y���f�[�^�̓ǂݍ���
	char* soil_filename  = argv[3];
	std::vector<double> velocity;
	std::vector<double> density;
	std::vector<double> thickness;
	// �����l���ǂݍ���
	// field��vector���o�^
	std::vector< std::vector<double>* > field;
	field.push_back( &thickness );
	field.push_back( &velocity  );
	field.push_back( &density );
    //�t�@�C�����J��
    ifstream istrm( soil_filename );
    if( !istrm.is_open() ){
        std::cerr << soil_filename << " not found." << std::endl;
	    return 1;
    }
    //�s���ƂɃf�[�^���i�[�����B�����l�����Ɛ[�x�������ʂ̃x�N�g���ɓ������B
    vector<string> strSoilProp;  //�����l�p�������z���FVs,Ro,Q
    vector<string> strSoilDepth; //�[�x���E�p�������z���FDepth
    string strbuf;               //�e���|�����o�b�t�@
    bool bProp = false;          //�����lpart�t���O
    bool bDepth = false;         //�[�xpart�t���O
    while( getline(istrm,strbuf)){
        //�w�b�_���s����#���Ȃ��Ȃ��܂ł͋��ǂ�
        //#�ł͂Ȃ��Ȃ��CbProp=false�̂Ƃ�true�ɂ���
        if( strbuf.compare(0,1,"#") != 0 && !bProp){
            bProp =true;
        }
        //bProp=true�ŕ����l���ǂݍ��݁C�[���̃f�[�^���o�Ă�����bDepth=true�ɂ���
        if(bProp && strbuf.compare(0,1,"#")==0 ) {
            bDepth = true;
        }
        //���ۂɃf�[�^���ǂݍ��ޕ���
        if(bProp && !bDepth ){
            //�����l�f�[�^�����̓ǂݍ���
            strSoilProp.push_back(strbuf);
        }else if( bProp && bDepth ){
            //�[�x�f�[�^�����̓ǂݍ���
            //�s���Ɂ��������Γǂݔ��΂�
            if(strbuf.compare(0,1,"#") != 0){
                strSoilDepth.push_back(strbuf);
            }
        }
        //���̒i�K�ł́C�s��������#�ȍ~�̃R�����g�͂��̂܂܎c���B
        //�f�[�^���z���Ɋi�[�����i�K�ŃR�����g�������B
    }
    // �ϑ��_���Ƃ̌v�Z
    //for( int k = stationnum-1 ; k < stationnum  ; k++ )
    //for( int k = 0 ; k < stationnum  ; k++ )
    for( unsigned int k = lStNum ; k < lEndNum ; k++){

        printf("CALC-No. = %d\n ", k+1 );//�v�Z�ϑ��_�ʒu

		double eelatitude  = atof(Elatitude[0].c_str());  // �f�w�����_�ܓx
		double eelongitude = atof(Elongitude[0].c_str()); // �f�w�����_�o�x
		double eedepth     = atof(Edepth[0].c_str());     // �f�w�����_�[�x
		double Samplefreq  = atof(ESamplefrq[0].c_str());
		double Sample      = atof(ESampletime[0].c_str());

		// �t�[���G�p���T�C�Y
		int oridatasize =  (int)(Samplefreq * Sample);
		int datasize;
		for(int nn = 0 ;  nn < 16 ; nn++ ){
			if( oridatasize > (int)pow(2.0, nn) ) datasize = (int)pow(2.0,(nn+1));
		}

		spec.SetRoPQ(atof(Em_dRo_pq[0].c_str()));
		spec.SetRoSB(atof(Em_dRo_sb[0].c_str()));
		spec.SetBetaPQ(atof(Em_dBeta_pq[0].c_str()));
		spec.SetBetaSB(atof(Em_dBeta_sb[0].c_str()));
		spec.SetM0PQ(atof(Em_dM0_pq[0].c_str()));
		spec.SetStress(atof(Em_dSig_pq[0].c_str()));
		spec.SetRadiation(atof(Em_dFrad[0].c_str()));
		spec.SetFMax(atof(Em_dFmax_pq[0].c_str()));
		spec.SetM(atof(Em_dM[0].c_str()));

		spec.m_cQval.Set(atof(Em_dDepend_1[0].c_str()), atof(Em_dDepend_2[0].c_str()),
                         atof(Em_dF[0].c_str()), atof(Em_dDepend_3[0].c_str()));

		//sslatitude = slatitude[k];
		//sslongitude = slongitude[k];
		//ssdepth = sdepth[k];
		sslatitude = atof(slatitude[k].c_str());
		sslongitude = atof(slongitude[k].c_str());
		ssdepth = atof(sdepth[k].c_str());

		Distance( eelongitude, eelatitude, eedepth ,sslongitude, sslatitude, ssdepth , dR_s);

		//------------------ �y���p�����[�^�̓ǂݍ��� ----------------------

		double dt = 1.0/(double)Samplefreq;

		//char* layer_filename  = argv[4];
		//std::vector<double> velocity;
		//std::vector<double> density;
		//std::vector<double> thickness;

		// �����l���ǂݍ���
		// field��vector���o�^
		//std::vector< std::vector<double>* > field;
		//field.push_back( &thickness );
		//field.push_back( &velocity  );
		//field.push_back( &density );

		//unsigned int numofLayer = read_csv( layer_filename, field );
        //<- 20060130 �w���Ƃɕ����l���ǂݍ��ރ��[�`�����ύX
        //�x�N�g���ϐ���������
        velocity.clear();
        density.clear();
        thickness.clear();
		unsigned int numofLayer = read_soil( strSoilProp, strSoilDepth[k], field );
        //->
		if( numofLayer == 0 ){
			//std::cerr << layer_filename << " not found." << std::endl;
			std::cerr << " layer not found." << std::endl;
			return 1;
		}

		//------------------ WAV�t�@�C���̓ǂݍ��� ----------------------------

		char* wav_filename  = argv[4];
		std::vector<double> SAngle;

		// �����l���ǂݍ���
		// field��vector���o�^
		std::vector< std::vector<double>* > wavfield;
		wavfield.push_back( &SAngle );

		unsigned int numofWav = read_csv( wav_filename, wavfield );

		if( numofWav == 0 ){
			std::cerr << wav_filename << " not found." << std::endl;
			return 1;
		}

		//------------------ �o�͔g�`�t�@�C���̐ݒ� --------------------------

		//char* output_filename = argv[6];
        string stempfilename;
        stempfilename  = sSASFileName[k].substr(1,sSASFileName[k].length()-1);
		char* output_filename ;//= (char*)stempfilename.c_str() ;
        output_filename = new char[stempfilename.length() + 1];
        strcpy(output_filename,stempfilename.c_str());
        //strlen(output_filename);
		//----------------------- �v�Z�����яo�� -----------------------------
		int Ssize = SAngle.size()-1;//���������_���ʑ��̐�

		spec.GetBasenWave(dR_s/1000 ,(int)Samplefreq, datasize, wave, wave2, SAngle, Ssize);

		// �n�k���Ղ̔g�`�t�@�C���o��
		double maxval = *std::max_element( wave2.begin(), wave2.end() );
		double minval = *std::min_element( wave2.begin(), wave2.end() );
		double absmax = 0;

		if(fabs( maxval ) > fabs( minval )) absmax = fabs( maxval );
		else absmax = fabs( minval );

		printf("SB-WAVE MAX(gal) = %f\n ",absmax );

		knet.set_sampling_frequency((int)Samplefreq);
		knet.set_duration(oridatasize/Samplefreq);
		knet.set_max_value( absmax );
		knet.set_data(wave2);
		ost1.open(output_filename);
		if(ost1.is_open()) ost1 << knet;
		ost1.close();

		/*------------------------------------------------------------*/

		// SH�g�̐������ˉ����I�u�W�F�N�g���쐬���ĉ������v�Z����
		CVerticalSHResponse vertShResp( velocity, density, thickness );
		vertShResp.CalcResponse( datasize, dt );

		// �n�k���Տ��̔g�`���ŉ��ʂɓ��˂����Ƃ��̉����g�`���o�͂���
		std::vector<std::complex<double> > inputWave( vertShResp.NumofSample() );
		if( !inputw( inputWave, dt, wave2 ) ) std::cerr << "inputwave error" << std::endl;
		CFastFourierTransform::Trans( inputWave );

		// ���g�������擾��Convolution���o��
		std::vector<std::complex<double> > response( vertShResp.NumofSample() );
		std::vector<std::complex<double> > outputWave( vertShResp.NumofSample() );

        //iLayer==0�̎������o�͂��Ȃ��̂ŁA���[�v�ł܂킷�K�v�Ȃ��̂ł́H�H
		for( unsigned int iLayer = 0; iLayer < numofLayer; iLayer++ ){
			// ���g�������擾
			vertShResp.GetSpectrum( iLayer, response, dt );
			// Convolution
			for( unsigned int j = 0; j < outputWave.size(); j++ ) {
                outputWave[j] = (inputWave[j] * response[j]);
            }
            // ���n��
            CFastFourierTransform::Inverse( outputWave );

			// �H�w�I���Ղ̔g�`�t�@�C���o��
			if( iLayer == 0 ){

				std::ostringstream filename_strm;
				filename_strm << "S" << output_filename;

				std::ostringstream station_code_strm;
				station_code_strm << "S" << output_filename;

				std::ofstream strm( filename_strm.str().c_str() );
				if( strm.is_open() ){
					output_in_knet( strm, outputWave, dt, station_code_strm.str(), oridatasize, Samplefreq );
					strm.close();
				}else{
                    std::cerr << "file out error." << std::endl;
                }
			}
		}
        delete output_filename;
	}
	return 0;
}

#include <cmath>
bool inputw( std::vector< std::complex<double> >& x, double dt, std::vector<double>& wave )
{
	int j = 0;
	for( unsigned int i = 0; i < x.size(); i++ ){
        //x[i].real( wave[i] );
        //x[i] =  complex<double>(wave[i],0) ;
        x[i] = wave[i];
	}
	return true;
}
