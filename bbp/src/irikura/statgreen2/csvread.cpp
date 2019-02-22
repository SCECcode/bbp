#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <vector>
#include <complex>
#include <algorithm>

#include "knetascii.h"
#include "csvread.h"

using namespace std;

unsigned int read_csv( const char* filename, std::vector< std::vector<double>* > field )
{
	std::ifstream istrm( filename );
	if( !istrm.is_open() ){
		return 0;
	}

	double valbuf;
	std::string strbuf;
	while( std::getline( istrm, strbuf ) ){
		std::istringstream strstrm( strbuf );
		for( std::vector< std::vector<double>* >::iterator ite = field.begin(); ite != field.end(); ite++ ){
			strstrm >> valbuf;
			(*ite)->push_back( valbuf );
		}
	}

	//���ׂẴt�B�[���h�������v�f���������Ȃ��ꍇ�A0��Ԃ�
	std::vector< std::vector<double>* >::iterator ite = field.begin();
	unsigned int size = (*ite)->size();
	ite++;
	while( ite != field.end() ){
		if( size != (*ite)->size() ) return 0;
		ite++;
	}

	//����I��
	return size;
}
//---------------------------------------------------------------------------
unsigned int read_fault( const char* filename, std::vector< std::vector<string>* > field )
{
	std::ifstream istrm( filename );
	if( !istrm.is_open() ){
		return 0;
	}

	string valbuf;
	string strbuf;
    vector< vector<string>* >::iterator ite = field.begin();
	while( std::getline( istrm, strbuf ) ){
        if(strbuf.compare(0,1,"#") == 0){
            //���ł͂��܂�s�͓ǂݔ�΂��B
            continue;
        }
		std::istringstream strstrm( strbuf );
        for(;;){
			strstrm >> valbuf;
            if(valbuf.compare(0,1,"#") == 0){
                //�����ł��烋�[�v�𔲂���B
                break;
            }
			(*ite)->push_back( valbuf ); ite++;

		}
	}

	//���ׂẴt�B�[���h�������v�f���������Ȃ��ꍇ�A0��Ԃ�
	//std::vector< std::vector<string>* >::iterator ite = field.begin();
	ite = field.begin();
	unsigned int size = (*ite)->size();
	ite++;
	while( ite != field.end() ){
		if( size != (*ite)->size() ) return 0;
		ite++;
	}

	//����I��
	return size;
}
//---------------------------------------------------------------------------
unsigned int read_sta( const char* filename, std::vector< std::vector<string>* > field )
{
	std::ifstream istrm( filename );
	if( !istrm.is_open() ){
		return 0;
	}

	string valbuf;
	std::string strbuf;
	while( std::getline( istrm, strbuf ) ){
		std::istringstream strstrm( strbuf );
		for( std::vector< std::vector<string>* >::iterator ite = field.begin(); ite != field.end(); ite++ ){
			strstrm >> valbuf;
			(*ite)->push_back( valbuf );
		}
	}

	//���ׂẴt�B�[���h�������v�f���������Ȃ��ꍇ�A0��Ԃ�
	std::vector< std::vector<string>* >::iterator ite = field.begin();
	unsigned int size = (*ite)->size();
	ite++;
	while( ite != field.end() ){
		if( size != (*ite)->size() ) return 0;
		ite++;
	}

	//����I��
	return size;
}
//---------------------------------------------------------------------------
unsigned int read_soil( std::vector<std::string> strProp , std::string strline , std::vector< std::vector<double>* > field )
{

	double valbuf;
    double dUlayer=0.0;
    std::string stemp,snum,svs,sro,sq;
	std::istringstream strstrm( strline );
	std::istringstream strstrmprop;
    std::string strbuf;
    int i=0;
    std::vector< std::vector<double>* >::iterator ite;
    strstrm >> stemp;//�A��
    strstrm >> stemp;//Xp
    strstrm >> stemp;//Yp
	while( !strstrm.eof() ){
		for( ite = field.begin(); ite != field.end(); i++){
            strstrm >> stemp;
            //#���ł��炻��ȍ~�̓R�����g�Ƃ��Ė߂�B
            //if(strcmp(stemp.substr(0,1).c_str(),"#") == 0){
            if(stemp.compare(0,1,"#") == 0){
                //�s���́����������烋�[�v�𔲂���
                ;
                break;
                //return 99999;
            }
            //�e�w�̕����l���e�ϐ��Ɋi�[����
            strstrmprop.clear();
            strstrmprop.str( strProp[i]);
            strstrmprop >> snum;
            strstrmprop >> svs;
            strstrmprop >> sro;
            strstrmprop >> sq;

            valbuf = atof(stemp.c_str());
            if(dUlayer == valbuf){
                //�[�x���O��Ɠ����ꍇ�C�Ώۂ̑w�͑��݂��Ȃ��B
                if(i == strProp.size() - 1)
                    break;
            }else{
			    //strstrm >> valbuf;
                //�w���Ƃ��ăf�[�^���i�[�CVS�C���x���i�[
    			(*ite)->push_back( valbuf - dUlayer );   ite++;
    			(*ite)->push_back( atof(svs.c_str()) );  ite++;
    			(*ite)->push_back( atof(sro.c_str()) );  ite++;
                //�[�x���X�V
                dUlayer = valbuf;
            }
		}
    }
    //�ŉ��w�̕�����������B�w����0�Ƃ���B
    strstrmprop.clear();
    strstrmprop.str( strProp[i]); //i�͍ő�̈����ɂȂ��Ă���
    strstrmprop >> snum;
    strstrmprop >> svs;
    strstrmprop >> sro;
    strstrmprop >> sq;
    ite = field.begin();
    (*ite)->push_back( 0.0 );   ite++;
    (*ite)->push_back( atof(svs.c_str()) );  ite++;
    (*ite)->push_back( atof(sro.c_str()) );  ite++;


	//���ׂẴt�B�[���h�������v�f���������Ȃ��ꍇ�A0��Ԃ�
	//std::vector< std::vector<double>* >::iterator ite = field.begin();
    ite = field.begin();
	unsigned int size = (*ite)->size();
	ite++;
	while( ite != field.end() ){
		if( size != (*ite)->size() ) return 0;
		ite++;
	}

	//����I��
	return size;
}

