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

	//すべてのフィールドが同じ要素数を持たない場合、0を返す
	std::vector< std::vector<double>* >::iterator ite = field.begin();
	unsigned int size = (*ite)->size();
	ite++;
	while( ite != field.end() ){
		if( size != (*ite)->size() ) return 0;
		ite++;
	}

	//正常終了
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
            //＃ではじまる行は読み飛ばす。
            continue;
        }
		std::istringstream strstrm( strbuf );
        for(;;){
			strstrm >> valbuf;
            if(valbuf.compare(0,1,"#") == 0){
                //＃がでたらループを抜ける。
                break;
            }
			(*ite)->push_back( valbuf ); ite++;

		}
	}

	//すべてのフィールドが同じ要素数を持たない場合、0を返す
	//std::vector< std::vector<string>* >::iterator ite = field.begin();
	ite = field.begin();
	unsigned int size = (*ite)->size();
	ite++;
	while( ite != field.end() ){
		if( size != (*ite)->size() ) return 0;
		ite++;
	}

	//正常終了
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

	//すべてのフィールドが同じ要素数を持たない場合、0を返す
	std::vector< std::vector<string>* >::iterator ite = field.begin();
	unsigned int size = (*ite)->size();
	ite++;
	while( ite != field.end() ){
		if( size != (*ite)->size() ) return 0;
		ite++;
	}

	//正常終了
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
    strstrm >> stemp;//連番
    strstrm >> stemp;//Xp
    strstrm >> stemp;//Yp
	while( !strstrm.eof() ){
		for( ite = field.begin(); ite != field.end(); i++){
            strstrm >> stemp;
            //#がでたらそれ以降はコメントとして戻る。
            //if(strcmp(stemp.substr(0,1).c_str(),"#") == 0){
            if(stemp.compare(0,1,"#") == 0){
                //行末の＃を見つけたらループを抜ける
                ;
                break;
                //return 99999;
            }
            //各層の物性値を各変数に格納する
            strstrmprop.clear();
            strstrmprop.str( strProp[i]);
            strstrmprop >> snum;
            strstrmprop >> svs;
            strstrmprop >> sro;
            strstrmprop >> sq;

            valbuf = atof(stemp.c_str());
            if(dUlayer == valbuf){
                //深度が前回と同じ場合，対象の層は存在しない。
                if(i == strProp.size() - 1)
                    break;
            }else{
			    //strstrm >> valbuf;
                //層厚としてデータを格納，VS，密度を格納
    			(*ite)->push_back( valbuf - dUlayer );   ite++;
    			(*ite)->push_back( atof(svs.c_str()) );  ite++;
    			(*ite)->push_back( atof(sro.c_str()) );  ite++;
                //深度を更新
                dUlayer = valbuf;
            }
		}
    }
    //最下層の物性を加える。層厚は0とする。
    strstrmprop.clear();
    strstrmprop.str( strProp[i]); //iは最大の引数になっている
    strstrmprop >> snum;
    strstrmprop >> svs;
    strstrmprop >> sro;
    strstrmprop >> sq;
    ite = field.begin();
    (*ite)->push_back( 0.0 );   ite++;
    (*ite)->push_back( atof(svs.c_str()) );  ite++;
    (*ite)->push_back( atof(sro.c_str()) );  ite++;


	//すべてのフィールドが同じ要素数を持たない場合、0を返す
	//std::vector< std::vector<double>* >::iterator ite = field.begin();
    ite = field.begin();
	unsigned int size = (*ite)->size();
	ite++;
	while( ite != field.end() ){
		if( size != (*ite)->size() ) return 0;
		ite++;
	}

	//正常終了
	return size;
}

