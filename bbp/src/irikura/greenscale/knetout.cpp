#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <complex>
#include <algorithm>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "knetascii.h"
#include "knetout.h"

void output_in_knet( std::ostream& strm, std::vector<double>& x, double dt,knetascii knet)
{

	std::vector<double> realx( x.size() );
	for( int i = 0; i < x.size(); i++ ) realx[i] = x[i];
	double maxval = *std::max_element( realx.begin(), realx.end() );
	double minval = *std::min_element( realx.begin(), realx.end() );
	double absmax = 0;
	if(fabs( maxval ) > fabs( minval )) absmax = fabs( maxval );
	else absmax = fabs( minval );

	knet.set_denominator( absmax );
	//knet.set_direction( "SH" );
	//knet.set_duration( realx.size() * dt );
	knet.set_max_value( absmax );
	knet.set_numerator( 1 << 23 );
	knet.set_data( realx );
	//knet.set_origin_time( "2001/01/01 00:00:00");
	//knet.set_record_time( "2001/01/01 00:00:00");
	//knet.set_sampling_frequency( knet.get_sampling_frequency() );
	//knet.set_site_code( "station" );

	strm << knet;
}
//---------------------------------------------------------------------------

void output_in_kneta( std::ostream& strm, std::vector<double>& x, double dt, int oridatasize, int Samplefrq )
{
	knetascii knet;
	knet.kind_of_data_ =  KNET_ACC;
	std::vector<double> realx( x.size() );
	for( int i = 0; i < x.size(); i++ ) realx[i] = x[i];
	realx.resize(oridatasize);
	double maxval = *std::max_element( realx.begin(), realx.end() );
	double minval = *std::min_element( realx.begin(), realx.end() );
	double absmax = 0;
	if(fabs( maxval ) > fabs( minval )) absmax = fabs( maxval );
	else absmax = fabs( minval );

	knet.set_data( realx );
	knet.set_denominator( absmax );
	knet.set_direction( "SH" );
	knet.set_duration( realx.size() * dt );
	knet.set_max_value( absmax );
	knet.set_numerator( 1 << 23 );
	knet.set_origin_time( "2001/01/01 00:00:00");
	knet.set_record_time( "2001/01/01 00:00:00");
	knet.set_sampling_frequency( Samplefrq );
//	knet.set_site_code( station_code );

	strm << knet;
}
//---------------------------------------------------------------------------
void output_in_knetv( std::ostream& strm, std::vector<double>& x, double dt, int oridatasize, int Samplefrq )
{
	knetascii knet;
	knet.kind_of_data_ =  KNET_VEL;
	std::vector<double> realx( x.size() );
	for( int i = 0; i < x.size(); i++ ) realx[i] = x[i];
	realx.resize(oridatasize);
	double maxval = *std::max_element( realx.begin(), realx.end() );
	double minval = *std::min_element( realx.begin(), realx.end() );
	double absmax = 0;
	if(fabs( maxval ) > fabs( minval )) absmax = fabs( maxval );
	else absmax = fabs( minval );

	knet.set_data( realx );
	knet.set_denominator( absmax );
	knet.set_direction( "SH" );
	knet.set_duration( realx.size() * dt );
	knet.set_max_value( absmax );
	knet.set_numerator( 1 << 23 );
	knet.set_origin_time( "2001/01/01 00:00:00");
	knet.set_record_time( "2001/01/01 00:00:00");
	knet.set_sampling_frequency( Samplefrq );
//	knet.set_site_code( station_code );

	strm << knet;
}
//---------------------------------------------------------------------------
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
unsigned int read_string( const char* filename, std::vector< std::vector<string>* > field )
{
	std::ifstream istrm( filename );
	if( !istrm.is_open() ){
		return 0;
	}

	string valbuf;
	std::string strbuf;
	while( std::getline( istrm, strbuf ) ){
        if(strbuf.compare(0,1,"#") == 0){
            continue;
        }
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
unsigned int read_char( const char* filename, std::vector< std::vector<char>* > field )
{
	std::ifstream istrm( filename );
	if( !istrm.is_open() ){
		return 0;
	}

	double valbuf;
	std::string strbuf;
	while( std::getline( istrm, strbuf ) ){
		std::istringstream strstrm( strbuf );
		for( std::vector< std::vector<char>* >::iterator ite = field.begin(); ite != field.end(); ite++ ){
			strstrm >> valbuf;
			(*ite)->push_back( valbuf );
		}
	}

	//すべてのフィールドが同じ要素数を持たない場合、0を返す
	std::vector< std::vector<char>* >::iterator ite = field.begin();
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
string& fnReplace(string& str , const string sb, const string sa)
{
    string::size_type n,nb=0;

    while((n = str.find(sb,nb)) != string::npos){
        str.replace(n,sb.size(),sa);
        nb = n + sa.size();
    }
    return str;
}
//---------------------------------------------------------------------------
string fnGetFileName(string str)
{
    int nDot;
    string strName;
    //'.'までの位置を探して，名前部分を取得する
    nDot = str.find('.');
    strName = str.substr(0,nDot);

    return strName;

}
//---------------------------------------------------------------------------
string fnGetExtName(string str)
{
    int nDot,nLeng;
    string strExt;
    //'.'までの位置を探して，拡張子部分を取得する。（.を含む）
    nDot = str.find('.');
    nLeng = str.length();
    strExt = str.substr(nDot,nLeng-nDot);

    return strExt;

}
//---------------------------------------------------------------------------

