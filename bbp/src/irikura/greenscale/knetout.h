#include <stdio.h>
#include <math.h>
#include <complex>
#include <vector>
#include <stdlib.h>
#include <string.h>
#include "knetascii.h"
//
void output_in_knet( std::ostream& strm, std::vector<double>& x, double dt ,knetascii knet);
void output_in_kneta( std::ostream& strm, std::vector<double>& x, double dt, int oridatasize, int Samplefrq );
void output_in_knetv( std::ostream& strm, std::vector<double>& x, double dt, int oridatasize, int Samplefrq );
//
unsigned int read_csv( const char* filename, std::vector< std::vector<double>* > field );
unsigned int read_string( const char* filename, std::vector< std::vector<std::string>* > field );
unsigned int read_fault( const char* filename, std::vector< std::vector<std::string>* > field );
unsigned int read_char( const char* filename, std::vector< std::vector<char>* > field );
//
string& fnReplace(string& str , const string sb, const string sa);
string fnGetFileName(string str);
string fnGetExtName(string str);

