#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <vector>
#include <complex>
#include <algorithm>

#include "knetascii.h"
#include "knetout.h"

using namespace std;

void output_in_knet( std::ostream& strm, std::vector< std::complex<double> >& x, double dt, std::string station_code, int oridatasize, double Samplefreq )
{
	knetascii knet;
	std::vector<double> realx( x.size() );
	for( int i = 0; i < realx.size() ; i++ ) realx[i] = (x[i].real()/(x.size())*(8388608/2000));
	realx.resize(oridatasize);
	double maxval = *std::max_element( realx.begin(), realx.end() );
	double minval = *std::min_element( realx.begin(), realx.end() );
	double absmax = 0;
	if(fabs( maxval ) > fabs( minval )) absmax = fabs( maxval );
	else absmax = fabs( minval );

	printf("1DAMP-WAVE MAX(gal) = %f\n", absmax/(8388608/2000));

	knet.set_data( realx );
	knet.set_denominator( absmax );
	knet.set_direction( "SH" );
	knet.set_duration( realx.size() * dt );
	knet.set_max_value( absmax/(8388608/2000));
	knet.set_numerator( 1 << 23 );
	knet.set_origin_time( "2001/01/01 00:00:00");
	knet.set_record_time( "2001/01/01 00:00:00");
	knet.set_sampling_frequency( Samplefreq );
	knet.set_site_code( station_code );

	strm << knet;
}

void output_in_plane( std::ostream& strm, std::vector< std::complex<double> >& x, double dt)
{
	for( int i = 0; i < x.size(); i++ ){
		strm << i * dt << " " << x[i].real() << std::endl;
	}
}
