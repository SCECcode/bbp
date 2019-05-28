#ifdef __cplusplus
extern "C" {
#endif

void output_in_knet( std::ostream& strm, std::vector< std::complex<double> >& x, double dt, std::string station_code, int oridatasize, double Samplefreq);
void output_in_plane( std::ostream& strm, std::vector< std::complex<double> >& x, double dt);

#ifdef __cplusplus
}
#endif
