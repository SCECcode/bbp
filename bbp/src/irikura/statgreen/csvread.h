#ifdef __cplusplus
extern "C" {
#endif

unsigned int read_csv( const char* filename, std::vector< std::vector<double>* > field );
unsigned int read_sta( const char* filename, std::vector< std::vector<string>* > field );
unsigned int read_fault( const char* filename, std::vector< std::vector<string>* > field );

unsigned int read_soil( std::vector<std::string> strProp ,
                        std::string strline ,
                        std::vector< std::vector<double>* > field );

#ifdef __cplusplus
}
#endif
