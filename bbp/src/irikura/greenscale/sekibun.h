#ifdef __cplusplus
extern "C" {
#endif

void Ssekibun(std::vector<double> *indata,
              std::vector<double> *outdata,
              int n,
              double s_time,
              double term,
              double gain);

#ifdef __cplusplus
}
#endif