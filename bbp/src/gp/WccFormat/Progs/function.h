//double pointer to handle realloc                                            
void wcc_resamp_arbdt(int param_string_len, char** param_string,
		      float** s1, struct statdata* head1);
void wcc_siteamp14(int param_string_len, char** param_string, float** s1, struct statdata* head1);
float wcc_getpeak(int param_string_len, char** param_string, float* seis, struct statdata* head1);
void wcc_tfilter (int param_string_len, char** param_string, float* s1, struct statdata* shead1);
void wcc_add(int param_string_len, char** param_string, float* s1, struct statdata* shead1, float* s2, struct statdata* shead2, float* p, struct statdata* shead3);
void integ_diff(int param_string_len, char** param_string, float* seis, struct statdata* shead);

void *check_malloc(size_t);
void *check_realloc(void *, size_t);
float *read_wccseis(char *, struct statdata *, float *, int);
void write_wccseis(char *, struct statdata *, float *, int);
FILE *fopfile(char*, char*);
int opfile_ro(char *);
int opfile(char *);
int croptrfile(char *);
int reed(int, void *, int);
int rite(int, void *, int);
void getheader(char *,struct statdata *);

void swap_in_place(int,char *);
void rotate(int, float *, float *, float *, float *, float *);
void forfft(struct complex *, int, int);
void invfft(struct complex *, int, int);
void dft(struct complex *, struct complex *, int, int);
void cfft_r(struct complex *, int, int);
void czero(struct complex *, int);
void makedir(char *);
int getname(char *, char *, int *);
void set_fullpath(char *, char *, char *);
void lp_filter(struct complex *, struct complex *, int, float *, float *, int);
void hp_filter(struct complex *, struct complex *, int, float *, float *, int);
