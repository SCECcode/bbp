void *check_malloc(size_t);
void *check_realloc(void *, size_t);
FILE *fopfile(char*, char*);
int opfile_ro(char *);
int opfile(char *);
int croptrfile(char *);
int reed(int, void *, int);
int rite(int, void *, int);

void fortran_rite(int ,int , ...);
void conv2vrup(struct velmodel *,struct velmodel *,float *,float *,float *,float *,float *);
void get_ard_srf(struct standrupformat *,int,int,float *,float *,float *,float *,float *,float *,float *,struct geoprojection *);
void get_radazi(float *,float *,float *,float *,float *,float *,float *,float *,float *,float *);
void set_ne(float *,float *,float *,float *,float *,float *);
void set_ll(float *,float *,float *,float *,float *,float *);

void get_gfpars(struct gfparam *);
void find_4gf(struct gfparam,struct gfheader *,float *,float *,float *,float *);
void read_4gf(char *,char *,float *,int,struct gfheader *,struct gfparam,float *,int *,float *,float *);
void read_gf(char *,char *,float *,int,struct gfheader *,struct gfparam);

void rand_init(float *,float *,long *,int,int,int,int,int,int);
double gaus_rand(float *,float *,long *);
double sfrand(long *);

void sum_4gf(float *,int,float *,struct gfheader *,int,int,float *,float *,float *,int,struct mechparam);
void timeshift_gf(float *,int,float *,struct gfheader *,int,float *,float *,int);
void mech_4gf(float *,float *,struct gfheader *,struct gfparam gfp,int,struct mechparam,float *,float *);

void read_beroza(struct beroza *,char *rfile,float *);
void get_brmpars(struct beroza *,int,int,float *,float *,float *,float *);
void beroza_stf(struct beroza *,int,int,float *,float *,float *,int,float *,float *);

void read_okumura(struct okumura *,char *rfile,float *);
void get_ormpars(struct okumura *,int,int,float *,float *,float *,float *);
void okumura_stf(struct okumura *,int,int,float *,float *,float *,int,float *);

void read_gene(struct gene *,char *rfile,float *);
void get_grmpars(struct gene *,int,int,float *,float *,float *,float *,float *);
void gene_stf(struct gene *,int,int,float *,float *,float *,int,float *);

void read_rob(struct rob *,char *rfile,float *);
void get_rrmpars(struct rob *,int,int,float *,float *,float *,float *,float *,float *);
void rob_stf(struct rob *,int,int,float *,float *,float *,int,float *,float *);
int gen_rob_stf(struct rob *,int,int,float *,int,float *,float *);

void get_srfpars(struct standrupformat *,int,int,float *,float *,float *,float *,float *,struct mechparam *);
void get_srfparsOLD(struct standrupformat *,int,int,float *,float *,float *,int,int,struct mechparam *);
void srf_stf(struct standrupformat *,int,int,float *,float *,float *,int,float *,struct mechparam,float *);

char *skipval(int,char *);
void do_cnvlv(float *,float *,int,float *,int);
void sum_nostf(float *,float *,float *,int);
void write_seis(char *,char *,char *,char *,float *,float *,int,float *);

void swap_in_place(int,char *);

double nt_tol(float,int);

void getname_gf(char *str,char *gfname,struct gfheader *gfh,struct gfparam gfpar);
int check_name(char *str,char *list,int n,int blen);

void get_rupt(struct velmodel *, float *, float *, float *, float *, float *, double *, double *, float *);
void zapit(float *, int);
void reset_wt_4gf(struct gfheader *, float *, float *, float *);
void taper_norm(float *, float *, int);
void norm(float *, float *, int);
void cfft_r(struct complex *, int, int);
void gcproj(float *, float *, float *, float *, float *, double *, double *, double *, double *, int);
void resample(float *, int, float *, int, int, int, float *, float *);
void get_indx(float *, float *, float *, struct sgtindex *, struct geoprojection *);
void get_master_list(struct sgtparams *, int, long long *, int *);
void forfft(struct complex *, int, int);
void invfft(struct complex *, int, int);
void mech_sgt(float *, float *, struct sgtheader *, struct sgtparams *, int, struct mechparam, float *);
void sum_sgt(float *, int, float *, struct sgtparams *, struct sgtheader *, int, float *, float *, struct mechparam);
void read_sgt(struct sgtfileparams *, struct sgtmaster *, struct sgtindex *, struct sgtheader *, float *);
void *sgt_subset(struct sgtfileparams *, struct sgtfileparams *, struct sgtmaster *, struct sgtindex *, int, long long *, char *);
void set_geoproj(struct sgtmaster *, struct geoprojection *);
void find_sgt(struct sgtparams *, struct sgtmaster *, struct sgtindex *, struct sgtindex *, struct sgtindex *, float *, float *);
void makedir(char *);
void *get_sgt_subset(struct sgtfileparams *, struct sgtmaster *, struct sgtindex *, struct sgtheader *, int, long long *, float *);
