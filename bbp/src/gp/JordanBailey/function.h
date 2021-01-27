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
void find_4gf_v2(struct gfparam,struct gfheader *,float *,float *,float *,float *,float *);
void read_4gf(char *,char *,float *,int,struct gfheader *,struct gfparam,float *,int *,float *,float *);
void read_4gf_v2(char *,char *,float *,int,struct gfheader *,struct gfparam,int *,float *,float *);
void read_gf(char *,char *,float *,int,struct gfheader *,struct gfparam);

void reset_wt_4gf(struct gfheader *,float *,float *,float *);
void reset_wt_4gf_v2(struct gfheader *,float *,float *,float *,float *);

void rand_init(float *,float *,long *,int,int,int,int,int,int);
double gaus_rand(float *,float *,long *);
double sfrand(long *);

void sum_4gf(float *,int,float *,struct gfheader *,int,int,float *,float *,float *,int,struct mechparam);
void timeshift_gf(float *,int,float *,struct gfheader *,int,float *,float *,int);
void mech_4gf(float *,float *,struct gfheader *,struct gfparam gfp,int,struct mechparam,float *,float *);

void sum_4gf_v2(float *,int,float *,struct gfheader *,int,int,float *,float *,float *,int,struct mechparam,float *);
void timeshift_gf_v2(float *,int,float *,struct gfheader *,int,float *,float *,int);
void mech_4gf_v2(float *,float *,struct gfheader *,struct gfparam gfp,int,struct mechparam,float *,float *);

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
void get_srfpars_v2(struct standrupformat *,int,int,float *,float *,struct mechparam *);
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
