void *check_malloc(size_t);
void *check_realloc(void *,size_t);
FILE *fopfile(char*, char*);
int opfile_ro(char *);
int opfile(char *);
int croptrfile(char *);
int reed(int, void *, int);
int rite(int, void *, int);

double frand(void);
double sfrand(long *);
double gaussian_rand(float *,float *,long *);

int gen_esg2006_stf(float *,float *,float *,int,float *,float *);
int gen_2tri_stf(float *,float *,float *,int,float *,float *);
int gen_ucsb_stf(float *,float *,float *,int,float *,float *);
int gen_ucsb2_stf(float *,float *,float *,int,float *,float *);
int gen_ucsbT_stf(float *,float *,float *,int,float *,float *);
int gen_ucsbvT_stf(float *,float *,float *,int,float *,float *);
int gen_brune_stf(float *,float *,float *,int,float *,float *);
int gen_cos_stf(float *,float *,float *,int,float *,float *);
int gen_seki_stf(float *,float *,float *,int,float *,float *);
int gen_ji2003_stf(float *,float *,float *,float *,int,float *);

void set_ne(float *,float *,float *,float *,float *,float *);
void set_ll(float *,float *,float *,float *,float *,float *);
void swap_in_place(int,char *);

void init_plane_srf(struct standrupformat *,float *,float *,int,int,float *,float *,float *,float *,float *,float *,float *,float *,float *);
/*
void load_slip_srf(struct standrupformat *,struct stfpar *,struct pointsource *);
void load_rupt_srf(struct standrupformat *,struct pointsource *,float *,float *);
*/

void read_srf(struct standrupformat *,char *,int);
void write_srf(struct standrupformat *,char *,int);
void write_srf1(struct standrupformat *,char *,int);
void write_srf2(struct standrupformat *,char *,int);
void copy_hcmnt(struct standrupformat *srf2,struct standrupformat *srf1);
void load_command_srf(struct standrupformat *,int,char **);
void load_seed_srf(struct standrupformat *srf,int,int);

void free_srf_stf(struct standrupformat *);

int write_xyz(char *,struct standrupformat *,char *,int,int,float *,float *,int,int,int,float *,float *,int,int);
void write_maxsvf(char *,struct standrupformat *,char *,int,float *,int);
void get_vmax2slip(char *,struct standrupformat *,char *,int);
int write_lld(char *,struct standrupformat *,int,float *,float *,char *,int);

void get_moment(struct standrupformat *,struct velmodel *,int);
void get_moment_rate(struct standrupformat *,struct velmodel *,int);
void read_velmodel(char *,struct velmodel *); 

void sum_srf(struct standrupformat *,struct standrupformat *,struct standrupformat *,float *);
void join_srf(struct standrupformat *,struct standrupformat *,struct standrupformat *);
void select_depths_srf(struct standrupformat *,struct standrupformat *,float *,float *);
void scale_srf(struct standrupformat *,struct standrupformat *,float *);

void replace_sdr(struct standrupformat *,struct standrupformat *,struct standrupformat *,int);
void read_Fvelmodel(char *, struct velmodel *);
void zapit(float *, int);
void gcproj(float *,float *,float *,float *,float *,double *,double *,double *,double *,int);
void gen_matrices(double *,double *,float *,float *,float *);
void dump_sliprate(char *,struct standrupformat *,int,float *,float *,float *);
