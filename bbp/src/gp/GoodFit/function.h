void *check_malloc(int);
void *check_realloc(void *, int);
FILE *fopfile(char*, char*);
int opfile_ro(char *);
int opfile(char *);
int croptrfile(char *);
int reed(int, void *, int);
int rite(int, void *, int);

void swap_in_place(int,char *);

void read_specfile(char *,float *,float *,int,int,int);
void read_bbpfile(char *,float *,float *,float *,int);
void read_bbpfile_3comp(char *,float *,float *,float *,float *,int);
