  /*
   *  libpar.h include file.
   *
   *  Provide function definitions/prototypes for the 
   *  routines found in libpar.a.  Note these routines may change if 
   *  a new release of libpar is received.  Unfortunately this file 
   *  and the release are independent.
   */
    
#ifndef _LIBPAR_H
#define _LIBPAR_H 1

#ifdef __cplusplus
extern "C" {
#endif
    
extern int	countarg(char *name, char *type);
extern void	endarg (void);
extern void	endpar (void);
extern int	getarg (char *name, char *type, void *ptr_to_some_type);
extern int	getpar (char *name, char *type, void *ptr_to_some_type);
extern int	lenarg (char *name);
extern int	mstpar (char *name, char *type, void *ptr_to_some_type);
extern void	setarg (char *list, char *subname);
extern int	setpar (int argc, char **argv);
extern char    *getspar(char *name, char *defvalue);
extern char    *mstspar(char *name);
extern char    *getsarg(char *name, char *defvalue);
extern int	getbpar(char *name, int defvalue);
extern int	getdpar(char *name, int defvalue);
extern float	getfpar(char *name, float defvalue);
extern double	getffpar(char *name, double defvalue);
extern int      getlocation(char *keyname, char *location, int fatal);

#ifdef __cplusplus
}
#endif

#endif	/* _LIBPAR_H */
