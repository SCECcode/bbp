/* 
 * NAME
 *	getpar
 *
 * DESCRIPTION
 *	retreive command line arguments.  Acquired from Caltech.
 *
 *	See getpar.3 for details.
 * copyright (c) Robert W. Clayton
 *		 Seismological Laboratory
 *		 Caltech
 *		 Pasadena, CA 91125
 *
 * Getpar routines:
 *
 * Externally visible routines:
 *
 *		setpar(argc,argv)
 *		getpar(name,type,valptr)
 *		mstpar(name,type,valptr)
 *		endpar()
 *
 * To get C-version:
 *		cc -c getpar.c
 *
 * To get F77-version:
 *		cp getpar.c fgetpar.c
 *		cc -c -DFORTRAN fgetpar.c
 *		rm fgetpar.c
 *
 * To get the environment processing stuff add the flag
 *-DENVIRONMENT to each of the cc's above.
 *
 *  Modification History:
 *  ---------------------
 *   3    94    Blair Zajac     Changed 
 *   7    91    Jeff Given      Added getspar/mstspar: new functions interface
 *                              for string parameters (C only)
 *   6    91    Jeff Given      Added escape for new line in par files
 *  04/01/91    Jeff Given      Added parameter substitution capability
 *
 *  06/08/90	Glenn Nelson	(via Richard Stead) added "ENDPAR" command
 *				line feature, added fix for booleans in
 *				parfiles, added better way of obtaining
 *				environment, corrected spelling
 *  06/08/90	Richard Stead	cleaned up some of lint's complaints
 *  05/25/90	Doug Neuhauser	Changed arglist elements from pointers to
 *				integer offsets into ARGBUF to allow
 *				arbitrary realloc calls.
 * ??/??/??	Rob Clayton	Original coding.
 *
 *SccsId: @(#)getpar.c	56.1 10/25/93
 */
#include	<stdio.h>
#include	<stdlib.h>
#include	<stdarg.h>
#include	<string.h>
#include	"libget.h"

#define MAXLINE		5120	/* max length of line in par file */
#define MAXNAME		256	/* max length of name */
#define MAXVALUE	5120	/* max length of value */
#define MAXFILENAME	256	/* max length of par file name */
#define MAXVECTOR	20	/* max # of elements for unspecified vectors */
#define GETPAR_ERROR	100	/* exit status for getpar error */
#define GETPAR_STOP	101	/* exit status for STOP or mstpar */
#define MAXPARLEVEL	8	/* max recursion level for par files */

#ifdef FORTRAN
#define GETPAR	getpar_
#define MSTPAR	mstpar_
#define ENDPAR	endpar_
#else
#define GETPAR	getpar
#define MSTPAR	mstpar
#define ENDPAR	endpar
#endif

#define INIT	 1	/* bits for FLAGS (ext_par.argflags) */
#define STOP	 2
#define LIST	 4
#define END_PAR	 8
#define VERBOSE	16

#define LISTINC		32	/* increment size for arglist */
#define BUFINC		1024	/* increment size for argbuf */

struct arglist		/* structure of list set up by setpar */
   {
	int argname_offset;
	int argval_offset;
	int hash;
   };
struct ext_par		/* global variables for getpar */
   {
	char *progname;
	int argflags;
	struct arglist *arglist;
	struct arglist *arghead;
	char *argbuf;
	int nlist;
	int nbuf;
	int listmax;
	int bufmax;
	FILE *listout;
  }	ext_par;


/* abbreviations: */
#define AL 		struct arglist
#define PROGNAME	ext_par.progname
#define FLAGS		ext_par.argflags
#define ARGLIST		ext_par.arglist
#define ARGHEAD		ext_par.arghead
#define ARGBUF		ext_par.argbuf
#define NLIST		ext_par.nlist
#define NBUF		ext_par.nbuf
#define LISTMAX		ext_par.listmax
#define BUFMAX		ext_par.bufmax
#define LISTFILE	ext_par.listout

static int gp_getvector(char *list, char *type, void *val);
static int gp_compute_hash(char *s);
static char *gp_fgets(char *line, int maxline, FILE *file);
static FILE *gp_create_dump(char *fname, char *filetype);
static void gp_add_entry(char *name, char *value);
static void gp_close_dump(FILE *file);
static void gp_do_par_file(char *fname, int level);
static void gp_subpar(char **apl, char **apv);
static void gp_getpar_err(char *subname, char *format, ...);
static void gp_do_environment(int ac, char **av);

int
#ifdef FORTRAN
setpar_()
#else
setpar(int ac, char **av)	/* set up arglist & process INPUT command */
#endif
   {
	register char *pl, *pn, *pv;
	char  t, name[MAXNAME], value[MAXVALUE];

	FILE *file;
	int i, addflags, nevlist, endsetpar = 0;
	struct arglist *alptr;
	
	char  *apl, *apv;

#ifdef FORTRAN
	int ac; char **av;
	extern int xargc; extern char **xargv;
	ac= xargc; av= xargv;
#endif

	if(av != (char **) NULL)
	{
		PROGNAME = *av;
	}
	else
	{
		PROGNAME = (char *) NULL;
		ac = 0;
	}
	
	FLAGS= INIT;
	LISTFILE= stderr;

	ARGLIST= NULL;
	ARGBUF = NULL;
	NLIST= NBUF= LISTMAX= BUFMAX= 0;
#ifdef ENVIRONMENT
	gp_do_environment(ac,av);
#endif
	nevlist= NLIST;
	while(--ac > 0 && endsetpar == 0)
	   {
		av++;
		pl= *av;
		while(*pl == ' ' || *pl == '\t') pl++;
		/* get name */
		pn= name;
		while(*pl != '=' && *pl != '\0') *pn++ = *pl++;
		*pn++ = '\0';
		/* get value */
		if(*pl == '=') pl++;
		*value = '\0';
		pv=value;

		if(*pl == '"' || *pl == '\'')
		{
			t= *pl++;
			while(*pl != '\0')
			{
				if(*pl == t)
				{
					if(pl[-1] != '\\') break;
					pv[-1]= t;
					pl++;
				}
				else
				{
					if(*pl == '$')
					{
						apl = pl;
						apv = pv;
						gp_subpar(&apl, &apv);
						pl = apl;
						pv = apv;
					}
					else *pv++ = *pl++;
				}
			}
		}
		else	
		{
			while(*pl) 
				if(*pl == '$')
				{
					apl = pl;
					apv = pv;
					gp_subpar(&apl, &apv);
					pl = apl;
					pv = apv;
				}
				else *pv++ = *pl++;
		}
		*pv= '\0';
		if(name[0] == '-') gp_add_entry("SWITCH",&name[1]);
		else		gp_add_entry(name,value);
		if(strcmp("par",name)==0) /* par file */
			gp_do_par_file(value,1);

	/* Added by Glenn Nelson (nelson@ollie.UCSC.EDU) to allow mixture
	   of getpar() and ordinary command line stuff. */
		if (strcmp("ENDPAR", name) == 0) endsetpar = 1;
	   }

	ARGHEAD= ARGLIST;

#ifdef ENVIRONMENT
	*value= '\0';
	if(GETPAR("NOENV","b",value)) ARGHEAD= ARGLIST+ nevlist;
#endif
	addflags= 0;
	*value= '\0';
	if(GETPAR("STOP","b",value)) addflags |= STOP;
	*value= '\0';
	if(GETPAR("VERBOSE","b",value)) addflags |= VERBOSE;
	*value= '\0';
#ifdef FORTRAN
	if(GETPAR("LIST","s",value, 0, 0, 0))
#else
	if(GETPAR("LIST","s",value))
#endif
	   {
		addflags |= LIST;
		LISTFILE =gp_create_dump(value,"list");
	   }
	*value= '\0';
#ifdef FORTRAN
	if(GETPAR("INPUT","s",value, 0, 0, 0))
#else
	if(GETPAR("INPUT","s",value))
#endif
	   {
		file =gp_create_dump(value,"list input");
		fprintf(file,"%s: getpar input listing\n",PROGNAME);
		for(i=0, alptr=ARGLIST; i<NLIST; i++, alptr++)
		   {
			fprintf(file,"%3d: %16s = %s\n",
				i,ARGBUF+alptr->argname_offset, 
				ARGBUF+alptr->argval_offset);
		   }
		gp_close_dump(file);
	   }
	FLAGS |= addflags;

	/* Added by Glenn Nelson (nelson@ollie.UCSC.EDU) to allow setpar()
	   to terminate before all command line args are exhausted. */
	return ac;
   }

/* add an entry to arglist, expanding memory if necessary */
void gp_add_entry(char *name, char *value)
   {
	struct arglist *alptr;
	int len;
	register char *ptr;

	/* check arglist memory */
	if(NLIST >= LISTMAX)
	   {
		LISTMAX += LISTINC;
		if(ARGLIST == NULL)
			ARGLIST= (AL *)malloc(LISTMAX * sizeof(AL));
		 else	ARGLIST= (AL *)realloc(ARGLIST,LISTMAX * sizeof(AL));
	   }
	/* check argbuf memory */
	len= strlen(name) + strlen(value) + 2; /* +2 for terminating nulls */
	while(NBUF+len >= BUFMAX)
	   {
		BUFMAX += BUFINC;
		if(ARGBUF == NULL)
			ARGBUF= (char *)malloc(BUFMAX);
		 else	ARGBUF= (char *)realloc(ARGBUF,BUFMAX);
	   }
	if(ARGBUF == NULL || ARGLIST == NULL)
		gp_getpar_err("setpar","cannot allocate memory");

	/* add name */
	alptr= ARGLIST + NLIST;
	alptr->hash= gp_compute_hash(name);
	alptr->argname_offset = NBUF;
	ptr= ARGBUF + NBUF;
	do *ptr++ = *name; while(*name++);

	/* add value */
	NBUF += len;
	alptr->argval_offset= ptr - ARGBUF;
	do *ptr++ = *value; while(*value++);
	NLIST++;
   }

#define BETTER_WAY	/* The environment is always available (as
			   suggested by Glenn Nelson (nelson@ollie.UCSC.EDU) */

void gp_do_environment(int ac, char **av)
   {
	char **ae;
	register char *pl, *pn, *pv;
	char name[MAXNAME], value[MAXVALUE], t;
#ifdef BETTER_WAY
	extern char     **environ;
#endif

	/* The environ pointer ae, is assumed to have a specific relation
	   to the arg pointer av. This may not be portable. */
#ifndef BETTER_WAY
	ae= av +(ac+1);
	if(ae == NULL) return;
#else
	ae = environ;
#endif

	while(*ae != NULL)
	   {
		pl= *ae++;
		while(*pl == ' ' || *pl == '\t') pl++;
		/* get name */
		pn= name;
		while(*pl != '=' && *pl != '\0') *pn++ = *pl++;
		*pn = '\0';
		if(strcmp("NOENV",pn) == 0) return;

		/* get value */
		if(*pl == '=') pl++;
		pv= value;
		if(*pl == '"' || *pl == '\'')
		   {
			t= *pl++;
			while(*pl != '\0')
			   {
				if(*pl == t)
				   {
					if(pl[-1] != '\\') break;
					pv[-1]= t;
					pl++;
				   }
				 else	*pv++ = *pl++;
			   }
		   }
		 else	while(*pl) *pv++ = *pl++;
		*pv= '\0';
		gp_add_entry(name,value);
	   }
   }

void ENDPAR(void)  /* free arglist & argbuf memory, & process STOP command */
   {
	if(ARGLIST != NULL) free(ARGLIST);
	if(ARGBUF  != NULL) free(ARGBUF);
	ARGBUF=  NULL;
	ARGLIST= NULL;
	if(FLAGS & STOP)
	   {
		fprintf(stderr,"%s[endpar]: stop due to STOP in input\n",
			PROGNAME);
		exit(GETPAR_STOP);
	   }
	FLAGS= END_PAR;	/* this stops further getpar calls */
   }

int
#ifdef FORTRAN
mstpar_(char *name, char *type, int *val, int dum1, int dum2, int lens)
/* dum1 & dum2 are extra args that fortran puts in */
#else
mstpar(char *name, char *type, void *val)
#endif
   {
	int cnt;
	char *typemess;
#ifndef	FORTRAN
	if( (cnt= GETPAR(name,type,val)) > 0) return(cnt);
#else
	if( (cnt= GETPAR(name,type,val,dum1,dum2,lens)) > 0) return(cnt);
#endif
	/* The following line corrects a common input error */
	if(type[1]=='v') { type[1]= type[0]; type[0]='v'; }

	switch(*type)
	   {
		case 'd': typemess= "an integer";	break;
		case 'f': typemess= "a float";		break;
		case 'F': typemess= "a double";		break;
		case 's': typemess= "a string";		break;
		case 'b': typemess= "a boolean";	break;
		case 'v': switch(type[1])
			   {
				case 'd': typemess= "an integer vector"; break;
				case 'f': typemess= "a float vector"; 	 break;
				case 'F': typemess= "a double vector";	 break;
				default : typemess= "unknown vector (error)";
					break;
			   }
			  break;
		default : typemess= "unknown (error)";	break;
	   }
	gp_getpar_err("mstpar","must specify value for '%s', expecting %s",
		name,typemess);
	return 0;
   }

int
#ifdef FORTRAN
getpar_(char *name, char *type, void* val, int dum1, int dum2, int lens)
/* dum1 & dum2 are extra args that fortran puts in */
#else
getpar(char *name, char *type, void *val)
#endif
   {
	register char *sptr;
	register struct arglist *alptr;
	double *dbl;
	float *flt;
	int *intgr;
	int h, hno, hyes, found;
#ifdef FORTRAN
	int termval;
#endif
	char line[MAXLINE], *str, *noname;
	if(FLAGS & END_PAR)
		gp_getpar_err("getpar","called after endpar");
	if( (FLAGS & INIT) == 0)
		gp_getpar_err("getpar","not initialized with setpar");
	if(FLAGS & VERBOSE)
		fprintf(stderr,"getpar: looking for %s\n",name);

	/* The following line corrects a common input error */
	if(type[1]=='v') { type[1]= type[0]; type[0]='v'; }

	found=0;

	if(NLIST <= 0) return(found);

	/*  
	 *  if list is NULL then return NOW; the following "for" loop
	 *  the pointer ARGLIST + (NLIST-1) gives ARGLIST -1 which
	 *  does not test to < ARGHEAD for some reason
         */

	if(*type == 'b') goto boolean;

	h= gp_compute_hash(name);

	/* search list backwards, stopping at first find */
	for(alptr= ARGLIST +(NLIST-1); alptr >= ARGHEAD; alptr--)
	   {
		if(alptr->hash != h) continue;
		if(strcmp(ARGBUF+alptr->argname_offset,name) != 0) continue;
		str= ARGBUF + alptr->argval_offset;
		switch(*type)
		   {
			case 'd':
			        intgr = (int *) val;
				*intgr= atoi(str);
				found=1;
				break;
			case 'f':
				flt= (float *) val;
				*flt= atof(str);
				found=1;
				break;
			case 'F':
				dbl= (double *) val;
				*dbl= atof(str);
				found=1;
				break;
			case 's':
                                sptr= (char *) val;
                                while(*str) *sptr++ = *str++;
                                *sptr= '\0';
#ifdef FORTRAN

/*  If we are in fortran, overwrite the null terminator and pad out 
 *  with blanks ;  if this is called internally from setpar_ lens
 *  will be 0 and the \0 will be left in place
 */

				while(sptr < (((char *)val)+lens)) 
					*sptr++ = ' ';
#endif
                                found=1;
                                break;
			case 'v':
				found= gp_getvector(str,type,val);
				break;
			default:
				gp_getpar_err("getpar",
					"unknown conversion type %s",type);
				break;
		   }
		break;
	   }
	goto list;

boolean:

	noname= line;
	sprintf(noname,"no%s",name);
	hno = gp_compute_hash(noname);
	hyes= gp_compute_hash(  name);
	found=0;
	/* search list backwards, stopping at first find */
	for(alptr= ARGLIST +(NLIST-1); alptr >= ARGHEAD; alptr--)
	   {
		if(alptr->hash != hno && alptr->hash != hyes) continue;
		if(strcmp(ARGBUF+alptr->argname_offset,  name)== 0)
		   {
			if(*(ARGBUF+alptr->argval_offset) == '\0')
			  *( (int *) val)= 1;
			else
			  *( (int *) val)= (int)atol(ARGBUF+alptr->argval_offset);
			found++;
			break;
		   }
		if(strcmp(ARGBUF+alptr->argname_offset,noname)== 0)
		   {	*( (int *) val)= 0; found++; break; }
	   }
   list:
	if(FLAGS & LIST)
	   {
		switch(*type)
		   {
			case 'd': sprintf(line,"(int) = %d",*( (int *) val));
				break;
			case 'f': flt= (float *)val;
				  sprintf(line,"(flt) = %14.6e",*flt); break;
			case 'F': dbl= (double *)val;
				  sprintf(line,"(dbl) = %14.6e",*dbl); break;
#ifdef FORTRAN
                        case 's':
				strcpy(line,"(str) = ");
				termval = (lens > MAXLINE-1-8) ? MAXLINE-1-8 : lens;
				strncpy(line + 8, (char *)val, termval);
				line[termval+8] = '\0';
				
				break;
#else
                        case 's': sprintf(line,"(str) = %s", (char *) val);
			        break;
#endif
			case 'b': sprintf(line,"(boo) = %d",*( (int *) val));
				break;
			case 'v': switch(type[1])
				   {
					/* should list these out */
					case 'd': sprintf(line,"(int vec)");
						break;
					case 'f': sprintf(line,"(flt vec)");
						break;
					case 'F': sprintf(line,"(dbl vec)");
						break;
					default : sprintf(line," vec type error");
						break;
				   }
				  break;
			default : sprintf(line," type error"); break;
		   }
		fprintf(LISTFILE,"%16s (%s) %s \n",name,
			(found ? "set":"def"),line);
	   }
	return(found);
   }

FILE *gp_create_dump(char *fname, char *filetype)
   {
	FILE *temp;

	if(*fname == '\0') return(stderr);
	if(strcmp(fname,"stderr") == 0) return(stderr);
	if(strcmp(fname,"stdout") == 0) return(stdout);
	if( (temp= fopen(fname,"w")) != NULL) return(temp);
	fprintf(stderr,"%s[setpar]: cannot create %s file %s\n",
		PROGNAME,filetype,fname);
	return(stderr);
   }

void gp_close_dump(FILE *file)
   {
	if(file == stderr || file == stdout) return;
	fclose(file);
   }

int gp_compute_hash(char *s)
   {
	register int h;
	h= s[0];
	if(s[1]) h |= (s[1])<<8;	else return(h);
	if(s[2]) h |= (s[2])<<16;	else return(h);
	if(s[3]) h |= (s[3])<<24;
	return(h);
   }

void gp_do_par_file(char *fname, int level)
   {
	register char *pl, *pn, *pv;
	char t, line[MAXLINE], name[MAXNAME], value[MAXVALUE];
	FILE *file;
	char *apl, *apv;
 
	if(level > MAXPARLEVEL)
		gp_getpar_err("setpar","%d (too many) recursive par file",level);
		
	if(*fname == '\0') return;
	
	if( (file=fopen(fname,"r"))==NULL)
		gp_getpar_err("setpar","cannot open par file %s",fname);

	while( gp_fgets(line,MAXLINE,file) != NULL )
	   {
		pl= line;
		/* loop over entries on each line */
	loop:	while(*pl==' ' || *pl=='\t') pl++;
		if(*pl=='\0'|| *pl=='\n') continue;
		if(*pl=='#') continue; /* comments on rest of line */

		/* get name */
		pn= name;
		while(*pl != '=' && *pl != '\0' && *pl != ' '
			&& *pl != '\n'		/* FIX by Glenn Nelson */
			&& *pl != '\t') *pn++ = *pl++;
		*pn = '\0';

		if(*pl == '=') pl++;

		/* get value */

		*value= '\0';
		pv= value;
		
		if(*pl == '"' || *pl == '\'')
		{
			t= *pl++;
			while(*pl != '\0' && *pl != '\n')
			{
				if(*pl == t)
				{
					if(pl[-1] != '\\') 
					{
						pl++;
						break;
					}
					pv[-1]= t;
					pl++;
				}
				else
				{
					if(*pl == '$')
					{
						apl = pl;
						apv = pv;
						gp_subpar(&apl, &apv);
						pl = apl;
						pv = apv;
					}
					else *pv++ = *pl++;
				}
			}
		}
		else   
		{
			while(*pl != '\0' && *pl != '\n'
			      && *pl != '\t' && *pl != ' ') 
				if(*pl == '$')
				{
					apl = pl;
					apv = pv;
					gp_subpar(&apl, &apv);
					pl = apl;
					pv = apv;
				}
				else *pv++ = *pl++;
		}
		*pv= '\0';

		gp_add_entry(name,value);
		if(strcmp("par",name) == 0)
			gp_do_par_file(value,level+1);
		goto loop;
	   }
	fclose(file);
   }

void gp_getpar_err(char *subname, char *format, ...)
   {
        va_list ap;
	va_start(ap, format);
	(void) fprintf(stderr,"\n***** ERROR in %s[%s] *****\n\t",
		       (PROGNAME == NULL ? "(unknown)" : PROGNAME),subname);
	(void) vfprintf(stderr, format, ap);
	va_end(ap);
	(void) fprintf(stderr,"\n");
	exit(GETPAR_ERROR);
   }
int gp_getvector(char *list, char *type, void *val)
   {
	register char *p;
	register int index, cnt;
	char *valptr;
	int limit;
	int ival, *iptr;
	float fval, *fptr;
	double dval, *dptr;

	limit= MAXVECTOR;
	if(type[2] == '(' || type[2] == '[') limit= (int)atol(&type[3]);
	if(limit <= 0)
		gp_getpar_err("getpar","bad limit=%d specified",limit);
	index= 0;
	p= list;
	while(*p != '\0'  && index < limit)
	   {
		cnt=1;
	 backup: /* return to here if we find a repetition factor */
		while(*p == ' ' || *p == '\t') p++;
		if(*p == '\0') return(index);
		valptr= p;
		while( *p != ',' && *p != '*' && *p != 'x' && *p != 'X' &&
			*p != '\0') p++;
		if(*p == '*' || *p == 'x' || *p == 'X')
		   {
			cnt= (int)atol(valptr);
			if(cnt <= 0)
				gp_getpar_err("getpar",
					"bad repetition factor=%d specified",
					 cnt);
			if(index+cnt > limit) cnt= limit - index;
			p++;
			goto backup;
		   }
		switch(type[1])
		   {
			case 'd':
				iptr= (int *) val;
				ival= (int)atol(valptr);
				while(cnt--) iptr[index++] = ival;
				break;
			case 'f':
				fptr= (float *) val;
				fval= atof(valptr);
				while(cnt--) fptr[index++] = fval;
				break;
			case 'F':
				dptr= (double *) val;
				dval= atof(valptr);
				while(cnt--) dptr[index++] = dval;
				break;
			default:
				gp_getpar_err("getpar",
					"bad vector type=%c specified",type[1]);
				break;
		   }
		if(*p != '\0') p++;
	   }
	return(index);
   }

/*  This allows parameter substitution */

void gp_subpar(char **apl, char **apv)
{
	register char *pl, *pv;
	char     subname[MAXNAME];
	int      valid = 0;
	char     *bpl, *bpv;
	
	if((*apl)[-1] == '\\')
	{
		(*apv)[-1] = (*apl)[0];
		(*apl)++;
		return;
	}

	pl = *apl;
	pl++;

	if(*pl == '(' )
	{
		pv=subname;
		pl++;
		while(*pl != ')' && *pl != '\0')
		{
			if(*pl == '$')
			{
				bpl = pl;
				bpv = pv;
				gp_subpar(&bpl, &bpv);
				pl = bpl;
				pv = bpv;
			}
			else *pv++ = *pl++;
		}
		*pv = '\0';
		if(*pl == ')' )
		{
			pl++;
			pv = *apv;
			valid = 1;
			ARGHEAD= ARGLIST;
#ifndef	FORTRAN
			if(GETPAR(subname, "s", pv))
#else
			if(GETPAR(subname, "s", pv, 0, 0, 0))
#endif
			{
				pv += strlen(pv);
			}
			*apv=pv;
			*apl=pl;
				
		}
	}
	if(valid)
	{
		return;
	}
	else
	{
		**apv = **apl;
		(*apv)++;
		(*apl)++;
	}
}

char *gp_fgets(char *line, int maxline, FILE *file)
{
	register char *p;
	register int i;
	int q;

        p = line;
	*p = '\0';

	if    ( (q = getc(file)) == EOF)     return((char *) NULL);
	else if( (unsigned char) q == '\n') return(line);

	*p = (unsigned char) q;
	p++;

	for (i=1; (q = getc(file)) != EOF && i < maxline;)
	{
		*p = (unsigned char) q;
		if(*p == '\0') return(line);
	        if(*p == '\n')
		{
			if( *(p-1) == '\\') 
			{
				p--;
				i--;
			}
			else break;
		}
		else 
		{
			p++;
			i++;
		}
	}
	*p = '\0';
	
	return line;
}
#ifndef FORTRAN
char *mstspar(char *parname)
{
	char parvalue[MAXVALUE];
	
	mstpar(parname, "s", parvalue);
	
	return(strdup(parvalue));
}
char *getspar(char *parname, char *defvalue)
{
	char parvalue[MAXVALUE];
	
	parvalue[0] = '\0';

	if(defvalue)
	{
		strcpy(parvalue, defvalue);
	}
	
	getpar(parname, "s", parvalue);

	if(! defvalue && parvalue[0] == '\0') return( (char *) 0 );
	
	else return(strdup(parvalue));

}
#endif
