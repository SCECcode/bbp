/* 
 * NAME
 *	getarg
 *
 * DESCRIPTION
 *	get subroutine arguments from a string.  Acquired from Caltech.
 *
 *	See getarg.3 for details.
 *
 * AUTHOR
 * copyright (c) Robert W. Clayton
 *		 Seismological Laboratory
 *		 Caltech
 *		 Pasadena, CA 91125
 *
 * 24 Jul 1991  J Given         Added function getsarg() for dynamic
 *                              string allocation
 *
 * 25 Apr 1991	Cynde K. Smith	Added vector of string capability for C
 *			       	version only.
 * 26 Apr 1991	Cynde K. Smith  Replaced calls to atof and atoi with
 *			       	strtod and strtol to allow error checking.
 * 09 May 1991  Cynde K. Smith  Added countarg() and cntarg_() subroutines.
 * 14 May 1991  Cynde K. Smith  Added lenarg() subroutine.
 * 
 * Getarg routines:
 *
 * Externally visable routines:
 *
 *		setarg(argc,argv)
 *		countarg(name,type)
 *		lenarg(name)
 *		getarg(name,type,valptr)
 *		endarg()
 *
 * To get C-version:
 *		cc -c getarg.c
 *
 * To get F77-version:
 *		cp getarg.c fgetarg.c
 *		cc -c -DFORTRAN fgetarg.c
 *		rm fgetarg.c
 *
 *SccsId: @(#)getarg.c	56.1 10/25/93
 */
#include	<stdio.h>
#include	<stdlib.h>
#include	<stdarg.h>
/*
#ifndef hpux 
#include	<floatingpoint.h>
#endif
*/
#include	<string.h>
#include	"libget.h"

#define MAXNAME		256	/* max length of name */
#define MAXVALUE	5120	/* max length of value */
#define MAXVECTOR	20	/* max # of elements for unspecified vectors */
#define GETARG_ERROR	-1	/* error status for getarg error */

#define INIT	 1	/* bits for FLAGS (ext_arg.argflags) */
#define END_PAR	 2

#define LISTINC		32	/* increment size for arglist */
#define BUFINC		1024	/* increment size for argbuf */

struct arglist		/* structure of list set up by setarg */
   {
	char *argname;
	char *argval;
	int hash;
   };
struct ext_arg		/* global variables for getarg */
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
   }	ext_arg;

/* abbreviations: */
#define AL 		struct arglist
#define PROGNAME	ext_arg.progname
#define FLAGS		ext_arg.argflags
#define ARGLIST		ext_arg.arglist
#define ARGHEAD		ext_arg.arghead
#define ARGBUF		ext_arg.argbuf
#define NLIST		ext_arg.nlist
#define NBUF		ext_arg.nbuf
#define LISTMAX		ext_arg.listmax
#define BUFMAX		ext_arg.bufmax

static int ga_add_entry(char *name, char *value);
static int ga_getarg_err(char *subname, char *format, ...);
static int ga_compute_hash(char *s);
static int ga_getvector(char *list, char *type, void *val);

void
#ifdef FORTRAN
setarg_(char *list, char *subname, int dum1, int dum2)
#else
setarg(char *list, char *subname)
#endif
   {
	register char *pl, *pn, *pv;
	register int i;
	int      lenlist;
	char t, name[MAXNAME], value[MAXVALUE];
	

	PROGNAME= subname;
	FLAGS= INIT;

	ARGLIST = NULL;
	ARGBUF = NULL;
	NLIST= NBUF= LISTMAX= BUFMAX= 0;

#ifdef 	FORTRAN
	lenlist = dum1;
#else
	lenlist = strlen(list);
#endif	
	if(lenlist == 0) return;
	if(list == NULL) return;

	pl= list;
	/* loop over entries on each line */

	for(i=0; i<lenlist && *pl != '\0'; i++)
	{
		while(*pl==' ' || *pl=='\t') pl++;
		if(*pl=='\0'|| *pl=='\n') continue;
		
		/* get name */
		pn= name;
		while(*pl != '=' && *pl != '\0' && *pl != ' '
		      && *pl != '\t') *pn++ = *pl++;
		*pn = '\0';
		if(*pl == '=') pl++;


		/* get value */

		*value= '\0';
		pv= value;
		if(*pl == '\"' || *pl == '\'')
		{
			t = *pl++;
			while(*pl != '\0')
			{
				if(*pl == t)
				{
					   if (pl[-1] != '\\' &&
					       (pl[1] == ' ' || pl[1] == '\0'))
					     {
					       pl++;
					       break;
					     }
				}
				*pv++ = *pl++;
			}
		}
		else	
		{
			while(*pl && *pl != ' ' && *pl != '\t') 
				*pv++ = *pl++;
		}
		
		*pv= '\0';
		if (ga_add_entry(name,value) == GETARG_ERROR)
		    return;
	}
}

/* add an entry to arglist, expanding memory if necessary */
int ga_add_entry(char *name, char *value)
   {
	struct arglist *alptr;
	int len;
	register char *ptr;

	/*fprintf(stderr,"getarg: adding %s (%s)\n",name,value);*/
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
	if(NBUF+len >= BUFMAX)
	   {
		BUFMAX += BUFINC;
		if(ARGBUF == NULL)
			ARGBUF= (char *)malloc(BUFMAX);
		 else	ARGBUF= (char *)realloc(ARGBUF,BUFMAX);
	   }
	if(ARGBUF == NULL || ARGLIST == NULL)
		return ga_getarg_err("setarg","cannot allocate memory");

	/* add name */
	alptr= ARGLIST + NLIST;
	alptr->hash= ga_compute_hash(name);
	ptr= alptr->argname= ARGBUF + NBUF;
	do *ptr++ = *name; while(*name++);

	/* add value */
	NBUF += len;
	alptr->argval= ptr;
	do *ptr++ = *value; while(*value++);
	NLIST++;
	return 0;
   }

void 
#ifdef FORTRAN
endarg_()
#else
endarg(void) /* free arglist & argbuf memory, & process STOP command */
#endif
   {
	if(ARGLIST != NULL) free(ARGLIST);
	if(ARGBUF  != NULL) free(ARGBUF);
	ARGBUF=  NULL;
	ARGLIST= NULL;
	FLAGS= END_PAR;	/* this stops further getarg calls */
   }


/* count the number of arguments for a particular parameter name */

int
#ifdef FORTRAN
cntarg_(char *name, char *type, int dum1, int dum2)
#else
countarg(char *name, char *type)
#endif
{
	int 	found, h;
	char 	*str, *ptr;
	register struct arglist *alptr;
	
	if (FLAGS & END_PAR)
	{
		return ga_getarg_err("countarg","called after endarg");
	}
	
	if ((FLAGS & INIT) == 0)
	{
		return ga_getarg_err("countarg","not initialized with setarg");
	}
	

	if (NLIST == 0 || ARGLIST == NULL) 
	{
		return (0);
	}
	

	found=0;

	h = ga_compute_hash(name);

	/*  
	 *  if list is NULL then return NOW; the following "for" loop
	 *  the pointer ARGLIST + (NLIST-1) gives ARGLIST -1 which
	 *  does not test to < ARGHEAD for some reason
         */
	
	if (NLIST <= 0)
	{ 
		return (0);
	}
	

	/* search list backwards, stopping at first find */
	for (alptr = ARGLIST +(NLIST-1); alptr >= ARGLIST; alptr--)
	{
		if (alptr->hash != h) 
			continue;
		if (strcmp(alptr->argname, name) != 0) 
			continue;
		str = alptr->argval;
		ptr = str;
		
		switch (*type)
		{
	       	case 'd':
		case 'f':
		case 'F':
			/*
			 * Count the number of commas to detemine list
			 * size.  If the str isn't NULL than there
			 * is at least one element and at least one
			 * more than there are commas.
			 */
			while ((ptr = strchr(ptr, ',')) != NULL)
			{
					found++;
					ptr++;
			}
			if (str != NULL)
				found++;
			break;
	       	case 's':
			while ((ptr = strchr(ptr, ',')) != NULL)
			{
				if (ptr[-1] != '\\')
					found++;
				ptr++;
			}
			if (str != NULL)
				found++;
			break;
	       	default:
		       	return ga_getarg_err("countarg",
		       		      "unknown conversion type %s",type);
		       	break;
		}
		break;
	}
	return (found);
}	

int lenarg(char *name)
{
	int 	h, len, new_len;
	char 	*str, *ptr1, *ptr2;
	register struct arglist *alptr;
	
	if (FLAGS & END_PAR)
	{
		return ga_getarg_err("lenarg","called after endarg");
	}
	
	if ((FLAGS & INIT) == 0)
	{
		return ga_getarg_err("lenarg","not initialized with setarg");
	}

	if (NLIST == 0 || ARGLIST == NULL) 
	{
		return (0);
	}
	

	h = ga_compute_hash(name);

	/*  
	 *  if list is NULL then return NOW; in the following "for"
	 *  loop the pointer ARGLIST + (NLIST-1) gives ARGLIST -1
	 *  which does not test to < ARGHEAD for some reason
         */
	
	if (NLIST <= 0)
	{ 
		return (0);
	}
	

	/* search list backwards, stopping at first find */
	for (alptr = ARGLIST +(NLIST-1); alptr >= ARGLIST; alptr--)
	{
		if (alptr->hash != h) 
			continue;
		if (strcmp(alptr->argname, name) != 0) 
			continue;
		str = alptr->argval;
		ptr1 = str;
		len = new_len = 0;
	       	while ((ptr2 = strchr(ptr1, ',')) != NULL)
	       	{
	       		new_len += ptr2 - ptr1;
			if (ptr2[-1] != '\\')
			{
				if (new_len > len)
					len = new_len;
				new_len = 0;
			}
			ptr2++;
			ptr1 = ptr2;
		}
			
		if (ptr1)
			new_len = strlen(ptr1);
		if (new_len > len)
			len = new_len;
		
		break;
	}
	return (len);
}	



int
#ifdef FORTRAN
fgtarg_(char *name, char *type, void *val, int dum1, int dum2, int lens)
/* dum1 & dum2 are extra args that fortran puts in */
#else
getarg(char *name, char *type, void *val)
#endif
   {
	register char *sptr;
	register struct arglist *alptr;
	double *dbl;
	float *flt;
	int *intr;
	int h, hno, hyes, found;
	char noname[MAXNAME+2], *str, *ptr;

	/*fprintf(stderr,"looking for %s, type=%s\n",name,type);*/
	if(FLAGS & END_PAR)
		return ga_getarg_err("getarg","called after endarg");
	if( (FLAGS & INIT) == 0)
		return ga_getarg_err("getarg","not initialized with setarg");
	if (val == NULL)
		return ga_getarg_err("getarg", "NULL pointer value");
	
	if(NLIST == 0 || ARGLIST == NULL) return(0);

	/* The following line corrects a common input error */
	if(type[1]=='v') { type[1]= type[0]; type[0]='v'; }

	found=0;

	if(*type == 'b') goto boolean;

	h= ga_compute_hash(name);

	/*  
	 *  if list is NULL then return NOW; the following "for" loop
	 *  the pointer ARGLIST + (NLIST-1) gives ARGLIST -1 which
	 *  does not test to < ARGHEAD for some reason
         */
	
	if(NLIST <= 0) return(found);


	/* search list backwards, stopping at first find */
	for(alptr= ARGLIST +(NLIST-1); alptr >= ARGLIST; alptr--)
	   {
		/*fprintf(stderr,"getarg: checking %s\n",alptr->argname);*/
		if(alptr->hash != h) continue;
		if(strcmp(alptr->argname,name) != 0) continue;
		str= alptr->argval;
		switch(*type)
		   {
			case 'd':
			        intr= (int *) val;
				*intr= (int) strtol(str, &ptr, 0);
				if (ptr == str)
					found = 0;
				else
					found = 1;
				break;
			case 'f':
				flt= (float *) val;
				*flt= (float) strtod(str, &ptr);
				if (ptr == str)
					found = 0;
				else
					found = 1;
				break;
			case 'F':
				dbl= (double *) val;
				*dbl= strtod(str, &ptr);
				if (ptr == str)
					found = 0;
				else
					found = 1;
				break;
			case 's':
                                sptr= (char *) val;
                                while(*str) *sptr++ = *str++;
#ifdef FORTRAN
				while(sptr < (((char *)val)+lens)) 
					*sptr++ = ' ';
#else
                                *sptr= '\0';
#endif
                                found=1;
                                break;
			case 'v':
/* 
 *  For now, it's an error to get a Fortran vector of strings,
 *  If ga_getvector returns an error, stop processing and
 *  return GETARG_ERROR (-1)
 */

#ifdef FORTRAN
				if (type[1] == 's')
					return ga_getarg_err("getarg",
							"fortran string vector");
				else
					if ((found = ga_getvector(str, type, val))
					    == GETARG_ERROR)
						return found;
#else
				if ((found = ga_getvector(str,type,val)) ==
				    GETARG_ERROR)
					return found;
#endif
				break;
			default:
				return ga_getarg_err("getarg",
					"unknown conversion type %s",type);
				break;
		   }
		break;
	   }
	return(found);
boolean:

	/*  
	 *  if list is NULL then return NOW; the following "for" loop
	 *  the pointer ARGLIST + (NLIST-1) gives ARGLIST -1 which
	 *  does not test to < ARGHEAD for some reason
         */
	
	if(NLIST <= 0) return(found);

	sprintf(noname,"no%s",name);
	hno = ga_compute_hash(noname);
	hyes= ga_compute_hash(  name);
	found=0;
	/* search list backwards, stopping at first find */
	for(alptr= ARGLIST +(NLIST-1); alptr >= ARGLIST; alptr--)
	   {
		if(alptr->hash != hno && alptr->hash != hyes) continue;
		if(strcmp(alptr->argname,  name)== 0)
		   {
			if(alptr->argval[0] == '\0')
			  *( (int *) val)= 1;
			else
			  *( (int *) val)= atol(alptr->argval);
			found++;
			break;
		   }
		if(strcmp(alptr->argname,noname)== 0)
		   {	*( (int *) val)= 0; found++; break; }
	   }
	return(found);
   }

int ga_compute_hash(char *s)
   {
	register int h;
	h= s[0];
	if(s[1]) h |= (s[1])<<8;	else return(h);
	if(s[2]) h |= (s[2])<<16;	else return(h);
	if(s[3]) h |= (s[3])<<24;
	return(h);
   }

int ga_getvector(char *list, char *type, void *val)
   {
	register char *p;
	register int index, cnt;
	char *valptr, **strptr, *svalptr, sval[MAXVALUE], *ptr;
	int limit;
	int ival, *iptr;
	float fval, *fptr;
	double dval, *dptr;

	limit= MAXVECTOR;
	if(type[2] == '(' || type[2] == '[') limit= atol(&type[3]);
	if(limit <= 0)
		return ga_getarg_err("getarg","bad limit=%d specified",limit);
	/*fprintf(stderr,"limit=%d\n",limit);*/
	index= 0;
	p= list;
	while(*p != '\0'  && index < limit)
	   {
		cnt=1;
	 backup: /* return to here if we find a repetition factor */
		while(*p == ' ' || *p == '\t') p++;
		if(*p == '\0') return(index);
		valptr= p;
	getvalue: /* return here if valid value in char*[] arg */
		while( *p != ',' && *p != '*' && *p != 'x' && *p != 'X' &&
			*p != '\0') p++;
		if (type[1] == 's' && ((*p == ',' && p[-1] == '\\')
				       || *p == '*' || *p == 'x' || *p == 'X'))
		{	p++; goto getvalue; }
		if((*p == '*' || *p == 'x' || *p == 'X') && 
		   type[1] != 's')
		   {
			cnt= atol(valptr);
			if(cnt <= 0)
				return ga_getarg_err("getarg",
					"bad repetition factor=%d specified",
					 cnt);
			if(index+cnt > limit) cnt= limit - index;
			p++;
			goto backup;
		   }
		/*fprintf(stderr,"index=%d cnt=%d p=%s$\n",index,cnt,p);*/
		switch(type[1])
		   {
			case 'd':
				iptr= (int *) val;
				ival= (int) strtol(valptr, &ptr, 0);
				if (ptr == valptr)
					index = 0;
				else
				{
					if (iptr+index == NULL)
						return ga_getarg_err("getarg",
			       		          "NULL vector ptr at index %d.", index);
					while(cnt--)
					       	iptr[index++] = ival;
				}
				break;
			case 'f':
				fptr= (float *) val;
				fval= (float) strtod(valptr, &ptr);
				if (ptr == valptr)
					index = 0;
				else
				{
					if (fptr+index == NULL)
						return ga_getarg_err("getarg",
				       	          "NULL vector ptr at index %d", index);
					while(cnt--) 
						fptr[index++] = fval;
				}
				break;
			case 'F':
				dptr= (double *) val;
				dval= strtod(valptr, &ptr);
				if (ptr == valptr)
					index = 0;
				else
				{
					if (dptr+index == NULL)
						return ga_getarg_err("getarg",
				       	          "NULL vector ptr at index %d", index);
					while(cnt--) 
						dptr[index++] = dval;
				}
				break;
			case 's':
				svalptr = sval;
				strptr = (char **) val;
				*svalptr = '\0';
				while (*valptr != '\'' && *valptr != '\0')
				{
					if (*valptr == '\\')
						valptr++;
				       	*svalptr++ = *valptr++;
				}
                                *svalptr= '\0';
				if (strptr+index == NULL)
					return ga_getarg_err("getarg",
		       	       		  "NULL vector ptr at index %d", index);
				strcpy(strptr[index++], sval);
				if (*p != '\0') p++;
				break;
			default:
				return ga_getarg_err("getarg",
					"bad vector type=%c specified",type[1]);
				break;
		   }
		/* if conversion couldn't be made return 0 */
		if (index == 0) 
			break;
		if (*p != '\0') 
			p++;
	   }
	return(index);
   }

int
ga_getarg_err(char *subname, char *format, ...)
   {
        va_list ap;
	va_start(ap, format);
	(void) fprintf(stderr,"\n***** ERROR in %s[%s] *****\n\t",
		(PROGNAME == NULL ? "(unknown)" : PROGNAME),subname);
	(void) vfprintf(stderr, format, ap);
	va_end(ap);
	(void) fprintf(stderr,"\n");
	return GETARG_ERROR;
   }
#ifndef FORTRAN
char *getsarg(char *parname, char *defvalue)
{
	char parvalue[1024];

	parvalue[0] = '\0';

	if(defvalue)
	{
		strcpy(parvalue, defvalue);
	}
	
	getarg(parname, "s", parvalue);

	if(! defvalue && parvalue[0] == '\0') return( (char *) 0 );

	else return(strdup(parvalue));
}
#endif
