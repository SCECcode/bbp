/*
 * Copyright 1992 Science Applications International Corporation.
 *
 * NAME
 *	getbpar()
 *      getdpar()
 *      getfpar()
 *      getffpar()
 * 
 * FILE 
 *	defpar.c
 *
 * SYNOPSIS
 *      par->duration = getffpar ("duration", 200.0);
 *	
 * DESCRIPTION
 *
 * DIAGNOSTICS
 *
 * FILES
 *
 * NOTES
 * 
 * SEE ALSO
 *
 * AUTHOR
 * 	Rick Jenkins  10/08/92
 *
 */

#include "libget.h"

int getbpar(char *parname, int defvalue)
{
	int parvalue = defvalue;
	
	getpar (parname, "b", &parvalue);

	return (parvalue);
}

int getdpar(char *parname, int defvalue)
{
	int parvalue = defvalue;
	
	getpar (parname, "d", &parvalue);

	return (parvalue);
}

float getfpar(char *parname, float defvalue)
{
	float parvalue = defvalue;
	
	getpar (parname, "f", &parvalue);

	return (parvalue);
}

double getffpar(char *parname, double defvalue)
{
	double parvalue = defvalue;
	
	getpar (parname, "F", &parvalue);

	return (parvalue);
}
