/*******************************************************************
*			sacio.c
* SAC I/O functions:
*	read_sachead	read SAC header
*	read_sac	read SAC binary data
*	read_sac2	read SAC data with cut option
*	write_sac	write SAC binary data
*	wrtsac2		write 2 1D arrays as XY SAC data
*	sachdr		creat new sac header
*	rdsac0_		fortran wraper for read_sac
*	wrtsac0_	fortran write 1D array as SAC binary data
*	wrtsac2_	fortran wraper for wrtsac2
*	wrtsac3_	wrtsac0 with component orientation cmpaz/cmpinc
*	swab4		reverse byte order for integer/float
*       trim            trim the space or tab in a string
*********************************************************************/


#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include "sac.h"

/***********************************************************

  read_sachead

  Description:	read binary SAC header from file.

  Author:	Lupei Zhu

  Arguments:	const char *name 	file name
		SACHEAD *hd		SAC header to be filled

  Return:	0 if success, -1 if failed

  Modify history:
	05/29/97	Lupei Zhu	Initial coding
************************************************************/

int	read_sachead(const char	*name,
		SACHEAD		*hd
	)
{
  FILE		*strm;

  if ((strm = fopen(name, "rb")) == NULL) {
     fprintf(stderr, "Unable to open %s\n",name);
     return -1;
  }

  if (fread(hd, sizeof(SACHEAD), 1, strm) != 1) {
     fprintf(stderr, "Error in reading SAC header %s\n",name);
     fclose(strm);
     return -1;
  }

#ifdef i386
  swab4((char *) hd, HD_SIZE);
#endif

  fclose(strm);
  return 0;

}


/***********************************************************

  read_sac

  Description:	read binary data from file. If succeed, it will return
		a float pointer to the read data array. The SAC header
		is also filled. A NULL pointer is returned if failed.

  Author:	Lupei Zhu

  Arguments:	const char *name 	file name
		SACHEAD *hd		SAC header to be filled

  Return:	float pointer to the data array, NULL if failed

  Modify history:
	09/20/93	Lupei Zhu	Initial coding
	12/05/96	Lupei Zhu	adding error handling
	12/06/96	Lupei Zhu	swap byte-order on PC
************************************************************/

float*	read_sac(const char	*name,
		SACHEAD		*hd
	)
{
  FILE		*strm;
  float		*ar;
  unsigned	sz;

  if ((strm = fopen(name, "rb")) == NULL) {
     fprintf(stderr, "Unable to open %s\n",name);
     return NULL;
  }

  if (fread(hd, sizeof(SACHEAD), 1, strm) != 1) {
     fprintf(stderr, "Error in reading SAC header %s\n",name);
     return NULL;
  }

#ifdef i386
  swab4((char *) hd, HD_SIZE);
#endif

  sz = hd->npts*sizeof(float);
  if ((ar = (float *) malloc(sz)) == NULL) {
     fprintf(stderr, "Error in allocating memory for reading %s\n",name);
     return NULL;
  }

  if (fread((char *) ar, sz, 1, strm) != 1) {
     fprintf(stderr, "Error in reading SAC data %s\n",name);
     return NULL;
  }

  fclose(strm);

#ifdef i386
  swab4((char *) ar, sz);
#endif

  return ar;

}



/***********************************************************

  write_sac

  Description:	write SAC binary data.

  Author:	Lupei Zhu

  Arguments:	const char *name 	file name
		SACHEAD hd		SAC header
		const float *ar		pointer to the data

  Return:	0 if succeed; -1 if failed

  Modify history:
	09/20/93	Lupei Zhu	Initial coding
	12/05/96	Lupei Zhu	adding error handling
	12/06/96	Lupei Zhu	swap byte-order on PC
************************************************************/

int	write_sac(const char	*name,
		SACHEAD		hd,
		const float	*ar
	)
{
  FILE		*strm;
  unsigned	sz;
  float		*data;
  int		error = 0;

  sz = hd.npts*sizeof(float);
  if (hd.iftype == IXY) sz *= 2;

  if ((data = (float *) malloc(sz)) == NULL) {
     fprintf(stderr, "Error in allocating memory for writing %s\n",name);
     error = 1;
  }

  if ( !error && memcpy(data, ar, sz) == NULL) {
     fprintf(stderr, "Error in copying data for writing %s\n",name);
     error = 1;
  }

#ifdef i386
  swab4((char *) data, sz);
  swab4((char *) &hd, HD_SIZE);
#endif

  if ( !error && (strm = fopen(name, "w")) == NULL ) {
     fprintf(stderr,"Error in opening file for writing %s\n",name);
     error = 1;
  }

  if ( !error && fwrite(&hd, sizeof(SACHEAD), 1, strm) != 1 ) {
     fprintf(stderr,"Error in writing SAC header for writing %s\n",name);
     error = 1;
  }

  if ( !error && fwrite(data, sz, 1, strm) != 1 ) {
     fprintf(stderr,"Error in writing SAC data for writing %s\n",name);
     error = 1;
  }

  free(data);
  fclose(strm);

  return (error==0) ? 0 : -1;

}


/*****************************************************

  swab4

  Description:	reverse byte order for float/integer

  Author:	Lupei Zhu

  Arguments:	char *pt	pointer to byte array
		int    n	number of bytes

  Return:	none

  Modify history:
	12/03/96	Lupei Zhu	Initial coding

************************************************************/

void	swab4(	char	*pt,
		int	n
	)
{
  int i;
  char temp;
  for(i=0;i<n;i+=4) {
    temp = pt[i+3];
    pt[i+3] = pt[i];
    pt[i] = temp;
    temp = pt[i+2];
    pt[i+2] = pt[i+1];
    pt[i+1] = temp;
  }
}

/***********************************************************

  sachdr

  Description:	creat a new SAC header

  Author:	Lupei Zhu

  Arguments:	float	dt		sampling interval
		int	ns		number of points
		float	b0		starting time

  Return:	SACHEAD

  Modify history:
	09/20/93	Lupei Zhu	Initial coding
************************************************************/

SACHEAD    sachdr( float dt, int ns, float b0)
{
  SACHEAD	hd = sac_null;
  hd.npts = ns;
  hd.delta = dt;
  hd.b = b0;
  hd.o = 0.;
  hd.e = b0+(ns-1)*hd.delta;
  hd.iztype = IO;
  hd.iftype = ITIME;
  hd.leven = TRUE;
  return hd;
}
 


/***********************************************************

  wrtsac2

  Description:	write 2 arrays into XY SAC data.

  Author:	Lupei Zhu

  Arguments:	const char *name	file name
		int	ns		number of points
		const float *x		x data array
		const float *y		y data array

  Return:	0 succeed, -1 fail

  Modify history:
	09/20/93	Lupei Zhu	Initial coding
	12/05/96	Lupei Zhu	adding error handling
	12/06/96	Lupei Zhu	swap byte-order on PC
************************************************************/

int	wrtsac2(const char	*name,
		int		n,
		const float	*x,
		const float	*y
	)
{
  SACHEAD	hd = sac_null;
  float		*ar;
  unsigned	sz;
  int		exit_code;

  hd.npts = n;
  hd.iftype = IXY;
  hd.leven = FALSE;

  sz = n*sizeof(float);

  if ( (ar = (float *) malloc(2*sz)) == NULL ) {
     fprintf(stderr, "error in allocating memory%s\n",name);
     return -1;
  }

  if (memcpy(ar, x, sz) == NULL) {
     fprintf(stderr, "error in copying data %s\n",name);
     free(ar);
     return -1;
  }
  if (memcpy(ar+sz, y, sz) == NULL) {
     fprintf(stderr, "error in copying data %s\n",name);
     free(ar);
     return -1;
  }

  exit_code = write_sac(name, hd, ar);
  
  free(ar);
  
  return exit_code;

}

/*for fortran--read evenly-spaced data */
void rdsac0_(const char *name2, float *dt, int *ns, float *b0, float *ar) {
   int i;
   SACHEAD hd;
   float *temp;
   char *name;

   name = trim(name2);
   temp = read_sac(name,&hd);
   *dt = hd.delta;
   *ns = hd.npts;
   *b0 = hd.b;
   for(i=0;i<*ns;i++) ar[i]=temp[i];
   free(temp);
   free(name);
}

/* for fortran--write evenly-spaced data */
void    wrtsac0_(const char *name2, float *dtp, int *nsp, float *b0p, 
float *distp, const float *ar) {
  SACHEAD hd;
  int ns, exit_code;
  float dt,b0,dist;
  char *name;

  ns = *nsp;
  dt = *dtp;
  b0 = *b0p;
  dist = *distp;
  name = trim(name2);
  hd = sachdr(dt,ns,b0);
  hd.dist = dist;
  exit_code=write_sac(name, hd, ar);
  free(name);
}

/* for fortran--write x-y data */
void    wrtsac2_(const char *name2, int n, const float *x, const float *y) {
  char *name;
  name = trim(name2);
  wrtsac2(name, n, x, y);
  free(name);
}

/* for fortran--write evenly-spaced data with comp orientation */
void    wrtsac3_(const char *name2, float dt, int ns, float b0, float dist, float cmpaz, float cmpinc, const float *ar) {
  SACHEAD hd;
  int exit_code;
  char *name;
  name = trim(name2);
  hd = sachdr(dt,ns,b0);
  hd.dist = dist;
  hd.cmpaz=cmpaz;
  hd.cmpinc=cmpinc;
  exit_code=write_sac(name, hd, ar);
  free(name);
}

/***********************************************************

  read_sac2

  Description:	read portion of data from file. If succeed, it will return
		a float pointer to the read data array. The SAC header
		is also filled. A NULL pointer is returned if failed.

  Author:	Lupei Zhu

  Arguments:	const char *name 	file name
		SACHEAD *hd		SAC header to be filled
		int	tmark,		time mark in sac header
		float	t1		begin time is tmark + t1
		int	npts		number of points

  Return:	float pointer to the data array, NULL if failed

  Modify history:
	11/08/00	Lupei Zhu	Initial coding
************************************************************/

float*	read_sac2(const char	*name,
		SACHEAD		*hd,
		int		tmark,
		float		t1,
		int		npts
	)
{
  FILE		*strm;
  int		i, nt1, nt2;
  float		*ar, *fpt;

  if ((strm = fopen(name, "rb")) == NULL) {
     fprintf(stderr, "Unable to open %s\n",name);
     return NULL;
  }

  if (fread(hd, sizeof(SACHEAD), 1, strm) != 1) {
     fprintf(stderr, "Error in reading SAC header %s\n",name);
     return NULL;
  }

#ifdef i386
  swab4((char *) hd, HD_SIZE);
#endif

  t1 += *( (float *) hd + 10 + tmark);
  nt1 = (int) rint((t1-hd->b)/hd->delta);
  nt2 = nt1+npts-1;
  if (nt1>=hd->npts-1 || nt2<=0) {
     fprintf(stderr,"data not in the specified window %s\n",name);
     return NULL;
  }
  if ((ar = (float *) malloc(npts*sizeof(float))) == NULL) {
     fprintf(stderr, "Error in allocating memory for reading %s\n",name);
     return NULL;
  }
  for(i=0;i<npts;++i) ar[i]=0.;

  if (nt1 > 0 ) {
     if (fseek(strm,nt1*sizeof(float),SEEK_CUR) < 0) {
	fprintf(stderr, "error in seek %s\n",name);
	return NULL;
     }
     fpt = ar;
  } else {
     fpt = ar-nt1;
     nt1 = 0;
  }
  if (nt2>hd->npts-1) nt2=hd->npts-1;
  i = nt2-nt1+1;
  if (fread((char *) fpt, sizeof(float), i, strm) != i) {
     fprintf(stderr, "Error in reading SAC data %s\n",name);
     return NULL;
  }
  fclose(strm);
#ifdef i386
  swab4((char *) fpt, i*sizeof(float));
#endif

  hd->npts = npts;
  hd->b = t1;
  hd->e = hd->b+npts*hd->delta;
  return ar;

}
/***********************************************************
  trim

  Description: Trim spaces and tabs from beginning and end

  Author:	Pengcheng Liu

  Arguments:	const char *name 	file name

  Return:	char pointer to the new string, NULL if failed

************************************************************/
char *trim(const char *name)
{
        int i,j, ns;
        char *s;
        
        ns= strlen (name);
        if(ns > 70) ns=70;
        s = (char *) malloc(ns*sizeof(char));
        s = strncpy(s, name,ns);

	i=0;
        while((i<ns-1)  && (s[i]<=' ')) {
                i++;
        }
        if(i>0) {
                for(j=0; j<ns-i; ++j) {
                        s[j]=s[j+i];
                }
        }

        ns=ns-i-1;
        for (j=0; j< ns; ++j) 
           if( s[j]<=' ')
               break;  
        s[j]='\0';
        return s;
}


void get_npts_(const char *name2, int *npt) {
    SACHEAD hd;
    int temp;
    char *name;
    name = trim(name2);
    temp = read_sachead(name,&hd);
    *npt=hd.npts;
    free(name);
}

void rsac1_(const char *name2, float *ar, int *ns, float *b0, float *dt, int *nl, int *ier) {
   int i;
   SACHEAD hd;
   float *temp;
   char *name;
   *ier = 1;
   name = trim(name2);
   temp = read_sac(name,&hd);
   *dt = hd.delta;
   *ns = hd.npts;
   *b0 = hd.b;
   if (*ns > *nl)  *ns = *nl;
   for(i=0;i<*ns;++i) ar[i]=temp[i];
   *ier = 0;
   free(temp);
   free(name);
}

/* for fortran--write evenly-spaced data */
void    wsac1_(const char *name2, float *ar, int *nsp, float *b0p, 
float *dtp,  int *ier) {
  SACHEAD hd;
  int ns, i, exit_code;
  float dt,b0,dist;
  char *name;

  *ier = 1;
  ns = *nsp;
  dt = *dtp;
  b0 = *b0p;
  dist = 1.0;
  name = trim(name2);   
  hd = sachdr(dt,ns,b0);
  hd.dist = dist;
  exit_code=write_sac(name, hd, ar); 
  *ier = 0;
  free(name);
}
