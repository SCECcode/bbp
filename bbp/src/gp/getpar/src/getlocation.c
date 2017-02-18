#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "libget.h"

#define LOC_FILE "/LOCATIONS"


static int getlin(FILE *fd, char *line)
{
  char c;
  int n;

  n=0;
  while( (c=getc(fd)) != EOF)
    {
      if(c == '\n')
	{
	  *line= '\0';
	  return(n);
	}
      n++;
      *line++ = c;
    }
  return(EOF);
}


int getlocation(char *keyname, char *location, int fatal)
{
  FILE *fd;
  int i, n;
  char line[128], name[128], loc[128];

  if( (fd= fopen(LOC_FILE,"r")) == NULL)
    {
      if(fatal)
	{
	  fprintf(stderr,"cannot open LOCATION file -fatal\n");
	  exit(-1);
	}
      else	return(-1);
    }
  while( (n=getlin(fd,line)) != EOF)
    {
      for(i=0; i<n; i++)
	{
	  if(line[i] == '#')	/* truncate line */
	    {
	      n= i;
	      line[n]= '\0';
	      break;
	    }
	}
      if(n == 0) continue;	/* blank line */
      n= sscanf(line,"%s %s",name,loc);
      if( n != 2 ) continue;	/* bad  syntax on line, ignore */
      if(strcmp(name,keyname) == 0)
	{
	  strcpy(location,loc);
	  fclose(fd);
	  return(1);
	}
    }
  if(fatal)
    {
      fprintf(stderr,"cannot find %s in LOCATION\n",keyname);
      fclose(fd);
      exit(-1);
    }
  fclose(fd);
  return(0);
}
