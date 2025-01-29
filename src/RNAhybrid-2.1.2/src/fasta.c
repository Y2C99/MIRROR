/* Copyright (C) 2004 Marc Rehmsmeier, Peter Steffen, Matthias Hoechsmann */

/* This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version. */

/* This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. */

/* You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA */

#include "fasta.h"


void nextAC(FILE *f, char *ac)
{
  fpos_t startloc;
  char s [MAXLINE];
  
  fgetpos(f,&startloc);
  while (fgets(s,MAXLINE,f)!=NULL && s[0]!='>');
  if (s==NULL) {
    fsetpos(f,&startloc);
    printf("Error in nextAC. Aborting\n");
    exit(2);
  }
  sscanf(s,">%s",ac);
}


void nextSQ(FILE *f, char *sq, int maxseqlength)
{
  fpos_t startloc;
  char s [MAXLINE];
  
  int offset = 0;
  fgetpos(f,&startloc);
  if (fgets(sq,MAXLINE,f)!=NULL) {
    while (sq[offset]!='\0')
      offset++;
    offset--; // remove newline
    if (offset > maxseqlength) {
      return;
    }
    fgetpos(f,&startloc);
    while (fgets(sq+offset,MAXLINE,f)!=NULL)
      if (sq[offset]=='>') {
	sq[offset]='\0';
	fsetpos(f,&startloc);
	break;
      }
      else {
	while (sq[offset]!='\0')
	  offset++;
	offset--; // remove newline
	if (offset > maxseqlength) {
	  return;
	}
	fgetpos(f,&startloc);
      }
    sq[offset]='\0';
    return;
  }
  fsetpos(f,&startloc);
  printf("Error in nextSQ. Aborting\n");
  exit(2);
}


int end(FILE *f)
{
  fpos_t startloc;
  char s [MAXLINE];
  
  fgetpos(f,&startloc);
  if (fgets(s,MAXLINE,f)==NULL)
    return 1;
  else {
    fsetpos(f,&startloc);
    return 0;
  }
}

