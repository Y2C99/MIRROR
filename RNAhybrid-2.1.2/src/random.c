/* Copyright (C) 2004 Marc Rehmsmeier, Peter Steffen, Matthias Hoechsmann */

/* This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version. */

/* This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. */

/* You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA */

#include "random.h"
#include "globals.h"
#include "mt19937-1.h"
#include <sys/types.h>
#include <stdio.h>
#include <string.h>
#include <math.h>



char *alphabet = "ACGT";


void read_dinucleotide_frequencies(float **freq_di, FILE *f)
{
  char s [MAXLINE];
  int i,j;

  while (fgets(s,MAXLINE,f)!=NULL) {
    if (s[0] != '#') {
      switch (toupper(s[0])) {
      case 'A':
	i=0;
	break;
      case 'C':
	i=1;
	break;
      case 'G':
	i=2;
	break;
      case 'T':
	i=3;
	break;
      case 'U':
	i=3;
	break;
      }
      switch (toupper(s[1])) {
      case 'A':
	j=0;
	break;
      case 'C':
	j=1;
	break;
      case 'G':
	j=2;
	break;
      case 'T':
	j=3;
	break;
      case 'U':
	j=3;
	break;
      }
      sscanf(s,"%*s %f",&(freq_di[i][j]));
    }
  }
}


char* random_sequence(int seqlen, float *fdf, float **fdf_di)
{
  int i, c;
  float s;

  char *sequence = (char *) calloc(seqlen+1,sizeof(char));

  for (i=0; i<seqlen; i++) {
    s = genrand();
    c = 0;
    if (i==0) {
      while (fdf[c]<s) c++;
    }
    else {
      while (fdf_di[strchr(alphabet,sequence[i-1])-alphabet][c]<s) c++;
    }
    sequence[i]=alphabet[c];
  }
  sequence[seqlen] = '\0';

  return sequence;
}


float normal_random_number()
{
  float u1, u2, v1,v2,s,z;

  do {

    u1 = genrand();
    u2 = genrand();

    v1 = 2*u1-1;
    v2 = 2*u2-1;

    s = pow(v1,2) + pow(v2,2);

  } while (s > 1);

  z = sqrt(-2*log(s)/s)*v1;

  return z;
}
