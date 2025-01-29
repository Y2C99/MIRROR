/* Copyright (C) 2004 Marc Rehmsmeier, Peter Steffen, Matthias Hoechsmann */

/* This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version. */

/* This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. */

/* You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA */

#include "input.h"
#include <string.h>


void convert_x(){
  int i;
  char c;
  for (i=0; i<=m; i++) {
    c=x[i];
    if      (c=='a' || c=='A') x[i]=A;
    else if (c=='c' || c=='C') x[i]=C;
    else if (c=='g' || c=='G') x[i]=G;
    else if (c=='u' || c=='U') x[i]=U;
    else if (c=='t' || c=='T') x[i]=U;
    else                       x[i]=N;
  }
}

void convert_y(){
  int i;
  char c;
  for (i=0; i<=n; i++) {
    c=y[i];
    if      (c=='a' || c=='A') y[i]=A;
    else if (c=='c' || c=='C') y[i]=C;
    else if (c=='g' || c=='G') y[i]=G;
    else if (c=='u' || c=='U') y[i]=U;
    else if (c=='t' || c=='T') y[i]=U;
    else                       y[i]=N;
  }
}


void remove_whitespace(char *sequence)
{
  int len = strlen(sequence);
  int i;

  for (i=0; i<len; i++)
    if (sequence[i]<33) {
      memmove(sequence+i,sequence+i+1,len-i-1);
      len--;
      i--;
    }

  sequence[len] = 0;
}

