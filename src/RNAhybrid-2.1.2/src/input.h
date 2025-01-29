/* Copyright (C) 2004 Marc Rehmsmeier, Peter Steffen, Matthias Hoechsmann */

/* This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version. */

/* This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. */

/* You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA */

#ifndef input_h
#define input_h


extern char *x;   /* input string */
extern int  m;    /* input length */
extern char *y;   /* input string */
extern int  n;    /* input length */


#define A 0
#define C 1
#define G 2
#define U 3
#define N 4
#define X 5 /* this is an additional letter for masking out hits */

#define ALPHASIZE 6


#define inpx(I) x[I]
#define inpy(I) y[I]


void convert_x();
void convert_y();

void remove_whitespace(char *sequence);

#endif
