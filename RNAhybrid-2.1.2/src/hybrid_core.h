/* Copyright (C) 2004 Marc Rehmsmeier, Peter Steffen, Matthias Hoechsmann */

/* This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version. */

/* This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. */

/* You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA */

#ifndef hybrid_core_h
#define hybrid_core_h


char *x;   /* input string */
int  m;    /* input length */
char *y;   /* input string */
int  n;    /* input length */


int iloop_upper_limit;
int bloop_upper_limit;


int helix_start, helix_end;


int gx;
char *t1, *t2, *t3, *t4;
char *r1, *r2, *r3, *r4;

int a1, a2, a3, a4, a5, a6;

void calc_unpaired_left_top(int i1, int j1, int i2, int j2);

void calc_unpaired_left_bot(int i1, int j1, int i2, int j2);

void calc_closed(int i1, int j1, int i2, int j2);

double calc_hybrid(int i1, int j1, int i2, int j2);

void mainloop(int bflag, int hit_number, int eflag, float energy_cutoff, int pflag, float pvalue_cutoff, int compact_output, char *target_ac, char *target_sq, char *query_ac, char *query_sq, int gflag, int plot_format, float xi, float theta);


#endif
