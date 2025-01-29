/* Copyright (C) 2004 Marc Rehmsmeier, Peter Steffen, Matthias Hoechsmann */

/* This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version. */

/* This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. */

/* You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA */

#ifndef energy_h
#define energy_h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "minmax.h"
#include "input.h"

double e;
double t;
double temp;
double r;

double mloop_close;
double free_base_penalty;
double helix_penalty;

double wkn;

double npp;
double pbp;

double il_asym_ar[16][16];

double stack_dg_ar  [ALPHASIZE][ALPHASIZE][ALPHASIZE][ALPHASIZE];
double tstackh_dg_ar[ALPHASIZE][ALPHASIZE][ALPHASIZE][ALPHASIZE];
double tstacki_dg_ar[ALPHASIZE][ALPHASIZE][ALPHASIZE][ALPHASIZE];
double hl_tetra_ar[ALPHASIZE][ALPHASIZE][ALPHASIZE][ALPHASIZE][ALPHASIZE][ALPHASIZE];

double hl_ent_ar[31];
double bl_ent_ar[31];
double il_ent_ar[31];

double dr_dangle_dg_ar[ALPHASIZE][ALPHASIZE][ALPHASIZE];
double dl_dangle_dg_ar[ALPHASIZE][ALPHASIZE][ALPHASIZE];

double int11_ar[ALPHASIZE][ALPHASIZE][ALPHASIZE][ALPHASIZE][ALPHASIZE][ALPHASIZE];
double int21_ar[ALPHASIZE][ALPHASIZE][ALPHASIZE][ALPHASIZE][ALPHASIZE][ALPHASIZE][ALPHASIZE];
double int22_ar[ALPHASIZE][ALPHASIZE][ALPHASIZE][ALPHASIZE][ALPHASIZE][ALPHASIZE][ALPHASIZE][ALPHASIZE];

char canPair[ALPHASIZE][ALPHASIZE];

void init_constants();

#define compl(I,J) canPair[I][J]

double sr_energy(int i, int j);

double dr_energy(int i, int j);

double dli_energy(int i, int j);

double dl_energy(int i, int j);

double dri_energy(int i, int j);
 
double dangles (int i, int j, int i2, int j2, int k, int l, int k2,int l2);

double sspenalty(int a);

double log_interp(int size);

double hl_ent(int size);

double bl_ent(int size);

#define il_ent(size) ((size)<=30 ? il_ent_ar[(size)] : il_ent_ar[30] + log_interp((size)))

int lengthOf(int i, int j);

#define il_asym(u,v) (il_asym_ar[u][v])

/* double il_asym (int sl, int sr); */

double top_stack(int lb, int rb);

double bot_stack(int lb, int rb);

double bl_stacking (int t, int b, int i, int j);

#define il_stack_open(i,j) tstacki_dg_ar[inpx(i)][inpx((i)+1)][inpy((j)+1)][inpy(j)]

#define il_stack_close(i,j) tstacki_dg_ar[inpy(j)][inpy((j)-1)][inpx((i)-1)][inpx(i)]

double int_special(int i, int j, int t, int b);

double do_il_special(int i, int j, int k, int l, int u, int v, double e);

/* double do_il(int i, int j, int k, int l, int u, int v, double e); */

#define do_il(i,j,k,l,u,v,e) ((e)+il_stack_close((l)+1,(v)+1)+il_ent((l)-(k)+(v)-(u))+il_asym((l)-(k),(v)-(u)))

void init_il_asym_ar();

void init_stack_dg_ar();

void init_tstackh_dg_ar();

void init_hl_tetra_ar();

void init_bl_ent_ar();

void init_il_ent_ar();

void init_tstacki_dg_ar();

void init_dr_dangle_dg_ar();

void init_dl_dangle_dg_ar();

void init_int11_ar();

void init_int21_ar();

void init_int22_ar();

void init_canPair();

void init_energies();




#endif
