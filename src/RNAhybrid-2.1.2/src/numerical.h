/* Copyright (C) 2004 Marc Rehmsmeier, Peter Steffen, Matthias Hoechsmann */

/* This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version. */

/* This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. */

/* You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA */

#ifndef numerical_h
#define numerical_h


int compare_floats (const void *a, const void *b);

void empirical_cumulative_distribution_function(float **xcdf, float **ycdf, float *bin_number, float *xs, int sample_size);

float mean(float *xs, int sample_size);

float variance(float *xs, int sample_size);

void linear_regression(float *slope, float *intercept, float *xs, float *ys, int sample_size);

void estimate_evd_parameters(int *used_sample_size, float *xi, float *theta, float *normalised_energies, int sample_size);

void mnbrak(float *ax, float *bx, float *cx, float *fa, float *fb, float *fc, float (*func)(float));

float brent(float ax, float bx, float cx, float (*)(float), float tol, float *xmin);

#endif
