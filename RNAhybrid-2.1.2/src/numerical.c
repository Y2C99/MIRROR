/* Copyright (C) 2004 Marc Rehmsmeier, Peter Steffen, Matthias Hoechsmann */

/* This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version. */

/* This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. */

/* You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA */

#include "globals.h"
#include "minmax.h"
#include "math.h"


int compare_floats (const void *a, const void *b)
{
  const float *da = (const float *) a;
  const float *db = (const float *) b;

  return (*da > *db) - (*da < *db);
}


void empirical_cumulative_distribution_function(float **xcdf, float **ycdf, float *bin_number, float *xs, int sample_size)
{
  float min_e;
  float max_e;

  float bin_width, bin;

  int k, i;

  qsort(xs,sample_size,sizeof(float),compare_floats);

  (*bin_number) = min(sample_size,BINNUMBER);

  *xcdf = (float *) calloc(BINNUMBER,sizeof(float));
  *ycdf = (float *) calloc(BINNUMBER+1,sizeof(float));


  min_e = 0.0;
  max_e = xs[sample_size-1];

  bin_width = fabs(max_e - min_e) / ((float) (*bin_number));

  k = 0;
  (*ycdf)[0] = 0.0;

  bin = 0.0;

  i = 0;

  for (k=1; k<=(*bin_number); k++) {
    bin = bin + bin_width;
    (*xcdf)[k-1] = bin;
    (*ycdf)[k] = (*ycdf)[k-1];
    while (i<sample_size && xs[i] <= bin) {
      (*ycdf)[k]++;
      i++;
    }
  }

  for (k=1;k<=(*bin_number);k++)
    (*ycdf)[k] /= sample_size;
}


float mean(float *xs, int sample_size)
{
  float sum = 0.0;

  int k;
  for (k=0; k<sample_size; k++)
    sum += xs[k];

  return (sum/sample_size);
}



float variance(float *xs, int sample_size)
{
  float m = mean(xs,sample_size);

  float sum = 0.0;
  int k;
  for (k=0; k<sample_size; k++)
    sum += pow(m-xs[k],2);

  return (sum / (sample_size-1));
}

void linear_regression(float *slope, float *intercept, float *xs, float *ys, int sample_size)
{
  float var, s, sx, sy, sxx, sxy, delta;
  int k;

  var = variance(ys,sample_size);
  s = 0.0;
  for (k=0; k<sample_size; k++)
    s += 1.0/var;
  sx = 0.0;
  for (k=0; k<sample_size; k++)
    sx += xs[k]/var;
  sy = 0.0;
  for (k=0; k<sample_size; k++)
    sy += ys[k]/var;
  sxx = 0.0;
  for (k=0; k<sample_size; k++)
    sxx += xs[k]*xs[k]/var;
  sxy = 0.0;
  for (k=0; k<sample_size; k++)
    sxy += xs[k]*ys[k]/var;
  delta = s*sxx-sx*sx;
  (*slope) = (s*sxy-sx*sy)/delta;
  (*intercept) = (sxx*sy-sx*sxy)/delta;
}


void estimate_evd_parameters(int *used_sample_size, float *xi, float *theta, float *normalised_energies, int sample_size)
{
  float *xcdf, *ycdf;
  float bin_number;
  int start, end;
  
  /*   calculate empirical cdf: */
  empirical_cumulative_distribution_function(&xcdf,&ycdf,&bin_number,normalised_energies,sample_size);

  /*   delete values below normalised energy of FITLOWERCUTOFFABSOLUTE: */
  start=0;

  while (start<bin_number && xcdf[start] < FITLOWERCUTOFFABSOLUTE)
    start++;

  /*   delete top FITUPPERCUTOFFPERCENT values: */
  end = (int) bin_number * (100.0 - FITUPPERCUTOFFPERCENT) / 100.0 - 1;

  (*used_sample_size) = end-start+1;

  if (end<=start) {
    (*xi) = 0.0;
    (*theta) = 0.0;
  }
  else {
    /*   do linear regression on log(-log) transformed cdf: */
    float slope, intercept;
    int k;
    for (k=start; k<=end; k++)
      ycdf[k+1] = log(-log(ycdf[k+1]));

    linear_regression(&slope,&intercept,xcdf+start,ycdf+start+1,end-start+1);

    (*theta) = -1.0/slope;
    (*xi) = intercept * (*theta);
  }

  free(xcdf);
  free(ycdf);
}


#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY 1.0e-20
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define SIGN(a,b) ((b) > 0.0 ? fabs(a) : -fabs(a))
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

void mnbrak(float *ax, float *bx, float *cx, float *fa, float *fb, float *fc, float (*func)(float))
{
	float ulim,u,r,q,fu,dum;

	*fa=(*func)(*ax);
	*fb=(*func)(*bx);
	if (*fb > *fa) {
		SHFT(dum,*ax,*bx,dum)
		SHFT(dum,*fb,*fa,dum)
	}
	*cx=(*bx)+GOLD*(*bx-*ax);
	*fc=(*func)(*cx);
	while (*fb > *fc) {
		r=(*bx-*ax)*(*fb-*fc);
		q=(*bx-*cx)*(*fb-*fa);
		u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/
			(2.0*SIGN(MAX(fabs(q-r),TINY),q-r));
		ulim=(*bx)+GLIMIT*(*cx-*bx);
		if ((*bx-u)*(u-*cx) > 0.0) {
			fu=(*func)(u);
			if (fu < *fc) {
				*ax=(*bx);
				*bx=u;
				*fa=(*fb);
				*fb=fu;
				return;
			} else if (fu > *fb) {
				*cx=u;
				*fc=fu;
				return;
			}
			u=(*cx)+GOLD*(*cx-*bx);
			fu=(*func)(u);
		} else if ((*cx-u)*(u-ulim) > 0.0) {
			fu=(*func)(u);
			if (fu < *fc) {
				SHFT(*bx,*cx,u,*cx+GOLD*(*cx-*bx))
				SHFT(*fb,*fc,fu,(*func)(u))
			}
		} else if ((u-ulim)*(ulim-*cx) >= 0.0) {
			u=ulim;
			fu=(*func)(u);
		} else {
			u=(*cx)+GOLD*(*cx-*bx);
			fu=(*func)(u);
		}
		SHFT(*ax,*bx,*cx,u)
		SHFT(*fa,*fb,*fc,fu)
	}
}

#undef GOLD
#undef GLIMIT
#undef TINY
#undef MAX
#undef SIGN
#undef SHFT


#define ITMAX 100
#define CGOLD 0.3819660
#define ZEPS 1.0e-10
#define SIGN(a,b) ((b) > 0.0 ? fabs(a) : -fabs(a))
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

float brent(float ax, float bx, float cx, float (*f)(float), float tol, float *xmin)
{
	int iter;
	float a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
	float e=0.0;

	a=((ax < cx) ? ax : cx);
	b=((ax > cx) ? ax : cx);
	x=w=v=bx;
	fw=fv=fx=(*f)(x);
	for (iter=1;iter<=ITMAX;iter++) {
		xm=0.5*(a+b);
		tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
		if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
			*xmin=x;
			return fx;
		}
		if (fabs(e) > tol1) {
			r=(x-w)*(fx-fv);
			q=(x-v)*(fx-fw);
			p=(x-v)*q-(x-w)*r;
			q=2.0*(q-r);
			if (q > 0.0) p = -p;
			q=fabs(q);
			etemp=e;
			e=d;
			if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
				d=CGOLD*(e=(x >= xm ? a-x : b-x));
			else {
				d=p/q;
				u=x+d;
				if (u-a < tol2 || b-u < tol2)
					d=SIGN(tol1,xm-x);
			}
		} else {
			d=CGOLD*(e=(x >= xm ? a-x : b-x));
		}
		u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
		fu=(*f)(u);
		if (fu <= fx) {
			if (u >= x) a=x; else b=x;
			SHFT(v,w,x,u)
			SHFT(fv,fw,fx,fu)
		} else {
			if (u < x) a=u; else b=u;
			if (fu <= fw || w == x) {
				v=w;
				w=u;
				fv=fw;
				fw=fu;
			} else if (fu <= fv || v == x || v == w) {
				v=u;
				fv=fu;
			}
		}
	}
	printf("Too many iterations in BRENT\n");
	*xmin=x;
	return fx;
}

#undef ITMAX
#undef CGOLD
#undef ZEPS
#undef SIGN

