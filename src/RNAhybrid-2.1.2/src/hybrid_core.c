/* Copyright (C) 2004 Marc Rehmsmeier, Peter Steffen, Matthias Hoechsmann */

/* This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version. */

/* This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. */

/* You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA */

/* compiled by the ADP compiler, version 0.8.342                                    */
/* source file: rnahybrid/HybridMFE3_simple_3tab.lhs                                */
/* command:                                                                         */
/* adpcompile -c rnahybrid/HybridMFE3_simple_3tab.lhs -al mfe enum -bt -bts -o rnahybrid/hybridBack4.c */
/* -------------------------------------------------------------------------------- */

#include "config.h"
#include "hybrid_core.h"
#include "plot.h"
#include "minmax.h"
#include "fasta.h"
#include "input.h"
#include "energy.h"
#include "globals.h"
#include <stdio.h>  
#include <string.h> 


static char const rcsid[] = "$Id: hybrid_core.c,v 1.6 2004/11/15 20:41:36 marc Exp $";


extern int iloop_upper_limit;
extern int bloop_upper_limit;

/* data structures                                                                  */
/* -------------------------------------------------------------------------------- */

struct str1 {
   double alg_mfe;
   struct str_Hybrid *alg_enum;
};

/* signature                                                                        */
/* -------------------------------------------------------------------------------- */

#define SIGID__NTID 1
#define SIGID_Ult 2
#define SIGID_Ulb 3
#define SIGID_Eds 4
#define SIGID_Edt 5
#define SIGID_Edb 6
#define SIGID_Sr 7
#define SIGID_Bt 8
#define SIGID_Bb 9
#define SIGID_Il 10
#define SIGID_El 11
#define SIGID_Nil 12


struct str_Hybrid {
   int utype;
   void *entry;
};

struct str_Hybrid *new_Hybrid(int u, void *entry)
{
   struct str_Hybrid *t;

   t=(struct str_Hybrid *) calloc(1, sizeof(struct str_Hybrid ));
   t->utype = u;
   t->entry = entry;
   return(t);
}

/* signature operators                                                              */
/* -------------------------------------------------------------------------------- */

/* operator _NTID                                                                   */
/* -------------------------------------------------------------------------------- */

struct str__NTID {
   struct str1 a1;
   struct str1 (*f1)(int , int , int , int );
   int i11, j11, i21, j21;
};

struct str_Hybrid *new__NTID(struct str1 (*f1)(int , int , int , int ), int i11, int j11, int i21, int j21)
{
   struct str__NTID *t;

   t=(struct str__NTID *) calloc(1, sizeof(struct str__NTID ));
   t->f1 = f1;
   t->i11 = i11;
   t->j11 = j11;
   t->i21 = i21;
   t->j21 = j21;
   return(new_Hybrid(SIGID__NTID, t));
}


/* operator Ult                                                                     */
/* -------------------------------------------------------------------------------- */

struct str_Ult {
   int a1;
   int a2;
   struct str1 a3;
   struct str1 (*f3)(int , int , int , int );
   int i13, j13, i23, j23;
};

struct str_Hybrid *new_Ult(int a1, int a2, struct str1 (*f3)(int , int , int , int ), int i13, int j13, int i23, int j23)
{
   struct str_Ult *t;

   t=(struct str_Ult *) calloc(1, sizeof(struct str_Ult ));
   t->a1 = a1;
   t->a2 = a2;
   t->f3 = f3;
   t->i13 = i13;
   t->j13 = j13;
   t->i23 = i23;
   t->j23 = j23;
   return(new_Hybrid(SIGID_Ult, t));
}


/* operator Ulb                                                                     */
/* -------------------------------------------------------------------------------- */

struct str_Ulb {
   int a1;
   int a2;
   struct str1 a3;
   struct str1 (*f3)(int , int , int , int );
   int i13, j13, i23, j23;
};

struct str_Hybrid *new_Ulb(int a1, int a2, struct str1 (*f3)(int , int , int , int ), int i13, int j13, int i23, int j23)
{
   struct str_Ulb *t;

   t=(struct str_Ulb *) calloc(1, sizeof(struct str_Ulb ));
   t->a1 = a1;
   t->a2 = a2;
   t->f3 = f3;
   t->i13 = i13;
   t->j13 = j13;
   t->i23 = i23;
   t->j23 = j23;
   return(new_Hybrid(SIGID_Ulb, t));
}


/* operator Eds                                                                     */
/* -------------------------------------------------------------------------------- */

struct str_Eds {
   int a1;
   int a2;
   struct str1 a3;
   struct str1 (*f3)(int , int , int , int );
   int i13, j13, i23, j23;
};

struct str_Hybrid *new_Eds(int a1, int a2, struct str1 (*f3)(int , int , int , int ), int i13, int j13, int i23, int j23)
{
   struct str_Eds *t;

   t=(struct str_Eds *) calloc(1, sizeof(struct str_Eds ));
   t->a1 = a1;
   t->a2 = a2;
   t->f3 = f3;
   t->i13 = i13;
   t->j13 = j13;
   t->i23 = i23;
   t->j23 = j23;
   return(new_Hybrid(SIGID_Eds, t));
}


/* operator Edt                                                                     */
/* -------------------------------------------------------------------------------- */

struct str_Edt {
   int a1;
   int a2;
   struct str1 a3;
   struct str1 (*f3)(int , int , int , int );
   int i13, j13, i23, j23;
};

struct str_Hybrid *new_Edt(int a1, int a2, struct str1 (*f3)(int , int , int , int ), int i13, int j13, int i23, int j23)
{
   struct str_Edt *t;

   t=(struct str_Edt *) calloc(1, sizeof(struct str_Edt ));
   t->a1 = a1;
   t->a2 = a2;
   t->f3 = f3;
   t->i13 = i13;
   t->j13 = j13;
   t->i23 = i23;
   t->j23 = j23;
   return(new_Hybrid(SIGID_Edt, t));
}


/* operator Edb                                                                     */
/* -------------------------------------------------------------------------------- */

struct str_Edb {
   int a1;
   int a2;
   struct str1 a3;
   struct str1 (*f3)(int , int , int , int );
   int i13, j13, i23, j23;
};

struct str_Hybrid *new_Edb(int a1, int a2, struct str1 (*f3)(int , int , int , int ), int i13, int j13, int i23, int j23)
{
   struct str_Edb *t;

   t=(struct str_Edb *) calloc(1, sizeof(struct str_Edb ));
   t->a1 = a1;
   t->a2 = a2;
   t->f3 = f3;
   t->i13 = i13;
   t->j13 = j13;
   t->i23 = i23;
   t->j23 = j23;
   return(new_Hybrid(SIGID_Edb, t));
}


/* operator Sr                                                                      */
/* -------------------------------------------------------------------------------- */

struct str_Sr {
   int a1;
   int a2;
   struct str1 a3;
   struct str1 (*f3)(int , int , int , int );
   int i13, j13, i23, j23;
};

struct str_Hybrid *new_Sr(int a1, int a2, struct str1 (*f3)(int , int , int , int ), int i13, int j13, int i23, int j23)
{
   struct str_Sr *t;

   t=(struct str_Sr *) calloc(1, sizeof(struct str_Sr ));
   t->a1 = a1;
   t->a2 = a2;
   t->f3 = f3;
   t->i13 = i13;
   t->j13 = j13;
   t->i23 = i23;
   t->j23 = j23;
   return(new_Hybrid(SIGID_Sr, t));
}


/* operator Bt                                                                      */
/* -------------------------------------------------------------------------------- */

struct str_Bt {
   int a1;
   int a2;
   int a3;
   int a4;
   int a5;
   struct str1 a6;
   struct str1 (*f6)(int , int , int , int );
   int i16, j16, i26, j26;
};

struct str_Hybrid *new_Bt(int a1, int a2, int a3, int a4, int a5, struct str1 (*f6)(int , int , int , int ), int i16, int j16, int i26, int j26)
{
   struct str_Bt *t;

   t=(struct str_Bt *) calloc(1, sizeof(struct str_Bt ));
   t->a1 = a1;
   t->a2 = a2;
   t->a3 = a3;
   t->a4 = a4;
   t->a5 = a5;
   t->f6 = f6;
   t->i16 = i16;
   t->j16 = j16;
   t->i26 = i26;
   t->j26 = j26;
   return(new_Hybrid(SIGID_Bt, t));
}


/* operator Bb                                                                      */
/* -------------------------------------------------------------------------------- */

struct str_Bb {
   int a1;
   int a2;
   int a3;
   int a4;
   int a5;
   struct str1 a6;
   struct str1 (*f6)(int , int , int , int );
   int i16, j16, i26, j26;
};

struct str_Hybrid *new_Bb(int a1, int a2, int a3, int a4, int a5, struct str1 (*f6)(int , int , int , int ), int i16, int j16, int i26, int j26)
{
   struct str_Bb *t;

   t=(struct str_Bb *) calloc(1, sizeof(struct str_Bb ));
   t->a1 = a1;
   t->a2 = a2;
   t->a3 = a3;
   t->a4 = a4;
   t->a5 = a5;
   t->f6 = f6;
   t->i16 = i16;
   t->j16 = j16;
   t->i26 = i26;
   t->j26 = j26;
   return(new_Hybrid(SIGID_Bb, t));
}


/* operator Il                                                                      */
/* -------------------------------------------------------------------------------- */

struct str_Il {
   int a1;
   int a2;
   int a3;
   int a4;
   int a5;
   int a6;
   struct str1 a7;
   struct str1 (*f7)(int , int , int , int );
   int i17, j17, i27, j27;
};

struct str_Hybrid *new_Il(int a1, int a2, int a3, int a4, int a5, int a6, struct str1 (*f7)(int , int , int , int ), int i17, int j17, int i27, int j27)
{
   struct str_Il *t;

   t=(struct str_Il *) calloc(1, sizeof(struct str_Il ));
   t->a1 = a1;
   t->a2 = a2;
   t->a3 = a3;
   t->a4 = a4;
   t->a5 = a5;
   t->a6 = a6;
   t->f7 = f7;
   t->i17 = i17;
   t->j17 = j17;
   t->i27 = i27;
   t->j27 = j27;
   return(new_Hybrid(SIGID_Il, t));
}

/* operator El                                                                      */
/* -------------------------------------------------------------------------------- */

struct str_El {
   int a1;
   int a2;
   int a3;
   int a4;
   int a5;
   int a6;
};

struct str_Hybrid *new_El(int a1, int a2, int a3, int a4, int a5, int a6)
{
   struct str_El *t;

   t=(struct str_El *) calloc(1, sizeof(struct str_El ));
   t->a1 = a1;
   t->a2 = a2;
   t->a3 = a3;
   t->a4 = a4;
   t->a5 = a5;
   t->a6 = a6;
   return(new_Hybrid(SIGID_El, t));
}


/* operator Nil                                                                     */
/* -------------------------------------------------------------------------------- */

struct str_Nil {
   int a1;
   int a2;
   int a3;
   int a4;
};

struct str_Hybrid *new_Nil(int a1, int a2, int a3, int a4)
{
   struct str_Nil *t;

   t=(struct str_Nil *) calloc(1, sizeof(struct str_Nil ));
   t->a1 = a1;
   t->a2 = a2;
   t->a3 = a3;
   t->a4 = a4;
   return(new_Hybrid(SIGID_Nil, t));
};


/* signature pretty printer                                                         */
/* -------------------------------------------------------------------------------- */

/* int gx; */
/* char *t1, *t2, *t3, *t4; */
/* char *r1, *r2, *r3, *r4; */

/* int a1, a2, a3, a4, a5, a6; */

char *sx(char *t, int i)
{
  if (x[i]==A) t[0]='A'; else 
  if (x[i]==C) t[0]='C'; else 
  if (x[i]==G) t[0]='G'; else 
  if (x[i]==U) t[0]='U'; else
               t[0]='N';
  t[1]=0;
  return(&t[1]);
}

char *shift(char *t, int u, int v)
{
  if ((v-u)==0) {t[0]=' '; t[1] = ' '; t[2] = 0;  return(&t[2]);}
  else {t[0]=' '; t[1] = 0;   return(&t[1]);}
}

char *one_space(char *t)
{
  t[0]=' '; t[1] = 0; return(&t[1]);
}


char *ssx(char *t, int i, int j)
{
  int k;
  for (k=i+1;k<=j;k++) t=sx(t,k); 
  return(&t[0]);
}


char *sy(char *t, int i)
{
  if (y[i]==A) t[0]='A'; else 
  if (y[i]==C) t[0]='C'; else 
  if (y[i]==G) t[0]='G'; else 
  if (y[i]==U) t[0]='U'; else
               t[0]='N';
  t[1]=0;
  return(&t[1]);
}

char *ssy(char *t, int i, int j)
{
  int k;
  for (k=i+1;k<=j;k++) t=sy(t,k);
  return(&t[0]);
}


char *blanks(char *t, int i, int j)
{
  int k,l;
  l=j-i;
  for (k=0;k<=l-1;k++) t[k]=' ';
  t[l]=0;
  return(&t[l]);
}

char *blanks2(char *t, int i, int j)
{
  int k,l;
  l=max(j-i-1,0);
  for (k=0;k<=l-1;k++) t[k]=' ';
  t[l]=0;
  return(&t[l]);
}

char *asym1(char *t, int l,int r,int u,int v)
{
  int k,ln;
  ln=max((v-u)-(r-l),0);
  for (k=0;k<=ln-1;k++) t[k]=' ';
  t[ln]=0;
  return(&t[ln]);
}

char *asym2(char *t, int l,int r,int u,int v)
{
  int k,ln;
  ln=max((r-l)-(v-u),0);
  for (k=0;k<=ln-1;k++) t[k]=' ';
  t[ln]=0;
  return(&t[ln]);
}


void pp_str_Hybrid(struct str1 l)
{
   struct str_Hybrid *c;

   c = l.alg_enum;

      if (c->utype == SIGID__NTID) {
            pp_str_Hybrid(((struct str__NTID *)(c->entry))->a1);
      } else 
      if (c->utype == SIGID_Ult) {
/* 	   gx++; */
           pp_str_Hybrid(((struct str_Ult *)(c->entry))->a3);
      } else 
      if (c->utype == SIGID_Ulb) {
            r1 = one_space(r1);
            r2 = one_space(r2);
            r3 = one_space(r3);
            a2 =((struct str_Ulb *)(c->entry))->a2; 
            r4 = sy(r4,a2);
            pp_str_Hybrid(((struct str_Ulb *)(c->entry))->a3);
      } else 
      if (c->utype == SIGID_Eds) {
	   a1=((struct str_Eds *)(c->entry))->a1; 
           a2=((struct str_Eds *)(c->entry))->a2;
            r1 = sx(r1,a1);
            r2 = one_space(r2);
            r3 = one_space(r3);
            r4 = sy(r4,a2);
            pp_str_Hybrid(((struct str_Eds *)(c->entry))->a3);
      } else 
      if (c->utype == SIGID_Edt) {
	   a1=((struct str_Edt *)(c->entry))->a1;
            r1 = sx(r1,a1);
            r2 = one_space(r2);
            r3 = one_space(r3);
            r4 = one_space(r4);
            pp_str_Hybrid(((struct str_Edt *)(c->entry))->a3);
      } else 
      if (c->utype == SIGID_Edb) {
	   a2=((struct str_Edb *)(c->entry))->a2;
            r1 = one_space(r1);
            r2 = one_space(r2);
            r3 = one_space(r3);
            r4 = sy(r4,a2);
            pp_str_Hybrid(((struct str_Edb *)(c->entry))->a3);
      } else 
      if (c->utype == SIGID_Sr) {
	   a1=((struct str_Sr *)(c->entry))->a1;
           a2=((struct str_Sr *)(c->entry))->a2;
            r1 = one_space(r1);
            r2 = sx(r2,a1);
            r3 = sy(r3,a2);
            r4 = one_space(r4);
            pp_str_Hybrid(((struct str_Sr *)(c->entry))->a3);
      } else 
      if (c->utype == SIGID_Bt) {
	   a1=((struct str_Bt *)(c->entry))->a1;
	   a2=((struct str_Bt *)(c->entry))->a2;
	   a3=((struct str_Bt *)(c->entry))->a3;
	   a4=((struct str_Bt *)(c->entry))->a4;
            r1 = one_space(r1);
            r1 = ssx(r1,a3,a4);
            r2 = sx(r2,a1);
            r2 = blanks(r2,a3,a4);
            r3 = sy(r3,a2);
            r3 = blanks(r3,a3,a4);
            r4 = one_space(r4);
            r4 = blanks(r4,a3,a4);
            pp_str_Hybrid(((struct str_Bt *)(c->entry))->a6);
      } else 
      if (c->utype == SIGID_Bb) {
	   a1=((struct str_Bb *)(c->entry))->a1;
	   a2=((struct str_Bb *)(c->entry))->a2;
	   a3=((struct str_Bb *)(c->entry))->a3;
	   a4=((struct str_Bb *)(c->entry))->a4;
	   a5=((struct str_Bb *)(c->entry))->a5;

            r1 = one_space(r1);
            r1 = blanks(r1,a4,a5);
            r2 = sx(r2,a1);
            r2 = blanks(r2,a4,a5);
            r3 = sy(r3,a2);
            r3 = blanks(r3,a4,a5);
            r4 = one_space(r4);
            r4 = ssy(r4,a4,a5);
            pp_str_Hybrid(((struct str_Bb *)(c->entry))->a6);
      } else 
      if (c->utype == SIGID_Il) {
	   a1=((struct str_Il *)(c->entry))->a1;
	   a2=((struct str_Il *)(c->entry))->a2;
	   a3=((struct str_Il *)(c->entry))->a3;
	   a4=((struct str_Il *)(c->entry))->a4;
	   a5=((struct str_Il *)(c->entry))->a5;
	   a6=((struct str_Il *)(c->entry))->a6;

            r1 = one_space(r1);
            r1 = ssx(r1,a3,a4);
            r1 = asym1(r1,a3,a4,a5,a6);

            r2 = sx(r2,a1);
            r2 = blanks(r2,a3,a4);
            r2 = asym1(r2,a3,a4,a5,a6);

            r3 = sy(r3,a2);
            r3 = blanks(r3,a5,a6);
            r3 = asym2(r3,a3,a4,a5,a6);

            r4 = one_space(r4);
            r4 = ssy(r4,a5,a6);
            r4 = asym2(r4,a3,a4,a5,a6);

            pp_str_Hybrid(((struct str_Il *)(c->entry))->a7);
      } else 
      if (c->utype == SIGID_El) {
	   a1=((struct str_El *)(c->entry))->a1;
	   a2=((struct str_El *)(c->entry))->a2;
	   a3=((struct str_El *)(c->entry))->a3;
	   a4=((struct str_El *)(c->entry))->a4;
	   a5=((struct str_El *)(c->entry))->a5;
	   a6=((struct str_El *)(c->entry))->a6;
 	if ((a4-a3) == 0) {
            r1 = one_space(r1);
            r1 = blanks(r1,a5,a6);

            r2 = sx(r2,a3);
            r2 = blanks(r2,a5,a6);

            r3 = sy(r3,a5);
            r3 = blanks(r3,a5,a6);

            r4 = one_space(r4);
            r4 = ssy(r4,a5,a6);
	} else {
            r1 = one_space(r1);
            r1 = sx(r1,a3+1);
            r1 = blanks2(r1,a5,a6);

            r2 = sx(r2,a3);
            r2 = one_space(r2);
            r2 = blanks2(r2,a5,a6);

            r3 = sy(r3,a5);
            r3 = one_space(r3);
            r3 = blanks2(r3,a5,a6);

            r4 = shift(r4,a5,a6);
            r4 = ssy(r4,a5,a6);
	}
      } else 
      if (c->utype == SIGID_Nil) {
	r1[0] = '\0';
	r2[0] = '\0';
	r3[0] = '\0';
	r4[0] = '\0';
      }
}


/* free structures                                                                  */
/* -------------------------------------------------------------------------------- */

void free_str_Hybrid(struct str_Hybrid *c)
{
   struct str1 l, l2;

   if (c != NULL) {
      if (c->utype == SIGID__NTID) {
         l = ((struct str__NTID *)(c->entry))->a1;
         free_str_Hybrid(l.alg_enum);
         free(c->entry);
      } else 
      if (c->utype == SIGID_Ult) {
         l = ((struct str_Ult *)(c->entry))->a3;
         free_str_Hybrid(l.alg_enum);
         free(c->entry);
      } else 
      if (c->utype == SIGID_Ulb) {
         l = ((struct str_Ulb *)(c->entry))->a3;
         free_str_Hybrid(l.alg_enum);
         free(c->entry);
      } else 
      if (c->utype == SIGID_Eds) {
         l = ((struct str_Eds *)(c->entry))->a3;
         free_str_Hybrid(l.alg_enum);
         free(c->entry);
      } else 
      if (c->utype == SIGID_Edt) {
         l = ((struct str_Edt *)(c->entry))->a3;
         free_str_Hybrid(l.alg_enum);
         free(c->entry);
      } else 
      if (c->utype == SIGID_Edb) {
         l = ((struct str_Edb *)(c->entry))->a3;
         free_str_Hybrid(l.alg_enum);
         free(c->entry);
      } else 
      if (c->utype == SIGID_Sr) {
         l = ((struct str_Sr *)(c->entry))->a3;
         free_str_Hybrid(l.alg_enum);
         free(c->entry);
      } else 
      if (c->utype == SIGID_Bt) {
         l = ((struct str_Bt *)(c->entry))->a6;
         free_str_Hybrid(l.alg_enum);
         free(c->entry);
      } else 
      if (c->utype == SIGID_Bb) {
         l = ((struct str_Bb *)(c->entry))->a6;
         free_str_Hybrid(l.alg_enum);
         free(c->entry);
      } else 
      if (c->utype == SIGID_Il) {
         l = ((struct str_Il *)(c->entry))->a7;
         free_str_Hybrid(l.alg_enum);
         free(c->entry);
      } else 
      if (c->utype == SIGID_El) {
         free(c->entry);
      } else 
      if (c->utype == SIGID_Nil) {
         free(c->entry);
      }
   }
   free(c);
}

/* structure builder                                                                */
/* -------------------------------------------------------------------------------- */

struct str1 build_str_Hybrid(struct str1 cl)
{
   struct str1 l;
   struct str_Hybrid *c;

   c = cl.alg_enum;
   if (c->utype == SIGID__NTID) {
     ((struct str__NTID *)(c->entry))->a1 = (*((struct str__NTID *)(c->entry))->f1)(((struct str__NTID *)(c->entry))->i11, ((struct str__NTID *)(c->entry))->j11, ((struct str__NTID *)(c->entry))->i21, ((struct str__NTID *)(c->entry))->j21);
   } else 
   if (c->utype == SIGID_Ult) {
      ((struct str_Ult *)(c->entry))->a3 = (*((struct str_Ult *)(c->entry))->f3)(((struct str_Ult *)(c->entry))->i13, ((struct str_Ult *)(c->entry))->j13, ((struct str_Ult *)(c->entry))->i23, ((struct str_Ult *)(c->entry))->j23);
   } else 
   if (c->utype == SIGID_Ulb) {
     ((struct str_Ulb *)(c->entry))->a3 = (*((struct str_Ulb *)(c->entry))->f3)(((struct str_Ulb *)(c->entry))->i13, ((struct str_Ulb *)(c->entry))->j13, ((struct str_Ulb *)(c->entry))->i23, ((struct str_Ulb *)(c->entry))->j23);
   } else 
   if (c->utype == SIGID_Eds) {
     ((struct str_Eds *)(c->entry))->a3 = (*((struct str_Eds *)(c->entry))->f3)(((struct str_Eds *)(c->entry))->i13, ((struct str_Eds *)(c->entry))->j13, ((struct str_Eds *)(c->entry))->i23, ((struct str_Eds *)(c->entry))->j23);
   } else 
   if (c->utype == SIGID_Edt) {
     ((struct str_Edt *)(c->entry))->a3 = (*((struct str_Edt *)(c->entry))->f3)(((struct str_Edt *)(c->entry))->i13, ((struct str_Edt *)(c->entry))->j13, ((struct str_Edt *)(c->entry))->i23, ((struct str_Edt *)(c->entry))->j23);
   } else 
   if (c->utype == SIGID_Edb) {
     ((struct str_Edb *)(c->entry))->a3 = (*((struct str_Edb *)(c->entry))->f3)(((struct str_Edb *)(c->entry))->i13, ((struct str_Edb *)(c->entry))->j13, ((struct str_Edb *)(c->entry))->i23, ((struct str_Edb *)(c->entry))->j23);
   } else 
   if (c->utype == SIGID_Sr) {
     ((struct str_Sr *)(c->entry))->a3 = (*((struct str_Sr *)(c->entry))->f3)(((struct str_Sr *)(c->entry))->i13, ((struct str_Sr *)(c->entry))->j13, ((struct str_Sr *)(c->entry))->i23, ((struct str_Sr *)(c->entry))->j23);
   } else 
   if (c->utype == SIGID_Bt) {
     ((struct str_Bt *)(c->entry))->a6 = (*((struct str_Bt *)(c->entry))->f6)(((struct str_Bt *)(c->entry))->i16, ((struct str_Bt *)(c->entry))->j16, ((struct str_Bt *)(c->entry))->i26, ((struct str_Bt *)(c->entry))->j26);
   } else 
   if (c->utype == SIGID_Bb) {
     ((struct str_Bb *)(c->entry))->a6 = (*((struct str_Bb *)(c->entry))->f6)(((struct str_Bb *)(c->entry))->i16, ((struct str_Bb *)(c->entry))->j16, ((struct str_Bb *)(c->entry))->i26, ((struct str_Bb *)(c->entry))->j26);
   } else 
   if (c->utype == SIGID_Il) {
     ((struct str_Il *)(c->entry))->a7 = (*((struct str_Il *)(c->entry))->f7)(((struct str_Il *)(c->entry))->i17, ((struct str_Il *)(c->entry))->j17, ((struct str_Il *)(c->entry))->i27, ((struct str_Il *)(c->entry))->j27);
   } else 
   if (c->utype == SIGID_El) {
   } else 
   if (c->utype == SIGID_Nil) {
   }
   return(cl);
}

/* table declarations                                                               */
/* -------------------------------------------------------------------------------- */

double **tbl_unpaired_left_top;
double **tbl_unpaired_left_bot;
double **tbl_closed;

/* forward declarations                                                             */
/* -------------------------------------------------------------------------------- */

double calc_hybrid(int i1, int j1, int i2, int j2);

/* table calculations                                                               */
/* -------------------------------------------------------------------------------- */

/* table calculation for production hybrid                                          */
/* -------------------------------------------------------------------------------- */

double calc_hybrid(int i1, int j1, int i2, int j2)
{
   double v1, v2, v3, v4, v5;

   /* ---------------------------------- start of --------------------------------- */
   /* --------------------- v1 = nil <<< tt(uregion, uregion) --------------------- */
   if (((j1-i1) >= 0) && ((j2-i2) >= 0)) {
      v1 = 0;
      /* No iteration neccessary! */
   }
   else {
      v1 = 65000;
   }
   /* --------------------- v1 = nil <<< tt(uregion, uregion) --------------------- */
   /* ---------------------------------- finished --------------------------------- */

   /* -------------------------- v2 = p unpaired_left_top ------------------------- */
   if (((j1-i1) >= 1) && ((j2-i2) >= 1)) {
      v2 = tbl_unpaired_left_top[i1][i2];
   }
   else {
      v2 = 65000;
   }
   /* ------------------------------- v3 = p closed ------------------------------- */
   if (((j1-i1) >= 1) && ((j2-i2) >= 1)) {
      v3 = tbl_closed[i1][i2];
   }
   else {
      v3 = 65000;
   }
   v4 = v2 < v3 ? v2 : v3;
   v5 = v1 < v4 ? v1 : v4;
   return(v5);
}

/* table calculation for production unpaired_left_top                               */
/* -------------------------------------------------------------------------------- */

void calc_unpaired_left_top(int i1, int j1, int i2, int j2)
{
   double v1, v2, v3;

   /* ---------------------------------- start of --------------------------------- */
   /* ---------- v1 = ult <<< (tt(lbase, empty)) ~~~ p unpaired_left_top ---------- */
   if (((j1-i1) >= 2) && ((j2-i2) >= 1)) {
      v1 = tbl_unpaired_left_top[i1+1][i2];
      /* No iteration neccessary! */
   }
   else {
      v1 = 65000;
   }
   /* ---------- v1 = ult <<< (tt(lbase, empty)) ~~~ p unpaired_left_top ---------- */
   /* ---------------------------------- finished --------------------------------- */

   /* -------------------------- v2 = p unpaired_left_bot ------------------------- */
   if (((j1-i1) >= 1) && ((j2-i2) >= 1) && ((i2 < helix_start) || (i2 >= helix_end))) {
      v2 = tbl_unpaired_left_bot[i1][i2];
   }
   else {
      v2 = 65000;
   }
   v3 = v1 < v2 ? v1 : v2;
   /* ------------------------- assign table entry result ------------------------- */
   if (((j1-i1) >= 1) && ((j2-i2) >= 1)) {
      tbl_unpaired_left_top[i1][i2] = v3;
   }
}

/* table calculation for production unpaired_left_bot                               */
/* -------------------------------------------------------------------------------- */

void calc_unpaired_left_bot(int i1, int j1, int i2, int j2)
{
   double v1, v2, v3, v4, v5, v6, v7;

   /* ---------------------------------- start of --------------------------------- */
   /* ---------- v1 = ulb <<< (tt(empty, lbase)) ~~~ p unpaired_left_bot ---------- */
   if (((j1-i1) >= 1) && ((j2-i2) >= 2) && ((i2 < helix_start-1) || (i2 >= helix_end))) {
      v1 = tbl_unpaired_left_bot[i1][i2+1];
      /* No iteration neccessary! */
   }
   else {
      v1 = 65000;
   }
   /* ---------- v1 = ulb <<< (tt(empty, lbase)) ~~~ p unpaired_left_bot ---------- */
   /* ---------------------------------- finished --------------------------------- */

   /* ---------------------------------- start of --------------------------------- */
   /* ---------------- v2 = eds <<< (tt(lbase, lbase)) ~~~ p closed --------------- */
   if (((j1-i1) >= 2) && ((j2-i2) >= 2) && compl(x[i1+2],y[i2+2]) && ((i2 < helix_start) || (i2 >= helix_end))) {
      v2 = (tbl_closed[i1+1][i2+1] + dl_energy((i1+1) + 1, (i2+1) + 1)) + dr_energy((i1+1) + 1, (i2+1) + 1);
      /* No iteration neccessary! */
   }
   else {
      v2 = 65000;
   }
   /* ---------------- v2 = eds <<< (tt(lbase, lbase)) ~~~ p closed --------------- */
   /* ---------------------------------- finished --------------------------------- */

   /* ---------------------------------- start of --------------------------------- */
   /* ---------------- v3 = edt <<< (tt(lbase, empty)) ~~~ p closed --------------- */
   if (((j1-i1) >= 2) && ((j2-i2) >= 1) && compl(x[i1+2],y[i2+1])) {
      v3 = tbl_closed[i1+1][i2] + dl_energy((i1+1) + 1, (i2) + 1);
      /* No iteration neccessary! */
   }
   else {
      v3 = 65000;
   }
   /* ---------------- v3 = edt <<< (tt(lbase, empty)) ~~~ p closed --------------- */
   /* ---------------------------------- finished --------------------------------- */

   /* ---------------------------------- start of --------------------------------- */
   /* ---------------- v4 = edb <<< (tt(empty, lbase)) ~~~ p closed --------------- */
   if (((j1-i1) >= 1) && ((j2-i2) >= 2) && compl(x[i1+1],y[i2+2]) && ((i2 < helix_start) || (i2 >= helix_end))) {
      v4 = tbl_closed[i1][i2+1] + dr_energy((i1) + 1, (i2+1) + 1);
      /* No iteration neccessary! */
   }
   else {
      v4 = 65000;
   }
   /* ---------------- v4 = edb <<< (tt(empty, lbase)) ~~~ p closed --------------- */
   /* ---------------------------------- finished --------------------------------- */

   v5 = v3 < v4 ? v3 : v4;
   v6 = v2 < v5 ? v2 : v5;
   v7 = v1 < v6 ? v1 : v6;
   /* ------------------------- assign table entry result ------------------------- */
   if (((j1-i1) >= 1) && ((j2-i2) >= 1)) {
      tbl_unpaired_left_bot[i1][i2] = v7;
   }
}

/* table calculation for production closed                                          */
/* -------------------------------------------------------------------------------- */

void calc_closed(int i1, int j1, int i2, int j2)
{
   double v1, v2, v3, v4, v5, v6, v7, v7b, v7c, v7d, v7e, v7f, v7g, v8, v9, v10, v11, v12;
   int k;
   int k2;
   int k3;
   int k4;

   /* ---------------------------------- start of --------------------------------- */
   /* - v1 = sr <<< (tt(lbase, lbase) `with` (pairingTTcross compl)) ~~~ p closed - */
   if (((j1-i1) >= 2) && ((j2-i2) >= 2)) {
      if (compl(x[i1+1], y[i2+1])) {
         v1 = sr_energy(i1+1, i2+1) + tbl_closed[i1+1][i2+1];
         /* No iteration neccessary! */
      }
      else {
         v1 = 65000;
      }
   }
   else {
      v1 = 65000;
   }
   /* - v1 = sr <<< (tt(lbase, lbase) `with` (pairingTTcross compl)) ~~~ p closed - */
   /* ---------------------------------- finished --------------------------------- */

   /* ---------------------------------- start of --------------------------------- */
   /*  v3 = bt <<< (tt(lbase, lbase) `with` (pairingTTcross compl)) ~~~ (tt(region, empty) `with` (sizeTT 1 15 0 0)) ~~~ p closed  */
   if (((j1-i1) >= 3) && ((j2-i2) >= 2) && compl(x[i1+1], y[i2+1]) && ((i2 < helix_start) || (i2 >= helix_end-1))) {
      v3 = 65000;
      for (k=i1+2; k<=min(i1+bloop_upper_limit+1, j1-1); k++) {
	if (x[k]==X)
	  break;
	  v2 = (tbl_closed[k][i2+1] + bl_stacking((k) - (i1+1), 0, i1+1, i2+1)) + bl_ent((k) - (i1+1));
	  /* No iteration neccessary! */
	v3 = v2 < v3 ? v2 : v3;
      }
   }
   else {
     v3 = 65000;
   }
   /*  v3 = bt <<< (tt(lbase, lbase) `with` (pairingTTcross compl)) ~~~ (tt(region, empty) `with` (sizeTT 1 15 0 0)) ~~~ p closed  */
   /* ---------------------------------- finished --------------------------------- */

   /* ---------------------------------- start of --------------------------------- */
   /*  v5 = bb <<< (tt(lbase, lbase) `with` (pairingTTcross compl)) ~~~ (tt(empty, region) `with` (sizeTT 0 0 1 15)) ~~~ p closed  */
   if (((j1-i1) >= 2) && ((j2-i2) >= 3) & compl(x[i1+1], y[i2+1])) {
      v5 = 65000;
      for (k2=i2+2; k2<=min(i2+bloop_upper_limit+1, j2-1); k2++) {
        if ((k2 > helix_start) && (k2 <= helix_end))
	  break;
	v4 = (tbl_closed[i1+1][k2] + bl_stacking(0, (k2) - (i2+1), i1+1, i2+1)) + bl_ent((k2) - (i2+1));
	  /* No iteration neccessary! */
	v5 = v4 < v5 ? v4 : v5;
      }
   }
   else {
     v5 = 65000;
   }
   /*  v5 = bb <<< (tt(lbase, lbase) `with` (pairingTTcross compl)) ~~~ (tt(empty, region) `with` (sizeTT 0 0 1 15)) ~~~ p closed  */
   /* ---------------------------------- finished --------------------------------- */

   /* ---------------------------------- start of --------------------------------- */
   /*  v7 = il <<< (tt(lbase, lbase) `with` (pairingTTcross compl)) ~~~ (tt(region, region) `with` (sizeTT 1 15 1 15)) ~~~ p closed  */
   
   if (((j1-i1) >= 3) && ((j2-i2) >= 3) && compl(x[i1+1], y[i2+1])) {
     v7 = 65000;
     /* special internal loops: */
     for (k3=i1+2; k3<=min(i1+min(3,iloop_upper_limit+1), j1-1); k3++) {
       if (x[k3]==X)
	 break;
       for (k4=i2+2; k4<=min(i2+min(3,iloop_upper_limit+1), j2-1); k4++) {
	 if ((k4 > helix_start) && (i2 < helix_end-1))
	   break;
	 if (compl(x[k3+1],y[k4+1])) {
	   v6 = do_il_special(i1+1, i2+1, i1+1, k3, i2+1, k4, tbl_closed[k3][k4]);
	   /* No iteration neccessary! */
	   v7 = v6 < v7 ? v6 : v7;
	 }
       }
     }
     v7g = 65000; 
     v7b = 65000;
     if ((x[k3]!=X)) { /*  && !((k4 > helix_start) && (i2 < helix_end-1))) { */
       /* normal internal loops: */
       for (k3=i1+2; k3<=min(i1+min(3,iloop_upper_limit+1), j1-1); k3++) {
	 if (x[k3]==X)
	   break;
	 for (k4=i2+4; k4<=min(i2+iloop_upper_limit+1, j2-1); k4++) {
	   if ((k4 > helix_start) && (i2 < helix_end-1))
	     break;
	   if (compl(x[k3+1],y[k4+1])) {
	     v6 = do_il(i1+1, i2+1, i1+1, k3, i2+1, k4, tbl_closed[k3][k4]);
	     /* No iteration neccessary! */
	     v7b = v6 < v7b ? v6 : v7b;
	   }
	 }
       }
       if (v7b < 65000)
	 v7b += il_stack_open(i1+1,i2+1);
       v7c = 65000;
       /* normal internal loops: */
       if ((x[k3]!=X)) { /*  && !((k4 > helix_start) && (i2 < helix_end-1))) { */
	 for (k3=i1+4; k3<=min(i1+iloop_upper_limit+1, j1-1); k3++) {
	   if (x[k3]==X)
	     break;
	   for (k4=i2+2; k4<=min(i2+min(3,iloop_upper_limit+1), j2-1); k4++) {
	     if ((k4 > helix_start) && (i2 < helix_end-1))
	       break;
	     if (compl(x[k3+1],y[k4+1])) {
	       v6 = do_il(i1+1, i2+1, i1+1, k3, i2+1, k4, tbl_closed[k3][k4]);
	       /* No iteration neccessary! */
	       v7c = v6 < v7c ? v6 : v7c;
	     }
	   }
	 }
	 if (v7c < 65000)
	   v7c += il_stack_open(i1+1,i2+1);
	 v7d = 65000;
	 /* normal internal loops: */
	 if ((x[k3]!=X)) { /*  && !((k4 > helix_start) && (i2 < helix_end-1))) { */
	   for (k3=i1+4; k3<=min(i1+iloop_upper_limit+1, j1-1); k3++) {
	     if (x[k3]==X)
	       break;
	     for (k4=i2+4; k4<=min(i2+iloop_upper_limit+1, j2-1); k4++) {
	       if ((k4 > helix_start) && (i2 < helix_end-1))
		 break;
	       if (compl(x[k3+1],y[k4+1])) {
		 v6 = do_il(i1+1, i2+1, i1+1, k3, i2+1, k4, tbl_closed[k3][k4]);
		 /* No iteration neccessary! */
		 v7d = v6 < v7d ? v6 : v7d;
	       }
	     }
	   }
	   if (v7d < 65000)
	     v7d += il_stack_open(i1+1,i2+1);
	 }

       }
       v7e = v7b < v7c ? v7b : v7c;
       v7f = v7d < v7e ? v7d : v7e;
       v7g = v7g < v7f ? v7g : v7f;
     }
     v7 = v7g < v7 ? v7g : v7;
   }
   else {
     v7 = 65000;
   }
   /*  v7 = il <<< (tt(lbase, lbase) `with` (pairingTTcross compl)) ~~~ (tt(region, region) `with` (sizeTT 1 15 1 15)) ~~~ p closed  */
   /* ---------------------------------- finished --------------------------------- */

   /* ---------------------------------- start of --------------------------------- */
   /*  v8 = el <<< (tt(lbase, lbase) `with` (pairingTTcross compl)) ~~~ (tt(uregion, uregion))  */
   if (((j1-i1) >= 1) && ((j2-i2) >= 1) && (j1==i1+1 || x[i1+2]!='X') && ((i2 >= helix_end-1) || (helix_end > j2))) {
      if (compl(x[i1+1], y[i2+1])) {
         v8 = ((((j1) - (i1+1)) > 0) ? dli_energy(i1+1, i2+1) : 0) + ((((j2) - (i2+1)) > 0) ? dri_energy(i1+1, i2+1) : 0);
         /* No iteration neccessary! */
      }
      else {
         v8 = 65000;
      }
   }
   else {
      v8 = 65000;
   }
   /*  v8 = el <<< (tt(lbase, lbase) `with` (pairingTTcross compl)) ~~~ (tt(uregion, uregion))  */
   /* ---------------------------------- finished --------------------------------- */


   v9 = v7 < v8 ? v7 : v8;
   v10 = v5 < v9 ? v5 : v9;
   v11 = v3 < v10 ? v3 : v10;
   v12 = v1 < v11 ? v1 : v11;
   /* ------------------------- assign table entry result ------------------------- */
   if (((j1-i1) >= 1) && ((j2-i2) >= 1)) {
      tbl_closed[i1][i2] = v12;
   }
}

/* forward declarations for backtracing functions                                   */
/* -------------------------------------------------------------------------------- */

struct str1 back_unpaired_left_bot(int i1, int j1, int i2, int j2);

struct str1 back_closed(int i1, int j1, int i2, int j2);

struct str1 back_hybrid(int i1, int j1, int i2, int j2);

/* backtracing code                                                                 */
/* -------------------------------------------------------------------------------- */

/* table calculation for production hybrid                                          */
/* -------------------------------------------------------------------------------- */

struct str1 back_hybrid(int i1, int j1, int i2, int j2)
{
   struct str1 v1, v2, v3, v4, v5;

   int k, best_k;

   /* ---------------------------------- start of --------------------------------- */
   /* --------------------- v1 = nil <<< tt(uregion, uregion) --------------------- */
   if (((j1-i1) >= 0) && ((j2-i2) >= 0)) {
      v1.alg_mfe = 0;
      v1.alg_enum = new_Nil(i1, j1, i2, j2);
      /* No iteration neccessary! */
   }
   else {
      v1.alg_mfe = 65000;
      v1.alg_enum = NULL;
   }
   /* --------------------- v1 = nil <<< tt(uregion, uregion) --------------------- */
   /* ---------------------------------- finished --------------------------------- */

   /* -------------------------- v2 = p unpaired_left_top ------------------------- */
   /* +------------------------------------------------------------------------------------ */
   /* Nonterminal unpaired_left_top is implemented as a tabulated                           */
   /* function which yields atomar results. Since we are in list context,                   */
   /* we need to wrap the result of unpaired_left_top into a single list element.           */
   /* +------------------------------------------------------------------------------------ */


/*    Find best unpaired_left_bot without backtracking recursion to avoid */
/*      stack overflow (hand made): */

   if (((j1-i1) >= 1) && ((j2-i2) >= 1) && ((i2 < helix_start) || (i2 >= helix_end))) {

     v2.alg_mfe = 65000;
     v2.alg_enum = NULL;

     best_k = -1;

     for (k=i1; k<j1; k++) {
       if (tbl_unpaired_left_bot[k][i2] < v2.alg_mfe) {
	 v2.alg_mfe = tbl_unpaired_left_bot[k][i2];
	 best_k = k;
       }
     }

     if (best_k != -1)
       v2.alg_enum = new__NTID(back_unpaired_left_bot, best_k, j1, i2, j2);

   }
   else {
     v2.alg_mfe = 65000;
     v2.alg_enum = NULL;
   }

/*    if (((j1-i1) >= 1) && ((j2-i2) >= 1)) { */
/*       v2.alg_mfe = tbl_unpaired_left_top[i1][i2]; */
/*       v2.alg_enum = new__NTID(back_unpaired_left_top, i1, j1, i2, j2); */
/*    } */
/*    else { */
/*       v2.alg_mfe = 65000; */
/*       v2.alg_enum = NULL; */
/*    } */




   /* ------------------------------- v3 = p closed ------------------------------- */
   /* +---------------------------------------------------------------------------- */
   /* Nonterminal closed is implemented as a tabulated                              */
   /* function which yields atomar results. Since we are in list context,           */
   /* we need to wrap the result of closed into a single list element.              */
   /* +---------------------------------------------------------------------------- */
   if (((j1-i1) >= 1) && ((j2-i2) >= 1)) {
      v3.alg_mfe = tbl_closed[i1][i2];
      v3.alg_enum = new__NTID(back_closed, i1, j1, i2, j2);
   }
   else {
      v3.alg_mfe = 65000;
      v3.alg_enum = NULL;
   }
   /* ---------------------------- v4 = minimum(v2, v3) --------------------------- */
   if (v2.alg_mfe < v3.alg_mfe) {
     gx = best_k + 1;
     v4 = v2;
     free_str_Hybrid(v3.alg_enum);
   }
   else {
     gx = 0 + 1;
     v4 = v3;
     free_str_Hybrid(v2.alg_enum);
   }
/*    v4 = v2.alg_mfe < v3.alg_mfe ? v2 : v3; */
   /* ---------------------------- v5 = minimum(v1, v4) --------------------------- */
   if (v1.alg_mfe < v4.alg_mfe) {
     gx = 0;
     v5 = v1;
     free_str_Hybrid(v4.alg_enum);
   }
   else {
     v5 = v4;
     free_str_Hybrid(v1.alg_enum);
   }
/*    v5 = v1.alg_mfe < v4.alg_mfe ? v1 : v4; */
   /* ------------------------- build candidate structures ------------------------ */

   return(build_str_Hybrid(v5));
}

/* table calculation for production unpaired_left_top                               */
/* -------------------------------------------------------------------------------- */

/* removed after hand made back_hybrid optimisation !*/



/* table calculation for production unpaired_left_bot                               */
/* -------------------------------------------------------------------------------- */

struct str1 back_unpaired_left_bot(int i1, int j1, int i2, int j2)
{
   struct str1 v1, v2, v3, v4, v5, v6, v7;

   /* ---------------------------------- start of --------------------------------- */
   /* ---------- v1 = ulb <<< (tt(empty, lbase)) ~~~ p unpaired_left_bot ---------- */
   if (((j1-i1) >= 1) && ((j2-i2) >= 2) && ((i2 < helix_start-1) || (i2 >= helix_end))) {
      v1.alg_mfe = tbl_unpaired_left_bot[i1][i2+1];
      v1.alg_enum = new_Ulb(i1, i2+1, back_unpaired_left_bot, i1, j1, i2+1, j2);
      /* No iteration neccessary! */
   }
   else {
      v1.alg_mfe = 65000;
      v1.alg_enum = NULL;
   }
   /* ---------- v1 = ulb <<< (tt(empty, lbase)) ~~~ p unpaired_left_bot ---------- */
   /* ---------------------------------- finished --------------------------------- */

   /* ---------------------------------- start of --------------------------------- */
   /* ---------------- v2 = eds <<< (tt(lbase, lbase)) ~~~ p closed --------------- */
   if (((j1-i1) >= 2) && ((j2-i2) >= 2) && compl(x[i1+2],y[i2+2]) && ((i2 < helix_start) || (i2 >= helix_end))) {
      v2.alg_mfe = (tbl_closed[i1+1][i2+1] + dl_energy((i1+1) + 1, (i2+1) + 1)) + dr_energy((i1+1) + 1, (i2+1) + 1);
      v2.alg_enum = new_Eds(i1+1, i2+1, back_closed, i1+1, j1, i2+1, j2);
      /* No iteration neccessary! */
   }
   else {
      v2.alg_mfe = 65000;
      v2.alg_enum = NULL;
   }
   /* ---------------- v2 = eds <<< (tt(lbase, lbase)) ~~~ p closed --------------- */
   /* ---------------------------------- finished --------------------------------- */

   /* ---------------------------------- start of --------------------------------- */
   /* ---------------- v3 = edt <<< (tt(lbase, empty)) ~~~ p closed --------------- */
   if (((j1-i1) >= 2) && ((j2-i2) >= 1) && compl(x[i1+2],y[i2+1])) {
      v3.alg_mfe = tbl_closed[i1+1][i2] + dl_energy((i1+1) + 1, (i2) + 1);
      v3.alg_enum = new_Edt(i1+1, i2, back_closed, i1+1, j1, i2, j2);
      /* No iteration neccessary! */
   }
   else {
      v3.alg_mfe = 65000;
      v3.alg_enum = NULL;
   }
   /* ---------------- v3 = edt <<< (tt(lbase, empty)) ~~~ p closed --------------- */
   /* ---------------------------------- finished --------------------------------- */

   /* ---------------------------------- start of --------------------------------- */
   /* ---------------- v4 = edb <<< (tt(empty, lbase)) ~~~ p closed --------------- */
   if (((j1-i1) >= 1) && ((j2-i2) >= 2) && compl(x[i1+1],y[i2+2]) && ((i2 < helix_start) || (i2 >= helix_end))) {
      v4.alg_mfe = tbl_closed[i1][i2+1] + dr_energy((i1) + 1, (i2+1) + 1);
      v4.alg_enum = new_Edb(i1, i2+1, back_closed, i1, j1, i2+1, j2);
      /* No iteration neccessary! */
   }
   else {
      v4.alg_mfe = 65000;
      v4.alg_enum = NULL;
   }
   /* ---------------- v4 = edb <<< (tt(empty, lbase)) ~~~ p closed --------------- */
   /* ---------------------------------- finished --------------------------------- */

   /* ---------------------------- v5 = minimum(v3, v4) --------------------------- */
   if (v3.alg_mfe < v4.alg_mfe) {
     v5 = v3;
     free_str_Hybrid(v4.alg_enum);
   }
   else {
     v5 = v4;
     free_str_Hybrid(v3.alg_enum);
   }
/*    v5 = v3.alg_mfe < v4.alg_mfe ? v3 : v4; */
   /* ---------------------------- v6 = minimum(v2, v5) --------------------------- */
   if (v2.alg_mfe < v5.alg_mfe) {
     v6 = v2;
     free_str_Hybrid(v5.alg_enum);
   }
   else {
     v6 = v5;
     free_str_Hybrid(v2.alg_enum);
   }
/*    v6 = v2.alg_mfe < v5.alg_mfe ? v2 : v5; */
   /* ---------------------------- v7 = minimum(v1, v6) --------------------------- */
   if (v1.alg_mfe < v6.alg_mfe) {
     v7 = v1;
     free_str_Hybrid(v6.alg_enum);
   }
   else {
     v7 = v6;
     free_str_Hybrid(v1.alg_enum);
   }
/*    v7 = v1.alg_mfe < v6.alg_mfe ? v1 : v6; */
   /* ------------------------- build candidate structures ------------------------ */

   return(build_str_Hybrid(v7));
}

/* table calculation for production closed                                          */
/* -------------------------------------------------------------------------------- */

struct str1 back_closed(int i1, int j1, int i2, int j2)
{
   struct str1 v1, v2, v3, v4, v5, v6, v7, v7b, v7c, v7d, v7e, v7f, v7g, v8, v9, v10, v11, v12;
   int k;
   int k2;
   int k3;
   int k4;

   /* ---------------------------------- start of --------------------------------- */
   /* - v1 = sr <<< (tt(lbase, lbase) `with` (pairingTTcross compl)) ~~~ p closed - */
   if (((j1-i1) >= 2) && ((j2-i2) >= 2)) {
      if (compl(x[i1+1], y[i2+1])) {
         v1.alg_mfe = sr_energy(i1+1, i2+1) + tbl_closed[i1+1][i2+1];
         v1.alg_enum = new_Sr(i1+1, i2+1, back_closed, i1+1, j1, i2+1, j2);
         /* No iteration neccessary! */
      }
      else {
         v1.alg_mfe = 65000;
         v1.alg_enum = NULL;
      }
   }
   else {
      v1.alg_mfe = 65000;
      v1.alg_enum = NULL;
   }
   /* - v1 = sr <<< (tt(lbase, lbase) `with` (pairingTTcross compl)) ~~~ p closed - */
   /* ---------------------------------- finished --------------------------------- */

   /* ---------------------------------- start of --------------------------------- */
   /*  v3 = bt <<< (tt(lbase, lbase) `with` (pairingTTcross compl)) ~~~ (tt(region, empty) `with` (sizeTT 1 15 0 0)) ~~~ p closed  */
   if (((j1-i1) >= 3) && ((j2-i2) >= 2) && compl(x[i1+1], y[i2+1]) && ((i2 < helix_start) || (i2 >= helix_end-1))) {
      v3.alg_mfe = 65000;
      v3.alg_enum = NULL;
      for (k=i1+2; k<=min(i1+bloop_upper_limit+1, j1-1); k++) {
	if (x[k]==X)
	  break;
            v2.alg_mfe = (tbl_closed[k][i2+1] + bl_stacking((k) - (i1+1), 0, i1+1, i2+1)) + bl_ent((k) - (i1+1));
            v2.alg_enum = new_Bt(i1+1, i2+1, i1+1, k, i2+1, back_closed, k, j1, i2+1, j2);
            /* No iteration neccessary! */
         /* ------------------------- v3 = minimum(v2, v3) ------------------------ */
         if (v2.alg_mfe < v3.alg_mfe) {
           free_str_Hybrid(v3.alg_enum);
           v3 = v2;
         } else free_str_Hybrid(v2.alg_enum);
      }
   }
   else {
      v3.alg_mfe = 65000;
      v3.alg_enum = NULL;
   }
   /*  v3 = bt <<< (tt(lbase, lbase) `with` (pairingTTcross compl)) ~~~ (tt(region, empty) `with` (sizeTT 1 15 0 0)) ~~~ p closed  */
   /* ---------------------------------- finished --------------------------------- */

   /* ---------------------------------- start of --------------------------------- */
   /*  v5 = bb <<< (tt(lbase, lbase) `with` (pairingTTcross compl)) ~~~ (tt(empty, region) `with` (sizeTT 0 0 1 15)) ~~~ p closed  */
   if (((j1-i1) >= 2) && ((j2-i2) >= 3) && compl(x[i1+1], y[i2+1])) {
      v5.alg_mfe = 65000;
      v5.alg_enum = NULL;
      for (k2=i2+2; k2<=min(i2+bloop_upper_limit+1, j2-1); k2++) {
        if ((k2 > helix_start) && (k2 <= helix_end))
	  break;
	v4.alg_mfe = (tbl_closed[i1+1][k2] + bl_stacking(0, (k2) - (i2+1), i1+1, i2+1)) + bl_ent((k2) - (i2+1));
	v4.alg_enum = new_Bb(i1+1, i2+1, i1+1, i2+1, k2, back_closed, i1+1, j1, k2, j2);
	/* No iteration neccessary! */
	/* ------------------------- v5 = minimum(v4, v5) ------------------------ */
	/* v5 = v4.alg_mfe < v5.alg_mfe ? v4 : v5; */
	if (v4.alg_mfe < v5.alg_mfe) {
	  free_str_Hybrid(v5.alg_enum);
	  v5 = v4;
	} else free_str_Hybrid(v4.alg_enum);
      }
   }
   else {
      v5.alg_mfe = 65000;
      v5.alg_enum = NULL;
   }
   /*  v5 = bb <<< (tt(lbase, lbase) `with` (pairingTTcross compl)) ~~~ (tt(empty, region) `with` (sizeTT 0 0 1 15)) ~~~ p closed  */
   /* ---------------------------------- finished --------------------------------- */

   /* ---------------------------------- start of --------------------------------- */
   /*  v7 = il <<< (tt(lbase, lbase) `with` (pairingTTcross compl)) ~~~ (tt(region, region) `with` (sizeTT 1 15 1 15)) ~~~ p closed  */
   if (((j1-i1) >= 3) && ((j2-i2) >= 3) && compl(x[i1+1], y[i2+1])) {

     v7.alg_mfe = 65000;
     v7.alg_enum = NULL;
     /* special internal loops: */
     for (k3=i1+2; k3<=min(i1+min(3,iloop_upper_limit+1), j1-1); k3++) {
       if (x[k3]==X)
	 break;
       for (k4=i2+2; k4<=min(i2+min(3,iloop_upper_limit+1), j2-1); k4++) {
	 if ((k4 > helix_start) && (i2 < helix_end-1))
	   break;
	 if (compl(x[k3+1],y[k4+1])) {
	   v6.alg_mfe = do_il_special(i1+1, i2+1, i1+1, k3, i2+1, k4, tbl_closed[k3][k4]);
	   v6.alg_enum = new_Il(i1+1, i2+1, i1+1, k3, i2+1, k4, back_closed, k3, j1, k4, j2);
	   /* No iteration neccessary! */
	   if (v6.alg_mfe < v7.alg_mfe) {
	     free_str_Hybrid(v7.alg_enum);
	     v7 = v6;
	   } else free_str_Hybrid(v6.alg_enum);
	 }
       }
     }
     v7g.alg_mfe = 65000;
     v7g.alg_enum = NULL;
     v7b.alg_mfe = 65000;
     v7b.alg_enum = NULL;
     if ((x[k3]!=X)) { /*  && !((k4 > helix_start) && (i2 < helix_end-1))) { */

       /* normal internal loops: */
       for (k3=i1+2; k3<=min(i1+min(3,iloop_upper_limit+1), j1-1); k3++) {
	 if (x[k3]==X)
	   break;
	 for (k4=i2+4; k4<=min(i2+iloop_upper_limit+1, j2-1); k4++) {
	   if ((k4 > helix_start) && (i2 < helix_end-1))
	     break;
	   if (compl(x[k3+1],y[k4+1])) {
	     v6.alg_mfe = do_il(i1+1, i2+1, i1+1, k3, i2+1, k4, tbl_closed[k3][k4]);
	     v6.alg_enum = new_Il(i1+1, i2+1, i1+1, k3, i2+1, k4, back_closed, k3, j1, k4, j2);
	     /* No iteration neccessary! */
	     if (v6.alg_mfe < v7b.alg_mfe) {
	       free_str_Hybrid(v7b.alg_enum);
	       v7b = v6;
	     } else free_str_Hybrid(v6.alg_enum);
	   }
	 }
       }
       if (v7b.alg_mfe < 65000)
	 v7b.alg_mfe += il_stack_open(i1+1,i2+1);
       
       v7c.alg_mfe = 65000;
       v7c.alg_enum = NULL;
       if ((x[k3]!=X)) { /*  && !((k4 > helix_start) && (i2 < helix_end-1))) { */

	 for (k3=i1+4; k3<=min(i1+iloop_upper_limit+1, j1-1); k3++) {
	   if (x[k3]==X)
	     break;
	   for (k4=i2+2; k4<=min(i2+min(3,iloop_upper_limit+1), j2-1); k4++) {
	     if ((k4 > helix_start) && (i2 < helix_end-1))
	       break;
	     if (compl(x[k3+1],y[k4+1])) {
	       v6.alg_mfe = do_il(i1+1, i2+1, i1+1, k3, i2+1, k4, tbl_closed[k3][k4]);
	       v6.alg_enum = new_Il(i1+1, i2+1, i1+1, k3, i2+1, k4, back_closed, k3, j1, k4, j2);
	       /* No iteration neccessary! */
	       if (v6.alg_mfe < v7c.alg_mfe) {
		 free_str_Hybrid(v7c.alg_enum);
		 v7c = v6;
	       } else free_str_Hybrid(v6.alg_enum);
	     }
	   }
	 }
	 if (v7c.alg_mfe < 65000)
	   v7c.alg_mfe += il_stack_open(i1+1,i2+1);
       }

       v7d.alg_mfe = 65000;
       v7d.alg_enum = NULL;
       if ((x[k3]!=X)) { /*  && !((k4 > helix_start) && (i2 < helix_end-1))) { */

	 /* normal internal loops: */
	 for (k3=i1+4; k3<=min(i1+iloop_upper_limit+1, j1-1); k3++) {
	   if (x[k3]==X)
	     break;
	   for (k4=i2+4; k4<=min(i2+iloop_upper_limit+1, j2-1); k4++) {
	     if ((k4 > helix_start) && (i2 < helix_end-1))
	       break;
	     if (compl(x[k3+1],y[k4+1])) {
	       v6.alg_mfe = do_il(i1+1, i2+1, i1+1, k3, i2+1, k4, tbl_closed[k3][k4]);
	       v6.alg_enum = new_Il(i1+1, i2+1, i1+1, k3, i2+1, k4, back_closed, k3, j1, k4, j2);
	       /* No iteration neccessary! */
	       if (v6.alg_mfe < v7d.alg_mfe) {
		 free_str_Hybrid(v7d.alg_enum);
		 v7d = v6;
	       } else free_str_Hybrid(v6.alg_enum);
	     }
	   }
	 }
	 if (v7d.alg_mfe < 65000)
	   v7d.alg_mfe += il_stack_open(i1+1,i2+1);

       }

       if (v7b.alg_mfe < v7c.alg_mfe) {
	 free_str_Hybrid(v7c.alg_enum);
	 v7e = v7b;
       }
       else {
	 free_str_Hybrid(v7b.alg_enum);
	 v7e = v7c;
       }

       if (v7e.alg_mfe < v7d.alg_mfe) {
	 free_str_Hybrid(v7d.alg_enum);
	 v7f = v7e;
       }
       else {
	 free_str_Hybrid(v7e.alg_enum);
	 v7f = v7d;
       }

       if (v7g.alg_mfe < v7f.alg_mfe) {
	 free_str_Hybrid(v7f.alg_enum);
       }
       else {
	 free_str_Hybrid(v7g.alg_enum);
	 v7g = v7f;
       }
     }
     if (v7.alg_mfe < v7g.alg_mfe) {
       free_str_Hybrid(v7g.alg_enum);
     }
     else {
       free_str_Hybrid(v7.alg_enum);
       v7 = v7g;
     }
   }
   else {
      v7.alg_mfe = 65000;
      v7.alg_enum = NULL;
   }
   /*  v7 = il <<< (tt(lbase, lbase) `with` (pairingTTcross compl)) ~~~ (tt(region, region) `with` (sizeTT 1 15 1 15)) ~~~ p closed  */
   /* ---------------------------------- finished --------------------------------- */

   /* ---------------------------------- start of --------------------------------- */
   /*  v8 = el <<< (tt(lbase, lbase) `with` (pairingTTcross compl)) ~~~ (tt(uregion, uregion))  */
   if (((j1-i1) >= 1) && ((j2-i2) >= 1)  && (j1==i1+1 || x[i1+2]!='X') && ((i2 >= helix_end-1) || (helix_end > j2))) {
      if (compl(x[i1+1], y[i2+1])) {
         v8.alg_mfe = ((((j1) - (i1+1)) > 0) ? dli_energy(i1+1, i2+1) : 0) + ((((j2) - (i2+1)) > 0) ? dri_energy(i1+1, i2+1) : 0);
         v8.alg_enum = new_El(i1+1, i2+1, i1+1, j1, i2+1, j2);
         /* No iteration neccessary! */
      }
      else {
         v8.alg_mfe = 65000;
         v8.alg_enum = NULL;
      }
   }
   else {
      v8.alg_mfe = 65000;
      v8.alg_enum = NULL;
   }
   /*  v8 = el <<< (tt(lbase, lbase) `with` (pairingTTcross compl)) ~~~ (tt(uregion, uregion))  */
   /* ---------------------------------- finished --------------------------------- */


   /* ---------------------------- v9 = minimum(v7, v8) --------------------------- */
   if (v7.alg_mfe < v8.alg_mfe) {
     v9 = v7;
     free_str_Hybrid(v8.alg_enum);
   }
   else {
     v9 = v8;
     free_str_Hybrid(v7.alg_enum);
   }
/*    v9 = v7.alg_mfe < v8.alg_mfe ? v7 : v8; */

   /* --------------------------- v10 = minimum(v5, v9) --------------------------- */
   if (v5.alg_mfe < v9.alg_mfe) {
     v10 = v5;
     free_str_Hybrid(v9.alg_enum);
   }
   else {
     v10 = v9;
     free_str_Hybrid(v5.alg_enum);
   }
/*    v10 = v5.alg_mfe < v9.alg_mfe ? v5 : v9; */

   /* --------------------------- v11 = minimum(v3, v10) -------------------------- */
   if (v3.alg_mfe < v10.alg_mfe) {
     v11 = v3;
     free_str_Hybrid(v10.alg_enum);
   }
   else {
     v11 = v10;
     free_str_Hybrid(v3.alg_enum);
   }
/*    v11 = v3.alg_mfe < v10.alg_mfe ? v3 : v10; */

   /* --------------------------- v12 = minimum(v1, v11) -------------------------- */
   if (v1.alg_mfe < v11.alg_mfe) {
     v12 = v1;
     free_str_Hybrid(v11.alg_enum);
   }
   else {
     v12 = v11;
     free_str_Hybrid(v1.alg_enum);
   }
/*    v12 = v1.alg_mfe < v11.alg_mfe ? v1 : v11; */

   /* ------------------------- build candidate structures ------------------------ */

   return(build_str_Hybrid(v12));
}

/* table memory allocation                                                          */
/* -------------------------------------------------------------------------------- */

void tableAlloc(int m, int n)
{
   int i, j, k, dim1, dim2;

   /* --- memory allocation for tbl_unpaired_left_bot, yield size: ((1,m),(1,n)) -- */
   dim1 = m-1;
   tbl_unpaired_left_bot=(double **) calloc(dim1+1, sizeof(double *));
   for (i=0; i<=dim1; i++) {
      dim2 = n-1;
      tbl_unpaired_left_bot[i]=(double *) calloc(dim2+1, sizeof(double));
   }
   /* -------- memory allocation for tbl_closed, yield size: ((1,m),(1,n)) -------- */
   dim1 = m-1;
   tbl_closed=(double **) calloc(dim1+1, sizeof(double *));
   for (i=0; i<=dim1; i++) {
      dim2 = n-1;
      tbl_closed[i]=(double *) calloc(dim2+1, sizeof(double));
   }
   /* --- memory allocation for tbl_unpaired_left_top, yield size: ((1,m),(1,n)) -- */
   dim1 = m-1;
   tbl_unpaired_left_top=(double **) calloc(dim1+1, sizeof(double *));
   for (i=0; i<=dim1; i++) {
      dim2 = n-1;
      tbl_unpaired_left_top[i]=(double *) calloc(dim2+1, sizeof(double));
   }
}

/* main dynamic programming loop                                                    */
/* -------------------------------------------------------------------------------- */


void mainloop(int bflag, int hit_number, int eflag, float energy_cutoff, int pflag, float pvalue_cutoff, int compact_output, char *target_ac, char *target_sq, char *query_ac, char *query_sq, int gflag, int plot_format, float xi, float theta)
{
   int i1, i2;
   double v2;
   struct str1 l;

   int
     no_nil_so_far,
     iterate,
     hit_count,
     k, k1, k2, sepPos,
     ali_length,
     hit_length,
     t_len;

   int target_has_left_dangle;
   int target_has_right_dangle;


   float
     normalised_energy,
     pvalue;

   char hitcount_str[50];
   char *filename;

   char *conc_seq;
   char *conc_struct;


   no_nil_so_far = 1;
   hit_count = 0;

   for (iterate=1; iterate && (!bflag || hit_count<hit_number); ) {

     for (i1=m; i1>=0; i1--) {
       for (i2=n; i2>=0; i2--) {
	 calc_unpaired_left_bot(i1, m, i2, n);
	 calc_closed           (i1, m, i2, n);
	 calc_unpaired_left_top(i1, m, i2, n);
       }
     }

     /* ----------------------------- show axiom: hybrid ---------------------------- */
     v2 = calc_hybrid(0, m, 0, n);

     l = back_hybrid(0, m, 0, n);

/*      gx = 1; */

     pp_str_Hybrid(l);

     normalised_energy = v2/log(m*n);
     pvalue = 1-exp(-exp(-(-normalised_energy-xi)/theta));

     if ((!eflag || v2 <= energy_cutoff) && (!pflag || pvalue <= pvalue_cutoff) && no_nil_so_far) {

       hit_count++;

       if (compact_output) { /* compact output */
	 printf("%s:%d:%s:%d:%.1f:%.6f:%d:%s:%s:%s:%s\n",target_ac,m,query_ac,n,v2,pvalue,gx,t1,t2,t3,t4);
       }
       else { /* verbose output */
	 printf("target: %s\n", target_ac);
         printf("length: %d\n", m);
	 printf("miRNA : %s\n", query_ac);
         printf("length: %d\n", n);
	 printf("\n");
	 printf("mfe: %.1f kcal/mol\n", v2);
	 printf("p-value: %.6f\n", pvalue);
	 printf("\n");
	 printf("position  %d\n", gx);
	 printf("target 5' %s 3'\n", t1);
	 printf("          %s   \n", t2);
	 printf("          %s   \n", t3);
	 printf("miRNA  3' %s 5'\n", t4);
	 printf("\n");
	 printf("\n");
       }

       if (gflag) { /* plot hybridisation */
	 ali_length = strlen(t1);

	 conc_seq    = (char *) calloc(4 * ali_length +3 + 1,sizeof(char));
	 conc_struct = (char *) calloc(4 * ali_length +3 + 1,sizeof(char));

	 /* make dot-parantheses notation: */

	 k1 = 0;
	 for (k2=0; k2<ali_length; k2++) {
	   if (t1[k2]!=' ') {
	     conc_seq[k1]=t1[k2];
	     conc_struct[k1]='.';
	     k1++;
	   }
	   else if (t2[k2]!=' ') {
	     conc_seq[k1]=t2[k2];
	     conc_struct[k1]='(';
	     k1++;
	   }
	 }
	 conc_seq[k1]='N';
	 conc_seq[k1+1]='N';
	 conc_seq[k1+2]='N';
	 conc_struct[k1]='.';
	 conc_struct[k1+1]='.';
	 conc_struct[k1+2]='.';
	 k1+=3;
	 sepPos = k1;
	 for (k2=ali_length-1; k2>=0; k2--) {
	   if (t4[k2]!=' ') {
	     conc_seq[k1]=t4[k2];
	     conc_struct[k1]='.';
	     k1++;
	   }
	   else if (t3[k2]!=' ') {
	     conc_seq[k1]=t3[k2];
	     conc_struct[k1]=')';
	     k1++;
	   }
	 }

	 conc_seq[k1]    = '\0';
	 conc_struct[k1] = '\0';

/* 	 printf("\n%s\n%s\n",conc_seq,conc_struct); */

	 sprintf(hitcount_str,"%d",hit_count);
	 filename = (char *) calloc(strlen(target_ac)+strlen(query_ac)+strlen(hitcount_str)+7, sizeof(char));

	 if (plot_format==PSPLOT || plot_format==ALLPLOT) {
	   strcpy(filename,target_ac);
	   strcat(filename,"_");
	   strcat(filename,query_ac);
	   strcat(filename,"_");
	   strcat(filename,hitcount_str);
	   strcat(filename,".ps");

#ifdef HAVE_LIBG2	 
	   hybridPlot(conc_seq,conc_struct,v2,sepPos,filename,0,PSPLOT);
#endif
	 }
	 if (plot_format==PNGPLOT || plot_format==ALLPLOT) {
	   strcpy(filename,target_ac);
	   strcat(filename,"_");
	   strcat(filename,query_ac);
	   strcat(filename,"_");
	   strcat(filename,hitcount_str);
	   strcat(filename,".png");

#ifdef HAVE_LIBG2	 
	   hybridPlot(conc_seq,conc_struct,v2,sepPos,filename,0,PNGPLOT);
#endif
	 }
	 if (plot_format==JPGPLOT || plot_format==ALLPLOT) {
	   strcpy(filename,target_ac);
	   strcat(filename,"_");
	   strcat(filename,query_ac);
	   strcat(filename,"_");
	   strcat(filename,hitcount_str);
	   strcat(filename,".jpg");

#ifdef HAVE_LIBG2	 
	   hybridPlot(conc_seq,conc_struct,v2,sepPos,filename,0,JPGPLOT);
#endif
	 }

	 free(filename);
	 free(conc_seq);
	 free(conc_struct);
       }

       t_len = strlen(t1);

       if (t_len==0)
	 no_nil_so_far = 0;

       k = 0;
       while (k < t_len && t1[k]==' ' && t2[k]==' ')
	 k++;
       target_has_left_dangle = k < t_len && t1[k] != ' ';

       k = t_len-1;
       while (k >= 0 && t1[k]==' ' && t2[k]==' ')
	 k--;
       target_has_right_dangle = k >=0 && t1[k] != ' ';

       /* count hit length: */
       hit_length = 0;

       for (k=0; k<t_len; k++)
	 if (t1[k]!=' ')
	   hit_length++;

       for (k=0; k<t_len; k++)
	 if (t2[k]!=' ')
	   hit_length++;

      hit_length -= target_has_left_dangle + target_has_right_dangle;

/*        printf("\n\nhitlength:%d\n\n",hit_length);  */

       /* mask hit: */
       for (k=0; k<hit_length; k++)
/* 	 x[gx+(t1[0]!=' ')+k]=X; */
	 x[gx+target_has_left_dangle+k]=X;


/*        for (k=1; k<=m; k++) */
/* 	 printf("%d",x[k]); */
/*        printf("\n\n"); */

       if (!eflag && !bflag)
	 iterate = 0;

     }
     else
       iterate = 0;


     /* free candidate */
     free_str_Hybrid(l.alg_enum);
     /* initialize strings */
     r1 = t1;
     r2 = t2;
     r3 = t3;
     r4 = t4;
   }
}



