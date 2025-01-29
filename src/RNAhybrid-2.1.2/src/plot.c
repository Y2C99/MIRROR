/* Copyright (C) 2004 Marc Rehmsmeier, Peter Steffen, Matthias Hoechsmann */

/* This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version. */

/* This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. */

/* You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA */

#include "config.h"
#include "plot.h"

#ifdef HAVE_LIBG2

#ifdef HAVE_LIBGD
#include <g2_gd.h>
#endif

void   loop(int i, int j, short *pair_table);

short *make_pair_table(const char *structure);

void *space(unsigned size);

void nrerror(const char message[]);

/* local variables for parsing routines */

float  *angle;
int    *loop_size, *stack_size;
int     lp, stk;

#ifndef PI
#define  PI       3.141592654
#endif
#define  PIHALF       PI/2.



void hybridPlot(char *seq, char *str, double energy, int sepPos, char *filename, int blackWhite, int format)
{
  double base_fontsize=12;
  float *X, *Y,min_X=0,max_X=0,min_Y=0,max_Y=0;
  unsigned int i;
  short *pair_table;
  int id_PS,id;
  int ps_color_black,ps_color_red,ps_color_blue;
  char energy_str[50];

#ifdef HAVE_LIBGD
  int id_PNG=0,id_JPG=0;
  int png_color_black,png_color_red,png_color_blue;
  int jpg_color_black,jpg_color_red,jpg_color_blue;
#endif

  unsigned int basenr_x=0, basenr_y=0;
  double xpos,ypos;
  char buf[2];

  buf[1]=0;
  
  //  assert(strlen(str) == strlen(seq));

  X = (float *) calloc(strlen(seq)+1,sizeof(float));
  Y = (float *) calloc(strlen(seq)+1,sizeof(float));

  pair_table = make_pair_table(str);
  i = simple_xy_coordinates(pair_table, X, Y);
  if(i!=strlen(str))
    printf("strange things happening in squigglePlot ...");
    
  // scale image
  //  for(i=0;i<strlen(str);i++)
  //    {
  //      X[i]*=options.scale;
  //      Y[i]*=options.scale;
  //    }  

  // calculate image dimensions
  for(i=0;i<strlen(str);i++)
    {
      min_X=min(min_X,X[i]);
      max_X=max(max_X,X[i]);
      min_Y=min(min_Y,Y[i]);
      max_Y=max(max_Y,Y[i]);
    }

  // add a border to image size
  min_X-=10;
  max_X+=10;
  min_Y-=10;
  max_Y+=10;

  // open device with respect to "format"
  id     = g2_open_vd();

  if(format==PSPLOT)
    {
      id_PS  = g2_open_EPSF(filename);
      g2_attach(id,id_PS);
    }

#ifdef HAVE_LIBGD
  if(format==PNGPLOT)
    {
      id_PNG=g2_open_gd(filename,(int)(max_X-min_X),(int)(max_Y-min_Y),g2_gd_png);
      g2_attach(id,id_PNG);
    }
  if(format==JPGPLOT)
    {
      id_JPG=g2_open_gd(filename,(int)(max_X-min_X),(int)(max_Y-min_Y),g2_gd_jpeg);
      g2_attach(id,id_JPG);
    }
#endif  
  
  g2_set_coordinate_system(id,-min_X,-min_Y,1,1);

  // drawing parameters
  g2_set_line_width(id,0.6);

  // mark 5' end
  g2_string(id,X[0]-20,Y[0],"5'");

  // print energy:
  sprintf(energy_str,"mfe: %.1f kcal/mol",energy);
  g2_string(id,min_X+40,min_Y+20,energy_str);

  // define colors
  if(blackWhite)
    {
      switch(format)
	{
	case PSPLOT:
	  ps_color_black=g2_ink(id_PS,0,0,0);
	  ps_color_red=g2_ink(id_PS,0,0,0);
	  ps_color_blue=g2_ink(id_PS,0,0,0);
	  break;
#ifdef HAVE_LIBGD
	case PNGPLOT:
	  png_color_black=g2_ink(id_PNG,0,0,0);
	  png_color_red=g2_ink(id_PNG,0,0,0);
	  png_color_blue=g2_ink(id_PNG,0,0,0);

	  ps_color_black=png_color_black;
	  ps_color_red=png_color_red;
	  ps_color_blue=png_color_blue;

	  break;
	case JPGPLOT:
	  jpg_color_black=g2_ink(id_JPG,0,0,0);
	  jpg_color_red=g2_ink(id_JPG,0,0,0);
	  jpg_color_blue=g2_ink(id_JPG,0,0,0);

	  ps_color_black=jpg_color_black;
	  ps_color_red=jpg_color_red;
	  ps_color_blue=jpg_color_blue;

	  break;  
#endif
	}
    }
  else
    {
      switch(format)
	{
	case PSPLOT:
	  ps_color_black=g2_ink(id_PS,0,0,0);
	  ps_color_red=g2_ink(id_PS,1,0,0);
	  ps_color_blue=g2_ink(id_PS,0,0.75,0);
	  break;
#ifdef HAVE_LIBGD
	case PNGPLOT:
	  png_color_black=g2_ink(id_PNG,0,0,0);
	  png_color_red=g2_ink(id_PNG,1,0,0);
	  png_color_blue=g2_ink(id_PNG,0,0.75,0);

	  ps_color_black=png_color_black;
	  ps_color_red=png_color_red;
	  ps_color_blue=png_color_blue;

	  break;
	case JPGPLOT:
	  jpg_color_black=g2_ink(id_JPG,0,0,0);
	  jpg_color_red=g2_ink(id_JPG,1,0,0);
	  jpg_color_blue=g2_ink(id_JPG,0,0.75,0);

	  ps_color_black=jpg_color_black;
	  ps_color_red=jpg_color_red;
	  ps_color_blue=jpg_color_blue;

	  break;
#endif
	}
    }

  // draw sequence
  g2_set_font_size(id,base_fontsize);
    for(i=0;i<strlen(str);i++)
   {
     if(i<sepPos)
       g2_pen(id,ps_color_red);
     else
       g2_pen(id,ps_color_blue);
         
     buf[0]=seq[i];
     xpos=X[i]-base_fontsize/2;
     ypos=Y[i]-4;
     if (i<=sepPos-4 || i > sepPos-1)
       g2_string(id,xpos,ypos,buf);
     
     // connection to next base
     if(i<strlen(str)-1)
       {
	 if(i<sepPos-4 || i>sepPos-1)
	   {
	     if(i<sepPos)
	       g2_pen(id,ps_color_red);
	     else
	       g2_pen(id,ps_color_blue);
	     g2_line(id,X[i],Y[i],X[i+1],Y[i+1]);
	   }

	 // draw circles at line endpoints
/* 	 if(i<sepPos) */
/* 	   g2_pen(id,ps_color_red); */
/* 	 else */
/* 	   g2_pen(id,ps_color_blue); */
/* 	 g2_filled_circle(id,X[i],Y[i],0.7);       // circles are drawn twice, but thats ok ... */
       }
   }
   
    // draw pairings
    // !!! pair_table indexing begins at 1 !!!
    for(i=0;i<strlen(str);i++)
      {
	if((unsigned short)pair_table[i+1]>i+1)
	  {
	    // pairs in both structures
	    g2_pen(id,ps_color_black);
	    g2_pen(id,ps_color_black);	   	    
	    g2_line(id,X[i],Y[i],X[pair_table[i+1]-1],Y[pair_table[i+1]-1]);
	  }
      }


    g2_flush(id);
    g2_close(id);

    free(pair_table);
    free(X);
    free(Y);
}



short *make_pair_table(const char *structure)
{
    /* returns array representation of structure.
       table[i] is 0 if unpaired or j if (i.j) pair.  */
   short i,j,hx;
   short length;
   short *stack;
   short *table;
   
   length = (short) strlen(structure);
   stack = (short *) space(sizeof(short)*(length+1));
   table = (short *) space(sizeof(short)*(length+2));
   table[0] = length;
   
   for (hx=0, i=1; i<=length; i++) {
      switch (structure[i-1]) {
       case '(': 
	 stack[hx++]=i;
	 break;
       case ')':
	 j = stack[--hx];
	 if (hx<0) {
	    fprintf(stderr, "%s\n", structure);
	    nrerror("unbalanced brackets in make_pair_table");
	 }
	 table[i]=j;
	 table[j]=i;
	 break;
       default:   /* unpaired base, usually '.' */
	 table[i]= 0;
	 break;
      }
   }
   if (hx!=0) {
      fprintf(stderr, "%s\n", structure);
      nrerror("unbalanced brackets in make_pair_table");
   }
   free(stack);
   return(table);
}


void *space(unsigned size)
{
  void *pointer;
  
  if ( (pointer = (void *) calloc(1, (size_t) size)) == NULL) {
      nrerror("SPACE allocation failure -> no memory");
  }
  return  pointer;
}


void nrerror(const char message[])       /* output message upon error */
{
  fprintf(stderr, "\n%s\n", message);
  exit(2);
}


int simple_xy_coordinates(short *pair_table, float *x, float *y)
{
  float INIT_ANGLE=0.;     /* initial bending angle */
  float INIT_X = 100.;     /* coordinate of first digit */
  float INIT_Y = 100.;     /* see above */
  float RADIUS =  15.;

  int i, length;
  float  alpha;

  length = pair_table[0];
  angle =      (float*) space( (length+5)*sizeof(float) );
  loop_size  =   (int*) space( 16+(length/5)*sizeof(int) );
  stack_size =   (int*) space( 16+(length/5)*sizeof(int) );
  lp = stk = 0;
  loop(0, length+1, pair_table);
  loop_size[lp] -= 2;     /* correct for cheating with function loop */
  
  alpha = INIT_ANGLE;
  x[0]  = INIT_X;
  y[0]  = INIT_Y;   

  for (i = 1; i <= length; i++) {
    x[i] = x[i-1]+RADIUS*cos(alpha);
    y[i] = y[i-1]+RADIUS*sin(alpha);
    alpha += PI-angle[i+1];
  }
  free(angle);
  free(loop_size);  
  free(stack_size); 

  return length;

}

/*---------------------------------------------------------------------------*/

void loop(int i, int j, short *pair_table)
             /* i, j are the positions AFTER the last pair of a stack; i.e
		i-1 and j+1 are paired. */
{
  int    count = 2;   /* counts the VERTICES of a loop polygon; that's
			   NOT necessarily the number of unpaired bases!
			   Upon entry the loop has already 2 vertices, namely
			   the pair i-1/j+1.  */

  int    r = 0, bubble = 0; /* bubble counts the unpaired digits in loops */

  int    i_old, partner, k, l, start_k, start_l, fill, ladder;
  int    begin, v, diff;
  float  polygon;

  short *remember;  

  remember = (short *) space((1+(j-i)/5)*2*sizeof(short));
    
  i_old = i-1, j++;         /* j has now been set to the partner of the
			       previous pair for correct while-loop
			       termination.  */
  while (i != j) {
    partner = pair_table[i];
    if ((!partner) || (i==0))
      i++, count++, bubble++;
    else {
      count += 2;
      k = i, l = partner;    /* beginning of stack */
      remember[++r] = k;
      remember[++r] = l;
      i = partner+1;         /* next i for the current loop */

      start_k = k, start_l = l;
      ladder = 0;
      do {
	k++, l--, ladder++;        /* go along the stack region */
      }
      while (pair_table[k] == l);

      fill = ladder-2;
      if (ladder >= 2) {
	angle[start_k+1+fill] += PIHALF;   /*  Loop entries and    */
	angle[start_l-1-fill] += PIHALF;   /*  exits get an        */
	angle[start_k]        += PIHALF;   /*  additional PI/2.    */
	angle[start_l]        += PIHALF;   /*  Why ? (exercise)    */
	if (ladder > 2) {
	  for (; fill >= 1; fill--) {
	    angle[start_k+fill] = PI;    /*  fill in the angles  */
	    angle[start_l-fill] = PI;    /*  for the backbone    */
	  }
	}
      }
      stack_size[++stk] = ladder;
      loop(k, l, pair_table);
    }
  }
  polygon = PI*(count-2)/(float)count; /* bending angle in loop polygon */
  remember[++r] = j;
  begin = i_old < 0 ? 0 : i_old;
  for (v = 1; v <= r; v++) {
    diff  = remember[v]-begin;
    for (fill = 0; fill <= diff; fill++)
      angle[begin+fill] += polygon;
    if (v > r)
      break;
    begin = remember[++v];
  }
  loop_size[++lp] = bubble;
  free(remember);
}

#endif
