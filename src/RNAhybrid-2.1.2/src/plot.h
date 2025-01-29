/* Copyright (C) 2004 Marc Rehmsmeier, Peter Steffen, Matthias Hoechsmann */

/* This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version. */

/* This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. */

/* You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA */

#ifndef plot_h
#define plot_h

#include "config.h"

#define  PSPLOT 0
#define PNGPLOT 1
#define JPGPLOT 2
#define ALLPLOT 3



#ifdef HAVE_LIBG2

#include <stdio.h>  
#include <g2.h>
#include <g2_PS.h>
#include "minmax.h"



void hybridPlot(char *seq, char *str, double energy, int sepPos, char *filename, int blackWhite, int format);


#endif

#endif
