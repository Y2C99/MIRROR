/* Copyright (C) 2004 Marc Rehmsmeier, Peter Steffen, Matthias Hoechsmann */

/* This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version. */

/* This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. */

/* You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA */

#ifndef globals_h
#define globals_h


#define MAXTARGET     2000
#define MAXQUERY        30
#define MAXLINE       1000
#define MAXSEQUENCE 100000

#define ILOOPUPPERLIMITDEFAULT 15
#define BLOOPUPPERLIMITDEFAULT 15

#define MEAN 500
#define STDDEV 300
#define SAMPLESIZE 5000

#define MIRNA_LENGTH_MEAN 22
#define MIRNA_LENGTH_STDDEV 0

#define MAXTARGETNUMBER 100000
#define MAXQUERYNUMBER 100000

#define MAXORTHOLOGNUMBER 10

#define BINNUMBER 500
#define FITLOWERCUTOFFABSOLUTE 2.0
#define FITUPPERCUTOFFPERCENT 1.0

#define SETNAME_3UTR_FLY   "3utr_fly"
#define SETNAME_3UTR_WORM  "3utr_worm"
#define SETNAME_3UTR_HUMAN "3utr_human"

#define XI_SLOPE_3UTR_FLY (-0.03144)
#define XI_INTERCEPT_3UTR_FLY (0.7201)
#define THETA_SLOPE_3UTR_FLY (-0.003429)
#define THETA_INTERCEPT_3UTR_FLY (0.01634)

#define XI_SLOPE_3UTR_WORM (-0.01074)
#define XI_INTERCEPT_3UTR_WORM (1.769)
#define THETA_SLOPE_3UTR_WORM (-0.001154)
#define THETA_INTERCEPT_3UTR_WORM (0.1419)

#define XI_SLOPE_3UTR_HUMAN (-0.01237)
#define XI_INTERCEPT_3UTR_HUMAN (1.901)
#define THETA_SLOPE_3UTR_HUMAN (-0.001678)
#define THETA_INTERCEPT_3UTR_HUMAN (0.1349)


#endif
