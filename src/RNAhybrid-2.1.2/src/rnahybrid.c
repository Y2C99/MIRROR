/* Copyright (C) 2004 Marc Rehmsmeier, Peter Steffen, Matthias Hoechsmann */

/* This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version. */

/* This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. */

/* You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA */

#include <stdio.h>  
#include <string.h> 
#include "config.h"
#include "plot.h"
#include "minmax.h"
#include "fasta.h"
#include "input.h"
#include "energy.h"
#include "globals.h"
#include "hybrid_core.h"


static char const rcsid[] = "$Id: rnahybrid.c,v 1.1.1.9.1.13 2004/11/15 20:41:38 marc Exp $";


float maximal_duplex_energy(char *seq);


int main(int argc, char **argv)                
{       
  int i;
 
  extern char *optarg;
  extern int optind;
  int c;
  int
    bflag = 0, // number of hits per target to show
    cflag = 0, // compact output
    dflag = 0, // parameters of extreme value distribution
               // for p-value calculation
    eflag = 0, // report all hits with energy better (ie. lower) than or equal to e
    fflag = 0, // force helix, arguments: from,to (positions in miRNA)
    gflag = 0, // make plots (graphics) of hybridisations
    hflag = 0, // help
    mflag = 0, // max target length
    nflag = 0, // max query length
    pflag = 0, // report all hits with p-values better than or equal to p
    sflag = 0, // data set flag, for organism and sequence specific settings
    qflag = 0, // query  input file
    tflag = 0, // target input file
    uflag = 0, // upper size for internal loops (per side)
    vflag = 0, // upper size for bulge loops
    errflag = 0;
  int
    maxtargetlength,
    maxquerylength;
  char
    *target_fn, // target filename
    *query_fn;  // query  filename
  float
    energy_cutoff = 0.0,
    pvalue_cutoff = 1.0,
    xi, theta,
    mde;

  int
    hit_number = 1;

  char target_ac[MAXLINE];
  char query_ac[MAXLINE];

  char *target_sq;
  char *query_sq;

  char *gflag_argument;
  char *plotfileextension;

  char *setname;

  int
    targetlength,
    querylength,
    plot_format;

  int
    fflag_start,
    fflag_end;

  float xi_slope, xi_intercept, theta_slope, theta_intercept;


  while ((c = getopt(argc,argv,"b:cd:e:f:g:hm:n:p:q:s:t:u:v:")) != EOF)
    switch(c) {
    case 'b':
      bflag=1;
      sscanf(optarg,"%d",&hit_number);
      break;
    case 'c':
      cflag=1;
      break;
    case 'd':
      dflag=1;
      sscanf(optarg,"%f,%f",&xi,&theta);
      break;
    case 'e':
      eflag=1;
      sscanf(optarg,"%f",&energy_cutoff);
      break;
    case 'f':
      fflag=1;
      sscanf(optarg,"%d,%d",&fflag_start,&fflag_end);
      break;
    case 'g':
      gflag=1;
      gflag_argument = optarg;
      break;
    case 'h':
      hflag=1;
      break;
    case 'm':
      mflag=1;
      sscanf(optarg,"%d",&maxtargetlength);
      break;
    case 'n':
      nflag=1;
      sscanf(optarg,"%d",&maxquerylength);
      break;
    case 'p':
      pflag=1;
      sscanf(optarg,"%f",&pvalue_cutoff);
      break;
    case 'q':
      qflag = 1;
      query_fn = optarg;
      break;
    case 's':
      sflag=1;
      setname = optarg;
      break;
    case 't':
      tflag = 1;
      target_fn = optarg;
      break;
    case 'u':
      uflag = 1;
      sscanf(optarg,"%d",&iloop_upper_limit);
      break;
    case 'v':
      vflag = 1;
      sscanf(optarg,"%d",&bloop_upper_limit);
      break;
    case '?':
      errflag = 1;
      break;
    }

  if (eflag && (energy_cutoff > 0)) {
    printf("Warning: energy cut-off positive. I'm converting it to minus the value you gave.\n");
    energy_cutoff = -energy_cutoff;
  }

  if (!uflag)
    iloop_upper_limit = ILOOPUPPERLIMITDEFAULT;
  if (!vflag)
    bloop_upper_limit = BLOOPUPPERLIMITDEFAULT;


#ifdef HAVE_LIBG2
  if (gflag) {
    if (strcmp(gflag_argument,"ps")==0)
      plot_format = PSPLOT;
#ifdef HAVE_LIBGD
    else if (strcmp(gflag_argument,"png")==0)
      plot_format = PNGPLOT;
    else if (strcmp(gflag_argument,"jpg")==0)
      plot_format = JPGPLOT;
    else if (strcmp(gflag_argument,"all")==0)
      plot_format = ALLPLOT;
#endif
    else {
      printf("\nWrong plot format.\n");
      errflag = 1;
    }
  }
#endif

  if (errflag || argc < 3 && !hflag) {
    printf("\nOption error. Type %s -h for usage.\n\n", argv[0]);
    exit(1);
  }
  else if (hflag) {
    printf("\nUsage: %s [options] [target sequence] [query sequence].\n\noptions:\n\n  -b <number of hits per target>\n  -c compact output\n  -d <xi>,<theta>\n  -f helix constraint\n  -h help\n  -m <max targetlength>\n  -n <max query length>\n  -u <max internal loop size (per side)>\n  -v <max bulge loop size>\n  -e <energy cut-off>\n  -p <p-value cut-off>\n  -s (3utr_fly|3utr_worm|3utr_human)\n  -g (ps|png|jpg|all)\n  -t <target file>\n  -q <query file>\n\nEither a target file has to be given (FASTA format)\nor one target sequence directly.\n\nEither a query file has to be given (FASTA format)\nor one query sequence directly.\n\nThe helix constraint format is \"from,to\", eg. -f 2,7 forces\nstructures to have a helix from position 2 to 7 with respect to the query.\n\n<xi> and <theta> are the position and shape parameters, respectively,\nof the extreme value distribution assumed for p-value calculation.\nIf omitted, they are estimated from the maximal duplex energy of the query.\nIn that case, a data set name has to be given with the -s flag.\n\n", argv[0]);

#ifndef HAVE_LIBG2
    printf("\nPS graphical output not supported.\n\n");
#endif

#ifndef HAVE_LIBGD
    printf("\nPNG and JPG graphical output not supported.\n\n");
#endif
    exit(0);
  }

  if (!dflag && !sflag) {
    printf("Error: without -d you have to give -s option. Aborting.\n");
    exit(2);
  }

  if (sflag) {
    if (strcmp(setname,SETNAME_3UTR_FLY) == 0) {
      xi_slope = XI_SLOPE_3UTR_FLY;
      xi_intercept = XI_INTERCEPT_3UTR_FLY;
      theta_slope = THETA_SLOPE_3UTR_FLY;
      theta_intercept = THETA_INTERCEPT_3UTR_FLY;
    }
    else if (strcmp(setname,SETNAME_3UTR_WORM) == 0) {
      xi_slope = XI_SLOPE_3UTR_WORM;
      xi_intercept = XI_INTERCEPT_3UTR_WORM;
      theta_slope = THETA_SLOPE_3UTR_WORM;
      theta_intercept = THETA_INTERCEPT_3UTR_WORM;
    }
    else if (strcmp(setname,SETNAME_3UTR_HUMAN) == 0) {
      xi_slope = XI_SLOPE_3UTR_HUMAN;
      xi_intercept = XI_INTERCEPT_3UTR_HUMAN;
      theta_slope = THETA_SLOPE_3UTR_HUMAN;
      theta_intercept = THETA_INTERCEPT_3UTR_HUMAN;
    }
    else {
      printf("Error: unknown set name. Aborting.\n");
      exit(2);
    }
  }

  if (!mflag)
    maxtargetlength = MAXTARGET;
  if (!nflag)
    maxquerylength = MAXQUERY;

  if (tflag && qflag) {
    target_sq = (char *) calloc(maxtargetlength+MAXLINE+1,sizeof(char));
    query_sq = (char *) calloc(maxquerylength+MAXLINE+1, sizeof(char));
    targetlength = maxtargetlength;
    querylength = maxquerylength;
  }

  if (tflag && !qflag) {
    target_sq = (char *) calloc(maxtargetlength+MAXLINE+1,sizeof(char));
    query_sq  = (char *) calloc(strlen(argv[argc-1])+1,sizeof(char));
    strcpy(query_ac,"command_line");
    strcpy(query_sq,argv[argc-1]);
    remove_whitespace(query_sq);
    targetlength = maxtargetlength;
    querylength = strlen(argv[argc-1]);
  }

  if (!tflag && qflag) {
    target_sq = (char *) calloc(strlen(argv[argc-1])+1,sizeof(char));
    query_sq = (char *) calloc(maxquerylength+MAXLINE+1, sizeof(char));
    strcpy(target_ac,"command_line");
    strcpy(target_sq,argv[argc-1]);
    remove_whitespace(target_sq);
    targetlength = strlen(argv[argc-1]);
    querylength = maxquerylength;
  }

  if (!tflag && !qflag) {
    target_sq = (char *) calloc(strlen(argv[argc-2])+1,sizeof(char));
    query_sq  = (char *) calloc(strlen(argv[argc-1])+1,sizeof(char));
    strcpy(target_ac,"command_line");
    strcpy(query_ac,"command_line");
    strcpy(target_sq,argv[argc-2]);
    remove_whitespace(target_sq);
    strcpy(query_sq,argv[argc-1]);
    remove_whitespace(query_sq);
    targetlength = strlen(argv[argc-2]);
    querylength = strlen(argv[argc-1]);
  }
                                             
  init_constants();                                            
  init_energies();

  tableAlloc(targetlength,querylength);

  /* allocate string space */
  r1 = t1 = (char *) calloc(2*max(targetlength,querylength), sizeof(char));
  r2 = t2 = (char *) calloc(2*max(targetlength,querylength), sizeof(char));
  r3 = t3 = (char *) calloc(2*max(targetlength,querylength), sizeof(char));
  r4 = t4 = (char *) calloc(2*max(targetlength,querylength), sizeof(char));

  x = (char *) calloc(targetlength+2, sizeof(char));     
  y = (char *) calloc( querylength+2, sizeof(char));     
 

  if (qflag) {

    FILE *query  = fopen( query_fn,"r");
    if (query==NULL) {
      printf("Error: Could not open query file. Aborting.\n");
      exit(2);
    }

    while (!end(query)) {

      nextAC(query,query_ac);
      nextSQ(query,query_sq,maxquerylength);

      remove_whitespace(query_sq);
    
      n = strlen(query_sq);

      if (!dflag) {
	// estimate evd parameters from maximal duplex energy
	mde = maximal_duplex_energy(query_sq);
        xi    =    xi_slope * mde +    xi_intercept;
        theta = theta_slope * mde + theta_intercept;
      }

      if (fflag) {
	helix_start = n - fflag_end;
	helix_end   = n - fflag_start + 1;
      }
      else {
	helix_start = 0;
	helix_end   = 0;
      }

      if (n > maxquerylength) {
	printf("query too long: %s\n", query_ac);
      }
      else {

	y[0]=' ';
	for (i=0;i<=n-1;i++) y[i+1] = query_sq[n-i-1];
	y[n+1]=0;
	convert_y();

	if (tflag) {
	  FILE *target = fopen(target_fn,"r");
	  if (target==NULL) {
	    printf("Error: Could not open target file. Aborting.\n");
	    exit(2);
	  }

	  while (!end(target)) {

	    nextAC(target,target_ac);
	    nextSQ(target,target_sq,maxtargetlength);

	    remove_whitespace(target_sq);

	    m = strlen(target_sq);

	    if (m > maxtargetlength) {
	      printf("target too long: %s\n", target_ac);
	    }
	    else {
	      strcpy(x, " ");                             
	      strcat(x, target_sq);
	      convert_x();

	      mainloop(bflag, hit_number,eflag,energy_cutoff,pflag,pvalue_cutoff,cflag,target_ac,target_sq,query_ac,query_sq,gflag,plot_format,xi,theta);
	    }
	    
	  }

	  fclose(target);
	}
	else {
	  m = strlen(target_sq);

	  strcpy(x, " ");                             
	  strcat(x, target_sq);
	  convert_x();

	  mainloop(bflag, hit_number,eflag,energy_cutoff,pflag,pvalue_cutoff,cflag,target_ac,target_sq,query_ac,query_sq,gflag,plot_format,xi,theta);
	}
      }
    }

    fclose(query);
  }
  else { /* !qflag */
    n = strlen(query_sq);

    if (!dflag) {
      // estimate evd parameters from maximal duplex energy
      mde = maximal_duplex_energy(query_sq);
      xi    =    xi_slope * mde +    xi_intercept;
      theta = theta_slope * mde + theta_intercept;
    }

    if (fflag) {
      helix_start = n - fflag_end;
      helix_end   = n - fflag_start + 1;
    }
    else {
      helix_start = 0;
      helix_end   = 0;
    }


    y[0]=' ';
    for (i=0;i<=n-1;i++) y[i+1] = query_sq[n-i-1];
    y[n+1]=0;
    convert_y();

    if (tflag) { /* !qflag && tflag */
      FILE *target = fopen(target_fn,"r");

      while (!end(target)) {

	nextAC(target,target_ac);
	nextSQ(target,target_sq,maxtargetlength);

	remove_whitespace(target_sq);

	m = strlen(target_sq);

	if (m > maxtargetlength) {
	  printf("target too long: %s\n", target_ac);
	}
	else {
	  strcpy(x, " ");                             
	  strcat(x, target_sq);
	  convert_x();

	  mainloop(bflag,hit_number,eflag,energy_cutoff,pflag,pvalue_cutoff,cflag,target_ac,target_sq,query_ac,query_sq,gflag,plot_format,xi,theta);
	}

      }

      fclose(target);
    }
    else { /* !qflag && !tflag */
      m = strlen(target_sq);

      strcpy(x, " ");                             
      strcat(x, target_sq);
      convert_x();

      mainloop(bflag, hit_number,eflag,energy_cutoff,pflag,pvalue_cutoff,cflag,target_ac,target_sq,query_ac,query_sq,gflag,plot_format,xi,theta);
    }
  }


   exit(0);                                    
}



float maximal_duplex_energy(char *seq_arg)
{
  int seq_len = strlen(seq_arg);
  char *seq = (char *) calloc(seq_len+1, sizeof(char));
  char *rc  = (char *) calloc(seq_len+1, sizeof(char));
  int i;
  char letter, compl_letter;
  float mde;
  char c;

  // convert sequence:
  for (i=0; i<seq_len; i++) {
    c=seq_arg[i];
    if      (c=='a' || c=='A') seq[i]=A;
    else if (c=='c' || c=='C') seq[i]=C;
    else if (c=='g' || c=='G') seq[i]=G;
    else if (c=='u' || c=='U') seq[i]=U;
    else if (c=='t' || c=='T') seq[i]=U;
    else                       seq[i]=N;
  }

  // make reverse complement:
  for (i=0; i<seq_len; i++) {
    c = seq[seq_len-i-1];
    if      (c==A) rc[i]=U;
    else if (c==C) rc[i]=G;
    else if (c==G) rc[i]=C;
    else if (c==U) rc[i]=A;
    else           rc[i]=N;
  }
  rc[seq_len] = '\0';


  // sum up stacked pair energies:
  mde = 0;
  for (i=0; i<seq_len-1; i++)
    mde += stack_dg_ar[rc[seq_len-i-1]][rc[seq_len-(i+1)-1]][seq[i+1]][seq[i]];

  free(seq);
  free(rc);

  return mde;
}

