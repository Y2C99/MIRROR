/* Copyright (C) 2004 Marc Rehmsmeier, Peter Steffen, Matthias Hoechsmann */

/* This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version. */

/* This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. */

/* You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA */

#include "hybrid_core.h"
#include "numerical.h"
#include "random.h"
#include "globals.h"
#include "minmax.h"
#include "mt19937-1.h"
#include <stdio.h>
#include <sys/types.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>


static char const rcsid[] = "$Id: rnaeffective.c,v 1.6 2004/11/16 15:41:12 marc Exp $";



int main(int argc, char **argv)
{
  extern char *optarg;
  extern int optind;
  int c;
  int
    dflag = 0, // file with dinucleotide parameters
    fflag = 0, // force helix, arguments: from,to (positions in miRNA)
    hflag = 0, // help
    kflag = 0, // number of sequences to generate
    lflag = 0, // length distribution parameters (mean, std deviation)
    mflag = 0, // max target length
    nflag = 0, // max query length
    qflag = 0, // query  input file
    sflag = 0, // make random sequences according to target file
    tflag = 0, // target input file
    uflag = 0, // upper size for internal loops (per side)
    vflag = 0, // upper size for bulge loops
    errflag = 0;

  int
    maxtargetlength,
    maxquerylength;
  char
    *query_fn,  // query filename
    *target_fn;

  char target_ac[MAXLINE];
  char query_ac[MAXLINE];

  char *target_sq;
  char *query_sq;

  char *freq_filename;

  int i1, i2;
  double v2, normalised_energy;

  int
    fflag_start,
    fflag_end;
  int
    targetlength,
    querylength;

  int mean, stddev, sample_size, used_sample_size;

  float pvalue;

  float *xi, *theta;

  FILE *f;

  int i,j,k,l,seqlen;

  char letter_one, letter_two;

  int index_one, index_two;

  int counter;

  float exponent, effective_number;

  float *xcdf, *ycdf, bin_number, error_sum, smallest_error_sum;

  

  float **normalised_energies;
  float *temp_energies;
  float *maximal_pvalues;
  float *joint_pvalues;

  float **freq_di = (float **) calloc(4,sizeof(float*));
  float **fdf_di = (float **) calloc(4,sizeof(float*));
  float *freq;
  float *fdf;

  for (i=0; i<4; i++) {
    freq_di[i] = (float *) calloc(4,sizeof(float));
    fdf_di[i]  = (float *) calloc(4,sizeof(float));
  }
  freq = (float *) calloc(4,sizeof(float));
  fdf  = (float *) calloc(4,sizeof(float));



  while ((c = getopt(argc,argv,"d:f:hk:l:m:n:q:st:u:v:")) != EOF)
    switch(c) {
    case 'd':
      dflag = 1;
      freq_filename = optarg;
      break;
    case 'f':
      fflag=1;
      sscanf(optarg,"%d,%d",&fflag_start,&fflag_end);
      break;
    case 'h':
      hflag=1;
      break;
    case 'k':
      kflag = 1;
      sscanf(optarg,"%d",&sample_size);
      break;
    case 'l':
      lflag = 1;
      sscanf(optarg,"%d,%d",&mean,&stddev);
      break;
    case 'm':
      mflag=1;
      sscanf(optarg,"%d",&maxtargetlength);
      break;
    case 'n':
      nflag=1;
      sscanf(optarg,"%d",&maxquerylength);
      break;
    case 'q':
      qflag = 1;
      query_fn = optarg;
      break;
    case 's':
      sflag = 1;
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

  if ((errflag || sflag && dflag || argc < 3) && !hflag) {
    printf("\nOption error. Type %s -h for usage.\n\n", argv[0]);
    exit(1);
  }
  else if (hflag) {
    printf("\nUsage: %s [options] [query sequence].\n\noptions:\n\n  -d <dinucleotide frequencies file>\n  -f helix constraint\n  -h help\n  -k <sample size>\n  -l <mean sequence length>,<std sequence length>\n  -m <max targetlength>\n  -n <max query length>\n  -u <max internal loop size (per side)>\n  -v <max bulge loop size>\n  -s randomise queries (only with -q, without -d)\n  -t <target file>\n  -q <query file>\n\nIf no query file is given, random sequences are generated according\nto the given dinucleotide distribution. Then <sample size> sequences\nare generated whose lengths are normally distributed with given mean\nand standard deviation. Default sample size is 5000, default mean\nand std are 22 and 0, respectively.\n\nIf a query file is given, and additionally the -s option, random sequences\nare generated according to the dinucleotide distribution of the query file.\n\nIf only a query file is given, it is used directly as a random database.\n\nThe query can also be given directly (makes only sense with the -s option).\n\nA target file has to be given (FASTA format).\n\n\nThe helix constraint format is \"from,to\", eg. -f 2,7 forces\nstructures to have a helix from position 2 to 7 with respect to the query.\n\n", argv[0]);
    exit(0);
  }

  if (!uflag)
    iloop_upper_limit = ILOOPUPPERLIMITDEFAULT;
  if (!vflag)
    bloop_upper_limit = BLOOPUPPERLIMITDEFAULT;

  if (!lflag) {
    mean = MIRNA_LENGTH_MEAN;
    stddev = MIRNA_LENGTH_STDDEV;
  }

  if (!kflag)
    sample_size = SAMPLESIZE;

  normalised_energies = (float **) calloc(MAXORTHOLOGNUMBER,sizeof(float *));
  xi = (float *) calloc(MAXORTHOLOGNUMBER,sizeof(float));
  theta = (float *) calloc(MAXORTHOLOGNUMBER,sizeof(float));
  
  if (qflag && !sflag) {
    for (i=0; i<MAXORTHOLOGNUMBER; i++)
      normalised_energies[i] = (float *) calloc(MAXQUERYNUMBER,sizeof(float));
    temp_energies = (float *) calloc(MAXQUERYNUMBER,sizeof(float));
    maximal_pvalues = (float *) calloc(MAXQUERYNUMBER,sizeof(float));
    joint_pvalues = (float *) calloc(MAXQUERYNUMBER,sizeof(float));
  }
  else {
    for (i=0; i<MAXORTHOLOGNUMBER; i++)
      normalised_energies[i] = (float *) calloc(sample_size,sizeof(float));
    temp_energies = (float *) calloc(sample_size,sizeof(float));
    maximal_pvalues = (float *) calloc(sample_size,sizeof(float));
    joint_pvalues = (float *) calloc(sample_size,sizeof(float));
  }

  if (!mflag)
    maxtargetlength = MAXTARGET;
  if (!nflag)
    maxquerylength = MAXQUERY;

  if (qflag) {
    query_sq = (char *) calloc(maxquerylength+MAXLINE+1, sizeof(char));
    querylength = maxquerylength;
  }
  else {
    query_sq  = (char *) calloc(strlen(argv[argc-1])+1,sizeof(char));
    strcpy(query_ac,"command_line");
    strcpy(query_sq,argv[argc-1]);
    remove_whitespace(query_sq);
    querylength = strlen(argv[argc-1]);
  }

  if (tflag) {
    target_sq = (char *) calloc(maxtargetlength+MAXLINE+1,sizeof(char));
    targetlength = maxtargetlength;
  }


  if (!tflag && argv[argc-1-(!qflag)][0] != '-') {
    /* target sequence given on command line */
    target_sq = (char *) calloc(strlen(argv[argc-1-(!qflag)])+1,sizeof(char));
    strcpy(target_ac,"command_line");
    strcpy(target_sq,argv[argc-1-(!qflag)]);
    remove_whitespace(target_sq);
    targetlength = strlen(target_sq);
  }


  if (sflag) {
    /* initialise dinucleotide frequencies: */
    for (i=0; i<4; i++)
      for (j=0; j<4; j++)
	freq_di[i][j] = 0.0;

    if (qflag) {
      /* open query file: */
      FILE *query = fopen(query_fn,"r");
      if (query==NULL) {
	printf("Error: Could not open query file. Aborting.\n");
	exit(2);
      }

      /* iterate over query sequences: */
      counter = 0;
      while (!end(query)) {

	nextAC(query,query_ac);
	nextSQ(query,query_sq,maxquerylength);

	remove_whitespace(query_sq);

	seqlen = strlen(query_sq);
	/* count dinucleotides: */
	for (k=0; k<seqlen-1; k++) {
	  letter_one = toupper(query_sq[k]);
	  letter_two = toupper(query_sq[k+1]);

	  if (letter_one == 'U')
	    letter_one = 'T';
	  if (letter_two == 'U')
	    letter_two = 'T';

	  index_one = (int) strchr(alphabet,letter_one);
	  index_two = (int) strchr(alphabet,letter_two);

	  if (index_one != 0 && index_two != 0) {
	    freq_di[index_one-(int) alphabet][index_two-(int) alphabet]++;
	    counter++;
	  }
	}
      }
      /* close query file: */
      fclose(query);
    } /* if qflag */
    else { /* !qflag */

      counter = 0;
      seqlen = strlen(query_sq);
      /* count dinucleotides: */
      for (k=0; k<seqlen-1; k++) {
	letter_one = toupper(query_sq[k]);
	letter_two = toupper(query_sq[k+1]);

	if (letter_one == 'U')
	  letter_one = 'T';
	if (letter_two == 'U')
	  letter_two = 'T';

	index_one = (int) strchr(alphabet,letter_one);
	index_two = (int) strchr(alphabet,letter_two);

	if (index_one != 0 && index_two != 0) {
	  freq_di[index_one-(int) alphabet][index_two-(int) alphabet]++;
	  counter++;
	}
      }

    }

    /* turn counts into relative frequencies: */
    for (i=0; i<4; i++)
      for (j=0; j<4; j++)
	freq_di[i][j] /= (float) counter;
    
  } /* if sflag */


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
 

  if (qflag && sflag || !qflag && argv[argc-1][0] != '-' && sflag)
    free(query_sq);


  if (dflag) {
    /*   read dinucleotide frequencies: */
    f = fopen(freq_filename,"r");
    read_dinucleotide_frequencies(freq_di,f);
    fclose(f);
  }


  /*   calculate single nucleotide frequencies: */
  for (i=0; i<strlen(alphabet); i++) {
    freq[i]=2*freq_di[i][i];
    for (j=0; j<strlen(alphabet); j++)
      if (i!=j)
	freq[i] += freq_di[i][j];
    for (j=0; j<strlen(alphabet); j++)
      if (i!=j)
	freq[i] += freq_di[j][i];
    freq[i] /= 2.0;
  }

  /*   make distribution functions: */
  fdf[0] = freq[0];
  for (i=1; i<strlen(alphabet); i++)
    fdf[i]=fdf[i-1]+freq[i];
  fdf[strlen(alphabet)-1]=1.0;
    
  for (i=0; i<strlen(alphabet); i++) {
    float row_sum = 0;
    for (j=0; j<strlen(alphabet); j++)
      row_sum+=freq_di[i][j];
    fdf_di[i][0]=freq_di[i][0]/row_sum;
    for (j=1; j<strlen(alphabet); j++)
      fdf_di[i][j]=fdf_di[i][j-1]+freq_di[i][j]/row_sum;
    fdf_di[i][strlen(alphabet)-1]=1.0;
  }


  sgenrand(time(NULL));

  /*  Search: */

  if (qflag && !sflag) {

    /* open query file: */
    FILE *query = fopen(query_fn,"r");
    if (query==NULL) {
      printf("Error: Could not open query file. Aborting.\n");
      exit(2);
    }

    /* iterate over query sequences: */
    counter = 0;
    while (!end(query)) {

      nextAC(query,query_ac);
      nextSQ(query,query_sq,maxquerylength);

      remove_whitespace(query_sq);

      seqlen = strlen(query_sq);

      n = strlen(query_sq);
    
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

      if (tflag) {

	FILE *target = fopen(target_fn,"r");
	if (target==NULL) {
	  printf("Error: Could not open target file. Aborting.\n");
	  exit(2);
	}

	k = 0;
	while (!end(target) && k<MAXTARGETNUMBER) {

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

	    for (i1=m; i1>=0; i1--) {
	      for (i2=n; i2>=0; i2--) {
		calc_unpaired_left_bot(i1, m, i2, n);
		calc_closed           (i1, m, i2, n);
		calc_unpaired_left_top(i1, m, i2, n);
	      }
	    }
	    v2 = calc_hybrid(0, m, 0, n);
	      
 	    normalised_energy = -v2/log(m*n);
 	    normalised_energies[k][l] = normalised_energy;

	    /* initialize strings */
	    r1 = t1;
	    r2 = t2;
	    r3 = t3;
	    r4 = t4;

	    k++;
	  }
	}
	fclose(target);
      }
      else { /* !tflag */
	printf("Shouldn't happen.\n");
	exit(1);
      }

      free(query_sq);

      l++;
    }
    for (i=0; i<k; i++) {
      for (j=0; j<l; j++)
	temp_energies[j] = normalised_energies[i][j];
      estimate_evd_parameters(&used_sample_size,&xi[i],&theta[i],temp_energies,l);
/*       printf("%f %f\n",xi[i],theta[i]); */
    }
    
    /* calculate p-values: */
    for (j=0; j<l; j++) {
      for (i=0; i<k; i++) {
	pvalue = 1-exp(-exp(-(normalised_energies[i][j]-xi[i])/theta[i]));
	if (i==0 || pvalue > maximal_pvalues[j])
	  maximal_pvalues[j] = pvalue;
      }
    }

/*     for (i=0; i<sample_size; i++) */
/*       printf("%f\n",maximal_pvalues[i]); */


    /* find effective number of orthologs: */
    effective_number = -1.0;
    for (exponent=1.0; exponent<(float) k+0.01; exponent=exponent+0.1) {

      /* joint pvalues: */
      for (j=0; j<l; j++)
	joint_pvalues[j] = pow(maximal_pvalues[j],exponent);

      /* empirical cdf: */
      empirical_cumulative_distribution_function(&xcdf,&ycdf,&bin_number,joint_pvalues,l);

      /* weighted squared error: */
      error_sum = 0.0;
      for (j=0; j<bin_number; j++)
	error_sum += 1.0/xcdf[j]*pow(ycdf[j+1]-xcdf[j],2.0);

      printf("%f %f\n",exponent,error_sum);

      /* check for new minimum: */
      if (effective_number == -1.0 || error_sum < smallest_error_sum) {
	smallest_error_sum = error_sum;
	effective_number = exponent;
      }
    }

    printf("Effective number: %f\n",effective_number);

  }
  else if (qflag && sflag || !qflag && sflag) {

    l = 0;
    while (l<sample_size) {
      do seqlen = (int) (normal_random_number()*stddev+mean); while (seqlen < 1 || seqlen > maxquerylength);
      query_sq = random_sequence(seqlen,fdf,fdf_di);

      n = strlen(query_sq);
    
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

      if (tflag) {

	FILE *target = fopen(target_fn,"r");
	if (target==NULL) {
	  printf("Error: Could not open target file. Aborting.\n");
	  exit(2);
	}

	k = 0;
	while (!end(target) && k<MAXTARGETNUMBER) {

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

	    for (i1=m; i1>=0; i1--) {
	      for (i2=n; i2>=0; i2--) {
		calc_unpaired_left_bot(i1, m, i2, n);
		calc_closed           (i1, m, i2, n);
		calc_unpaired_left_top(i1, m, i2, n);
	      }
	    }
	    v2 = calc_hybrid(0, m, 0, n);
	      
 	    normalised_energy = -v2/log(m*n);
 	    normalised_energies[k][l] = normalised_energy;

	    /* initialize strings */
	    r1 = t1;
	    r2 = t2;
	    r3 = t3;
	    r4 = t4;

	    k++;
	  }
	}
	fclose(target);
      }
      else { /* !tflag */
	printf("Shouldn't happen.\n");
	exit(1);
      }

      free(query_sq);

      l++;
    }
    for (i=0; i<k; i++) {
      for (j=0; j<l; j++)
	temp_energies[j] = normalised_energies[i][j];
      estimate_evd_parameters(&used_sample_size,&xi[i],&theta[i],temp_energies,l);
/*       printf("%f %f\n",xi[i],theta[i]); */
    }
    
    /* calculate p-values: */
    for (j=0; j<l; j++) {
      for (i=0; i<k; i++) {
	pvalue = 1-exp(-exp(-(normalised_energies[i][j]-xi[i])/theta[i]));
	if (i==0 || pvalue > maximal_pvalues[j])
	  maximal_pvalues[j] = pvalue;
      }
    }

/*     for (i=0; i<sample_size; i++) */
/*       printf("%f\n",maximal_pvalues[i]); */


    /* find effective number of orthologs: */
    effective_number = -1.0;
    for (exponent=1.0; exponent<(float) k+0.01; exponent=exponent+0.1) {

      /* joint pvalues: */
      for (j=0; j<l; j++)
	joint_pvalues[j] = pow(maximal_pvalues[j],exponent);

      /* empirical cdf: */
      empirical_cumulative_distribution_function(&xcdf,&ycdf,&bin_number,joint_pvalues,l);

      /* weighted squared error: */
      error_sum = 0.0;
      for (j=0; j<bin_number; j++)
	error_sum += 1.0/xcdf[j]*pow(ycdf[j+1]-xcdf[j],2.0);

      printf("%f %f\n",exponent,error_sum);

      /* check for new minimum: */
      if (effective_number == -1.0 || error_sum < smallest_error_sum) {
	smallest_error_sum = error_sum;
	effective_number = exponent;
      }
    }

    printf("Effective number: %f\n",effective_number);


  }
  else { /* !qflag && !sflag */
    printf("Refusing to base analysis on one query sequence. Aborting.\n");
    exit(2);
  }





}
