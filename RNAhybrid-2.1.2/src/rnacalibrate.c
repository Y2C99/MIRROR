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


static char const rcsid[] = "$Id: rnacalibrate.c,v 1.11 2004/11/16 15:42:33 marc Exp $";


struct output {
  char query_ac[MAXLINE];
  int used_sample_size;
  float location, scale;
};


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

  int mean, stddev, sample_size;

  struct output output;

  FILE *f;

  int i,j,k,seqlen;

  char letter_one, letter_two;

  int index_one, index_two;

  int counter;

  

  float *normalised_energies;

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
    printf("\nUsage: %s [options] [target sequence] [query sequence].\n\noptions:\n\n  -d <dinucleotide frequencies file>\n  -f helix constraint\n  -h help\n  -k <sample size>\n  -l <mean sequence length>,<std sequence length>\n  -m <max targetlength>\n  -n <max query length>\n  -u <max internal loop size (per side)>\n  -v <max bulge loop size>\n  -s randomise targets (only with -t, without -d)\n  -t <target file>\n  -q <query file>\n\nIf no target file is given, random sequences are generated according\nto the given dinucleotide distribution. Then <sample size> sequences\nare generated whose lengths are normally distributed with given mean\nand standard deviation. Default sample size is 5000, default mean\nand std are 500 and 300, respectively.\n\nIf a target file is given, and additionally the -s option, random sequences are generated according to the dinucleotide distribution of the target file.\n\nIf only a target file is given, it is used directly as a random database.\n\nThe target can also be given directly (makes only sense with the -s option).\n\nEither a query file has to be given (FASTA format)\nor one query sequence directly.\n\nThe helix constraint format is \"from,to\", eg. -f 2,7 forces\nstructures to have a helix from position 2 to 7 with respect to the query.\n\n", argv[0]);
    exit(0);
  }

  if (!uflag)
    iloop_upper_limit = ILOOPUPPERLIMITDEFAULT;
  if (!vflag)
    bloop_upper_limit = BLOOPUPPERLIMITDEFAULT;

  if (!lflag) {
    mean = MEAN;
    stddev = STDDEV;
  }

  if (!kflag)
    sample_size = SAMPLESIZE;

  if (tflag && !sflag)
    normalised_energies = (float *) calloc(MAXTARGETNUMBER,sizeof(float));
  else
    normalised_energies = (float *) calloc(sample_size,sizeof(float));

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

  if (tflag)
    target_sq = (char *) calloc(maxtargetlength+MAXLINE+1,sizeof(char));


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

    if (tflag) {
      /* open target file: */
      FILE *target = fopen(target_fn,"r");
      if (target==NULL) {
	printf("Error: Could not open target file. Aborting.\n");
	exit(2);
      }

      /* iterate over target sequences: */
      counter = 0;
      while (!end(target)) {

	nextAC(target,target_ac);
	nextSQ(target,target_sq,maxtargetlength);

	remove_whitespace(target_sq);

	seqlen = strlen(target_sq);
	/* count dinucleotides: */
	for (k=0; k<seqlen-1; k++) {
	  letter_one = toupper(target_sq[k]);
	  letter_two = toupper(target_sq[k+1]);

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
      /* close target file: */
      fclose(target);
    } /* if tflag */
    else { /* !tflag */

      counter = 0;
      seqlen = strlen(target_sq);
      /* count dinucleotides: */
      for (k=0; k<seqlen-1; k++) {
	letter_one = toupper(target_sq[k]);
	letter_two = toupper(target_sq[k+1]);

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


  targetlength = maxtargetlength;

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
 

  if (tflag && sflag || !tflag && argv[argc-1-(!qflag)][0] != '-' && sflag)
    free(target_sq);



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

	if (tflag && !sflag) {
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
	      normalised_energies[k] = normalised_energy;

/* 	      free(target_sq); */

	      /* initialize strings */
	      r1 = t1;
	      r2 = t2;
	      r3 = t3;
	      r4 = t4;

	      k++;
	    }
	  }
	  strcpy(output.query_ac,query_ac);
	  estimate_evd_parameters(&output.used_sample_size,&output.location,&output.scale,normalised_energies,k);
	  printf("%s %d %f %f\n",output.query_ac,output.used_sample_size,output.location,output.scale);
/* 	  printf("%s %f %f\n",query_ac,xi,theta); */

	  fclose(target);
	} /* if tflag && !sflag */
	else if (!tflag && !sflag && !dflag) {
	  m = strlen(target_sq);

	  k = 0;
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
	    normalised_energies[k] = normalised_energy;

	    free(target_sq);

	    /* initialize strings */
	    r1 = t1;
	    r2 = t2;
	    r3 = t3;
	    r4 = t4;
	
	    k++;
	  }

	  strcpy(output.query_ac,query_ac);
	  estimate_evd_parameters(&output.used_sample_size,&output.location,&output.scale,normalised_energies,k);
	  printf("%s %d %f %f\n",output.query_ac,output.used_sample_size,output.location,output.scale);

/* 	  printf("%s %f %f\n",query_ac,xi,theta); */

	} /* if !tflag && !sflag %% !dflag*/
	else { /* !tflag or sflag or dflag*/
	  k = 0;
	  while (k<sample_size) {
	    do seqlen = (int) (normal_random_number()*stddev+mean); while (seqlen < 1 || seqlen > maxtargetlength);
	    target_sq = random_sequence(seqlen,fdf,fdf_di);
	    
	    m = strlen(target_sq);

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
	    normalised_energies[k] = normalised_energy;

/* 	    printf("%f\n",normalised_energy); */

	    free(target_sq);

	    /* initialize strings */
	    r1 = t1;
	    r2 = t2;
	    r3 = t3;
	    r4 = t4;

	    k++;
	  }
	  strcpy(output.query_ac,query_ac);
	  estimate_evd_parameters(&output.used_sample_size,&output.location,&output.scale,normalised_energies,k);
	  printf("%s %d %f %f\n",output.query_ac,output.used_sample_size,output.location,output.scale);
	}
      }

    }
    fclose(query);
  } /* if qflag */
  else { /* !qflag */
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

    if (tflag && !sflag) {
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
	  normalised_energies[k] = normalised_energy;

	  /* initialize strings */
	  r1 = t1;
	  r2 = t2;
	  r3 = t3;
	  r4 = t4;

	  k++;
	}
      }
      strcpy(output.query_ac,query_ac);
      estimate_evd_parameters(&output.used_sample_size,&output.location,&output.scale,normalised_energies,k);
      printf("%s %d %f %f\n",output.query_ac,output.used_sample_size,output.location,output.scale);

      fclose(target);
    } /* if tflag && !sflag */
    else if (!tflag && !sflag && !dflag) {
      m = strlen(target_sq);

      k = 0;
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
	normalised_energies[k] = normalised_energy;

	free(target_sq);

	/* initialize strings */
	r1 = t1;
	r2 = t2;
	r3 = t3;
	r4 = t4;

	k++;
      }
      strcpy(output.query_ac,query_ac);
      estimate_evd_parameters(&output.used_sample_size,&output.location,&output.scale,normalised_energies,k);
      printf("%s %d %f %f\n",output.query_ac,output.used_sample_size,output.location,output.scale);
    } /* if !tflag && !sflag && !dflag */
    else { /* !tflag or sflag or dflag */
      k = 0;
      while (k<sample_size) {
	do seqlen = (int) (normal_random_number()*stddev+mean); while (seqlen < 1 || seqlen > maxtargetlength);
	target_sq = random_sequence(seqlen,fdf,fdf_di);
	
	m = strlen(target_sq);

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
	normalised_energies[k] = normalised_energy;

	free(target_sq);

	/* initialize strings */
	r1 = t1;
	r2 = t2;
	r3 = t3;
	r4 = t4;

	k++;
      }
      strcpy(output.query_ac,query_ac);
      estimate_evd_parameters(&output.used_sample_size,&output.location,&output.scale,normalised_energies,k);
      printf("%s %d %f %f\n",output.query_ac,output.used_sample_size,output.location,output.scale);
    } /* !tflag or sflag */
  }

  if (tflag && !sflag)
    free(target_sq);

  for (i=0; i<4; i++) {
    free(freq_di[i]);
    free(fdf_di[i]);
  }
  free(freq_di);
  free(fdf_di);
  free(freq);
  free(fdf);

  free(normalised_energies);

  free(query_sq);

  return 0;
}

 	

