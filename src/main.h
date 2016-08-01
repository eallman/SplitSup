#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <memory.h>
#include <assert.h>
#include <sys/resource.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

typedef char naym[11];

// for choosing SW analysis
extern int SW;

// for choosing rank condition
extern int rank;

// for DNA sequence data
extern int nsite, nseq, num_unique, include_gaps, num_no_gaps;
extern int blocksize;
extern int **ppBase, **ppBase_unique, *site_counter;
extern int **ppBase_full;

// for scores computataions
extern double fnorm2;
extern FILE *scoresfile;

// for binary encoding
extern int nints, nfit, num_unique_bin;
extern unsigned int *split_taxa, *split_taxaR, *split_taxaC;
extern unsigned int *pArr_u_bin, **ppBase_u_bin;

extern void free_Memory_bin();
extern int Data_Bin();
extern void getSVDscore(unsigned int *split_taxa, int sze);
