#include "las2.c"
#include "svdutil.c"
#include "svdlib.c"
#include "main.h"

/*********************************************************/
/****** Functions to compute the flattened matrix ********/
/*********************************************************/


/** Memory allocation for all global variable which are pointers
    that are used in binary encoding                             **/

int MemAlloc_bin() {
  int i;

  // split_taxa memory allocation in main() since global

  /* allocate space for binary encoding of splits */  
  split_taxaR = (unsigned int*)malloc((nints)*sizeof(unsigned int));
  if (split_taxaR==NULL)
    {
      printf("Can't memalloc split_taxaR.\n");
      return 0;
    }

  split_taxaC = (unsigned int*)malloc((nints)*sizeof(unsigned int));
  if (split_taxaC==NULL)
    {
      printf("Can't memalloc split_taxaC.\n");
      return 0;
    }

  /** Allocate array for unique site patterns                 **/
  /** Note: we allocate as contiguous block for use in qsort  **/
  /**       Sorting will also use 2 extra ints in each column **/

  pArr_u_bin = (unsigned int*)malloc(num_unique*(nints+2)*sizeof(unsigned int));
  if (pArr_u_bin==NULL)
    {
      printf("Can't memalloc pArr_u_bin.\n");
      return 0;
    }
  
  /** now allocate for pointers to patterns(rows of array) **/
  ppBase_u_bin=(unsigned int**)malloc(num_unique*sizeof(unsigned int*));
  if (ppBase_u_bin==NULL)
    {
      printf("Can't memalloc ppBase_u_bin.\n");
      return 0;
    }
  /** complete allocation by pointing ppBase_u_bin to pArr_u_bin **/
  for (i=0; i<num_unique; i++)
    {
      ppBase_u_bin[i]=pArr_u_bin+i*(nints+2);
    }

  return 1;

}

/** Free memory for all global variable which are pointers
    that are used in binary encoding                             **/

void free_Memory_bin() {
  int i;

  if (num_unique!=0){
    for (i=num_unique-1; i<0; i--){
      free(ppBase_u_bin[i]);
    }
  }
  free(ppBase_u_bin);

  free(pArr_u_bin);
  
  free(split_taxaR);
  free(split_taxaC);

}

/*************************************************************************/
/**  Compare binary vectors encoding patterns, using a mask for split   **/
/**  This has both row and column versions                              **/
/*************************************************************************/

int bin_compC(const void *n1, const void *n2)
{  
  int i;
  unsigned int a,b;

  const unsigned int *num1 =(const unsigned int *)n1;
  const unsigned int *num2 =(const unsigned int *)n2;

  for (i=0; i<nints; i++)
    {
      a=(num1[i] & split_taxaC[i]);
      b=(num2[i] & split_taxaC[i]);
      if ( a<b ) return -1;
      else if ( a>b ) return  1;
    }
  return 0;
}

int bin_compR(const void *n1, const void *n2)
{
  int i;
  unsigned int a,b;

  const unsigned int *num1 =(const unsigned int *)n1;
  const unsigned int *num2 =(const unsigned int *)n2;

  for (i=0; i<nints; i++)
    {
      a=(num1[i] & split_taxaR[i]);
      b=(num2[i] & split_taxaR[i]);
      if ( a<b ) return -1;
      else if ( a>b ) return  1;
    }
  return 0;
}


/** Create binary encoding of patterns for use in flattening **/
int Data_Bin()
{
  /** Note: All sites with states other than ATCG are removed **/

  int  badBase, i, j, jj, k, m;

  // memory allocation for binary encoding
  i = MemAlloc_bin();
  if (i==0) {
    printf("Failure to allocate memory.  Exiting ....\n");
    exit(1);
  }


  /** Encode unique patterns with ONLY ATCG in binary format **/
  jj=0;     // initialize pattern number in ppBase_u_bin

  for (j=0; j<num_unique; j++)  // loop on unique sites                                                       
    {
      badBase=0;  // initialize flag,  1 will indicate non-ATCG
      for (i=0; i<nseq; i++)  // loop on sequences 
        {
          if (ppBase_unique[i][j]<4) // make sure we have a ATCG 
            {
              k=i/nfit;      // compute word to code this taxon in
              m=i-k*nfit;    // and location in word, sort of

              if ( m == 0 )  /** if beginning new word **/
                {
		  /*... just store */
                  ppBase_u_bin[jj][k]=ppBase_unique[i][j]; 
                }
              else
                {
		  /*... otherwise, shift base and combine */
		  ppBase_u_bin[jj][k] ^= ( ppBase_unique[i][j] << (2*m) ); 
                }
            }
          else badBase=1; // a non-ATCG popped up --- 
	                  // this site will be ignored eventually                  
        }

      if ( badBase==0 ) // if this site had only ATCG                                                       
        {
	  // set high bit so all patterns gives non-zero encoding    
          ppBase_u_bin[jj][nints-1] ^= (1 << 30); 
          ppBase_u_bin[jj][nints]=site_counter[j]; // copy counts                                            
          jj++; // don't overwrite binary encoded pattern                                                    
        }
    }
  num_unique_bin=jj; // number of unique patterns in binary encoding                                         

  return 0;
}


/** Read splits of interest from file, compute svd and scores **/

void getSVDscore(unsigned int *split_taxa, int sze)
{  

  int j, jj;
  int num_unique_rows, num_unique_cols;
  long int* firstnz; 

  // initialize score
  double score=0;

  SMat S = NULL;
  SVDRec R = NULL;

  /* initialize split_taxa vector */
  for (j=0; j<nints; j++) {
    split_taxaR[j] = 0;
    split_taxaC[j] = 0;
  }

  // transfer split 
  for (j=0; j<nints; j++) split_taxaR[j] = split_taxa[j];   
  // complementary mask 
  for (j=0; j<nints; j++) split_taxaC[j] =~ split_taxaR[j]; 

  /* sort on rows, to get row indices, 
     then on columns, to get sparse encoding of flattening */

  /* sort for rows */
  qsort(*ppBase_u_bin, num_unique_bin, (nints+2)*sizeof(unsigned int), bin_compR);
  
  /* Enter row numbers into array */
  jj=0;  // initialize row number
  ppBase_u_bin[0][nints+1]=0;  // this is always in first row
  for (j=1; j<num_unique_bin; j++)
    {	   	  
      if (bin_compR(ppBase_u_bin[j],ppBase_u_bin[j-1])==1) jj++;
      ppBase_u_bin[j][nints+1]=jj; // record row number
    }
  num_unique_rows=jj+1;

  /* sort for columns */
  qsort(*ppBase_u_bin, num_unique_bin, (nints+2)*sizeof(unsigned int), bin_compC);
  
  /* allocate space for sparse matrix encoding */
  long int* rowindex = (long int*) calloc(num_unique_bin,sizeof(long int));
  double* values  = (double*) calloc(num_unique_bin,sizeof(double));

  /*  create vector encoding column numbers */
  jj=0; // first column number
  // copy value of first non-zero for SVD
  values[0] = (double)ppBase_u_bin[0][nints]/(double)num_no_gaps;

  // copy row index of first entry for SVD
  rowindex[0]=ppBase_u_bin[0][nints+1]; 
  ppBase_u_bin[0][nints+1]=0;  // first entry is always in first column

  /* Note: we're reusing this space to store column info as we remove row info */
  for (j=1; j<num_unique_bin; j++)
    { 
      // copy non-zeros for SVD
      values[j] = (double)(ppBase_u_bin[j][nints])/(double)num_no_gaps;
      rowindex[j]=ppBase_u_bin[j][nints+1]; // copy row indices for SVD
      if (bin_compC(ppBase_u_bin[j],ppBase_u_bin[j-1])==1)  // if new column...
	{ 
          jj++;      
	  ppBase_u_bin[jj][nints+1]=j; // save column coding here temporarily 
	}
    }
  num_unique_cols=jj+1;

  /* allocate final space for sparse matrix encoding */
  firstnz = (long int*) calloc( num_unique_cols+1 , sizeof(long int));

  for (j=0; j<num_unique_cols; j++)
    {
      firstnz[j]=ppBase_u_bin[j][nints+1];
    }

  /* copy over column coding for SVD */
  firstnz[j]=num_unique_bin;/* sparse format requires this too */
  
  /* Create Smat object to use in SVD routines */
  S = svdNewSMat(num_unique_rows,num_unique_cols,num_unique_bin);
  R = svdNewSVDRec();
  
  S->pointr = firstnz;
  S->rowind = rowindex;
  S->value = values;
  
  /* Compute SVD */
  extern long SVDVerbosity;
  SVDVerbosity=0; // print nothing from SVD routine

  R = svdLAS2A(S,rank);

  /** Scores computations **/
  for (j=0; j<rank; j++)   score = score+pow(R->S[j],2);
  score = sqrt(1-score/fnorm2);
  
  if (SW==0)
    {
      // no SW
      fprintf(scoresfile, "%1.15f\n", score);
    }
  else 
    {
      // SW analysis
      fprintf(scoresfile, "%1.15f\n",score);
    }
  
  fflush(0);
  
  free(firstnz);    
  free(rowindex);   
  free(values);
  
  S->pointr = NULL;
  S->rowind = NULL;
  S->value = NULL;
  
  svdFreeSMat(S);
  svdFreeSVDRec(R);
  
  num_unique_rows=0;
  num_unique_cols=0;

}
