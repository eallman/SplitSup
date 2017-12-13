/***************************************************************/
/********** SplitSup version 1.02                     **********/
/********** Written by L. S. Kubatko, E. S. Allman    **********/
/********** and J. Rhodes, July 2016                  **********/
/***************************************************************/

#include "main.h"

// set value to 1 if want timing information printed out
int showTimeElapsed = 1;

/************************/
/** global variables   **/
/************************/

/** general **/
naym *psname;
int nsite, nseq, num_unique, include_gaps, num_no_gaps, num_const;
int SW, rank;  

// parameters needed for SW analysis
int  blocksize;                        /* Window size for sliding window. */ 
int  slidesize;                        /* Offset for each iteration of sliding window.  */
int  min_sites;                        /* Minimum number of sites need for computing score */

int **ppBase_full, **ppBase, **ppBase_unique, *site_counter;
double fnorm2;

/** file variables:
    infile:      input DNA sequence file in phylip format
    splitsfile:  input file with splits
    scoresfile:  output file for scores                   **/

FILE *infile;
FILE *splitsfile;
FILE *scoresfile;

/** variables for binary encoding **/
int nints, nfit, num_unique_bin;
unsigned int *split_taxa;
unsigned int **ppBase_u_bin, *counts, *split_taxaR, *split_taxaC;
unsigned int *pArr_u_bin;

void parse_cmdline(int argc, char *argv[])
{
  int m, n,                              /* Loop counters. */
    l,                                   /* String length. */
    x,                                   /* Exit code. */
    ch;                                  /* Character buffer. */

  char s[256];                           /* String buffer. */
 
  switch (argc)
    {
    case 1:
      // Not performing sliding window analysis
      SW=0;
      printf("\n SplitSup version 1.02\n\n\n"); 
      break;

    case 2:
      // Not performing sliding window analysis
      // but setting user-defined rank condition
      SW=0;
      printf("\n SplitSup version 1.02\n\n\n"); 

      n = 1;

      switch( (int)argv[n][0] )            /* Check for option character. */
	{
	case '-': x = 0;                   /* Bail out if 1. */
	  l = strlen( argv[n] );

	  // exit if only a dash and nothing else
	  if (l==1) {
	    printf( " Illegal option \n\n" );
	    exit(1);
	  }

	  for( m = 1; m < l; ++m ) /* Scan through options. */
	    {
	      ch = (int)argv[n][m];
	      switch( ch )
		{
		  // number of singular values to compute
		case 'r':
		case 'R':
		  if( m + 1 >= l )
		    {
		      puts( " You did not indicate a user-specified rank!\n" );
		      exit( 1 );
		    }
		  else
		    {
		      strcpy( s, &argv[n][m+1] );
		      rank = atoi(s);
		      printf(" Testing for rank %d matrix flattenings.\n",rank);
		    }
		  x = 1;
		  break;

		default:  printf( " Illegal option code = %c\n\n", ch );
		  x = 1;      /* Not legal option. */
		  exit( 1 );
		  break;
		}
	      if( x == 1 )
		{
		  break;
		}
	    }
	  break;

	default:
	  /* Not option -- text. */
	  printf( " %s is not an option\n\n", argv[n] );
	  exit( 1 );
	  break;
	}

      break;

    case 5:
      // Sliding Window analysis

      for( n = 1; n < argc; n++ )              /* Scan through args. */
  	{
  	  switch( (int)argv[n][0] )            /* Check for option character. */
  	    {
  	    case '-': x = 0;                   /* Bail out if 1. */
  	      l = strlen( argv[n] );
  	      for( m = 1; m < l; ++m )         /* Scan through options. */
  		{
  		  ch = (int)argv[n][m];
  		  switch( ch )
  		    {
  		      // Sliding window size
  		    case 'w':
  		    case 'W':
  		      if( m + 1 >= l )
  			{
  			  puts( " You did not indicate a window size!\n" );
  			  exit( 1 );
  			}
  		      else
  			{
  			  strcpy( s, &argv[n][m+1] );
  			  blocksize = atoi(s);
  			}
  		      x = 1;
  		      break;

  		      // Offset size
  		    case 'o':
  		    case 'O':
  		      if( m + 1 >= l )
  			{
  			  puts( " You did not indicate a slide offset size!\n" );
  			  exit( 1 );
  			}
  		      else
  			{
  			  strcpy( s, &argv[n][m+1] );
  			  slidesize = atoi(s);
  			}
  		      x = 1;
  		      break;

  		      // Mininum of sites needed in window to compute score
  		    case 'm':
  		    case 'M':
  		      if( m + 1 >= l )
  			{
  			  puts( " You did not indicate the mininum number of gapless sites needed!\n" );
  			  exit( 1 );
  			}
  		      else
  			{
  			  strcpy( s, &argv[n][m+1] );
  			  min_sites = atoi(s);
  			}
  		      x = 1;
  		      break;

  		      // Choose Sliding Window analysis
  		    case 's':
  		    case 'S':
  		      x = 1;
  		      break;

  		    default:  printf( "  Illegal option code = %c\n\n", ch );
  		      x = 1;      /* Not legal option. */
  		      exit( 1 );
  		      break;
  		    }
  		  if( x == 1 )
  		    {
  		      break;
  		    }
  		}
  	      break;

  	    default:
  	      /* Not option -- text. */
  	      printf( "%s is not an option\n", argv[n] );
  	      exit( 1 );
  	      break;
  	    }

  	}

      //  Bail out if user-specified window size is smaller than min_sites
      if (min_sites > blocksize) {
	printf( "\n The minimum number of sites (%d) can not be larger than the window size (%d).\n",min_sites,blocksize);
	printf( " Exiting.\n\n");
	exit( 1 );
	break;
      }

      // Correctly gave parameters for SW analysis.  Set toggle.
      SW=1;
      printf("\n SplitSup version 1.02\n\n\n"); 
      printf(" Sliding window analysis with parameters:\n");
      printf("     window size is %d\n",blocksize);
      printf("     slide is %d\n",slidesize);
      printf("     %d gapless sites will be required for the window to be used\n\n\n",min_sites);
      break;

    case 6:
      // Sliding Window analysis and user-defined rank condition

      for( n = 1; n < argc; n++ )              /* Scan through args. */
  	{
  	  switch( (int)argv[n][0] )            /* Check for option character. */
  	    {
  	    case '-': x = 0;                   /* Bail out if 1. */
  	      l = strlen( argv[n] );
  	      for( m = 1; m < l; ++m )         /* Scan through options. */
  		{
  		  ch = (int)argv[n][m];
  		  switch( ch )
  		    {
  		      // Sliding window size
  		    case 'w':
  		    case 'W':
  		      if( m + 1 >= l )
  			{
  			  puts( " You did not indicate a window size!\n" );
  			  exit( 1 );
  			}
  		      else
  			{
  			  strcpy( s, &argv[n][m+1] );
  			  blocksize = atoi(s);
  			}
  		      x = 1;
  		      break;

  		      // Offset size
  		    case 'o':
  		    case 'O':
  		      if( m + 1 >= l )
  			{
  			  puts( " You did not indicate a slide offset size!\n" );
  			  exit( 1 );
  			}
  		      else
  			{
  			  strcpy( s, &argv[n][m+1] );
  			  slidesize = atoi(s);
  			}
  		      x = 1;
  		      break;

  		      // Mininum of sites needed in window to compute score
  		    case 'm':
  		    case 'M':
  		      if( m + 1 >= l )
  			{
  			  puts( " You did not indicate the mininum number of gapless sites needed!\n" );
  			  exit( 1 );
  			}
  		      else
  			{
  			  strcpy( s, &argv[n][m+1] );
  			  min_sites = atoi(s);
  			}
  		      x = 1;
  		      break;

  		      // Choose Sliding Window analysis
  		    case 's':
  		    case 'S':
  		      x = 1;
  		      break;

		    case 'r':
		    case 'R':
		      if( m + 1 >= l )
			{
			  puts( " You did not indicate a user-specified rank!\n" );
			  exit( 1 );
			}
		      else
			{
			  strcpy( s, &argv[n][m+1] );
			  rank = atoi(s);

			}
		      x = 1;
		      break;

  		    default:  printf( "  Illegal option code = %c\n\n", ch );
  		      x = 1;      /* Not legal option. */
  		      exit( 1 );
  		      break;
  		    }
  		  if( x == 1 )
  		    {
  		      break;
  		    }
  		}
  	      break;

  	    default:
  	      /* Not option -- text. */
  	      printf( "%s is not an option\n", argv[n] );
  	      exit( 1 );
  	      break;
  	    }

  	}

      //  Bail out if user-specified window size is smaller than min_sites
      if (min_sites > blocksize) {
	printf( "\n The minimum number of sites (%d) can not be larger than the window size (%d).\n",min_sites,blocksize);
	printf( " Exiting.\n\n");
	exit( 1 );
	break;
      }

      // Correctly gave parameters for SW analysis.  Set toggle.
      SW=1;
      printf("\n SplitSup version 1.02\n\n\n"); 
      printf(" Sliding window analysis with parameters:\n");
      printf("     window size is %d\n",blocksize);
      printf("     slide is %d\n",slidesize);
      printf("     %d gapless sites will be required for the window to be used\n\n\n",min_sites);
      if (rank!=4) printf(" Testing for rank %d matrix flattenings.\n",rank);
      break;

    default:

      printf( "   usage: %s [-rRANK]\n", argv[0]);
      printf("       OR\n");
      printf( "   usage: %s -s -wWINDOW_SIZE -oOFFSET_SIZE -mMIN_NUMBER_SITES [-rRANK]\n", argv[0] );

      exit( 1 );
      break;

    }
}


// Memory allocation for all global variables which are pointers 

int MemAlloc() {

  int i;

  ppBase_full = (int**)malloc(nseq*sizeof(int*));
  if (ppBase_full==NULL)
    {
      printf("     Can't memalloc ppBase_full.\n");
      return 0;
    }

  for (i=0; i<nseq; i++) {
    ppBase_full[i]=(int*)malloc((nsite+1)*sizeof(int));
    if (ppBase_full[i]==NULL) 
      {
	printf("     Can't memalloc ppBase_full[%d].\n",i);
	return 0;
      }
  }

  ppBase = (int**)malloc(nseq*sizeof(int*));
  if (ppBase==NULL)
    {
      printf("     Can't memalloc ppBase.\n");
      return 0;
    }
  
  for (i=0; i<nseq; i++)
    {
      ppBase[i]=(int*)malloc((blocksize+1)*sizeof(int));
      if (ppBase[i]==NULL)
	{
	  printf("     Can't memalloc ppBase[%d].\n",i);
	  return 0;
	}
    }
 
  ppBase_unique = (int**)malloc(nseq*sizeof(int*));
  if (ppBase_unique==NULL)
    {
      printf("     Can't memalloc ppBase_unique.\n");
      return 0;
    }
  
  for (i=0; i<nseq; i++)
    {
      ppBase_unique[i]=(int*)malloc((blocksize+5)*sizeof(int));
      if (ppBase_unique[i]==NULL)
	{
	  printf("     Can't memalloc ppBase_unique[%d].\n",i);
	  return 0;
	}
    }

  site_counter = (int*)malloc((blocksize+1)*sizeof(int));
  if (site_counter==NULL)
    {
      printf("     Can't memalloc site_counter.\n");
      return 0;
    }
    
  return 1;
}


// Transform character sequences (A,G,C,T) into numerical values (0,1,2,3) 

void transfer(char **cbase,int **nbase,int pres_loc, int trans_length) {

  int i, j;

  for (i=0; i<nseq; i++)  {

    for (j=0; j<trans_length; j++) {

      switch (cbase[i][j]) {
	    
      case 'A': case 'a': nbase[i][j+pres_loc] = 0;  break;
      case 'G': case 'g': nbase[i][j+pres_loc] = 1;  break;
      case 'C': case 'c': nbase[i][j+pres_loc] = 2;  break;
      case 'T': case 't': case 'U': case 'u': nbase[i][j+pres_loc] = 3;  break;
      case 'M': case 'm': nbase[i][j+pres_loc] = 5;  break;
      case 'R': case 'r': nbase[i][j+pres_loc] = 6;  break;
      case 'W': case 'w': nbase[i][j+pres_loc] = 7;  break;
      case 'S': case 's': nbase[i][j+pres_loc] = 8;  break;
      case 'Y': case 'y': nbase[i][j+pres_loc] = 9;  break;
      case 'K': case 'k': nbase[i][j+pres_loc] = 10;  break;
      case 'B': case 'b': nbase[i][j+pres_loc] = 11;  break;
      case 'D': case 'd': nbase[i][j+pres_loc] = 12;  break;
      case 'H': case 'h': nbase[i][j+pres_loc] = 13;  break;
      case 'V': case 'v': nbase[i][j+pres_loc] = 14;  break;
      case 'N': case 'n': nbase[i][j+pres_loc] = 15;  break;
      case '.':nbase[i][j+pres_loc] = nbase[0][j+pres_loc]; break;
      case '-':nbase[i][j+pres_loc]=4; break;
      case '?':nbase[i][j+pres_loc]=15; break;
      default: break;

      }
	  
    }
  }
}


// With the number of sequences (nseq), number of sites (nsite),    
// get  DNA sequences of taxa from datafile.                        
// Then transform DNA sequences into numerical sequences (ppBase)   

naym* getSequenceData(void) {
  
  int i, k;
  char **ppDNASeq;
  int charLen;
  
  charLen=(nsite+1)*sizeof(char);
  ppDNASeq=(char**)malloc(nseq*sizeof(char*));
  if(ppDNASeq==NULL){
    return NULL;
  }
  
  for (i=0; i<nseq; i++)
    {
      ppDNASeq[i]=(char*)malloc(charLen);
      if(ppDNASeq[i]==NULL){
	return NULL;
      }
    }
  
  psname = (naym*)malloc(nseq*sizeof(naym));
  
  if(psname==NULL){
    return NULL;
  }
  
  for (i=0; i<nseq; i++)
    {
      fscanf(infile, "%s %s", psname[i],ppDNASeq[i]);

      // No need to pad or load names, but done for using with other code.
      for (k=0; k<11; k++) if (psname[i][k] == 0) psname[i][k] = ' ';
 
    }
  
  transfer(ppDNASeq,ppBase_full,0,nsite);
  
  fclose(infile);
  
  for (i=0; i<nseq; i++) {
    free(ppDNASeq[i]);
  }
  
  free(ppDNASeq);
  
  return psname;
  
}


// Search all unique sites and count them - counts  are  
// stored in the array site_counter.                     

void unique_sites(void){

  int i,k,s,m=0,j=0,jj=0;
  int start;

  start=0;
  num_unique=0;
  num_no_gaps=0;
  num_const=0;
  for (i=0; i<blocksize; i++) site_counter[i]=0;

  // Find the first site with no gaps if include_gaps==0 

  if (include_gaps==0) {

    for (j=0; j<blocksize; j++) {      
      i=0;
      while (i<nseq && ppBase[i][j]!=4 && ppBase[i][j]!=15) { i++; }
      if (i==nseq) {
	start=j;
	break;
      }
      
    }

  }


  if (j==blocksize && start==0) start=blocksize;

  // First remove sites with all gaps if include_gaps==1 

  if (include_gaps==1) {

    for (j=0; j<blocksize; j++) {
      i=0; 
      while (i<nseq && ppBase[i][j]==4) { i++; }
      if (i==nseq) ppBase[0][j] = 99;
    }
  }

  for (j=start; j<blocksize; j++) {
    if ((ppBase[0][j]!=99 && include_gaps==1) || (ppBase[0][j]!=99 && include_gaps==0 && ppBase[0][j]!=4 && ppBase[0][j]!=15)) {
      site_counter[m]=1;
      for (i=0; i<nseq; i++) {
	ppBase_unique[i][m] = ppBase[i][j];
      }
      ppBase[0][j]=99;

      for (k=j+1; k<blocksize; k++) {
	i=0;
	while (i< nseq && ppBase_unique[i][m] == ppBase[i][k]) {
	  i++;
	}

	if (i==nseq) {
	  site_counter[m]++;
	  ppBase[0][k]=99;
	}

	if (i<nseq && (ppBase[i][k]==4 || ppBase[i][k]==15) && include_gaps==0) {
	  ppBase[0][k]=99;
	}

	if (i<nseq && include_gaps==0) {
	  for (s=0; s<nseq; s++) {
	    if (ppBase[s][k]==4 || ppBase[s][k]==15) ppBase[0][k]=99;
	  }
	}
      }
      m++;
    }
  }

  num_unique=m;

  for (i=0; i<num_unique; i++) 
    {
      num_no_gaps += site_counter[i];
      for (jj=0; jj<nseq; jj++) if (ppBase_unique[jj][i] != ppBase_unique[0][i]) break;
      if (jj==nseq) num_const += site_counter[i];
    }

}

//  Free pointers
void free_Memory(){
  int i;

  free(split_taxa);

  if (SW==1) {
    for (i=0; i<nseq; i++) {
      free(ppBase_full[i]);
    }
    free(ppBase_full);
  }

  for (i=0; i<nseq; i++) {
    free(ppBase[i]);
  }
  free(ppBase);

  for (i=0; i<nseq; i++) {
    free(ppBase_unique[i]);
  }
  free(ppBase_unique);

  for (i=nseq-1; i<0; i--)
    {
      free(psname[i]);
    }
  free(psname);

  free(site_counter);

}



/*************   The Main Program  *********************/


int main( int argc, char *argv[])
{

  time_t start,end;
  start=clock();  

  int  num_blocks; 

  naym *sname;
  int i, jj;

  int nsplits, split_size, taxon;
  int split_num;  
  int k, kk;
  
  int nn, mm, pp;

  include_gaps=0;
  fnorm2 = 0.0;


  /* Boolean for choosing sliding window analysis or not.
     SW: 0 for no SW analysis  (default)
     1 for choosing SW analyis                        */
  SW = -10;


  // set default value of rank = 4 for GM model with no mixtures
  rank = 4;

  /** process command line arguments **/
  parse_cmdline(argc, argv );


  // Read data from infile

  if ((infile=fopen("infile","r"))==NULL)
    {
      printf("\tCan't open infile.  Exiting.\n");
      exit(1);;
    }
  fscanf(infile, "%d %d", &nseq, &nsite);


  // set blocksize to be nsite for non-SW analysis; 

  if (SW==0)  
    {  
      blocksize=nsite;
    }


  // Allocate memory 
  i=MemAlloc();
  if (i==0) {
    printf(" Failure to allocate memory.  Exiting ....\n");
    exit(1);
  }

  // Call the function getSequenceData which will read sequences  
  // from infile and convert them to numeric using the function   
  // transform.  Returns NULL if a problem occurs.                
 
  sname = getSequenceData();
  if (sname==NULL) {
    printf(" Can't allocate memory in getSequenceData.  Exiting.\n");
    exit(1);
  }

  // open files for reading splits and for writing scores                                                                           
  if ((splitsfile=fopen("splitsfile","r"))==NULL)
    {
      printf(" Can't open splitsfile.  Exiting.\n");
      exit(1);
    }

  if ((scoresfile=fopen("scoresfile","w"))==NULL)
    {
      printf(" Can't open scoresfile.  Exiting.\n");
      exit(1);
    }


  // probably could be changed to 64-bit integer ?

  // compute quantities needed to store splits in binary format                                                                     
  nfit = 32/2;             // number of bases we can encode in one 32-bit integer                                                   
  nints = (nseq/nfit) + 1; // number of integers needed for encoding patterns                                                       
                           //       (plus at least 2 extra bits)   

  // memory allocation for split_taxa
  split_taxa = (unsigned int*)malloc((nints)*sizeof(unsigned int));
  if (split_taxa==NULL)
    {
      printf("     Can't memalloc split_taxa.\n");
      return 0;
    }

  switch (SW)
    {
    case 0:

      // user did not choose SW analysis

      // assign variables needed for analysis
      slidesize=blocksize;
      min_sites=1;

      ppBase=ppBase_full;

      // Find all unique site patterns and store them in ppBase_unique 
      unique_sites();

      fflush(0);

      printf(" The sequence data has %d gapless sites out of %d sites in the alignment.\n",num_no_gaps,nsite);

      // Start computations of scores

      // Get sum of squared counts, i.e. the square of the Froebenius norm 
      for (i=0; i<num_unique; i++) fnorm2 = fnorm2 + pow((double)site_counter[i]/(double)num_no_gaps,2);

      Data_Bin();

      fscanf(splitsfile, "%d", &nsplits);
      fprintf(scoresfile, "%d\n", nsplits);

      printf(" Analyzing the splits provided in the file splitsfile ....\n\n");

      // loop on number splits in file
      for (split_num=0; split_num<nsplits; split_num++)
	{
	  // get size of split from file
	  fscanf(splitsfile, "%d", &split_size);
	  fprintf(scoresfile, "%d ", split_size);

	  // zero out contents of split_taxa 
	  for (jj=0; jj<nints; jj++) split_taxa[jj] = 0;

	  // get taxa in split from file and create mask for split_taxa
	  for (i = 0; i < split_size; i++)
	    {
	      fscanf(splitsfile, "%d\n", &taxon);

	      kk=(taxon-1)/nfit;
	      k=taxon-1-kk*nfit;

	      // flip 2 bits for mask for this taxon
	      split_taxa[kk] ^= (3 << 2*k);

	    }

	  // get SVD score for current split
	  getSVDscore(split_taxa, split_size);
	}

      free_Memory_bin();

      break;

    case 1:   // Sliding window analysis

      // Compute the number of windows minus one
      num_blocks = (nsite-blocksize)/slidesize;

      fscanf(splitsfile, "%d", &nsplits);

      if (nsplits!=1) {
	printf(" The sliding window option only handles one split at a time.\n Please modify splitsfile to include only a single split.\n Exiting....\n\n");
	exit(1);
      }

      fscanf(splitsfile, "%d", &split_size);

      // zero out contents of split_taxa                                                                                              
      for (jj=0; jj<nints; jj++) split_taxa[jj] = 0;


      // get taxa in split from file and create mask for split_taxa                                                                   
      for (i = 0; i < split_size; i++)
	{
	  fscanf(splitsfile, "%d\n", &taxon);

    	  kk=(taxon-1)/nfit;
	  k=taxon-1-kk*nfit;

	  // flip 2 bits for mask for this taxon                                                                                      
	  split_taxa[kk] ^= (3 << 2*k);

	}


      for (nn=1; nn<num_blocks+2; nn++) {

	// transfer data for nnth block to ppBase; then proceed as usual
	for (mm=0; mm<nseq; mm++) {

	  for (pp=(nn-1)*slidesize; pp<(nn-1)*slidesize+blocksize; pp++) {
	    ppBase[mm][pp-((nn-1)*slidesize)] = ppBase_full[mm][pp];

	  }
	}

	//  Find all unique site patterns and store them in ppBase_unique  
	unique_sites();

        if (num_no_gaps >= min_sites) {

	  // Start computations of scores
	  fnorm2 = 0.0;
	  // Get sum of squared counts, i.e. the square of the Froebenius norm 
	  for (i=0; i<num_unique; i++) fnorm2 = fnorm2 + pow((double)site_counter[i]/(double)num_no_gaps,2);

	  Data_Bin();

	  // get SVD score for current split

	  fprintf(scoresfile,"%d %d %d ",1+(nn-1)*slidesize,num_no_gaps,num_const);
	  getSVDscore(split_taxa, split_size);

	  free_Memory_bin();
	}
      }

      break;

    default:
      printf(" %s not called correctly.  Exiting ....\n", argv[0]);
      exit(1);
    }

  fclose(splitsfile);
  fclose(scoresfile);

  printf(" Analysis complete. Results have been written to scoresfile.\n\n\n");

  free_Memory();

  if (showTimeElapsed==1) {
    end=clock();
    double cc = CLOCKS_PER_SEC;
    double sec =(end-start)/cc;
    printf(" Number of seconds elapsed is %0.2f\n\n",sec);   
  }       
 
  exit(0);   

}
