// PURPOSE: Perform k-way clustering analysis to identify HLA loci
// associated with variation in HIV RNA levels, using the description
// length as a relative measure of merit.

// USAGE: program_name <locus> <min partition> <max partition> <Lmin>, where

// <locus> is an integer specifies at which locus to partition the
// data, beginning with 1, and ending with the number of loci defined
// in the data header file, and <min partition> and <max partition>
// are integers that specify a subset of partitions for analysis -
// ordering is arbitrary, though consistent for any given data header
// file.  This was implemented to facilitate parallel execution for
// loci having large allelic diversity.  Default valus are the first
// locus and all possible partitions.  <Lmin> is a floating-point
// value for a description length cutoff; for partitions that yield
// description lengths above this value, no results will be reported;
// default is L0, corresponding to the description length for the
// 'null hypothesis' of all observations in one group.

// NB: We are not now using recursion to compute the description
// lengths for more than two groups but you can simulate this easily
// enough by recompiling with a greater value of maxClass

// AUTHORS: Tom Kepler and Peter Hraber (11/2001 - 3/2003)

// specify header file containing data
#include "patients.h" // concatenated supertype genotypes
//#include "HLAAsup.h"
//#include "HLABsup.h"
//#include "n379.h"
//#include "n479.h"

#define	maxClass  2
#define VERBOSE   1
#define DEBUG     1

#include <math.h>
#include <stdlib.h> // for qsort - used to Winsorize mean
#include <stdio.h>  // for fopen, sprintf, etc.
#include <string.h> // for memcpy
#include "etc.h" // define the data file above this, because it sets #defines used therein

int
main(int argc, char** argv)
{
  extern int X[nObs][nCov];  // stores allele identifiers
#if TRANSFORMED
  extern double Y[nObs];     // stores transformed viral load levels
#else
  double Y[nObs];        // stores transformed viral load levels
  extern double y[nObs]; // stores raw viral load levels
#endif
  double Lmin, L;
  moments* g;
  unsigned long m, p, q, mStar;
  int idMap[nObs];           // in indexAlleles and at end of main
  int a[nObs*(nCov/nLoci)];  // in indexAlleles and at end of main
  // better explain use of nCov/nLoci
  int iClass[nObs];
  int nAlleles;
  int lStar;
  int ind[nObs];
  int is_mapped;
  int isEmpty;
  int is_lower = 1;
  int i, j, k, l;
  int nSplits;
  int nClass = maxClass;
  int n[nClass];
  double IcostA,IcostB;

  // data in file patients.h have already been log_2-transformed
#if !TRANSFORMED
  transformData(y,Y,nObs,0.); /* use theta=0. for log-transform */
#endif

  g = malloc(sizeof(moments)); /* should error trap here */
  g->n = nObs;
#if CENSORED
  g = WDMoments(Y,g,THRESH);
#else
  g = Moments(Y,g);
#endif
  Lmin = L0(g);

  /* 'PARSE' COMMAND-LINE ARGUMENTS */
  l=0;
  if (argc>=2) l=atoi(argv[1])-1; /* set locus from command line */

  if (l >= nLoci) { 
    printf("ERROR: cannot analyze locus %d with only %d loci in data!\n", 
	   l+1, nLoci); 
    return(1);
  }

  nAlleles = indexAlleles(l,nObs,idMap,a);

#if VERBOSE
  printf("L0=%f\tmean=%f\tvar=%f\tnAlleles=%d\tnObs=%d\n",
	 Lmin,g->mean,g->var,nAlleles,nObs);
#endif
  free(g);
  g=0;
  p=0;
  q=pow(nClass,nAlleles)-1;

  if (argc >= 3) {
    p=atol(argv[2]); 
    q=p; 
  }

  if (argc >= 4)
    q=atol(argv[3]);

  if (p > q) { // swap them
    m=q;
    q=p;
    p=m;
  }

  if (q > pow(nClass,nAlleles)-1)
    q=pow(nClass,nAlleles)-1;

  IcostB=0;

  IcostA=nAlleles*log2(nClass);
  if (argc>=5) Lmin=atof(argv[4]); 

  nSplits = 0;

  for (m=p; m<=q; m++) { // iterate through partitions

#if DEBUG
    printf("%lu\t",m);
#endif
    for (k=0; k < nAlleles; k++) {
      ind[k] = (m / (int)pow(nClass, k)) % nClass; // map alleles to classes
#if DEBUG
      printf(" %d",ind[k]);
#endif
    }
#if DEBUG
    printf("\n");
#endif


// HAPLOID data
#if (nCov == nLoci)
    for (k=0; k < nClass; k++) {
      n[k] = 0;
      isEmpty = 1;

      // assign cases to classes and count cases per class
      for (i=0; i < nObs; i++) {
	if (ind[a[i]] == k) {
	  iClass[i] = k;
	  n[k]++;
	  isEmpty=0;
	}
      }
      if (n[k] == 1) /* this will yield zero variance */
	isEmpty = 1;

      if (isEmpty) break;
    }
#elif (maxClass > 2)
    // DIPLOID DATA, k > 2
    for (k=nClass-1; k >= 0; k--) {
        // assign cases to classes and count cases per class
        for (i=0; i < nObs; i++)
	    if (ind[a[i]] == k || ind[a[i+nObs]] == k)
	        iClass[i] = k;
    }
    isEmpty=0;
    for (k=0; k < nClass; k++)
        n[k] = 0;

    // count cases per class
    for (i=0; i < nObs; i++)
        n[iClass[i]]++;

    for (k=0; k < nClass; k++)
        if (n[k] <= 1) 
            isEmpty = 1;
#else
    // DIPLOID DATA, k==2
    isEmpty=0;

    n[0] = 0;
    //    n[1] = 0;
    /* assign classes and counts per class */
    for (i=0; i < nObs; i++) {
      if (ind[a[i]] == 0 || ind[a[i+nObs]] == 0) {
	iClass[i] = 0;
	n[0]++;
      }
      else {
	iClass[i] = 1;
	//	n[1]++;
      }
    } /* for i=1..nObs */
    n[1]=nObs-n[0];

    if (n[0]<=1 || n[1]<=1) isEmpty=1; /* n<=1 -> var=0 -> log(0) -> NaN */
#endif
    //    printf ("nObs=%d\tn[0]=%d\tn[1]=%d\n",nObs,n[0],n[1]);
    if (isEmpty) continue; /* skip partitions having empty bins*/

    nSplits++;

#if (maxClass == 2)
    L = computeL2(iClass,Y,nObs);
#else
    L = computeL(iClass,Y,nObs,nClass);
#endif

    if (L <= Lmin) { /* keep parameters assoc. with best L */
      Lmin = L;
      mStar = m;
      lStar = l;
      is_lower = 1;

#if VERBOSE
      printf("L'=%g\tl'=%d\tm'=%lu\n",Lmin,lStar+1,mStar);
#endif
    }
  } // end of partition iteration loop

  if (is_lower) {
    IcostB = log2(nSplits);

#if DEBUG
    printf("%lu\t",mStar);
#endif
    for (k=0; k < nAlleles; k++) {
      ind[k] = (mStar / (int)pow(nClass, k)) % nClass;
#if DEBUG
      printf(" %d",ind[k]);
#endif
    }
#if DEBUG
    printf("\n");
#endif

#if (nCov == nLoci)
    // re-assign classes and counts per class
    for (k=0; k < nClass; k++)
      for (i=0; i < nObs; i++)
	if (ind[a[i]] == k)
	  iClass[i] = k;
#elif (maxClass > 2)
    // DIPLOID DATA, k > 2
    for (k=nClass-1; k >= 0; k--) {
        // assign cases to classes and count cases per class
        for (i=0; i < nObs; i++)
	    if (ind[a[i]] == k || ind[a[i+nObs]] == k)
	        iClass[i] = k;
    }
    isEmpty=0;
    for (k=0; k < nClass; k++)
        n[k] = 0;

    // count cases per class
    for (i=0; i < nObs; i++)
        n[iClass[i]]++;

    for (k=0; k < nClass; k++)
        if (n[k] <= 1) 
            isEmpty = 1;
#else
    // DIPLOID DATA, k==2
    isEmpty=0;
    n[0] = 0;
    //    n[1] = 0;
    /* assign classes and counts per class */
    for (i=0; i < nObs; i++) {
      if (ind[a[i]] == 0 || ind[a[i+nObs]] == 0) {
	iClass[i] = 0;
	n[0]++;
      }
      else {
	iClass[i] = 1;
	//n[1]++;
      }
    } /* for i=1..nObs */
    n[1]=nObs-n[0];
#endif

    printf ("L*=%g\tL'=%f  L\"=%f\tl*=%d\tm*=%lu\n",Lmin,Lmin+IcostA,Lmin+IcostB,lStar+1,mStar);

#if (maxClass == 2)
    showL2(iClass,Y,nObs);
#else
    showL(iClass,Y,nObs,nClass);
#endif
    showClasses(ind,idMap,nAlleles,nClass);
    showCases(iClass,Y,nObs,l,mStar,nClass);
  }
  return(0);
} /* MAIN */
/***************************************************************************/
