/* main.c -- 
 *    USAGE: 
 *  PURPOSE: Perform all possible (legitimate) first-order splits.
 *           Each case has 6 (10?) different loci with 2 alleles at each.
 *           a[j][l][i] is the jth allele at locus l for case i
 *           At any split, we choose l = L and have a list of alleles ind[k].
 *           If a[j][L][i] = ind[k] for either j and any k, i belongs in 
 *           primary group.  If not, i belongs in secondary group.
 *           AUTHORS: Tom Kepler and Peter Hraber (11/30/2001)
 */

#define	maxClass 2
#define VERBOSE  1
#define DEBUG    0

//double log2;

#include <math.h>
#include <stdio.h>  // for printf
#include <stdlib.h> // for qsort - used to Winsorize mean
#include <string.h> // for memcpy

//#include "patient.h"
//#include "HLAAsup.h"
//#include "HLABsup.h"
//#include "mytest.h"
#include "n379.h"
#include "etc.h" // this must follow the data-file inclusion

int
main(int argc, char** argv)
{
  extern int X[nObs][nCov];
  extern double y[nObs]; // raw data
  int nClass = maxClass;
  double thresh;
  double Y[nObs]; // transformed data 
  moments* g;
  int n[nClass];

  transformData(y,Y,nObs,0.); /* use theta=0. for log-transform */

  g = malloc(sizeof(moments)); /* should error trap here */
  g->n = nObs;

  g = Moments(Y,g);
  printf(" g->mean=%f\t g->var=%f\tL=%f\n",g->mean,g->var, L0(g));

  // use estimators adjusted to censor below value of thresh
  thresh=1.;
  g = WDMoments(Y,g,log2(thresh));
  printf("g->Wmean=%f\tg->Dvar=%f\tL=%f\t(thresh=%f.)\n",g->mean,g->var,L0(g),thresh);

  thresh=300.;
  g = WDMoments(Y,g,log2(thresh));
  printf("g->Wmean=%f\tg->Dvar=%f\tL=%f\t(thresh=%f.)\n",g->mean,g->var,L0(g),thresh);

  thresh=400.;
  g = WDMoments(Y,g,log2(thresh));
  printf("g->Wmean=%f\tg->Dvar=%f\tL=%f\t(thresh=%f.)\n",g->mean,g->var,L0(g),thresh);

  thresh=500.;
  g = WDMoments(Y,g,log2(thresh));
  printf("g->Wmean=%f\tg->Dvar=%f\tL=%f\t(thresh=%f.)\n",g->mean,g->var,L0(g),thresh);

  thresh=1000.;
  g = WDMoments(Y,g,log2(thresh));
  printf("g->Wmean=%f\tg->Dvar=%f\tL=%f\t(thresh=%f.)\n",g->mean,g->var,L0(g),thresh);

  free(g);
  g=0;
  return(0);
} /* MAIN */
/***************************************************************************/
