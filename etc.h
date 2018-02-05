/***************************************************************************/
#define PI 3.141592653

typedef struct {
  double mean;
  double var;
  int n;
} moments;
/***************************************************************************/
static
int intcompare(const void *p1, const void *p2) {
  int i = *((int *)p1);
  int j = *((int *)p2);

  if (i > j)
    return (1);
  if (i < j)
    return (-1);
  return (0);
}
/***************************************************************************/
static
int doublecompare(const void *p1, const void *p2) {
  double i = *((double *)p1);
  double j = *((double *)p2);

  if (i > j)
    return (1);
  if (i < j)
    return (-1);
  return (0);
}
/***************************************************************************/
double
Mean(double* z, int n) {
  int i;
  double sum = 0.;

  for (i=0; i < n; i++)
    sum += z[i];

  return (sum /= (double)n);
}
/***************************************************************************/
double
Median(double* z, int n) {

  qsort(z,n,sizeof(double),doublecompare);
  // better test/check this!
  if (n % 2) // n odd
    return z[(int)(n/2)];
  else       // n even
    return (z[(int)(n/2)]+z[(int)(n/2)-1])/2.;
}
/***************************************************************************/
double Var(double* z, int n, double mean) {
  int i;
  double var = 0.;

  for (i=0; i < n; i++)
    var += pow(z[i],2);

  return (var - (double)n*pow(mean,2))/((double)n-1.);
}
/***************************************************************************/
moments*
Moments(double* z, moments* m) {

  int i;
  double sum = 0.;
  double var = 0.;

  for (i=0; i < m->n; i++) {
    sum += z[i];
    var += z[i]*z[i];
  }

  m->mean = sum / (double)m->n;
  m->var = (var - (double)m->n*pow(m->mean,2))/((double)m->n - 1.);

  return m;
}
/***************************************************************************/
moments*
WDMoments(double* z, moments* m, double theta) {

  /* Computes unbiased estimators for censored data;
   * theta is the lower observation threshold, or truncation point.
   *
   * When no observations are censored (<= theta), simply returns
   * the sample mean and variance.
   *
   * Winsorized mean estimator yields bias-corrected mean for j missing 
   * values from a (truncated/censored) normal distribution.
   * The Winsorized mean replaces the j-most outlying values with values
   * at j+1 (for left tail) or n-j (for right tail).
   * REF: Norman L. Johnson, Samuel Kotz, and N. Balakrishnan. (1994) 
   * Continuous Univariate Distributions. vol. I, second ed. 
   * Wiley Intersciences, NY, pp. 124, 146ff, & 156ff.
   * 
   * Variance is estimated using (the square of) Healy's adjustment
   * to Downton's linear estimator of the standard deviation.  
   * REFS:
   * MJR Healy. (1978). A mean difference estimator of standard deviation
   * in symmetrically censored normal samples. Biometrika 65(3):643-646.
   * MJR Healy. (1982). A linear difference estimator of standard deviation
   * in symmetrically censored normal samples. Applied Statistics 31(2):174-5.
   *
   * The value returned is a linear estimate of the standard deviation
   * from a sample of N values of X symmetrically trimmed to retain K values.
   */

  int i;
  int j=0;
  int k;
  int q;
  //  int rank[m->n];
  double b,d,p;
  double kprime;
  double sum=0.;
  double std=0.;
  double w[m->n];

  m->mean = 9E+99;
  m->var  = 9E+99;

  for (i=0; i < m->n; i++)
    if (z[i] <= theta)
      j++;

  --j; // include one observation at threshold value

  k = m->n - (2*j); /* number of observations remaining */
  p = k/(double)m->n; /* uncensored proportion of sample */

  if (p >= 0.98)
    m = Moments(z,m);

  else {
    memcpy(w,z,m->n*sizeof(double));
    qsort(w,m->n,sizeof(double),doublecompare); /* sort observations */

    kprime=(double)k*(k-0.5);

    sum += (double)j*w[j]; // for strict Winsorization

    for (i=j; i < m->n - j; i++) {
      sum += w[i];
      std += (2.*(i+1.-j) - k - 1.)*w[i]/kprime; // Healy '78
    }

    sum += (double)j*w[m->n-j];

    if (p > 0.5) {
      m->mean = sum / (double)m->n;
#if CALCDVAR
      p = 0.74 - p;
      b = exp(((((17.5265*p-3.91154)*p+1.32458)*p+0.403049)*p+1.72534)*p+1.06683); // Healy '82
#endif
    }
    else {
      if (p > 0.2)
	m->mean = sum / (double)m->n;
      else
      	m->mean = Median(w,m->n);

#if CALCDVAR
      b = exp(-1.01295*log(p)+0.837567);
#endif
    }

    //    printf("b=%f\n",b);

#if CALCDVAR
    m->var = pow(b*std, 2);
#endif
  }
  // pth - 01/28/03
#if !CALCDVAR
  m->var = Var(z,m->n,m->mean); // this works best for computing variance
#endif

  return m;
}
/***************************************************************************/
void
transformData(double* y, double* Y, int n, double lambda) {
  int i;

  if (lambda == 0.)
    for (i = 0; i < n; i++)
      Y[i]=log2(y[i]);
  else
    for (i = 0; i < n; i++)
      Y[i]=(pow(y[i],lambda)-1)/(lambda);
}
/***************************************************************************/
void
atransData(double* y, double* Y, int n) {
  int i;

  for (i=0; i<n; i++)
    Y[i]=atan(y[i]*PI);
}
/***************************************************************************/
double calcL(moments** u, int k) {

  double L=0.;
  int i;

  for (i=0; i < k; i++) {

    L += (u[i]->n)/2.;

    /* error variance */
    L += (u[i]->n - 1.)*log2(u[i]->var)/2.;

    /* variance of means */
    L += log2(u[i]->n * pow(u[i]->mean,2))/2.;

    L += log2(u[i]->n);
  }

  return L;
}
/***************************************************************************/
double calcLa(moments** u) {

  double L;

  L = u[0]->n + u[1]->n;

  /* error variance */
  L += (u[0]->n - 1.)*log2(u[0]->var);
  L += (u[1]->n - 1.)*log2(u[1]->var);

  /* variance of means */
  L += log2(u[0]->n * pow(u[0]->mean,2));
  L += log2(u[1]->n * pow(u[1]->mean,2));

  L /= 2.;

  L += log2(u[0]->n);
  L += log2(u[1]->n);

  return L;
}
/***************************************************************************/
double calcLb(moments** u) {

  double L;

  L = u[0]->n + u[1]->n;

  /* error variance */
  L += (u[0]->n - 1.)*log2(u[0]->var);
  L += (u[1]->n - 1.)*log2(u[1]->var);

  /* variance of means */
  L += log2((pow(pow(u[0]->mean,2) - pow(u[1]->mean,2),2))*(u[0]->n * u[1]->n)/4.);

  L += log2(u[0]->n);
  L += log2(u[1]->n);

  L /= 2.;

  L += log2(u[0]->n + u[1]->n);

  return L;
}
/***************************************************************************/
int
indexAlleles(int l, int n, int* id_map, int* a) {

  int i,j;
  int is_mapped,indx;
  int nAlleles = 0;

  /* map arbitrary locus identifiers to continuous integer indices */

  for (i=0; i < n; i++) { /* i indexes observations */
    for (j=0; j < nCov/nLoci; j++) { /* j indexes diploid loci */

      is_mapped = 0;
      /* scan id_map to see whether orig_id has already been indexed */
      for (indx = 0; indx <= nAlleles; indx++) {
	if (id_map[indx] == X[i][l*(nCov/nLoci)+j]) {
	  a[(j*nObs)+i] = indx; // okay for supertype genotypes
	  is_mapped = 1;
	  break;
	}
      }
      if (!is_mapped) {
	a[(j*nObs)+i] = nAlleles;
	id_map[nAlleles] = X[i][l*(nCov/nLoci)+j];
	nAlleles++;
      }
    } /* j=0,nCov */
  } /* i=1,n */

#if DEBUG
  printf("Indexed Alleles:\n");
  for (i = 0; i<nAlleles; i++) {
    printf("\t%d\t%d\n",i,id_map[i]);
  }
#endif

  return nAlleles;
}
/****************************************************************************/
double L0(moments* m) {
  /* should do error checking, to avoid invalid operations */
  return (m->n + (m->n - 1.)*log2(m->var) + log2(m->n*m->mean*m->mean))/2. + log2(m->n);
}
/***************************************************************************/
double
computeL2(int* iClass, double* z, int nobs) {

  int i,k;
  int n[2]; // we know nClass=2
  moments* u[2];
  double obs[2][nobs];
  double La;
#if CALCL2B
  double Lb;
#endif
#if CENSORED
  double theta=THRESH;
#endif

  /* initialize */
  for (k=0; k < 2; k++) {
    u[k] = malloc(sizeof(moments));
    u[k]->n = 0;
  }

  /* construct observation vectors */
  for (i=0; i < nobs; i++)
    obs[iClass[i]][u[iClass[i]]->n++] = z[i];

  /* compute sample means (mean[i]) and error variances (var[i]) */
  for (k=0; k < 2; k++) {
#if CENSORED
    u[k] = WDMoments(obs[k], u[k], theta);
#else
    u[k] = Moments(obs[k], u[k]);
#endif
    if (u[k]->var <= 0.) {

      for (i=0; i < 2; i++) {
	free(u[i]);
	u[i] = 0;
      }

      return 9E+99;
    }
  }

  /* compute L */
  La = calcLa(u);

#if CALCL2B
  if (pow((u[0]->mean*u[0]->mean - u[1]->mean*u[1]->mean),2)*(u[0]->n*u[1]->n)/4. > 0.)
    Lb = calcLb(u);
  else
    Lb = 9E+99;
#endif

#if (CALCL2B && DEBUG)
  printf ("%7.2f  %7.2f  0:( %3d %7.4f %7.4f ) 1:( %3d %7.4f %7.4f )\n",
	  La,Lb,u[0]->n,u[0]->mean,u[0]->var,u[1]->n,u[1]->mean,u[1]->var);
#elif DEBUG
  printf ("%7.2f  0:( %3d %7.4f %7.4f ) 1:( %3d %7.4f %7.4f )\n",
	  La,u[0]->n,u[0]->mean,u[0]->var,u[1]->n,u[1]->mean,u[1]->var);
#endif

  for (k=0; k < 2; k++) {
    free(u[k]);
    u[k] = 0;
  }

#if CALCL2B
  return (La < Lb) ? La : Lb;
#else
  return La;
#endif
}
/***************************************************************************/
double computeL(int* iClass, double* z, int nobs, int nClass) {

  int i,k;
  int n[maxClass];
  moments* u[maxClass];
  double obs[maxClass][nobs];
  double L;
#if CENSORED
  double theta=THRESH;
#endif
  /* initialize */
  for (i=0; i < nClass; i++) {
    u[i] = malloc(sizeof(moments));
    u[i]->n = 0;
  }

  /* construct observation vectors */
  for (i=0; i < nobs; i++)
    obs[iClass[i]][u[iClass[i]]->n++] = z[i];

  /* compute sample means (mean[i]) and error variances (var[i]) */
  for (i=0; i < nClass; i++) {
#if CENSORED
    u[i] = WDMoments(obs[i], u[i], theta);
#else
    u[i] = Moments(obs[i], u[i]);
#endif
    if (u[i]->var <= 0.) {
      for (k=0; k < nClass; k++) {
	free(u[k]);
	u[k] = 0;
      }
      return 9E+99;
    }
  }
  /* compute L */
  L = calcL(u, nClass);
#if DEBUG
  printf ("%7.2f\n", L);
  for (i=0; i < nClass; i++)
    printf("\t%d: ( %3d %7.4f %7.4f )\n",i,u[i]->n,u[i]->mean,u[i]->var);
#endif
  for (i=0; i < nClass; i++) {
    free(u[i]);
    u[i] = 0;
  }
  return L;
}
/***************************************************************************/
void
showL(int* iClass,double* z,int nobs,int nClass) {

  int i;
  int n[maxClass];
  moments* u[maxClass];
  double obs[maxClass][nobs];
  double L;
#if CENSORED
  double theta=THRESH;
#endif
  /* initialize */
  for (i=0; i < nClass; i++) {
    u[i] = malloc(sizeof(moments));
    u[i]->n = 0;
  }

  /* construct observation vectors */
  for (i=0; i < nobs; i++)
    obs[iClass[i]][u[iClass[i]]->n++] = z[i];

  /* compute sample means (mean[i]) and error variances (var[i]) */
  for (i=0; i < nClass; i++) {
#if CENSORED
    u[i] = WDMoments(obs[i], u[i], theta);
#else
    u[i] = Moments(obs[i], u[i]);
#endif
    if (u[i]->var <= 0.) {

      for (i=0; i < nClass; i++) {
	free(u[i]);
	u[i] = 0;
      }
    }
  }

  /* compute L */
  L = calcL(u, nClass);

  printf ("L=%7.2f\n", L);
  for (i=0; i < nClass; i++)
    printf ("\t%d:( %3d %7.4f %7.4f )\n", i, u[i]->n, u[i]->mean, u[i]->var);

  for (i=0; i < nClass; i++) {
    free(u[i]);
    u[i] = 0;
  }
}
/***************************************************************************/
double
showL2(int* iClass, double* z, int nobs) {
#if VERBOSE
  int i,k;
  int n[2]; // we know nClass=2
  moments* u[2];
  double obs[2][nobs];
  double La;
#if CALCL2B
  double Lb;
#endif
#if CENSORED
  double theta=THRESH;
#endif

  /* initialize */
  for (k=0; k < 2; k++) {
    u[k] = malloc(sizeof(moments));
    u[k]->n = 0;
  }

  /* construct observation vectors */
  for (i=0; i < nobs; i++)
    obs[iClass[i]][u[iClass[i]]->n++] = z[i];

  /* compute sample means (mean[i]) and error variances (var[i]) */
  for (k=0; k < 2; k++) {
#if CENSORED
    u[k] = WDMoments(obs[k], u[k], theta);
#else
    u[k] = Moments(obs[k], u[k]);
#endif
    if (u[k]->var <= 0.) {
      for (i=0; i < 2; i++) {
	free(u[i]);
	u[i] = 0;
      }
      return 9E+99;
    }
  }

  /* compute L */
  La = calcLa(u);

#if CALCL2B
  if (pow((u[0]->mean*u[0]->mean - u[1]->mean*u[1]->mean),2)*(u[0]->n*u[1]->n)/4. > 0.)
    Lb = calcLb(u);
  else
    Lb = 9E+99;

  printf ("%7.2f  %7.2f  0:( %3d %7.4f %7.4f ) 1:( %3d %7.4f %7.4f )\n",
	  La,Lb,u[0]->n,u[0]->mean,u[0]->var,u[1]->n,u[1]->mean,u[1]->var);
#else
  printf ("%7.2f  0:( %3d %7.4f %7.4f ) 1:( %3d %7.4f %7.4f )\n",
	  La,u[0]->n,u[0]->mean,u[0]->var,u[1]->n,u[1]->mean,u[1]->var);
#endif


  for (k=0; k < 2; k++) {
    free(u[k]);
    u[k] = 0;
  }

#if CALCL2B
  return (La < Lb) ? La : Lb;
#else
  return La;
#endif
#endif
}
/***************************************************************************/
void
showCases(int* iClass,double* Y,int n,int l,long m,int nClass) {
#if DEBUG
  int i;
  char* fname;
  FILE* fstream;

  fname = malloc(80*sizeof(char));

  sprintf(fname,"cases-n%d-l%d-k%d-m%ld.out",nObs,l+1,nClass,m);

  fstream = fopen(fname,"w");
  for (i=0; i < n; i++)
    fprintf (fstream, "%d\t%f\n",iClass[i],Y[i]);

  fclose(fstream);
  free(fname);
#endif
}
/***************************************************************************/
void
showClasses(int* ind, int* idMap, int nAlleles, int nClass) {
#if VERBOSE
  int k,i;

  for (k=0; k < nClass; k++) {
      printf ("%d:( ",k);

      for (i=0; i < nAlleles; i++)
	if (ind[i]==k)
	  printf("%d ",idMap[i]);

      printf (")\n");

  }
#endif
}
/***************************************************************************/
int
split(unsigned long m,int nClass,int nAlleles,int* ind,int* a,int* iClass) {

  int i,k;
  int isEmpty; 
  int n[maxClass];

#if DEBUG
  printf("%lu\t",m);
#endif
  for (k=0; k < nAlleles; k++) {
    ind[k] = (m / (int)pow(nClass, k)) % nClass;
#if DEBUG
    printf(" %d",ind[k]);
#endif
  }
#if DEBUG
  printf("\n");
#endif


// PTH 01/28/03 - this appears to work properly only for HAPLOID data
#if (nCov == nLoci)
  for (k=0; k < nClass; k++) {
    n[k] = 0;
    isEmpty = 1;

    /* assign classes and counts per class */
    for (i=0; i < nObs; i++) {
      if (ind[a[i]] == k) {
 	iClass[i] = k;
 	n[k]++;
 	isEmpty=0;
      }
    }
    if (n[k] == 1) // skip partitions that yield zero variance
      isEmpty = 1;
    if (isEmpty) break;
  }
#elif (maxClass > 2)
  // DIPLOID data, k > 2
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
  // DIPLOID DATA (k==2)
  n[0] = 0;
  /* n[1] = 0; */
  /* assign classes and counts per class */
  for (i=0; i < nObs; i++) {
    if (ind[a[i]] == 0 || ind[a[i+nObs]] == 0) {
      iClass[i] = 0;
      n[0]++;
    }
    else {
      iClass[i] = 1;
    }
  } /* for i=1..nObs */
  n[1]=nObs-n[0];

  if (n[0]<=1 || n[1]<=1) isEmpty=1; /* n<=1 -> var=0 -> log(0) -> NaN */
#endif

  return (isEmpty); /* skip partitions having empty bins*/
}
/***************************************************************************/
