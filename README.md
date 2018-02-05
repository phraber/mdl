## MDL

### Application of the minimum description length (MDL) principle to statistical genetics

### Code and data to support analysis of HLA allele associations with HIV setpoint.

This directory contains data and analytic source code to accompany two manuscripts.  Please cite the first one if you use it or a derivative; please cite the second one if you use the data for your own analyses.

[P.T. Hraber, B.T. Korber, S. Wolinksy, H. Erlich, E. Trachtenberg, and T.B. Kepler.  2003.  HLA and HIV infection progression: Application of the minimum description length principle to statistical genetics.  *Lecture Notes in Bioinformatics* 4345: 1-12, 2006.](https://link.springer.com/content/pdf/10.1007%2F11946465_1.pdf)

[E.A. Trachtenberg, B.T. Korber, C. Sollars, T.B. Kepler, P.T. Hraber, E. Hayes, R. Funkhouser, M. Fugate, J. Theiler, M. Hsu, K. Kunstman, S. Wu, J. Phair, H.A. Erlich, and S. Wolinsky.  2003.  HLA-A and HLA-B supertype alleles predict human immunodeficiency virus disease progression rate.  *Nature Medicine* 9(7):928-35. PMID:12819779 DOI:10.1038/nm893](https://www.nature.com/articles/nm893)

Before you compile, please examine the file main.c; it contains several important #include and #define statements near the top, before the beginning of the main function.

Data are provided as header (.h) files.  To analyze any given set of data, uncomment the corresponding #include line.  Only one such line may be defined (uncommented), or the compiler will be unhappy.  Please do not change the #define statements within a data/header file or, if you must, please do so only at your own risk.  Comments about the data format are contained in each header file.

To set the value of k (number of groups), modify this line (in main.c):
```c
#define	maxClass  2
```

To generate prolix output, ensure that both of these are true:
```c
#define VERBOSE   1
#define DEBUG     1
```

Setting the first to 1 and the second to 0 is good for most purposes.

To compile, simply type: gcc main.c -o mdl -lm

You may opt to rename the output file and/or set more compiler flags.

This source code is provided without warrantee.  Feedback, however,
especially of the positive variety, is most welcome, as is funding.
The code is admittedly ugly and suboptimal in places.  I do apologize.

See LICENSE.TXT for copyright details.

Peter T. Hraber
17 March 2003

PS: Only the last solution is reported as best, in case of ties.

## MANIFEST

### Source code
+ main.c: primary flow-of-control for analysis routine
+ etc.h: auxiliary analysis routines
+ mytest.c: simple file that computes summary stats for testing purposes.

### Data files
+ n479.h: 2-digit allele codes for Class II & Class I loci
+ n379.h: 2-digit allele codes for Class II & Class I loci, Caucasian subset 
+ HLAAsup.h: HLA-A supertypes (n=399, N=4)
+ HLABsup.h: HLA-B supertypes (n=352, N=5)

The above data were analyzed in the first paper cited.

+ patients.h: HLA "supertype genotype" alleles (n=293, N=15) as analyzed in the second paper cited above
+ mytest.h: a hypothetical set of data, for testing purposes

## That's all!
