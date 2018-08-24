#ifndef FB
#define FB 1

/* 
 * #include <stdio.h>
 * #include <stdlib.h>
 */

// compute forward and backward variables, and xi
void forwardbackwardC(int *hom, int *ns, int *nc, int *nt, int *ntimes, int *bt, int *et, 
					 double *init, double *trdens, double *dens, 
					 double *alpha, double *beta, double *sca, double *xi);

#endif
