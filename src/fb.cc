#include "fb.h"

extern "C" {

// 	# Forward-Backward algorithm (used in Baum-Welch)
// 	# Returns alpha, beta, and full data likelihood
// 	
// 	# A = T*K*K matrix with transition probabilities, from row to column!!!!!!!
// 	# B = T*K matrix with elements ab_{ij} = P(y_i|s_j)
// 	# init = K vector with initial probabilities
// 
// 	# NOTE: to prevent underflow, alpha and beta are scaled, using sca
// 	# NOTE: xi[t,i,j] = P(S[t] = j & S[t+1] = i) !!!NOTE the order of i and j!!!
	

// inputs are:
// 0) hom: whether the transition probs are homogeneous or not
// a) ns: the number of states
// b) nc: the number of cases
// c) nt: the number of rows of data
// d) ntimes: rows of data of individual cases
// 1) init: ns vector or ns by nc matrix
// 2) trdens: ns by ns matrix or nt by ns by ns array
// 3) dens: nt by ns matrix
// 
// outputs are:
// 1) alpha: nt by ns matrix
// 2) beta: nt by ns matrix
// 3) xi: nt by ns by ns array
// 4) sca: nt vector

// gamma is computed as alpha*beta/sca in R (no loop needed)

void forwardbackward(int *hom, int *ns, int *nc, int *nt, int *ntimes, int *bt, int *et,
					 double *init, double *trdens, double *dens, 
					 double *alpha, double *beta, double *sca, double *xi) {
		
	// compute forward variables
	// loop over cases
	for(int cas=0; cas<nc[0]; cas++) {
		
		// compute alpha1 for this case
		double sca1=0.0;
		for(int i=0; i<ns[0]; i++) {
			alpha[(bt[cas]-1)*ns[0]+i] = init[cas*ns[0]+i]*dens[(bt[cas]-1)*ns[0]+i];
			// compute scale for alpha1
			sca1 += alpha[(bt[cas]-1)*ns[0]+i];
		}
		sca[(bt[cas]-1)] = 1/sca1;
		// scale alpha1 and copy to alpha
		for(int i=0; i<ns[0]; i++) {
			alpha[(bt[cas]-1)*ns[0]+i] *= sca[(bt[cas]-1)];
		}
		
		// compute alphat without intervening matrices
		if(ntimes[cas]>0) {
			for(int t=bt[cas]; t<et[cas]; t++) {				
				// compute alphat
				for(int j=0; j<ns[0]; j++) {
					for(int i=0; i<ns[0]; i++) {
						if(hom[0]==1) alpha[t*ns[0]+j] += trdens[i*ns[0]+j]*alpha[(t-1)*ns[0]+i];
						else alpha[t*ns[0]+j] += trdens[(i*ns[0]+j)*nt[0]+t-1]*alpha[(t-1)*ns[0]+i];
					}
					alpha[t*ns[0]+j] *= dens[t*ns[0]+j];
					sca[t] += alpha[t*ns[0]+j];
				}
				sca[t] = 1/(sca[t]);				
				// scale alphat values
				for(int i=0; i<ns[0]; i++) {
					alpha[t*ns[0]+i] *= sca[t];
				}
			}
		}	
		
		// compute backward variables
		// initial backward vars for t=T
		for(int i=0; i<ns[0]; i++) {
			beta[(et[cas]-1)*ns[0]+i] = sca[et[cas]-1];
		}
		
 		// compute beta t-1
  		if(ntimes[cas]>1) {			
			// loop from T-1 to 1 for each case
 			for(int t=(et[cas]-1); t>=bt[cas]; t--) {				
				for(int i=0; i<ns[0]; i++) {
					for(int j=0; j<ns[0]; j++) {
						if(hom[0]==1) beta[(t-1)*ns[0]+i] += trdens[i*ns[0]+j]*dens[t*ns[0]+j]*beta[t*ns[0]+j];
						else beta[(t-1)*ns[0]+i] += trdens[(i*ns[0]+j)*nt[0]+t-1]*dens[t*ns[0]+j]*beta[t*ns[0]+j];
					}
					beta[(t-1)*ns[0]+i] *= sca[t-1];
				}
 			}
			
			// compute xi 
			for(int t=bt[cas]; t<et[cas]; t++) {
				for(int i=0; i<ns[0]; i++) {
					for(int j=0; j<ns[0]; j++) {
						if(hom[0]==1) xi[(i*ns[0]+j)*nt[0]+t-1] = alpha[(t-1)*ns[0]+i]*trdens[i*ns[0]+j]*dens[t*ns[0]+j]*beta[t*ns[0]+j];
						else xi[(i*ns[0]+j)*nt[0]+t-1] = alpha[(t-1)*ns[0]+i]*trdens[(i*ns[0]+j)*nt[0]+t-1]*dens[t*ns[0]+j]*beta[t*ns[0]+j];
					}
				}
			}			
 		}
						
	} // end cases
	
}

} // end extern "C"
