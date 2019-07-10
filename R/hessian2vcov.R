# 
# Ingmar Visser, Oct 10, 2018
# 
# Compute a corrected covariance using linear constraint matrix and hessian
# Note: the constraints should only be the linear equality constraints, not the 
# inequality constraints!!!
#

hessian2vcov <- function(hessian,lincon=NULL) {
	np <- dim(hessian)[1]
	if(!(dim(hessian)[1]==dim(hessian)[2])) stop("'hessian' should be a square matrix")
	if(!is.null(lincon)) { # deal with linear constraints
		nc <- ncol(lincon)
		if(np!=nc) stop("Nr of columns in linear constraint matrix not compatible with 'hessian' dimension")
		if(nrow(lincon)>0) {
			A <- lincon
			d <- hessian+t(lincon)%*%lincon
			di <- try(solve(d),silent=TRUE)
			if(class(di)=="try-error") {
				warning("Hessian singular, ses could not be computed.") 
				vcov <- 0 
			} else {
				ada <- A%*%di%*%t(A)
				adai <- try(solve(ada),silent=TRUE)
				if(class(adai)=="try-error") {
					warning("Near-singular hessian, ses may be bad.\n")
					diag(ada) <- diag(ada)*1.000001
					adai <- try(solve(ada))
					if(class(adai)=="try-error") {
						warning("Corrected hessian also singular, ses computed without contraints.\n")
					} else {
						vcov <- di-di%*%t(A)%*%adai%*%A%*%di
					}
				} else {
					vcov <- di-di%*%t(A)%*%adai%*%A%*%di
				} 
			}
		} else { # linear constraint matrix present but 0 rows
			vcov <- solve(hessian)
		}
	} else { # no linear constraints
		vcov <- solve(hessian)
	}
	vcov
}
