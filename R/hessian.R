# Ingmar Visser, 17 oktober 2018
# 
# Description
# The function hessian computes the hessian of a (dep)mix object at the current 
# parameter values; it has a method argument; the only method currently implemented
# is the finite differences method. 
#
# Details
# 
# The function optionally accepts arguments related to linear constraints, which are 
# neccessary to have when the hessian is used to compute the variance-covariance 
# matrix of the parameters. 
# 
# The function also checks whether parameters are estimated on the boundary and 
# leaves them out of the process when this is the case. 
# 
# Fixed parameters are similarly ignored when computing the hessian
# 
# Value
#
# The function returns a named list with the following elements:
# 
# hessian: the hessian of the parameters
# 
# elements: vector of length npar(object) indicating for each parameter whether it 
# is 'inc'luded, 'fix'ed, or estimated on the boundary, 'bnd'; the dimension of the hessian 
# is thus the number of non-fixed parameters minus the number of boundary parameters. 
#

setMethod("hessian", "mix",
    function(object, 
	tolerance=1e-6, 	
	method="finiteDifferences", ...) {

	if(is.nan(logLik(object))) stop("Log likelihood is 'NaN'; cannot compute hessian. ")
			
	# get the full set of parameters
	allpars <- getpars(object)
	
	# get parameter boundaries from the model object
	constraints <- getConstraints(object)
	par.u <- constraints$par.u
	par.l <- constraints$par.l
	
	# get fixed parameters
	fixed <- getpars(object,"fixed")
	
	# return vector with specification of 'inc'luded, 'fix'ed and 'bnd'ary parameters
	elements <- rep("inc",npar(object))
	
	# identify parameters that are on their boundary
	low <- which(sapply(as.numeric(allpars-par.l),all.equal,tolerance=tolerance,0)==TRUE)
	up <- which(sapply(as.numeric(allpars-par.u),all.equal,tolerance=tolerance,0)==TRUE)
	bnd <- union(low, up)
	
	# identify parameters that are fixed or on the boundary
	if(length(which(fixed)>0)) elements[which(fixed)] <- "fix"
	if(length(bnd)>0) elements[bnd] <- "bnd"
	
	# get the reduced set of parameters, ie the ones that the hessian will be computed for
	# only non-fixed parameters
	pars <- allpars[which(elements=="inc")]	
		
	# make loglike function that only depends on pars
	logl <- function(pars) {
		allpars[which(elements=="inc")] <- pars
		object <- setpars(object,allpars)
		ans <- -as.numeric(logLik(object))
		if(is.na(ans)) ans = 1000000 # remove magic number here!!!
		ans
	}
	
	fdh <- fdHess(pars,logl)
		
	return(list(hessian=fdh$Hessian,elements=elements))
}
)


