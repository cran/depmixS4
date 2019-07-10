#
# Ingmar Visser, 19-10-2018
#
# Description
#
# Computes standard errors for (dep)mix model parameters. 
#
# Details
#
# Standard errors are computed through the variance-covariance matrix
# which in turn is computed using the hessian and the linear constraints
# of the model. See ?vcov and ?hessian for more details on this. 
#
# Value
#
# A data.frame with columns: pars=parameter values, constr=whether the 
# parameter is 'inc'luded, 'fix'ed or estimated on the bound, 'bnd', 
# and the column ses with the standard errors. 
#

setMethod("standardError", "mix",
    function(object, fixed=NULL, equal=NULL, 
	conrows=NULL, conrows.upper=NULL, conrows.lower=NULL, 
	tolerance=1e-6, 	
	method="finiteDifferences", ...) {
		
	vc <- vcov(object,fixed=fixed,equal=equal,
		conrows=conrows,conrows.upper=conrows.upper,conrows.lower=conrows.lower,
		tolerance=tolerance,method=method, ...)
	
	ses <- sqrt(diag(vc$vcov))
	
	pars <- getpars(object)
	
	elements <- vc$elements
	
	parsinc <- pars[which(elements=="inc")]
		
	ret <- data.frame(par=pars, constr=elements, se=NA)
	
	ret$se[which(elements=="inc")] <- ses
	
	return(ret)
}
)