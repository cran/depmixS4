#
# Ingmar Visser, 19-10-2018
#
# Description
# 
# vcov computes the variance-covariance matrix of depmix model parameters
# 
# Details
# 
# The variance-covariance matrix of the parameters is computed using the hessian
# and the linear constraints in the model. First, the hessian is computed through
# the 'hessian' method for (dep)mix model. Second, the hessian and linear constraint
# matrix are passed to the 'hessian2vcov' function. 
# 
# Value
#
# The function returns a named list with the following elements:
#
# vcov: a square matrix of dimension the 'inc'luded parameters.
#
# elements: a vector of length npar(object) indicating 
# for each parameter whether it is 'inc'luded, 'fix'ed, or estimated on the 
# boundary, 'bnd'; the dimension of the vcov matrix is thus the number of 
# non-fixed parameters minus the number of boundary parameters. 
# 
# lincon: the linear constraint matrix used to compute the variance-covariance 
# matrix from the hessian; it only contains the parts of the linear constraint 
# matrix that relate to equality constraints; moreover, the columns related 
# to 'fix'ed and boundary ('bnd') parameters are left out. 
#

setMethod("vcov", "mix",
    function(object, fixed=NULL, equal=NULL, 
	conrows=NULL, conrows.upper=NULL, conrows.lower=NULL, 
	tolerance=1e-6, 	
	method="finiteDifferences", ...) {
	
	if(is.nan(logLik(object))) stop("Log likelihood is 'NaN'; cannot compute variance-covariance. ")
	
	# check for presence of constraints
	fi <- !is.null(fixed)
	cr <- !is.null(conrows)
	eq <- !is.null(equal)
	
	constr <- any(c(fi,cr,eq))
	
	# determine which parameters are fixed
	if(fi) {
		if(length(fixed)!=npar(object)) stop("'fixed' does not have correct length")
	} else {
		if(eq) {
			if(length(equal)!=npar(object)) stop("'equal' does not have correct length")
			fixed <- !pa2conr(equal)$free
		} else {
			fixed <- getpars(object,"fixed")
		}
	}
	
	# set those fixed parameters in the appropriate submodels
	object <- setpars(object,fixed,which="fixed")	
	
	# get standard constraints from (sub)models
	constraints <- getConstraints(object)

	lincon <- constraints$lincon
	lin.u <- constraints$lin.u
	lin.l <- constraints$lin.l
	par.u <- constraints$par.u
	par.l <- constraints$par.l
	
	if(class(object)=="depmix.fitted"|class(object)=="mix.fitted") {
		lincon <- object@conMat
		lin.u <- object@lin.upper
		lin.l <- object@lin.lower
	} 
	
	# incorporate equality constraints provided with the hessian function, if any
	if(eq) {
		if(length(equal)!=npar(object)) stop("'equal' does not have correct length")
		equal <- pa2conr(equal)$conr
		lincon <- rbind(lincon,equal)
		lin.u <- c(lin.u,rep(0,nrow(equal)))
		lin.l <- c(lin.l,rep(0,nrow(equal)))				
	}
	
	# incorporate general linear constraints, if any, via argument conrows
	if(cr) {
		if(ncol(conrows)!=npar(object)) stop("'conrows' does not have the right dimensions")
		lincon <- rbind(lincon,conrows)
		if(any(conrows.upper==0)) {
			lin.u <- c(lin.u,rep(0,nrow(conrows)))
		} else {
			if(length(conrows.upper)!=nrow(conrows)) stop("'conrows.upper does not have correct length")
			lin.u <- c(lin.u,conrows.upper)
		}
		if(any(conrows.lower==0)) {
			lin.l <- c(lin.l,rep(0,nrow(conrows)))
		} else {
			if(length(conrows.lower)!=nrow(conrows)) stop("'conrows.lower does not have correct length")
			lin.l <- c(lin.l,conrows.lower)
		}
	}
	
	# get the full set of parameters
	allpars <- getpars(object)

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

	# select only those columns of the constraint matrix that correspond to non-fixed parameters
	lincon <- lincon[,which(elements=="inc"),drop=FALSE]

	# remove redundant rows in lincon (all zeroes)
	allzero <- which(apply(lincon,1,function(y) all(y==0)))
	if(length(allzero)>0) {
		lincon <- lincon[-allzero,,drop=FALSE]
		lin.u <- lin.u[-allzero]
		lin.l <- lin.l[-allzero]
	}
	
	# remove rows of lincon with inequality constraints
	dflu <- lin.u-lin.l
	ineq <- which(dflu!=0)
	if(length(ineq)>0) {
		lincon <- lincon[-ineq,,drop=FALSE]
		lin.u <- lin.u[-ineq]
		lin.l <- lin.l[-ineq]
	}
	
	hs <- hessian(object)
	
	vc <- hessian2vcov(hs$hessian,lincon)
	
	return(list(vcov=vc,elements=elements,lincon=lincon))
	
}
)