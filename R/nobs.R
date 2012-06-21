setMethod("nobs", signature(object="mix"),
	function(object, ...) {
		nt <- sum(object@ntimes)
		n <- sum(!apply(object@dens,1,function(x) any(is.na(x))))
		if(n!=nt) warning("missing values detected; nobs is number of cases without any missing values")
		return(n)
	}
)
