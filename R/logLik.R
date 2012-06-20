# depends on getpars and nobs
setMethod("logLik",signature(object="depmix"),
	#function(object,method="lystig") { 
	function(object,method="fb") { #4/5/2012: set to fb as this is now in C
		if(method=="fb") ll <- fb(init=object@init,A=object@trDens,B=object@dens,ntimes=object@ntimes,stationary=object@stationary)$logLike
		if(method=="lystig") ll <- lystig(init=object@init,A=object@trDens,B=object@dens,ntimes=object@ntimes,stationary=object@stationary)$logLike
		attr(ll, "df") <- freepars(object)
		attr(ll, "nobs") <- nobs(object)
		class(ll) <- "logLik"
		ll
	}
)

# depends on getpars and nobs
setMethod("logLik",signature(object="mix"),
	#function(object,method="lystig") { 
	function(object,method="fb") { 
		if(method=="fb") ll <- fb(init=object@init,A=matrix(0,1,1),B=object@dens,ntimes=object@ntimes,stationary=TRUE)$logLike
		if(method=="lystig") ll <- lystig(init=object@init,A=matrix(0,1,1),B=object@dens,ntimes=object@ntimes,stationary=TRUE)$logLike
		attr(ll, "df") <- freepars(object)
		attr(ll, "nobs") <- nobs(object)
		class(ll) <- "logLik"
		ll
	}
)
