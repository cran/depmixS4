setMethod(
	"multistart", 
	signature(object="mix"),
	function(object, nstart=10, initIters=10, ...) {
		llbest <- as.numeric(logLik(object))
		# below will fail if llbest = NA, so in that case replace with largest number
		if(is.na(llbest)) llbest <- .Machine$double.xmax
		bestmodel <- object
		nfailed <- 0
		for(i in 1:nstart) {
			fmod <- try(fit(object, emcontrol=em.control(maxit=initIters)), silent=TRUE)
			if(inherits(fmod,"try-error")) {
				nfailed <- nfailed + 1
			} else {
				if(logLik(fmod) > llbest) {
					llbest <- logLik(fmod)
					bestmodel <- fmod
				}
			}
		}
		bestmodel <- fit(bestmodel, emcontrol=em.control(random.start=FALSE))
		if(nfailed > 0) {
			warning(nfailed,
				"out of",
				nstart, 
				"attempts failed; result is based on ", 
				nstart-nfailed, 
				"starting values.\n"
				)
		}
		return(bestmodel)
	}
)