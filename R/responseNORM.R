
setClass("NORMresponse",contains="GLMresponse")

# method 'fit'
# use: in EM (M step)
# returns: (fitted) response with (new) estimates of parameters

setMethod("fit","NORMresponse",
	function(object,w) {
		if(missing(w)) w <- NULL
		pars <- object@parameters
		if(!is.null(w)) {
			fit <- lm.wfit(x=object@x,y=object@y,w=w)
		} else {
			fit <- lm.fit(x=object@x,y=object@y)
		}
		pars$coefficients <- fit$coefficients
		if(!is.null(w)) {
			pars$sd <- sqrt(sum(w*fit$residuals^2/sum(w)))
		} else {
			pars$sd <- sd(fit$residuals)
		}
		object <- setpars(object,unlist(pars))
		object
	}
)

setMethod("logDens","NORMresponse",
	function(object) {
		dnorm(x=object@y,mean=predict(object),sd=object@parameters$sd,log=TRUE)
	}
)

setMethod("dens","NORMresponse",
	function(object,log=FALSE) {
		dnorm(x=object@y,mean=predict(object),sd=object@parameters$sd,log=log)
	}
)

setMethod("predict","NORMresponse",
	function(object) {
		object@x%*%object@parameters$coefficients
	}
)

setMethod("simulate",signature(object="NORMresponse"),
  function(object,nsim=1,seed=NULL,times) {
    if(!is.null(seed)) set.seed(seed)
    if(missing(times)) {
      # draw in one go
      mu <- predict(object)
    } else {
      mu <- predict(object)[times]
    }  
    nt <- length(mu)
    sd <- object@parameters$sd
    response <- rnorm(nt*nsim,mean=mu,sd=sd)
    #if(nsim > 1) response <- matrix(response,ncol=nsim)
    response <- as.matrix(response)
    return(response)
  }
)
