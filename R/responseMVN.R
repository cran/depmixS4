# Class 'MVNresponse' (multivariate normal response model)

cov2par <- function(x) {
	if(NROW(x) != NCOL(x)) stop("cov2par requires a square matrix") 
	x[lower.tri(x,diag=TRUE)]
}

par2cov <- function(x) {
	npar <- length(x)
	dim <- (sqrt(8*npar + 1) - 1)/2
	if(abs(dim - round(dim)) >= .Machine$double.eps^0.5) stop("number of parameters not suitable for par2cov")
	cov <- matrix(0.0,ncol=dim,nrow=dim)
	cov[lower.tri(cov,diag=TRUE)] <- x
	cov[upper.tri(cov)] <- t(cov)[upper.tri(cov)]
	cov
}

setClass("MVNresponse",
  representation(formula="formula"),
  contains="response"
)

setMethod("fit","MVNresponse",
	function(object,w) {
    if(missing(w)) w <- NULL
		pars <- object@parameters
		nas <- any(is.na(object@y))
		# if(NCOL(object@x) == 1 && all(object@x == 1)) {
		#   # intercept only model, use standard MVN estimation for EM
		#   if(!is.null(w)) {
		#     if(nas) {
		#       y <- object@y[rowSums(is.na(object@y))==0,,drop=FALSE]
		#       w <- w[rowSums(is.na(object@y))==0]
		#       mean <- apply(y,2,weighted.mean,w=w)
		#       cov <- cov.wt(y,wt=w,method = "ML")$cov
		#     } else {
		#       mean <- apply(object@y,2,weighted.mean,w=w)
		#       cov <- cov.wt(object@y,wt=w,method = "ML")$cov
		#     }
		#   } else {
		#     if(nas) {
		#       y <- object@y[rowSums(is.na(object@y))==0,,drop=FALSE]
		#       mean <- colMeans(y)
		#       cov <- cov.wt(y,method="ML")$cov
		#     } else {
		#       mean <- colMeans(object@y)
		#       cov <- cov.wt(object@y,method="ML")$cov
		#     }
		#   }
		#   object@parameters$coefficients <- mean
		#   object@parameters$Sigma <- cov2par(cov)
		#   return(object)
		# }
		if(!is.null(w)) {
		  if(nas) {
		    x <- object@x[rowSums(is.na(object@y))==0,,drop=FALSE]
		    y <- object@y[rowSums(is.na(object@y))==0,,drop=FALSE]
		    w <- w[rowSums(is.na(object@y))==0]
		    fit <- lm.wfit(x=x,y=y,w=w)
		  } else {
		    fit <- lm.wfit(x=object@x,y=object@y,w=w) #else fit <- lm.fit(x=object@x,y=object@y)
		  }
		} else {
		  if(nas) {
		    x <- object@x[rowSums(is.na(object@y))==0,,drop=FALSE]
		    y <- object@y[rowSums(is.na(object@y))==0,,drop=FALSE]
		    fit <- lm.fit(x=x,y=y)
		  } else {
		    fit <-  lm.fit(x=object@x,y=object@y)
		  }
		}
		object@parameters$coefficients <- fit$coefficients
		if(!is.null(w)) {
		  object@parameters$Sigma <- cov2par(cov.wt(x=fit$residuals,wt=w,method="ML")$cov)
		} else {
		  object@parameters$Sigma <- cov2par(cov.wt(x=fit$residuals,method="ML")$cov)
		}
		return(object)
	}
)


dm_dmvnorm <- function(x,mean,sigma,log=FALSE,logdet,invSigma) {
  # update 29/08/2018: using new function definition from mvtnorm 1.0-8
  # taken from mvtnorm package
  # using Cholesky decomposition
  # allows passing of logdet (sigma) and invsigma to save 
  # computation when called repeated times with same sigma
  #
  # mean can be a vector of length ncol(x) or a matrix of dim(x)  
  if(is.vector(x)) {
    x <- matrix(x, ncol = length(x))
  }
  p <- ncol(x)
  if(!missing(invSigma)) {
    warning("Use of logdet and invsigma are depreciated")
    # do the old style computation
    # as far as I'm aware, depmixS4 doesn't pass invSigma or logdet to this function
    # so this shouldn't really be necessary
    if(missing(mean)) {
      mean <- matrix(0, ncol = p)
    }
    if(is.vector(mean)) {
      mean <- matrix(mean, ncol = p)
    }
    # check consistency
    if (NCOL(x) != NCOL(invSigma)) {
      stop("x and sigma have non-conforming size")
    }
    if (NROW(invSigma) != NCOL(invSigma)) {
      stop("sigma must be a square matrix")
    }
    if (NCOL(invSigma) != NCOL(mean)) {
      stop("mean and sigma have non-conforming size")
    }
    if(missing(logdet)) {
      ev <- eigen(sigma, symmetric = TRUE, only.values = TRUE)$values
      if(!all(ev >= 0)) return(rep(NaN,nrow(x))) else logdet <- sum(log(ev))
    }
    if(NROW(mean) == NROW(x)) {
      # varying means
      
      # from "mahalanobis":    
      x <- x - mean
      distval <- rowSums((x %*% invSigma) * x)
      #names(retval) <- rownames(x)
      #retval
    } else {
      # constant mean
      if (length(mean) != NROW(invSigma)) {
        stop("mean and sigma have non-conforming size")
      }
      distval <- mahalanobis(x, center = mean, cov = invSigma, inverted=TRUE)
    }	
    logretval <- -(ncol(x) * log(2 * pi) + logdet + distval)/2
  } else {
    ## use Cholesky decomposition
    if(missing(mean)) {
      mean <- rep(0,p)
    }
    if(missing(sigma)) {
      sigma <- diag(p)
    }
    if(!missing(mean)) {
      ##if(!is.null(dim(mean))) 
        ##dim(mean) <- NULL
      if(NCOL(mean) != NCOL(sigma))
        stop("mean and sigma have non-conforming size")
    }
    if (!missing(sigma)) {
      if(p != ncol(sigma)) 
        stop("x and sigma have non-conforming size")
      if(!isSymmetric(sigma, tol = sqrt(.Machine$double.eps), 
                       check.attributes = FALSE)) 
        stop("sigma must be a symmetric matrix")
    }
    dec <- tryCatch(chol(sigma), error = function(e) e)
    if (inherits(dec, "error")) {
      if(all(dim(x) == dim(mean))) {
        x.is.mu <- rowSums(x != mean) == 0
      } else {
        if(length(mean) == p) {
          x.is.mu <- colSums(t(x) != mean) == 0
        }
      }
      #x.is.mu <- colSums(t(x) != mean) == 0
      logretval <- rep.int(-Inf, nrow(x))
      logretval[x.is.mu] <- Inf
    }
    else {
      if(all(dim(x) == dim(mean))) {
        tmp <- backsolve(dec, t(x - mean), transpose = TRUE)
      } else {
        tmp <- backsolve(dec, t(x) - mean, transpose = TRUE)
      }
      rss <- colSums(tmp^2)
      logretval <- -sum(log(diag(dec))) - 0.5 * p * log(2 * pi) - 0.5 * rss
    }
  }
  if(log) {
    return(logretval)
  } else {
    return(exp(logretval))
  }
}

setMethod("logDens","MVNresponse",
	function(object,...) {
		dm_dmvnorm(x=object@y,mean=predict(object),sigma=par2cov(object@parameters$Sigma),log=TRUE,...)
	}
)

setMethod("dens","MVNresponse",
	function(object,log=FALSE,...) {
		dm_dmvnorm(x=object@y,mean=predict(object),sigma=par2cov(object@parameters$Sigma),log=log,...)
	}
)


setMethod("predict","MVNresponse",
	function(object) {
		object@x%*%object@parameters$coefficients
	}
)

setMethod("simulate",signature(object="MVNresponse"),
	function(object,nsim=1,seed=NULL,times) {
		if(!is.null(seed)) set.seed(seed)
		if(missing(times)) {
			# draw in one go
			mu <- predict(object)
		} else {
			mu <- predict(object)[times,]
		}
		nt <- nrow(mu)
		if(nrow(object@parameters$coefficients==1)) response <- mvrnorm(nt*nsim,mu=mu[1,],Sigma=par2cov(object@parameters$Sigma))
		else {
			response <- matrix(0,nrow(mu),ncol(mu))
			for(i in 1:nrow(mu)) {
				response[i,] <- response <- mvrnorm(1,mu=mu[i,],Sigma=par2cov(object@parameters$Sigma))
			}
		}
		return(response)
	}
)

setMethod("MVNresponse",
	signature(formula="formula"),
	function(formula,data,pstart=NULL,fixed=NULL,na.action="na.pass",...) {
		call <- match.call()
		mf <- match.call(expand.dots = FALSE)
		m <- match(c("formula", "data"), names(mf), 0)
		mf <- mf[c(1, m)]
		mf$drop.unused.levels <- TRUE
		mf$na.action <- na.action
		mf[[1]] <- as.name("model.frame")
		mf <- eval(mf, parent.frame())
		x <- model.matrix(attr(mf, "terms"),mf)
		y <- model.response(mf)
		if(!is.matrix(y)) y <- matrix(y,ncol=1)
		
		nas <- any(is.na(y))
		if(nas) {
		  y_missing <- rowSums(is.na(y))
		  x_missing <- rowSums(is.na(x))
		  if(sum(x_missing[y_missing > 0]) > 0) {
		    stop("'depmixS4' does not currently handle covariates with missing data at places where the response is not missing.")
		  }
		}
		#if(any(is.na(x))) 
		parameters <- list()
		parameters$coefficients <- matrix(0.0,ncol=ncol(y),nrow=ncol(x))
		parameters$Sigma <- cov2par(diag(ncol(y)))
		npar <- length(unlist(parameters))
		parlow.coeff=rep(-Inf,length(unlist(parameters$coefficients)))
		parup.coeff=rep(Inf,length(unlist(parameters$coefficients)))
		parup.cov <- rep(Inf,length(unlist(parameters$Sigma)))
		mcov <- matrix(-Inf,ncol(y),ncol(y))
		diag(mcov) <- .Machine$double.eps
		parlow.cov <- cov2par(mcov)
		constr <- list(	
			parup = c(parup.coeff,parup.cov),
			parlow = c(parlow.coeff,parlow.cov)
		)
		if(is.null(fixed)) fixed <- as.logical(rep(0,npar))
		if(!is.null(pstart)) {
			if(length(pstart)!=npar) stop("length of 'pstart' must be",npar)
			parameters$coefficients <- matrix(pstart[1:(ncol(x)*ncol(y))],ncol(x),byrow=T)
			parameters$Sigma <- as.numeric(pstart[(length(parameters$coefficients)+1):length(pstart)])			
		}
		mod <- new("MVNresponse",formula=formula,parameters=parameters,fixed=fixed,x=x,y=y,npar=npar,constr=constr)
		mod
	}
)

setMethod("show","MVNresponse",
	function(object) {
		cat("Multivariate Normal Model, formula: ",sep="")
		print(object@formula)
		cat("Coefficients: \n")
		print(object@parameters$coefficients)
		cat("Sigma: \n")
		print(par2cov(object@parameters$Sigma))
	}
)

setMethod("setpars","MVNresponse",
	function(object, values, which="pars", prob=FALSE, ...) {
		npar <- npar(object)
		if(length(values)!=npar) stop("length of 'values' must be",npar)
		# determine whether parameters or fixed constraints are being set
		nms <- names(object@parameters$coefficients)
		if(length(values) == 0) return(object) # nothing to set;
		switch(which,
			"pars" = {
				object@parameters$coefficients <- matrix(values[1:length(object@parameters$coefficients)],ncol(object@x))
				st <- length(object@parameters$coefficients)+1
				object@parameters$Sigma <- as.numeric(values[st:(st+length(object@parameters$Sigma)-1)])
			},
			"fixed" = {
				object@fixed <- as.logical(values)
			}
		)
		names(object@parameters$coefficients) <- nms
		return(object)
	}
)

setMethod("getpars","MVNresponse",
	function(object,which="pars",...) {
		switch(which,
			"pars" = {
				parameters <- numeric()
				parameters <- unlist(object@parameters)
				pars <- parameters
			},
			"fixed" = {
				pars <- object@fixed
			}
		)
		return(pars)
	}
)
