
# 
# response, transition and prior models for DEPMIX models
# 23-3-2008
# Maarten Speekenbrink & Ingmar Visser

# 
# RESPONSE CLASS
# 

setClass("response",
	representation(parameters="list",
		fixed="logical",
		y = "matrix",
		x = "matrix",
		npar = "numeric" # this is not really needed as it simply is length(unlist(parameters))
	)
)

#
# RESPONSE CLASS METHODS
# 

setMethod("npar","response",
	function(object) {
		return(object@npar)
	}
)

setMethod("getdf","response",
	function(object) {
		return(sum(!object@fixed))
	}
)

# 
# SPECIFIC INSTANCES: GLM 
# 

setClass("GLMresponse",
	representation(formula="formula",
		family="ANY"
	),
	prototype(
		formula=.~.,
		family=gaussian()
	),
	contains="response"
)

# 
# TODO: review print and summary methods
# 

setMethod("GLMresponse",
	signature(formula="formula"),
	function(formula,data=NULL,family=gaussian(),pstart=NULL,fixed=NULL,prob=TRUE, ...) {
		call <- match.call()
		mf <- match.call(expand.dots = FALSE)
		m <- match(c("formula", "data"), names(mf), 0)
		mf <- mf[c(1, m)]
		mf$drop.unused.levels <- TRUE
		mf[[1]] <- as.name("model.frame")
		mf <- eval(mf, parent.frame())
		x <- model.matrix(attr(mf, "terms"),mf)
		y <- model.response(mf)
		if(!is.matrix(y)) y <- matrix(y,ncol=1)
		parameters <- list()
		parameters$coefficients <- vector("numeric",length=ncol(x))
		if(family$family=="gaussian") {
			parameters$sd <- 1
		}
		if(family$family=="binomial") {
			# FIX ME
		  y <- model.response(mf)
			switch(is(y)[1],
        factor = {
          y <- as.matrix(as.numeric(as.numeric(y)==1))
          n <- matrix(1,nrow=nrow(y))
        },
        matrix = {
          if(ncol(y) == 2) {
            n <- rowSums(y)
            y <- as.matrix(y[,1])
          } else {
            stop("model response not valid for binomial model")
          }
        },
        numeric = {
          if(sum(y %in% c(0,1)) != length(y)) stop("model response not valid for binomial model")
          n <- matrix(1,nrow=length(y))
          y <- as.matrix(y)
        },
        stop("model response not valid for binomial model")
			 # assume 1 success, rest not
			 #y <- as.numeric(as.numeric(y)==1)
		  )
	  }
		if(family$family=="multinomial") {
			if(is.factor(y)) y <- model.matrix(~y-1)
			if(is.numeric(y)) y <- model.matrix(~factor(y)-1)
			parameters$coefficients <- matrix(0,ncol=ncol(y),nrow=ncol(x))
			if(is.null(fixed)) {
				fixed <- parameters$coefficients
				fixed[,family$base] <- 1 
				fixed <- c(as.logical(t(fixed)))
			}
		}
		npar <- length(unlist(parameters))
		if(is.null(fixed)) fixed <- as.logical(rep(0,npar))
		if(!is.null(pstart)) {
			if(length(pstart)!=npar) stop("length of 'pstart' must be",npar)
			if(family$family=="multinomial") {
				if(family$link=="identity") parameters$coefficients[1,] <- family$linkfun(pstart[1:ncol(parameters$coefficients)])
				else {
					if(prob) parameters$coefficients[1,] <- family$linkfun(pstart[1:ncol(parameters$coefficients)],base=family$base)
					else parameters$coefficients[1,] <- pstart[1:ncol(parameters$coefficients)]
				}
				pstart <- matrix(pstart,ncol(x),byrow=TRUE)
				if(ncol(x)>1) parameters$coefficients[2:ncol(x),] <- pstart[2:ncol(x),]
			} else {
				parameters$coefficients <- family$linkfun(as.numeric(pstart[1:length(parameters$coefficients)]))
			}
			if(length(unlist(parameters))>length(parameters$coefficients)) {
				if(family$family=="gaussian") parameters$sd <- as.numeric(pstart[(length(parameters$coefficients)+1)])
			}
		}
		mod <- switch(family$family,
			gaussian = new("NORMresponse",formula=formula,family=family,parameters=parameters,fixed=fixed,x=x,y=y,npar=npar),
			binomial = new("BINOMresponse",formula=formula,family=family,parameters=parameters,fixed=fixed,x=x,y=y,npar=npar,n=n),
			multinomial = new("MULTINOMresponse",formula=formula,family=family,parameters=parameters,fixed=fixed,x=x,y=y,npar=npar),
			new("GLMresponse",formula=formula,family=family,parameters=parameters,fixed=fixed,x=x,y=y,npar=npar)
		)
		mod
	}
)

# 
# SPECIFIC INSTANCES: MULTINOMIAL 
# 
# for the transition models and the prior (y is missing, ie there is no
# response, and nstates must be provided as the number of categories
# neccessary in the mulinomial model)
# 

setClass("transInit",contains="GLMresponse")

# FIX ME: data is a necessary argument to determine the dimension of x, even when there
# are no covariates (and there are by definition no responses ...)
setMethod("transInit",
	signature(formula="formula"),
	function(formula,nstates,data=NULL,family=multinomial(),pstart=NULL,fixed=NULL,prob=TRUE, ...) {
		call <- match.call()
		mf <- match.call(expand.dots = FALSE)
		m <- match(c("formula", "data"), names(mf), 0)
		mf <- mf[c(1, m)]
		mf$drop.unused.levels <- TRUE
		mf[[1]] <- as.name("model.frame")
		mf <- eval(mf, parent.frame())
		x <- model.matrix(attr(mf, "terms"),mf)
		y <- matrix(1,ncol=1) # y is not needed in the transition and init models		    
		parameters <- list()
		if(is.null(nstates)) stop("'nstates' must be provided in call to trinModel")
		if(family$family=="multinomial") {
			parameters$coefficients <- matrix(0,ncol=nstates,nrow=ncol(x))
			if(is.null(fixed)) {
				fixed <- parameters$coefficients
				fixed[,family$base] <- 1 
				fixed <- c(as.logical(t(fixed)))
			}
		}
		npar <- length(unlist(parameters))
		if(is.null(fixed)) fixed <- rep(0,npar)
		if(!is.null(pstart)) {
			if(length(pstart)!=npar) stop("length of 'pstart' must be ",npar)
			if(family$family=="multinomial") {
				if(prob) {
					if(family$link=="identity") {
						parameters$coefficients[1,] <- family$linkfun(pstart[1:ncol(parameters$coefficients)])
					} else {
						parameters$coefficients[1,] <- family$linkfun(pstart[1:ncol(parameters$coefficients)],base=family$base)
					}
				} else {
					parameters$coefficients[1,] <- pstart[1:ncol(parameters$coefficients)]
				}
				pstart <- matrix(pstart,,ncol(x),byrow=TRUE)
				if(ncol(x)>1) parameters$coefficients[2:ncol(x),] <- pstart[2:ncol(x),]
			} else {
				if(family$link=="identity") parameters$coefficients <- family$linkfun(pstart[1:length(parameters$coefficients)])
				else parameters$coefficients <- family$linkfun(pstart[1:length(parameters$coefficients)],base=family$base)
			}
		}
		mod <- switch(family$family,
			multinomial = new("transInit",formula=formula,family=family,parameters=parameters,fixed=fixed,x=x,y=y,npar=npar),
			new("transInit",formula=formula,family=family,parameters=parameters,fixed=fixed,x=x,y=y,npar=npar)
		)
		mod
	}
)

setMethod("show","GLMresponse",
	function(object) {
		cat("Model of type ", object@family$family, ", formula: ",sep="")
		print(object@formula)
		cat("Coefficients: \n")
		print(object@parameters$coefficients)
		if(object@family$family=="gaussian") {
			cat("sd ",object@parameters$sd,"\n")
		}	
	}
)

setMethod("setpars","GLMresponse",
	function(object,values,which="pars",...) {
		npar <- npar(object)
		if(length(values)!=npar) stop("length of 'values' must be",npar)
		# determine whether parameters or fixed constraints are being set
		switch(which,
			"pars"= {
				if(object@family$family=="multinomial") {
					object@parameters$coefficients[1,] <- values[1:ncol(object@parameters$coefficients)]
					values <- matrix(values,,ncol(object@x),byrow=TRUE)
					if(ncol(object@x)>1) object@parameters$coefficients[2:ncol(object@x),] <- values[2:ncol(object@x),]
				} else {
					object@parameters$coefficients <- values[1:length(object@parameters$coefficients)]
				}
				if(length(unlist(object@parameters))>length(object@parameters$coefficients)) {
					if(object@family$family=="gaussian") object@parameters$sd <- values[(length(object@parameters$coefficients)+1)]
				}
			},
			"fixed" = {
				object@fixed <- as.logical(values)
			}
		)
		return(object)
	}
)

setMethod("getpars","GLMresponse",
	function(object,which="pars",...) {
		switch(which,
			"pars" = {
				parameters <- numeric()
				if(object@family$family=="multinomial") {
					# coefficient is usually a matrix here 		
					parameters <- c(t(object@parameters$coefficients)) # Why transpose?
				} else {
					parameters <- unlist(object@parameters)
				}
				pars <- parameters
			},
			"fixed" = {
				pars <- object@fixed
			}
		)
		return(pars)
	}
)

# Class 'rGLM' (extends 'rModel')
# use: when method 'fit' can use glm.fit

setClass("BINOMresponse",representation(n="matrix"),contains="GLMresponse")
setClass("MULTINOMresponse",contains="GLMresponse")
setClass("NORMresponse",contains="GLMresponse")

# Class 'MVNresponse' (multivariate normal response model)
setClass("MVNresponse",contains="response")

# method 'fit'
# use: in EM (M step)
# returns: (fitted) response with (new) estimates of parameters

setMethod("fit","GLMresponse",
	function(object,w) {
		pars <- object@parameters
		fit <- glm.fit(x=object@x,y=object@y,weights=w,family=object@family)
		pars$coefficients <- fit$coefficients
		object <- setpars(object,unlist(pars))
		object
	}
)

setMethod("fit","NORMresponse",
	function(object,w) {
		pars <- object@parameters
		fit <- lm.wfit(x=object@x,y=object@y,w=w)
		pars$coefficients <- fit$coefficients
		pars$sd <- sqrt(sum(w*fit$residuals^2/sum(w)))
		object <- setpars(object,unlist(pars))
		object
	}
)

setMethod("fit","MVNresponse",
	function(object,w) {
		pars <- object@parameters
		fit <- lm.wfit(x=object@x,y=object@y,w=w)
		object@parameters$coefficients <- fit$coefficients
		object@parameters$Sigma <- cov.wt(x=fit$residuals,wt=w)["cov"]
		object <- setpars(object,unlist(pars))
		object
	}
)

setMethod("fit","MULTINOMresponse",
	function(object,w) {
# 		require(nnet)
		pars <- object@parameters
		base <- object@family$base # delete me
		y <- object@y[,-base]
		x <- object@x
		mask <- matrix(1,nrow=nrow(pars$coefficients),ncol=ncol(pars$coefficients))
		mask[,base] <- 0
		fit <- nnet.default(x,y,weights=w,size=0,entropy=TRUE,skip=TRUE,mask=mask,rang=0,trace=FALSE)
		pars$coefficients <- matrix(fit$wts,ncol=ncol(pars$coefficients),nrow=nrow(pars$coefficients),byrow=TRUE)
		object <- setpars(object,unlist(pars))
		object
	}
)

# method 'logDens'
# use: instead of density slot in rModel
# returns: matrix with log(p(y|x,parameters))
setMethod("logDens","BINOMresponse",
	function(object) {
		dbinom(x=object@y,size=object@n,prob=predict(object),log=TRUE)
	}
)

setMethod("dens","BINOMresponse",
	function(object,log=FALSE) {
		dbinom(x=object@y,size=object@n,prob=predict(object),log=log)
	}
)

setMethod("logDens","MULTINOMresponse",
	function(object) {
		log(rowSums(object@y*predict(object)))
	}
)

setMethod("dens","MULTINOMresponse",
	function(object,log=FALSE) {
		if(log) log(rowSums(object@y*predict(object)))
		else rowSums(object@y*predict(object))
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

setMethod("logDens","transInit",
	function(object) {
		log(predict(object))
	}
)

setMethod("dens","transInit",
	function(object,log=FALSE) {
		if(log) log(predict(object))
		else predict(object)
	}
)

setMethod("logLik","GLMresponse",
	function(object) {
		sum(logDens(object))
	}
)

# method "predict"

setMethod("predict","GLMresponse",
	function(object) {
		object@family$linkinv(object@x%*%object@parameters$coefficients)
	}
)

setMethod("predict","MULTINOMresponse",
	function(object) {
		if(object@family$link=="identity") object@x%*%object@parameters$coefficients
		else {
			object@family$linkinv(object@x%*%object@parameters$coefficients,base=object@family$base)
		}
	}
)

setMethod("predict","NORMresponse",
	function(object) {
		object@x%*%object@parameters$coefficients
	}
)

setMethod("predict","transInit",
	function(object) {
		object@family$linkinv(object@x%*%object@parameters$coefficients,base=object@family$base)
	}
)

# 
# MULTINOMIAL LINK FUNCTION
# 

mlogit <- function(base=1) {
	# 	matrix formulation is possibly very inefficient?!?!?
	# 	moreover it does not admit of bases being different from 1??!?!?
	linkfun <- function(p,base) {
		lfun <- function(p,base) {
			p <- p/sum(p)
			beta <- numeric(length(p))
			if(any(p==1)) beta[which(p==1)]=Inf
			else beta[-base] <- log(p[-base]/p[base])
			return(beta)
		}
		if(is.matrix(p)) {
			beta <- t(apply(p,1,lfun,base=base))
		} else {
			beta <- lfun(p,base)
		}
		return(beta)
	}
	linkinv <- function(eta,base) {
		linv <- function(eta,base) {
			pp <- numeric(length(eta))
			if(any(is.infinite(eta))) {
				pp[which(is.infinite(eta))] <- 1
			} else {
				expb <- exp(eta)
				sumb <- sum(expb)
				pp[base] <- 1/sumb
				pp[-base] <- expb[-base]/sumb
			}
			return(pp)
		}
		if(is.matrix(eta)) {
			if(ncol(eta)==1) {
				pp <- as.matrix(apply(eta,1,linv,base=base)) # fixes problem with column matrix eta
			} else pp <- t(apply(eta,1,linv,base=base)) 	
		} else {
			pp <- linv(eta,base)
		}
		return(pp)
	}
	mu.eta <- function(eta) {
		if(length(eta)==1) return(eta-eta^2)
		if(is.vector(eta)) return(diag(eta)-outer(eta,eta))
	}
	valideta <- function(eta) {
		TRUE # fix me
	}
	
	name <- "mlogit"
	structure(list(linkfun=linkfun,
			linkinv=linkinv,
			mu.eta=mu.eta,
			valideta=valideta,
			name=name,
			base=base),
		class="link-glm")
}

# 
# MAKE FAMILY OBJECT
# 

multinomial <- function(link="mlogit",base=1) {
	# adapted from gaussian()
	linktemp <- substitute(link)
	if (!is.character(linktemp)) {
		linktemp <- deparse(linktemp)
		if (linktemp == "link") {
			warning("use of multinomial(link=link) is deprecated\n",
				domain = NA)
			linktemp <- eval(link)
			if (!is.character(linktemp) || length(linktemp) !=1)
			stop("'link' is invalid", domain = NA)
		}
	}
	okLinks <- c("mlogit")
	if (linktemp %in% okLinks) {
		if(linktemp == "mlogit") stats <- mlogit() else stats <- make.link(linktemp)
	} else {
		if (is.character(link)) {
			stats <- make.link(link)
			linktemp <- link
		} else {
			if (inherits(link, "link-glm")) {
				stats <- link
				if (!is.null(stats$name))
				linktemp <- stats$name
			} else {
				stop(gettextf("link \"%s\" not available for multinomial family; available links are %s",
						linktemp, paste(sQuote(okLinks), collapse = ", ")),
					domain = NA)
			}
		}
	}
	variance <- function(mu) {
		n <- length(mu)
		v <- diag(n)*outer(mu,1-mu) - (1-diag(n))*outer(mu,-mu)
	}
	validmu <- function(mu) {
		all(mu > 0) && all(mu < 1)
	}
	dev.resids <- function(y,mu,wt) {
		
	}
	initialize <- expression()
	structure(list(family = "multinomial", link = linktemp, linkfun = stats$linkfun,
			linkinv = stats$linkinv, variance = variance, dev.resids = dev.resids,
			mu.eta = stats$mu.eta, initialize = initialize, validmu = validmu, valideta = stats$valideta, base=base),
			class = "family")
}

setMethod("fit","transInit",
	function(object,w,ntimes) {
		pars <- object@parameters
		if(missing(w)) w <- NULL
		oldfit <- function() {
			tol <- 1e-5 # TODO: check global options
			pars <- object@parameters
			b <- pars$coefficients
			base <- object@family$base
			if(is.matrix(w)) nan <- which(is.na(rowSums(w))) else nan <- which(is.na(w))
			#vgam(cbind(w[,-base],w[,base]) ~ ) # what is this?
			y <- as.vector(t(object@family$linkinv(w[-c(nan,ntimes),-base],base=object@family$base)))
			x <- object@x[-c(nan,ntimes),]			
			if(!is.matrix(x)) x <- matrix(x,ncol=ncol(object@x))
			nt <- nrow(x)			
			Z <- matrix(ncol=length(b))
			Z <- vector()
			for(i in 1:nt) Z <- rbind(Z,t(bdiag(rep(list(x[i,]),ncol(w)-1))))			
			mu <- object@family$linkinv(x%*%b,base=base)			
			mt <- as.numeric(t(mu[,-base]))
			Dl <- Sigmal <- Wl <- list()			
			converge <- FALSE
			while(!converge) {
				b.old <- b
				for(i in 1:nt) {
					Dl[[i]] <- object@family$mu.eta(mu[i,-base])
					Sigmal[[i]] <- object@family$variance(mu[i,-base])
					Wl[[i]] <- Dl[[i]]%*%solve(Sigmal[[i]])%*%t(Dl[[i]]) # TODO: 
				}
				Sigma <- bdiag(Sigmal)
				D <- bdiag(Dl)
				W <- bdiag(Wl)
				
				b[,-base] <- as.numeric(b[,-base]) + solve(t(Z)%*%W%*%Z)%*%(t(Z)%*%D%*%solve(Sigma)%*%(y-mt))
				if(abs(sum(b-b.old)) < tol) converge <- TRUE
				mu <- object@family$linkinv(x%*%b,base=base)
				mt <- as.numeric(t(mu[,-base]))
			}
			pars$coefficients <- t(b) # TODO: setpars gets matrix in wrong order!!! Fix this in setpars.
			pars
		}
		
		vglmfit <- function() {		
			base <- object@family$base
			w <- cbind(w[,-base],w[,base])
			x <- slot(object,"x")
			fam <- slot(object,"family")
			fit <- vglm(w~x,fam)
			pars$coefficients[,-base] <- t(slot(fit,coefficients))  # TODO: setpars gets matrix in wrong order!!! Fix this in setpars.
			pars
		}
		
		nnetfit <- function() {
# 			require(nnet)
			pars <- object@parameters
			base <- object@family$base # delete me
			#y <- object@y[,-base]
			y <- object@y
			x <- object@x
			if(is.matrix(y)) na <- unlist(apply(y,2,function(x) which(is.na(x)))) else na <- which(is.na(y))
			if(is.matrix(x)) na <- c(na,unlist(apply(x,2,function(x) which(is.na(x))))) else na <- c(na,which(is.na(x)))
			if(!is.null(w)) na <- c(na,which(is.na(w)))
			y <- as.matrix(y)
			x <- as.matrix(x)
			na <- unique(na)
			x <- x[-na,]
			y <- y[-na,]
			y <- round(y) # delete me
			if(!is.null(w)) w <- w[-na]
			#mask <- matrix(1,nrow=nrow(pars$coefficients),ncol=ncol(pars$coefficients))
			#mask[,base] <- 0
			if(!is.null(w)) fit <- multinom(y~x-1,weights=w,trace=FALSE) else fit <- multinom(y~x-1,weights=w,trace=FALSE)
			ids <- vector(,length=ncol(y))
			ids[base] <- 1
			ids[-base] <- 2:ncol(y)
			pars$coefficients <- t(matrix(fit$wts,ncol=ncol(y))[-1,ids])
			object <- setpars(object,unlist(pars))
			#object
			pars
		}
		
# 		require(nnet)
		pars <- object@parameters
		base <- object@family$base # delete me
		#y <- object@y[,-base]
		y <- object@y
		x <- object@x
		if(is.matrix(y)) na <- unlist(apply(y,2,function(x) which(is.na(x)))) else na <- which(is.na(y))
		if(is.matrix(x)) na <- c(na,unlist(apply(x,2,function(x) which(is.na(x))))) else na <- c(na,which(is.na(x)))
		if(!is.null(w)) na <- c(na,which(is.na(w)))
		y <- as.matrix(y)
		x <- as.matrix(x)
		na <- unique(na)
		if(length(na)>0) {
			x <- x[-na,]
			y <- y[-na,]
			#y <- round(y) # delete me
			if(!is.null(w)) w <- w[-na]
		}
		#mask <- matrix(1,nrow=nrow(pars$coefficients),ncol=ncol(pars$coefficients))
		#mask[,base] <- 0
		if(!is.null(w)) fit <- multinom(y~x-1,weights=w,trace=FALSE) else fit <- multinom(y~x-1,trace=FALSE)
		ids <- vector(,length=ncol(y))
		ids[base] <- 1
		ids[-base] <- 2:ncol(y)
		pars$coefficients <- t(matrix(fit$wts,ncol=ncol(y))[-1,ids]) # why do we need to transpose?
		object <- setpars(object,unlist(pars))
		object
	}
)

