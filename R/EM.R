# 
# Maarten Speekenbrink 23-3-2008
# 

em <- function(object,...) {
	if(!is(object,"mix")) stop("object is not of class '(dep)mix'")
	call <- match.call()
	if(is(object,"depmix")) {
		call[[1]] <- as.name("em.depmix")
	} else {
		call[[1]] <- as.name("em.mix")
	}
	object <- eval(call, parent.frame())
	object
}

# em for lca and mixture models
em.mix <- function(object,maxit=100,tol=1e-8,crit=c("relative","absolute"),random.start=FALSE,verbose=FALSE,...) {
	
	if(!is(object,"mix")) stop("object is not of class 'mix'")
	
	crit <- match.arg(crit)
	
	ns <- nstates(object)
	ntimes <- ntimes(object)
	lt <- length(ntimes)
	et <- cumsum(ntimes)
	bt <- c(1,et[-lt]+1)
	
	converge <- FALSE
	j <- 0
	
	# compute responsibilities
	B <- apply(object@dens,c(1,3),prod)
	gamma <- object@init*B
	LL <- sum(log(rowSums(gamma)))
	# normalize
	gamma <- gamma/rowSums(gamma)
	
	if(random.start) {
		nr <- sum(ntimes(object))
		gamma <- matrix(runif(nr*ns,min=.0001,max=.9999),nr=nr,nc=ns)
		gamma <- gamma/rowSums(gamma)
	} 
	
	LL.old <- LL + 1
	
	while(j <= maxit & !converge) {
		
		# maximization
		
		# should become object@prior <- fit(object@prior)
		object@prior@y <- gamma[bt,,drop=FALSE]
		object@prior <- fit(object@prior, w=NULL,ntimes=NULL)
		object@init <- dens(object@prior)
		
		for(i in 1:ns) {
			for(k in 1:nresp(object)) {
				object@response[[i]][[k]] <- fit(object@response[[i]][[k]],w=gamma[,i])
				# update dens slot of the model
				object@dens[,k,i] <- dens(object@response[[i]][[k]])
			}
		}
		
		# expectation
		B <- apply(object@dens,c(1,3),prod)
		gamma <- object@init*B
		LL <- sum(log(rowSums(gamma)))
		# normalize
		gamma <- gamma/rowSums(gamma)
		
		# print stuff
		if(verbose&((j%%5)==0)) {
			cat("iteration",j,"logLik:",LL,"\n")
		}
		
		if(LL >= LL.old) {
		  if((crit == "absolute" &&  LL - LL.old < tol) || (crit == "relative" && (LL.old - LL)/LL.old  < tol)) {
			  cat("iteration",j,"logLik:",LL,"\n")
			  converge <- TRUE
			}
		} else {
		  # this should not really happen...
		  if(j > 0) warning("likelihood decreased on iteration",j)
		}

		LL.old <- LL
		j <- j+1

	}

	class(object) <- "mix.fitted"

	if(converge) {
		object@message <- switch(crit,
			relative = "Log likelihood converged to within tol. (relative change)",
			absolute = "Log likelihood converged to within tol. (absolute change)"
		)
	} else object@message <- "'maxit' iterations reached in EM without convergence."

	# no constraints in EM, except for the standard constraints ...
	# which are produced by the following (only necessary for getting df right in logLik and such)
	constraints <- getConstraints(object)
	object@conMat <- constraints$lincon
	object@lin.lower <- constraints$lin.l
	object@lin.upper <- constraints$lin.u
	
	object
	
}

# em for hidden markov models
em.depmix <- function(object,maxit=100,tol=1e-8,crit=c("relative","absolute"),random.start=FALSE,verbose=FALSE,...) {
	
	if(!is(object,"depmix")) stop("object is not of class '(dep)mix'")
	crit <- match.arg(crit)
	
	ns <- nstates(object)
	
	ntimes <- ntimes(object)
	lt <- length(ntimes)
	et <- cumsum(ntimes)
	bt <- c(1,et[-lt]+1)
	
	converge <- FALSE
	j <- 0
	
	# A <- object@trDens
	# B <- object@dens
	# init <- object@init
	
	# initial expectation
	fbo <- fb(init=object@init,A=object@trDens,B=object@dens,ntimes=ntimes(object),stationary=object@stationary)
	LL <- fbo$logLike
	LL.old <- LL + 1
	
	if(random.start) {
		nr <- sum(ntimes(object))
		fbo$gamma <- matrix(runif(nr*ns,min=.0001,max=.9999),nr=nr,nc=ns)
		fbo$gamma <- fbo$gamma/rowSums(fbo$gamma)
	}
	
	while(j <= maxit & !converge) {
		
		# maximization
				
		# should become object@prior <- fit(object@prior)
		object@prior@y <- fbo$gamma[bt,,drop=FALSE]
		object@prior <- fit(object@prior, w=NULL,ntimes=NULL)
		object@init <- dens(object@prior)
				
		trm <- matrix(0,ns,ns)
		for(i in 1:ns) {
			if(!object@stationary) {
				object@transition[[i]]@y <- fbo$xi[,,i]/fbo$gamma[,i]
				object@transition[[i]] <- fit(object@transition[[i]],w=as.matrix(fbo$gamma[,i]),ntimes=ntimes(object)) # check this
			} else {
				for(k in 1:ns) {
					trm[i,k] <- sum(fbo$xi[-c(et),k,i])/sum(fbo$gamma[-c(et),i])
				}
				# FIX THIS; it will only work with specific trinModels??
				object@transition[[i]]@parameters$coefficients <- switch(object@transition[[i]]@family$link,
					identity = object@transition[[i]]@family$linkfun(trm[i,]),
					mlogit = object@transition[[i]]@family$linkfun(trm[i,],base=object@transition[[i]]@family$base),
					object@transition[[i]]@family$linkfun(trm[i,])
				)
			}
			# update trDens slot of the model
			object@trDens[,,i] <- dens(object@transition[[i]])
		}
		
		for(i in 1:ns) {
			for(k in 1:nresp(object)) {
				object@response[[i]][[k]] <- fit(object@response[[i]][[k]],w=fbo$gamma[,i])
				# update dens slot of the model
				object@dens[,k,i] <- dens(object@response[[i]][[k]])
			}
		}
		
		# expectation
		fbo <- fb(init=object@init,A=object@trDens,B=object@dens,ntimes=ntimes(object),stationary=object@stationary)
		LL <- fbo$logLike
				
		if(verbose&((j%%5)==0)) cat("iteration",j,"logLik:",LL,"\n")
		
		if( (LL >= LL.old)) {
		  if((crit == "absolute" &&  LL - LL.old < tol) || (crit == "relative" && (LL.old - LL)/LL.old  < tol)) {
			  cat("iteration",j,"logLik:",LL,"\n")
			  converge <- TRUE
			}
		} else {
		  # this should not really happen...
		  if(j > 0) warning("likelihood decreased on iteration",j)
		}
		
		LL.old <- LL
		j <- j+1
		
	}
		
	class(object) <- "depmix.fitted"
	
	if(converge) {
		object@message <- switch(crit,
			relative = "Log likelihood converged to within tol. (relative change)",
			absolute = "Log likelihood converged to within tol. (absolute change)"
		)
	} else object@message <- "'maxit' iterations reached in EM without convergence."
	
	# no constraints in EM, except for the standard constraints ...
	# which are produced by the following (only necessary for getting df right in logLik and such)
	constraints <- getConstraints(object)
	object@conMat <- constraints$lincon
	object@lin.lower <- constraints$lin.l
	object@lin.upper <- constraints$lin.u
	
	object
}
