
# 
# Ingmar Visser, 11-6-2008
# 

# 
# DEPMIX CLASS BELOW THE MIX CLASS
# 

# 
# Class definition, accessor functions, print and summary methods
# 

# 
# MIX CLASS
# 

setClass("mix",
	representation(response="list", # response models
		prior="ANY", # the prior model (multinomial logistic)
		dens="array", # response densities (B)
		init="array", # usually called pi 
		nstates="numeric",
		nresp="numeric",
		ntimes="numeric",
		npars="numeric" # number of parameters
	)
)

# accessor functions
setMethod("npar","mix",
	function(object) return(object@npars)
)

setMethod("ntimes","mix",
	function(object) return(object@ntimes)
)

setMethod("nstates","mix",
	function(object) return(object@nstates)
)

setMethod("nresp","mix",
	function(object) return(object@nresp)
)


# 
# PRINT method
# 

setMethod("show","mix",
	function(object) {
		cat("Initial state probabilties model \n")
		print(object@prior)
		cat("\n")
		for(i in 1:object@nstates) {
			cat("Response model(s) for state", i,"\n\n")
			for(j in 1:object@nresp) {
				cat("Response model for response",j,"\n")
				print(object@response[[i]][[j]])
				cat("\n")
			}
			cat("\n")
		}
	}
)

# 
# SUMMARY method: to do
# 


# 
# Ingmar Visser, 23-3-2008
# 

# 
# Class definition, accessor functions, print and summary methods
# 

# 
# DEPMIX CLASS
# 

setClass("depmix",
	representation(transition="list", # transition models (multinomial logistic)
		trDens="array", # transition densities (A)
		stationary="logical"
	),
	contains="mix"
)

# 
# PRINT method
# 

setMethod("show","depmix",
	function(object) {
		cat("Initial state probabilties model \n")
		print(object@prior)
		cat("\n")
		for(i in 1:object@nstates) {
			cat("Transition model for state (component)", i,"\n")
			print(object@transition[[i]])
			cat("\n")
		}
		cat("\n")
		for(i in 1:object@nstates) {
			cat("Response model(s) for state", i,"\n\n")
			for(j in 1:object@nresp) {
				cat("Response model for response",j,"\n")
				print(object@response[[i]][[j]])
				cat("\n")
			}
			cat("\n")
		}
	}
)

setMethod("is.stationary",signature(object="depmix"),
  function(object) {
		return(object@stationary)
	}
)

setMethod("simulate",signature(object="depmix"),
  function(object,nsim=1,seed=NULL,...) {
    if(!is.null(seed)) set.seed(seed)
    ntim <- ntimes(object)
   	nt <- sum(ntim)
  	lt <- length(ntim)
  	et <- cumsum(ntim)
  	bt <- c(1,et[-lt]+1)

  	nr <- nresp(object)
  	ns <- nstates(object)

    # simulate state sequences first, then observations

    # random generation is slow when done separately for each t, so first draw
    #   variates for all t, and then determine state sequences iteratively
  	states <- array(,dim=c(nt,nsim))
  	states[bt,] <- simulate(object@prior,n=nsim,is.prior=T)
  	sims <- array(,dim=c(nt,ns,nsim))
  	for(i in 1:ns) {
      if(is.stationary(object)) {
        # TODO: this is a temporary fix!!! 
        sims[,i,] <- simulate(object@transition[[i]],nsim=nsim,times=rep(1,nt))
      } else {
        sims[,i,] <- simulate(object@transition[[i]],nsim=nsim)
      }
 	  }
 	  # track states
  	for(case in 1:lt) {
      for(i in (bt[case]+1):et[case]) {
        states[i,] <- sims[cbind(i,states[i-1,],1:nsim)]
      }
    }

    states <- as.vector(states)
    responses <- list(length=nr)
    #responses <- array(,dim=c(nt,nr,nsim))
    for(i in 1:nr) {
      tmp <- matrix(,nrow=nt*nsim,ncol=NCOL(object@response[[1]][[i]]@y))
      for(j in 1:ns) {
        tmp[states==j,] <- simulate(object@response[[j]][[i]],nsim=nsim)[states==j,]
      }
      responses[[i]] <- tmp
    }

    # generate new depmix.sim object
    class(object) <- "depmix.sim"
    object@states <- as.matrix(states)

    object@prior@x <- as.matrix(apply(object@prior@x,2,rep,nsim))
    for(j in 1:ns) {
      if(!is.stationary(object)) object@transition[[j]]@x <- as.matrix(apply(object@transition[[j]]@x,2,rep,nsim))
      for(i in 1:nr) {
        object@response[[j]][[i]]@y <- as.matrix(responses[[i]])
        object@response[[j]][[i]]@x <- as.matrix(apply(object@response[[j]][[i]]@x,2,rep,nsim))
      }
    }
    object@ntimes <- rep(object@ntimes,nsim)

   	# make appropriate array for transition densities
  	nt <- sum(object@ntimes)
  	if(is.stationary(object)) trDens <- array(0,c(1,ns,ns)) else trDens <- array(0,c(nt,ns,ns))

  	# make appropriate array for response densities
  	dns <- array(,c(nt,nr,ns))

  	# compute observation and transition densities
  	for(i in 1:ns) {
  		for(j in 1:nr) {
  			dns[,j,i] <- dens(object@response[[i]][[j]]) # remove this response as an argument from the call to setpars
  		}
  		trDens[,,i] <- dens(object@transition[[i]])
  	}

  	# compute initial state probabilties
  	object@init <- dens(object@prior)
    object@trDens <- trDens
    object@dens <- dns

    return(object)
  }
)



# 
# SUMMARY method: to do
# 



