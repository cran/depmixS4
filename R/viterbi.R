# 
# Maarten Speekenbrink, 23-3-2008
# 

viterbi <- function(object) {
	# returns the most likely state sequence
	nt <- sum(object@ntimes)
	lt <- length(object@ntimes)
	et <- cumsum(object@ntimes)
	bt <- c(1,et[-lt]+1)
		
	ns <- object@nstates
	
	delta <- psi <- matrix(nrow=nt,ncol=ns)
	state <- vector(length=nt)
	
	prior <- object@init
	
    A <- object@trDens
	B <- apply(log(object@dens),c(1,3),sum)
	
	for(case in 1:lt) {
		# initialization
		delta[bt[case],] <- - (log(prior[case,]) + B[bt[case],])
		psi[bt[case],] <- 0
		# recursion
		if(object@ntimes[case]>1) {
			for(i in ((bt[case]+1):et[case])) {
				for(j in 1:ns) {
					if(!object@stationary) {
						delta[i,j] <- min(delta[i-1,] - log(A[i,,j])) - B[i,j]
						k <- which.min(delta[i-1,] - log(A[i,,j]))
					} else {
						delta[i,j] <- min(delta[i-1,] - log(A[,,j])) - B[i,j]
						k <- which.min(delta[i-1,] - log(A[,,j]))
					}
					if(length(k) == 0) k <- 0
					psi[i,j] <- k
				}
			}
		}
		
		# trace maximum likely state
		state[et[case]] <- which.min(delta[et[case],])
		
		# this doesn't need a for loop does it???? FIX ME	  
		if(object@ntimes[case]>1) {
			for(i in (et[case]-1):bt[case]) {
				state[i] <- psi[i+1,state[i+1]]
			}
		}
	}
	
	delta <- data.frame(state,delta) 
	return(delta)
}

viterbi2 <- function(object) {
	# returns the most likely state sequence
	nt <- sum(object@ntimes)
	lt <- length(object@ntimes)
	et <- cumsum(object@ntimes)
	bt <- c(1,et[-lt]+1)
		
	ns <- object@nstates
	
	delta <- psi <- matrix(nrow=nt,ncol=ns)
	state <- vector(length=nt)
	
	prior <- object@init
	
	A <- object@trDens
	B <- apply((object@dens),c(1,3),prod)
	
	for(case in 1:lt) {
		# initialization
		delta[bt[case],] <- prior[case,]*B[bt[case],]
		delta[bt[case],] <- delta[bt[case],]/(sum(delta[bt[case],]))
		psi[bt[case],] <- 0
		# recursion
		if(object@ntimes[case]>1) {
			for(i in ((bt[case]+1):et[case])) {
				for(j in 1:ns) {
					if(!object@stationary) {
						delta[i,j] <- max(delta[i-1,]*(A[i,,j]))*B[i,j]
						k <- which.max(delta[i-1,]*A[i,,j])
					} else {
						delta[i,j] <- max(delta[i-1,]*A[,,j])*B[i,j]
						k <- which.max(delta[i-1,]*A[,,j])
					}
					if(length(k) == 0) k <- 0
					psi[i,j] <- k
				}
				delta[i,] <- delta[i,]/(sum(delta[i,]))
			}
		}
		
		# trace maximum likely state
		state[et[case]] <- which.max(delta[et[case],])
		
		# this doesn't need a for loop does it???? FIX ME	  
		if(object@ntimes[case]>1) {
			for(i in (et[case]-1):bt[case]) {
				state[i] <- psi[i+1,state[i+1]]
			}
		}
	}
	
	delta <- data.frame(state,delta) 	
	return(delta)
}


viterbi.fb <- function(A,B,prior) {
    # returns the most likely state sequence
    nt <- nrow(B)
    ns <- ncol(A)
    delta <- psi <- matrix(nrow=nt,ncol=ns)
    state <- vector(length=nt)
    # initialization
    delta[1,] <- - (log(prior) + log(B[1,]))
    psi[1,] <- 0
    # recursion
    for(i in 2:nt) {
        for(j in 1:ns) {
            delta[i,j] <- min(delta[i-1,] - log(A[,j])) - log(B[i,j])
            k <- which.min(delta[i-1,] - log(A[,j]))
            if(length(k) == 0) k <- 0
            psi[i,j] <- k
        }
    }
    #trace maximum likely state
    state[nt] <- which.min(delta[nt,])
    for(i in (nt-1):1) {
        state[i] <- psi[i+1,state[i+1]]
    }
    return(state)
}