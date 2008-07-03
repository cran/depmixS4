makeResponseModels <-
function(response,data=NULL,nstates,family,values=NULL,prob=TRUE,...) {
	
	resp <- response
	response <- list()
	nresppars <- 0
		
	# univariate response data
	if(class(resp)=="formula") {
		nresp <- 1
		for(i in 1:nstates) {
			response[[i]] <- list()
			response[[i]][[1]] <- GLMresponse(resp,data=data,family=family)
			nresppars <- nresppars + npar(response[[i]][[1]])
		}
	}
	
	# multivariate response data
	if(is.list(resp)) {
		nresp <- length(resp)
		for(i in 1:nstates) {
			response[[i]] <- list()
			for(j in 1:nresp) {
				response[[i]][[j]] <- GLMresponse(resp[[j]],data=data,family=family[[j]])
				nresppars <- nresppars + npar(response[[i]][[j]])
			}
		}
	}
	
	# set the starting values, if any
	if(!is.null(values)) {
		if(!(length(values)==nresppars)) stop(paste("'respstart' has incorrect length, it should be", nresppars, "\n"))
		for(i in 1:nstates) {
			for(j in 1:nresp) {
				bp <- npar(response[[i]][[j]])
				response[[i]][[j]] <- setpars(response[[i]][[j]],val=values[1:bp],prob=prob)
				bp <- bp+1
				values <- values[bp:length(values)]
			}
		}
	}
	
	return(response)
}

