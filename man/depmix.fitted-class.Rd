\name{depmix.fitted-class}

\docType{class}

\alias{depmix.fitted}
\alias{depmix.fitted-class}

\alias{depmix.fitted.classLik}
\alias{depmix.fitted.classLik-class}

\title{Class "depmix.fitted" (and  "depmix.fitted.classLik")}

\description{A fitted \code{\link{depmix}} model.}

\section{Slots}{
	
	A \code{depmix.fitted} object is a \code{depmix} object with three 
	additional slots, here is the complete list:
	
	\describe{

		\item{\code{response}:}{List of list of \code{response} objects.}

		\item{\code{transition}}{List of \code{transInit} objects.}
		
		\item{\code{prior}:}{\code{transInit} object.}
		
		\item{\code{dens}:}{Array of dimension sum(ntimes)*nresp*nstates
			providing the densities of the observed responses for each state.}
		
		\item{\code{trDens}:}{Array of dimension \code{sum(ntimes)}*nstates
			providing the probability of a state transition depending on the
			predictors.}
		
		\item{\code{init}:}{Array of dimension \code{length(ntimes)}*nstates with 
			the current predictions for the initial state probabilities.}
		
		\item{\code{stationary}:}{Logical indicating whether the transitions are
			time-dependent or not; for internal use.}
		
		\item{\code{ntimes}:}{A vector containing the lengths of independent time
			series; if data is provided, sum(ntimes) must be equal to
			nrow(data).}
	
		\item{\code{nstates}:}{The number of states of the model.}
	
		\item{\code{nresp}:}{The number of independent responses.}
	
		\item{\code{npars}:}{The total number of parameters of the model. This is not
			the degrees of freedom, ie there are redundancies in the
			parameters, in particular in the multinomial models for the
			transitions and prior.}
			
		\item{\code{message}:}{This provides some information on convergence,
			either from the EM algorithm or from Rdonlp2.}
		
		\item{\code{conMat}:}{The linear constraint matrix, which has zero rows
			if there were no constraints.}
		
		\item{\code{lin.lower}}{The lower bounds on the linear constraints.}

		\item{\code{lin.upper}}{The upper bounds on the linear constraints.}		
		
		\item{\code{posterior}:}{Posterior (Viterbi) state sequence.}
	}
}

\section{Details}{
	
	The print function shows some convergence information, and the summary
	method shows the parameter estimates.
	
}

\section{Extends}{
	\code{depmix.fitted} extends the \code{"\linkS4class{depmix}"} class directly. \code{depmix.fitted.classLik} is
	similar to \code{depmix.fitted}, the only difference being that the model is fitted
	by maximising the classification likelihood.
}

\author{Ingmar Visser & Maarten Speekenbrink}

\keyword{classes}
