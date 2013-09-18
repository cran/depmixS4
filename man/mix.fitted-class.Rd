\name{mix.fitted-class}

\docType{class}

\alias{mix.fitted-class}
\alias{mix.fitted.classLik-class}

\title{Class "mix.fitted" (and "mix.fitted.classLik")}

\description{A fitted \code{\link{mix}} model.}

\section{Slots}{
	
	A \code{mix.fitted} object is a \code{mix} object with three 
	additional slots, here is the complete list:
	
	\describe{

		\item{\code{response}:}{List of list of \code{response} objects.}

		\item{\code{prior}:}{\code{transInit} object.}
		
		\item{\code{dens}:}{Array of dimension sum(ntimes)*nresp*nstates
			providing the densities of the observed responses for each state.}
				
		\item{\code{init}:}{Array of dimension \code{length(ntimes)}*nstates with 
			the current predictions for the initial state probabilities.}
				
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
	Class \code{"\linkS4class{mix}"} directly. \code{mix.fitted.classLik} is
	similar to \code{mix.fitted}, the only difference being that the model is fitted
	by maximising the classification likelihood.
}

\author{Ingmar Visser & Maarten Speekenbrink}

\examples{
	showClass("mix.fitted")
}

\keyword{classes}
