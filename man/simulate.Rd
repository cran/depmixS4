\name{simulate}

\docType{method}

\alias{simulate}

\alias{simulate,response-method}
\alias{simulate,mix-method}
\alias{simulate,depmix-method}
\alias{simulate,GLMresponse-method}
\alias{simulate,transInit-method}

\alias{simulate,BINOMresponse-method}
\alias{simulate,GAMMAresponse-method}
\alias{simulate,MULTINOMresponse-method}
\alias{simulate,MVNresponse-method}
\alias{simulate,NORMresponse-method}
\alias{simulate,POISSONresponse-method}

\title{Methods to simulate from (dep-)mix models}

\description{

	Random draws from \code{(dep-)mix} objects.

}

\usage{

  \S4method{simulate}{depmix}(object, nsim=1, seed=NULL, ...)
  
  \S4method{simulate}{mix}(object, nsim=1, seed=NULL, ...)

  \S4method{simulate}{response}(object, nsim=1, seed=NULL, times, ...)
  
  \S4method{simulate}{GLMresponse}(object, nsim=1, seed=NULL, times, ...)
  
  \S4method{simulate}{transInit}(object, nsim=1, seed=NULL, times, is.prior=FALSE, ...)
}

\arguments{
	
	\item{object}{Object to generate random draws.  An object of class
	\code{mix}, \code{depmix}, \code{response} or \code{transInit}}

	\item{nsim}{The number of draws (one draw simulates a data set of the
	size that is defined by ntimes); defaults to 1.}
	
	\item{seed}{Set the seed.}
	
	\item{times}{(optional) An indicator vector indicating for which times in the
	complete series to generate the data.  For internal use.}
	
	\item{is.prior}{For \code{transInit} objects, indicates whether it is a prior 
	(init) model, or not (i.e., it is a transition model)}
	
	\item{...}{Not used currently.}
}

\details{

	For a \code{depmix} model, simulate generates \code{nsim} random state
	sequences, each of the same length as the observation sequence in the
	\code{depmix} model (i.e., \code{sum(ntimes(object))}.  The state
	sequences are then used to generate \code{nsim} observation sequences
	of thee same length.

	For a \code{mix} model, simulate generates \code{nsim} random class 
	assignments for each case. Those assigments are then used to generate
	observation/response values from the appropriate distributions. 

	Setting the \code{times} option selects the time points in the total
	state/observation sequence (i.e., counting is continued over ntimes).
	Direct calls of simulate with the \code{times} option are not recommended.

}

\value{

	For a \code{depmix} object, a new object of class \code{depmix.sim}. 
	
	For a \code{transInit} object, a state sequence.
	
	For a \code{response} object, an observation sequence.

}

\examples{

y <- rnorm(1000)
respst <- c(0,1,2,1)
trst <- c(0.9,0.1,0.1,0.9)

df <- data.frame(y=y)

mod <- depmix(y~1,data=df,respst=respst,trst=trst,inst=c(0.5,0.5),nti=1000,nst=2)

mod <- simulate(mod)

}
\author{Maarten Speekenbrink}

\keyword{methods}
