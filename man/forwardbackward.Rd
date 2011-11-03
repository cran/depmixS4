\name{forwardbackward}

\alias{forwardbackward}

\alias{forwardbackward,mix-method}

\alias{forwardbackward,depmix-method}

\title{ Forward and backward variables }

\description{

	Compute the forward and backward variables of a \code{depmix} object.

}

\usage{

	\S4method{forwardbackward}{mix}(object, return.all=TRUE, useC=TRUE, ...)

}

\arguments{
  
	\item{object}{A depmix object.}
  	
	\item{return.all}{If FALSE, only gamma and xi and the log likelihood
	are returned (which are the only variables needed in using EM).}
  	
	\item{useC}{Flag used to set whether the C-code is used to compute the
	forward, backward, gamma and xi variables or not; the R-code is basically
	obsolete (but retained for now for debugging purposes).}
	
	\item{\dots}{Not currently used.}

}

\value{
	
	\code{forwardbackward} returns a list of the following (the variables 
	are named after the notation from Rabiner, 1989):

	\item{alpha}{The forward variables.}
	
	\item{beta}{The backward variables.}
	
	\item{gamma}{The smoothed transition probabilities.}
	
	\item{xi}{The smoothed state probabilities.}
	
	\item{sca}{The scale factors (called lambda in Rabiner).}
	
	\item{logLike}{The log likelihood (which is \code{sum(log(sca))}.}
	
	If return.all=FALSE, only gamma, xi and the log likelihood are returned.
}

\references{ 	

	Lawrence R. Rabiner (1989).  A tutorial on hidden Markov models and
	selected applications in speech recognition.  \emph{Proceedings of
	IEEE}, 77-2, p.  267-295. 

}

\author{ Maarten Speekenbrink & Ingmar Visser }

\examples{
	# add some later
}

\keyword{ methods }

