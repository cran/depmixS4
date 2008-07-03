\name{forwardbackward}

\alias{forwardbackward}

\alias{forwardbackward,depmix-method}

\title{ Forward and backward variables }

\description{

	Compute the forward and backward variables of a \code{depmix} object.

}

\usage{

	\S4method{forwardbackward}{depmix}(object, return.all=TRUE, ...)

}

\arguments{
  
	\item{object}{A depmix object.}
  	
	\item{return.all}{If FALSE, only gamma and xi and the log likelihood
	are returned (which are the only variables needed in using EM).}
  	
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

