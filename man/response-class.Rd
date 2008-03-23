\name{response-class}

\docType{class}

\alias{response-class}

\title{Class "response"}

\description{

A generic \code{\link{response}} model for \code{\link{depmix}} models.

}

\arguments{
	\item{object}{An object of class "response".}
}

\section{Slots}{

	\describe{
    	\item{\code{parameters}:}{A (named) list of parameters.}

		\item{\code{fixed}:}{A logical vector indicating which parameters
			are fixed.}

		\item{\code{y}:}{A matrix with the actual response; possibly
			multivariate.}

		\item{\code{x}:}{A design matrix; possibly only an intercept term.}

		\item{\code{npar}:}{The number of parameters.}
	}
	
}

\details{

	This class offers a framework from which to build specific response
	models such as glm based responses or multinomial responses.  For
	extensibility, objects with class \code{response} should have at least
	methods: \code{dens} to return the \code{dens}'ity of responses, and
	\code{getpars} and \code{setpars} methods to get and set parameters.
	
}

\section{Accessor Functions}{
	The following functions should be used for accessing the corresponding
	slots:
	\describe{
		\item{\code{npar}:}{The number of parameters of the model.}
		\item{\code{getdf}:}{The number of non-fixed parameters.}
	}
}

\author{Ingmar Visser}

\keyword{classes}
