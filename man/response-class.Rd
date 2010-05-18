\name{response-class}

\docType{class}

\alias{response-class}

\title{Class "response"}

\description{

A generic \code{response} model for \code{\link{depmix}} models.

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
		
		\item{\code{constr}:}{Information on constraints.}
	}
	
}

\details{

	This class offers a framework from which to build specific response
	models such as glm based responses or multinomial responses.  For
	extensibility, objects with class \code{response} should have at least
	methods: \code{dens} to return the \code{dens}'ity of responses, and
	\code{getpars} and \code{setpars} methods to get and set parameters.
	
	The \code{constr} slot is used for information on constraints that are
	inherently part of a model; the only such constraints which are currently 
	used are 1) the sum constraint in multinomial models with identity link, 
	and 2) a lower bound of zero of sd parameters in gaussian distributions. 
	Such constraints are only used in fitting models with general optimization
	routines such as \code{Rsolnp}; the EM algorithm automagically respects the
	sum constraint. 
	\describe{
		\item{\code{lin}:}{Derivative of linear constraint.}
		\item{\code{linup}:}{Upper bounds for linear constraints.}
		\item{\code{linlow}:}{Lower bounds for linear constraints.}
		\item{\code{parup}:}{Upper bounds on parameters.}
		\item{\code{parlow}:}{Lower bounds on parameters.}
	}
	
}

\section{Accessor Functions}{
	The following functions should be used for accessing the corresponding
	slots:
	\describe{
		\item{\code{npar}:}{The number of parameters of the model.}
		\item{\code{getdf}:}{The number of non-fixed parameters.}
	}
}

\author{Ingmar Visser & Maarten Speekenbrink}

\keyword{classes}
