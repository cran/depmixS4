\name{response}

\docType{method}

\alias{response}
\alias{transInit}
\alias{GLMresponse}

\alias{transInit,formula-method}
\alias{GLMresponse,formula-method}

\alias{getdf}
\alias{getdf,response-method}

\alias{show,GLMresponse-method}

\title{Methods for creating depmix response models}

\description{

Create \code{response} objects for \code{\link{depmix}} models using
formulae and family objects.

}

\usage{
	
	GLMresponse(formula, data=NULL, family=gaussian(), pstart=NULL, 
		fixed=NULL, prob=TRUE, ...)
	
	transInit(formula, nstates, data=NULL, family=multinomial(),
		pstart=NULL, fixed=NULL, prob=TRUE, ...)
	
	\S4method{getdf}{response}(object)
}

\arguments{
	\item{formula}{A model \code{\link{formula}}.}
	\item{data}{An optional data.frame to interpret the variables from the 
		formula argument in.}	
	\item{family}{A family object; currently only multinomial() and
		gaussian() are allowed options.}
	\item{pstart}{Starting values for the coefficients and other
		parameters, ie the standard deviation for the gaussian() family.}
	\item{fixed}{Logical vector indicating which paramters are to be fixed.}
	\item{prob}{Logical indicating whether the starting values for
		multinomial() family models are probabilities or logistic
		parameters (see details).}
	\item{nstates}{The number of states of the model.}
	\item{object}{Object of class response.}
	\item{...}{Not used currently.}
}

\details{
	
	Both methods use the familiar formula interface from \code{\link{glm}}
	to specify how responses (or transition or prior parameters) depend on
	covariates/predictors.
	
	The \code{GLMresponse} model is an interface to the \code{glm}
	functions of which the functionality is leant: predict, fit and
	density functions.
	
	The \code{transInit} response model provides functionality for
	multinomial responses that are currently fit using nnet (this may
	change in the future but this should not affect the interface of this
	function).  Note that the transInit model actually lacks a reponse, ie
	the y-slot is empty, at the time of construction, as the transitions
	are not observed.
	
}

\value{
	
	\code{GLMresponse} and \code{transInit} return objects of class
	\code{GLMresponse} and \code{transInit} respectively; both classes
	extend the \code{\link{response-class}}. 
	
	\code{getdf} returns the number of free parameters of a
	\code{response} model.
	
}

\author{Ingmar Visser & Maarten Speekenbrink}

\keyword{methods}
