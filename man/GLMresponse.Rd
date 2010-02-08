\name{GLMresponse}

\docType{method}

\alias{GLMresponse}
\alias{GLMresponse,formula-method}

\alias{getdf}
\alias{getdf,response-method}

\alias{show,GLMresponse-method}

\title{Methods for creating depmix response models}

\description{

	Create \code{GLMresponse} objects for \code{\link{depmix}} models using
	formulae and family objects.

}

\usage{
	
	GLMresponse(formula, data=NULL, family=gaussian(), pstart=NULL, 
		fixed=NULL, prob=TRUE, ...)
	
	\S4method{getdf}{response}(object)
}

\arguments{
	\item{formula}{A model \code{\link{formula}}.}
	\item{data}{An optional data.frame to interpret the variables from the 
		formula argument in.}	
	\item{family}{A family object; }
	\item{pstart}{Starting values for the coefficients and other
		parameters, e.g. the standard deviation for the gaussian() family.}
	\item{fixed}{Logical vector indicating which paramters are to be fixed.}
	\item{prob}{Logical indicating whether the starting values for
		multinomial() family models are probabilities or logistic
		parameters (see details).}
	\item{object}{Object of class response.}
	\item{...}{Not used currently.}
}

\details{
	
	\code{GLMresponse} is the default driver for specifying response
	distributions of \code{depmix} models.  It uses the familiar formula
	interface from \code{\link{glm}} to specify how responses depend on
	covariates/predictors.

	Currently available options for the family argument are
	\code{binomial}, \code{gaussian}, \code{poisson}, \code{Gamma}, and
	\code{multinomial}.  Except for the latter option, the
	\code{GLMresponse} model is an interface to the \code{glm} functions of
	which the functionality is used: predict, fit and density functions.
	
	The \code{multinomial} model takes as link functions \code{mlogit}, the
	default, and then uses functionality from the \code{nnet} package to
	fit multinomial logistic models; it also takes \code{identity} as a
	link function.  The latter is typically faster and is hence preferred
	when no covariates are present.  Note that the sum constraint inherent
	in such a model is automatically respected in the EM algorithm, but
	must be supplied by the user in case of using Rdonlp optimization.
	
	See the \code{\link{responses}} help page for examples. 
	
}

\value{
	
	\code{GLMresponse} returns an object of class \code{GLMresponse} which
	extends the \code{\link{response-class}}.
	
	\code{getdf} returns the number of free parameters of a
	\code{response} model.
	
}

\seealso{
	
	\code{\link{makeDepmix}} has an example of specifying a model with a 
	multivariate normal response and an example of how to add a user-defined
	response model, in particular an ex-gauss distribution used for the
	\code{speed} data. 
}

\author{Ingmar Visser & Maarten Speekenbrink}

\keyword{methods}
