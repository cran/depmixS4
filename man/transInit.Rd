\name{transInit}

\docType{method}

\alias{transInit}

\alias{transInit,formula-method}

\alias{getdf,transInit-method}

\title{Methods for creating depmix transition and initial probability models}

\description{

Create \code{transInit} objects for \code{\link{depmix}} models using
formulae and family objects.

}

\usage{
	
	transInit(formula, nstates, data=NULL, family=multinomial(),
		pstart=NULL, fixed=NULL, prob=TRUE, ...)
		
		\S4method{getdf}{transInit}(object)

}

\arguments{
	\item{formula}{A model \code{\link{formula}}.}
	\item{data}{An optional data.frame to interpret the variables from the 
		formula argument in.}
	\item{family}{A family object; see details.}
	\item{pstart}{Starting values for the coefficients.}
	\item{fixed}{Logical vector indicating which paramters are to be fixed.}
	\item{prob}{Logical indicating whether the starting values for
		multinomial() family models are probabilities or logistic
		parameters (see details).}
	\item{nstates}{The number of states of the model.}
	\item{object}{An object of class \code{transInit}.}
	\item{...}{Not used currently.}
}

\details{
	
	The \code{transInit} model provides functionality for the multinomial
	probabilities of the transitions between states, as well as for the
	prior or initial probabilities.  These probabilities may depend on
	(time-varying) covariates.  The model can be used with link functions
	\code{mlogit} and \code{identity}; the latter is the default when no
	covariates are.  With the \code{mlogit} link function, the transition
	probabilities are modeled as baseline logistic multinomials (see
	Agresti, 2002, p.  272 ff.).
	
	Start values for the parameters may be provided using the \code{pstart}
	argument; these can be provided as probabilities, the default option,
	or as baseline logistic parameters, use the \code{prob} argument to
	specify the chosen option.  The default baseline category is set to 1,
	which can be modified through calling, say, family=multinomial(base=2).
	
	Note that the transInit model extends the \code{\link{response-class}},
	but that it actually lacks a reponse, i.e. the y-slot is empty, at the
	time of construction, as the transitions are not observed.
	
	\code{getdf} returns the number of free parameters of a transInit model.
	
}

\value{
	
	\code{transInit} return objects of class \code{transInit}; this class
	extends the \code{\link{response-class}}.
	
}

\references{

	Agresti, A. (2002).  \emph{Categorical Data Analysis}.  Wiley series in
	probability and mathematical statistics.  Wiley-Interscience, Hoboken,
	NJ, 2 edition.

}

\author{Ingmar Visser & Maarten Speekenbrink}

\keyword{methods}
