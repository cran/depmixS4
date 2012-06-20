\name{depmix-methods}

\docType{methods}

\alias{depmix-methods}

\alias{logLik}
\alias{logLik,depmix-method}
\alias{logLik,mix-method}

\alias{nobs}
\alias{nobs,depmix-method}
\alias{nobs,mix-method}

\alias{npar}
\alias{npar,depmix-method}
\alias{npar,mix-method}

\alias{freepars}
\alias{freepars,depmix-method}
\alias{freepars,mix-method}
\alias{freepars,depmix.fitted-method}
\alias{freepars,mix.fitted-method}

\alias{setpars}
\alias{setpars,depmix-method}
\alias{setpars,mix-method}

\alias{getpars}
\alias{getpars,depmix-method}
\alias{getpars,mix-method}

\title{ 'depmix' and 'mix' methods. }

\description{

Various methods for \code{depmix} and \code{mix} objects. 

}

\usage{
	
	\S4method{logLik}{depmix}(object,method="fb")
	\S4method{logLik}{mix}(object,method="fb")
	
	\S4method{nobs}{depmix}(object, ...)
	\S4method{nobs}{mix}(object, ...)
	
	\S4method{npar}{depmix}(object)
	\S4method{npar}{mix}(object)
	
	\S4method{freepars}{depmix}(object)
	\S4method{freepars}{mix}(object)
	
	\S4method{setpars}{depmix}(object,values, which="pars",...)
	\S4method{setpars}{mix}(object,values, which="pars",...)
	
	\S4method{getpars}{depmix}(object,which="pars",...)	
	\S4method{getpars}{mix}(object,which="pars",...)	
}
\arguments{

  \item{object}{A \code{depmix} or \code{mix} object.}
    
  \item{values}{To be used in \code{setpars} to set new parameter values; see
  the example.}

  \item{method}{The log likelihood can be computed by either the forward
  backward algorithm from Rabiner, 1989, or by the method of Lystig and
  Hughes, 2002.  The latter is the default as it is faster because in the
  forward backward routine the state and transition smoothed probabilities
  are also computed which are not neccessary for the log likelihood.  Those
  smoothed variables, and the forward and backward variables are accessible
  through the \code{\link{forwardbackward}} function.} 

	\item{which}{The default "pars" returns a vector of all parameters of a
	\code{depmix} object; the alternative value "fixed" returns a logical
	vector of the same length indicating which parameters are fixed.  The
	setpars functions sets parameters to new values; \code{setpars} also
	recomputes the dens, trans and init slots of \code{depmix} objects.  Note
	that the \code{getpars} and \code{setpars} functions for \code{depmix}
	objects simply call the functions of the same name for the response and
	transition models.}
	  
  \item{...}{Not used currently.}

}

\value{
  
	\item{logLik}{returns a \code{logLik} object with attributes \code{df} and \code{nobs}.}
	
	\item{nobs}{returns the number of observations (used in computing the BIC).}

	\item{npar}{returns the number of paramters of a model.}

	\item{freepars}{returns the number of non-redundant parameters.}

	\item{setpars}{returns a \code{(dep-)mix} object with new parameter values.}

	\item{getpars}{returns a vector with the current parameter values.}

}

\author{Ingmar Visser}

\examples{

# create a 2 state model with one continuous and one binary response
data(speed)
mod <- depmix(list(rt~1,corr~1),data=speed,nstates=2,family=list(gaussian(),multinomial()))

# get the loglikelihood of the model
logLik(mod)

# to see the ordering of parameters to use in setpars
mod <- setpars(mod, value=1:npar(mod))
mod

# to see which parameters are fixed (by default only baseline parameters in
# the multinomial logistic models for the transition models and the initial
# state probabilities model
mod <- setpars(mod, getpars(mod,which="fixed"))
mod

}

\keyword{methods}

