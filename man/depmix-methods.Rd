\name{depmix-methods}

\docType{methods}

\alias{depmix-methods}

\alias{logLik}
\alias{logLik,depmix-method}
\alias{logLik,mix-method}

\alias{logLik,depmix.fitted.classLik-method}
\alias{logLik,mix.fitted.classLik-method}

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

\alias{getmodel}
\alias{getmodel,depmix-method}
\alias{getmodel,mix-method}


\title{ 'depmix' and 'mix' methods. }

\description{

Various methods for \code{depmix} and \code{mix} objects. 

}

\usage{
	
	\S4method{logLik}{depmix}(object,method=c("fb","lystig","classification"),na.allow=TRUE)
	\S4method{logLik}{mix}(object,method=c("fb","lystig","classification"),na.allow=TRUE)
	
	\S4method{logLik}{depmix.fitted.classLik}(object,method=c("classification","fb","lystig"),na.allow=TRUE)
	\S4method{logLik}{mix.fitted.classLik}(object,method=c("classification","fb","lystig"),na.allow=TRUE)
	
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
	
	\S4method{getmodel}{depmix}(object,which="response",state=1,number=1)	
	\S4method{getmodel}{mix}(object,which="response",state=1,number=1)
}
\arguments{

  \item{object}{A \code{depmix} or \code{mix} object.}
    
  \item{values}{To be used in \code{setpars} to set new parameter values; see
  the example.}

  \item{method}{The log likelihood can be computed by either the forward
  backward algorithm (Rabiner, 1989), or by the method of Lystig and
  Hughes, 2002. The former is the default and implemented in a fast
  C routine. The forward-backward routine also computes the state and transition 
  smoothed probabilities, which are not directly neccessary for the log likelihood.  
  Those smoothed variables, and the forward and backward variables are accessible
  through the \code{\link{forwardbackward}} function. When method="classification",
  the classification likelihood is computed, which is the likelihood of the data 
  assuming the state sequence is known and equal to the maximum a posteriori state 
  sequence. The MAP state sequence is available through the \code{\link{viterbi}}
  function. The classification likelihood is comuted by default when calling the 
  logLik method on an a model fitted by maximising the classification likelihood.}
  
  \item{na.allow}{Allow missing observations? When set to FALSE,
  the logLik method will return NA in the presence of missing observations. 
  When set to TRUE, missing values will be ignored when computing the likelihood.
  When observations are partly missing (when a multivariate observation has missing
  values on only some of its dimensionis), this may give unexpected results.}

	\item{which}{\code{getpars} function: The default "pars" returns a vector
	of all parameters of a \code{depmix} object; the alternative value
	"fixed" returns a logical vector of the same length indicating which
	parameters are fixed.  The setpars functions sets parameters to new
	values; \code{setpars} also recomputes the dens, trans and init slots of
	\code{depmix} objects.  Note that the \code{getpars} and \code{setpars}
	functions for \code{depmix} objects simply call the functions of the same
	name for the response and transition models.
	
	\code{getmodel} function: possible values are "response" (the default),
	"prior" and "transition" to return the corresponding submodels.  }
	
	\item{state}{In \code{getmodel} this determines the submodel to be
	returned (together with number in the case of response models); when
	which="transition", \code{getmodel} returns the transition submodel for
	state=\code{state}.}

	\item{number}{In \code{getmodel} this determines the "response" model to
	be returned from state \code{state}.}
	  
  \item{...}{Not used currently.}

}

\value{
  
	\item{logLik}{returns a \code{logLik} object with attributes \code{df} and \code{nobs}.}
	
	\item{nobs}{returns the number of observations (used in computing the BIC).}

	\item{npar}{returns the number of paramaters of a model.}

	\item{freepars}{returns the number of non-redundant parameters.}

	\item{setpars}{returns a \code{(dep-)mix} object with new parameter values.}

	\item{getpars}{returns a vector with the current parameter values.}
	
	\item{getmodel}{returns a submodel of a \code{(dep-)mix} model; the prior
	model, one of the transition models (determined by argument state) or one
	of the response models (determined by arguments state and number).}

}

\author{Ingmar Visser & Maarten Speekenbrink}

\examples{

# create a 2 state model with one continuous and one binary response
data(speed)
mod <- depmix(list(rt~1,corr~1),data=speed,nstates=2,family=list(gaussian(),multinomial()))

getmodel(mod,"response",2,1)

getmodel(mod,"prior")

# get the loglikelihood of the model
logLik(mod)

# to see the ordering of parameters to use in setpars
mod <- setpars(mod, value=1:npar(mod))
mod

# to see which parameters are fixed (by default only baseline parameters in
# the multinomial logistic models for the transition models and the initial
# state probabilities model)
mod <- setpars(mod, getpars(mod,which="fixed"))
mod

}

\keyword{methods}

