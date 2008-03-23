
\name{depmix}

\docType{methods}

\alias{depmixS4}
\alias{depmix}
\alias{depmix,ANY-method}

\alias{logLik}
\alias{logLik,depmix-method}

\alias{AIC}
\alias{AIC,depmix-method}

\alias{BIC}
\alias{BIC,depmix-method}

\alias{nobs}
\alias{nobs,depmix-method}

\alias{npar}
\alias{npar,depmix-method}

\alias{freepars}
\alias{freepars,depmix-method}

\alias{setpars}
\alias{setpars,depmix-method}

\alias{getpars}
\alias{getpars,depmix-method}

\alias{show,depmix-method}
\alias{summary,depmix-method}

\title{Dependent Mixture Model Specifiction}

\description{
		
	depmixS4 is a framework for specifying and fitting dependent mixture
	models, otherwise known as hidden or latent Markov models.
	Optimization is done with the EM algorithm or optionally with Rdonlp2
	when constraints on parameters are imposed.  Models for different
	distributions can easily be added.
	
	\code{depmix} creates an object of class \code{depmix}, a
	dependent mixture model, otherwise known as hidden Markov 
	model.
	
}

\usage{
	
	depmix(response, data=NULL, nstates, transition=~1, family=gaussian(), 
		prior=~1, initdata=NULL, respstart=NULL, trstart=NULL, instart=NULL,
		ntimes=NULL,...)
	
	\S4method{summary}{depmix}(object)
	\S4method{logLik}{depmix}(object)
	\S4method{AIC}{depmix}(object)
	\S4method{BIC}{depmix}(object)
	\S4method{nobs}{depmix}(object)
	
	\S4method{npar}{depmix}(object)
	\S4method{freepars}{depmix}(object)
	\S4method{setpars}{depmix}(object,which="pars",...)
	\S4method{getpars}{depmix}(object,which="pars",...)	
	
}

\arguments{
	\item{object}{An object with class \code{depmix}.}
	
	\item{response}{The response to be modeled; either a formula or a list 
		of formulae in the multivariate case. See details.}
	
	\item{data}{An optional data.frame to interpret the variables in
		response and transition.}
	
	\item{nstates}{The number of states of the model.}

	\item{transition}{A one-sided formula specifying the model for the 
		transitions. See details.}
	
	\item{family}{A family argument for the response. This must be a list
		of family's if the response is multivariate.}
	
	\item{prior}{A one-sided formula specifying the density for the prior
		or initial state probabilities.}
	
	\item{initdata}{An optional data.frame to interpret the variables
		occuring in prior. The number of rows of this data.frame must be
		equal to the number of cases being modeled. See details.}
	
	\item{respstart}{Starting values for the parameters of the response
		models.}
	
	\item{trstart}{Starting values for the parameters of the transition
		models.}
	
	\item{instart}{Starting values for the parameters of the prior or
		initial state probability model.}

	\item{ntimes}{A vector specifying the lengths of individual, ie
		independent, time series. If not specified, the responses are
		assumed to form a single time series. If the data argument has an
		attribute ntimes, then this is used.}
	
	\item{which}{The default "pars" returns a vector of all parameters of a
		\code{depmix} object; the alternative value "fixed" return a
		logical vector of the same length indicating which parameters are
		fixed. The setpars functions sets parameters (or the logical fixed 
		vector) to new values; \code{setpars} also recomputes the dens,
		trans and init slots of \code{depmix} objects. Note that the 
		\code{getpars} and \code{setpars} functions for \code{depmix} 
		objects simply call the functions of the same name for the response
		and transition models.}
		
	\item{...}{Not used currently.}
	
}

\details{ 

	The function \code{depmix} creates an S4 object of class \code{depmix},
	which needs to be fitted using \code{\link[depmixS4]{depmix.fit}} to
	optimize the parameters.
	
	The response model(s) are created by call(s) to
	\code{\link[depmixS4]{response}} providing the family and optional
	predictors.  If response is a list of formulae, the response's are
	assumed to be independent conditional on the latent state.
	
	The transitions are modeled as a multinomial logistic model for each
	state.  Hence, the transition matrix can be modeled as time-dependent,
	depending on predictors.  The prior density is also modeled as a
	multinomial logistic.  Both are created by calls to
	\code{\link[depmixS4]{transInit}}.
	
	Starting values may be provided by the respective arguments.  The order
	in which parameters must be provided can be easily studied by using the
	\code{setpars} function (see example).
	
	Linear constraints on parameters can be provided as argument to the
	\code{\link[depmixS4]{depmix.fit}} function.
	
}

\value{
	
\code{depmix} returns an object of class \code{depmix} which has the
following slots:

	\item{response}{A list of a list of response models; the first
		index runs over states; the second index runs over the independent 
		responses in case a multivariate response is provided.}
	
	\item{transition}{A list of \code{transInit} models, ie multinomial
		logistic models with length the number of states.}

	\item{prior}{A multinomial logistic model for the initial state
		probabilities.}

	\item{dens,trDens,init}{See depmix-class help for details. For internal
		use.}
	
	\item{stationary}{Logical indicating whether the transitions are
		time-dependent or not; for internal use.}
	
	\item{ntimes}{A vector containing the lengths of independent time
		series; if data is provided, sum(ntimes) must be equal to
		nrow(data).}
	
	\item{nstates}{The number of states of the model.}
	
	\item{nresp}{The number of independent responses.}

	\item{npars}{The total number of parameters of the model. This is not
		the degrees of freedom, ie there are redundancies in the
		parameters, in particular in the multinomial models for the
		transitions and prior.}

\code{logLik}, \code{AIC}, and \code{BIC} return the respective values
associated with the current parameter values.  \code{nobs} returns the
number of observations, ie \code{sum(ntimes)} that is used in computing the
\code{BIC}.  \code{npar} returns the number of paramters of a model;
\code{freepars} returns the number of non-fixed parameters.

}

\author{Ingmar Visser \email{i.visser@uva.nl}} 

\seealso{

	\code{\link[depmixS4]{depmix.fit}}, \code{\link{transInit}},
	\code{\link{response}}.

}

\references{
	
	On hidden Markov models: Lawrence R. Rabiner (1989).  A tutorial on
	hidden Markov models and selected applications in speech recognition.
	\emph{Proceedings of IEEE}, 77-2, p.  267-295.
	
	On latent class models: A. L. McCutcheon (1987).  \emph{Latent class
	analysis}.  Sage Publications.
	
}

\examples{

# create a 2 state model with one continuous and one binary response
data(speed)
mod <- depmix(list(rt~1,corr~1),data=speed,nstates=2,family=list(gaussian(),multinomial()))
# print the model, formulae and parameter values
mod

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


