
\name{mix}

\docType{methods}

\alias{mix}
\alias{mix,ANY-method}

\alias{show,mix-method}
\alias{summary,mix-method}

\title{ Mixture Model Specifiction }

\description{
	
	\code{mix} creates an object of class \code{mix}, an (independent)
	mixture model (as a limit case of dependent mixture models in which all
	observed time series are of length 1), otherwise known latent class or
	mixture model.  For a short description of the package see
	\code{\link{depmixS4}}.
	
}

\usage{
	
	mix(response, data=NULL, nstates, family=gaussian(), 
		prior=~1, initdata=NULL, respstart=NULL, instart=NULL,...)	
	
}

\arguments{
	
	\item{response}{The response to be modeled; either a formula or a list 
		of formulae in the multivariate case; this interfaces to the glm
		distributions. See 'Details'.}
	
	\item{data}{An optional data.frame to interpret the variables in
		the response and transition arguments.}
	
	\item{nstates}{The number of states of the model.}
	
	\item{family}{A family argument for the response. This must be a list
		of family's if the response is multivariate.}
	
	\item{prior}{A one-sided formula specifying the density for the prior
		or initial state probabilities.}
	
	\item{initdata}{An optional data.frame to interpret the variables
		occuring in prior. The number of rows of this data.frame must be
		equal to the number of cases being modeled. See 'Details'.}
	
	\item{respstart}{Starting values for the parameters of the response
		models.}
		
	\item{instart}{Starting values for the parameters of the prior or
		initial state probability model.}
		
	\item{...}{Not used currently.}
		
}

\details{ 

	The function \code{mix} creates an S4 object of class \code{mix},
	which needs to be fitted using \code{\link{fit}} to optimize the
	parameters.
	
	The response model(s) are by default created by call(s) to
	\code{GLMresponse} using the \code{formula} and the \code{family}
	arguments, the latter specifying the error distribution.  See
	\code{\link{GLMresponse}} for possible values of the \code{family}
	argument for \code{glm}-type responses (ie a subset of the \code{glm}
	family options, and the multinomial).  Alternative response
	distributions are specified by using the \code{\link{makeDepmix}}
	function.  Its help page has examples of specifying a model with a
	multivariate normal response, as well as an example of adding a
	user-defined response model, in this case for the ex-gauss
	distribution.
	
	If \code{response} is a list of formulae, the \code{response}'s are
	assumed to be independent conditional on the latent state.
	
	The prior density is modeled as a multinomial logistic.  This model is
	created by a call to \code{\link{transInit}}.
	
	Starting values may be provided by the respective arguments.  The order
	in which parameters must be provided can be easily studied by using the
	\code{\link{setpars}} and \code{\link{getpars}} functions.
	
	Linear constraints on parameters can be provided as argument to the
	\code{\link{fit}} function.
	
	The print function prints the formulae for the response and prior
	models along with their parameter values.
	
}

\value{
	
\code{mix} returns an object of class \code{mix} which has the
following slots:

	\item{response}{A list of a list of response models; the first
		index runs over states; the second index runs over the independent 
		responses in case a multivariate response is provided.}
	
	\item{prior}{A multinomial logistic model for the initial state
		probabilities.}

	\item{dens,init}{See \code{\link{mix-class}} help for details.  For
	internal use.}
		
	\item{ntimes}{A vector made by \code{rep(1,nrow(data))}; for internal
	use only.}
	
	\item{nstates}{The number of states of the model.}
	
	\item{nresp}{The number of independent responses.}

	\item{npars}{The total number of parameters of the model.  Note: this
	is \emph{not} the degrees of freedom because there are redundancies in the
	parameters, in particular in the multinomial models for the transitions
	and prior probabilities.}

}

\author{Ingmar Visser}

\seealso{

	\code{\link{fit}}, \code{\link{transInit}}, \code{\link{GLMresponse}},
	\code{\link{depmix-methods}} for accessor functions to \code{depmix}
	objects.

}

\references{
	
	A. L. McCutcheon (1987).  \emph{Latent class analysis}.  Sage
	Publications.
	
}

\examples{

# four binary items on the balance scale task
data(balance)

# define a latent class model
instart=c(0.5,0.5)
set.seed(1)
respstart=runif(16)
# note that ntimes argument is used to make this a mixture model
mod <- mix(list(d1~1,d2~1,d3~1,d4~1), data=balance, nstates=2,
	family=list(multinomial(),multinomial(),multinomial(),multinomial()),
	respstart=respstart,instart=instart)
# to see the model
mod

}

\keyword{methods}


