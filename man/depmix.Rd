
\name{depmix}

\docType{methods}

\alias{depmix}
\alias{depmix,ANY-method}

\alias{show}
\alias{summary}

\alias{show,depmix-method}
\alias{summary,depmix-method}

\title{Dependent Mixture Model Specifiction}

\description{
	
	\code{depmix} creates an object of class \code{depmix}, a dependent
	mixture model, otherwise known as hidden Markov model.  For a short
	description of the package see \code{\link{depmixS4}}. See the vignette 
	for an introduction to hidden Markov models and the package. 
	
}

\usage{
	
	depmix(response, data=NULL, nstates, transition=~1, family=gaussian(), 
		prior=~1, initdata=NULL, respstart=NULL, trstart=NULL, instart=NULL,
		ntimes=NULL,...)
	
}

\arguments{
	
	\item{response}{The response to be modeled; either a formula or a list 
		of formulae (in the multivariate case); this interfaces to the glm
		and other distributions. See 'Details'.}
	
	\item{data}{An optional \code{data.frame} to interpret the variables in
		the response and transition arguments.}
	
	\item{nstates}{The number of states of the model.}

	\item{transition}{A one-sided formula specifying the model for the 
		transitions. See 'Details'.}
	
	\item{family}{A family argument for the response. This must be a list
		of \code{family}'s if the response is multivariate.}
	
	\item{prior}{A one-sided formula specifying the density for the prior
		or initial state probabilities.}
	
	\item{initdata}{An optional data.frame to interpret the variables
		occuring in \code{prior}.  The number of rows of this data.frame
		must be equal to the number of cases being modeled,
		\code{length(ntimes)}.  See 'Details'.}
	
	\item{respstart}{Starting values for the parameters of the response
		models.}
	
	\item{trstart}{Starting values for the parameters of the transition
		models.}
	
	\item{instart}{Starting values for the parameters of the prior or
		initial state probability model.}

	\item{ntimes}{A vector specifying the lengths of individual, i.e.
		independent, time series. If not specified, the responses are
		assumed to form a single time series, i.e. \code{ntimes=nrow(data)}.
		If the \code{data} argument has an attribute \code{ntimes}, then this 
		is used. The first example in \code{\link{fit}} uses this argument.}
		
	\item{...}{Not used currently.}
		
}

\details{ 

	The function \code{depmix} creates an S4 object of class \code{depmix},
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
	
	The transitions are modeled as a multinomial logistic model for each
	state.  Hence, the transition matrix can be modeled using time-varying
	covariates.  The prior density is also modeled as a multinomial
	logistic.  Both of these models are created by calls to
	\code{\link{transInit}}.
	
	Starting values for the initial, transition, and response models may be
	provided by their respective arguments.  NB: note that the starting
	values for the initial and transition models as well as of the
	multinomial logit response models are interpreted as probabilities, and
	internally converted to multinomial logit parameters.  The order in
	which parameters must be provided can be easily studied by using the
	\code{\link{setpars}} and \code{\link{getpars}} functions.
	
	Linear constraints on parameters can be provided as argument to the
	\code{\link{fit}} function.
	
	The print function prints the formulae for the response, transition and
	prior models along with their parameter values.
	
	Missing values are allowed in the data, but missing values in the 
	covariates lead to errors. 
	
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

	\item{dens,trDens,init}{See \code{\link{depmix-class}} help for
		details.  For internal use.}
	
	\item{stationary}{Logical indicating whether the transitions are
		time-dependent or not; for internal use.}
	
	\item{ntimes}{A vector containing the lengths of independent time
		series.}
	
	\item{nstates}{The number of states of the model.}
	
	\item{nresp}{The number of independent responses.}

	\item{npars}{The total number of parameters of the model.  Note: this
	is \emph{not} the degrees of freedom because there are redundancies in the
	parameters, in particular in the multinomial models for the transitions
	and prior probabilities.}

}

\note{
	Models are not fitted; the return value of \code{depmix} is a model
	specification without optimized parameter values. Use the \code{\link{fit}}
	function to optimize parameters, and to specify additional constraints.
}

\author{Ingmar Visser & Maarten Speekenbrink}

\seealso{

	\code{\link{fit}}, \code{\link{transInit}}, \code{\link{GLMresponse}},
	\code{\link{depmix-methods}} for accessor functions to \code{depmix}
	objects. 	
	
	For full control see the \code{\link{makeDepmix}} help page and its
	example section for the possibility to add user-defined response
	distributions. 

}

\references{
	
	Ingmar Visser and Maarten Speekenbrink (2010).  depmixS4: An R Package for
	Hidden Markov Models. \emph{Journal of Statistical Software, 36(7)}, p. 1-21.

	Lawrence R. Rabiner (1989).  A tutorial on hidden Markov models and
	selected applications in speech recognition.  \emph{Proceedings of
	IEEE}, 77-2, p.  267-295.
	
}

\examples{

# create a 2 state model with one continuous and one binary response
# ntimes is used to specify the lengths of 3 separate series
data(speed)	
mod <- depmix(list(rt~1,corr~1),data=speed,nstates=2,
    family=list(gaussian(),multinomial("identity")),ntimes=c(168,134,137))
# print the model, formulae and parameter values
mod
set.seed(1)
# fit the model by calling fit
fm <- fit(mod)

# Volatility of S & P 500 returns
# (thanks to Chen Haibo for providing this example)

data(sp500)

# fit some models
msp <- depmix(logret~1,nstates=2,data=sp500)
set.seed(1)
fmsp <- fit(msp)	

# plot posterior state sequence for the 2-state model
plot(ts(posterior(fmsp)[,2], start=c(1950,2),deltat=1/12),ylab="probability",
main="Posterior probability of state 1 (volatile, negative markets).",
frame=FALSE)

\dontrun{

# this creates data with a single change point with Poisson data
set.seed(3)
y1 <- rpois(50,1)
y2 <- rpois(50,2)
ydf <- data.frame(y=c(y1,y2))

# fit models with 1 to 3 states
m1 <- depmix(y~1,ns=1,family=poisson(),data=ydf)
set.seed(1)
fm1 <- fit(m1)
m2 <- depmix(y~1,ns=2,family=poisson(),data=ydf)
set.seed(1)
fm2 <- fit(m2)
m3 <- depmix(y~1,ns=3,family=poisson(),data=ydf)
set.seed(1)
fm3 <- fit(m3,em=em.control(maxit=500))

# plot the BICs to select the proper model
plot(1:3,c(BIC(fm1),BIC(fm2),BIC(fm3)),ty="b")

}

\dontrun{
# similar to the binomial model, data may also be entered in 
# multi-column format where the n for each row can be different
dt <- data.frame(y1=c(0,1,1,2,4,5),y2=c(1,0,1,0,1,0),y3=c(4,4,3,2,1,1))
# specify a mixture model ...
m2 <- mix(cbind(y1,y2,y3)~1,data=dt,ns=2,family=multinomial("identity"))
set.seed(1)
fm2 <- fit(m2)
# ... or dependent mixture model
dm2 <- depmix(cbind(y1,y2,y3)~1,data=dt,ns=2,family=multinomial("identity"))
set.seed(1)
fdm2 <- fit(dm2)
}



}

\keyword{methods}


