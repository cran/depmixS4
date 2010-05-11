\name{fit}

\docType{method}

\alias{fit}
\alias{depmix.fit}
\alias{fit,depmix-method}

\alias{show,depmix.fitted-method}
\alias{summary,depmix.fitted-method}

\alias{mix.fit}
\alias{fit,mix-method}

\alias{show,mix.fitted-method}
\alias{summary,mix.fitted-method}

\title{Fit 'depmix' or 'mix' models}

\description{

	\code{fit} optimizes parameters of \code{\link{depmix}} or
	\code{\link{mix}} models, optionally subject to general linear
	(in)equality constraints.

}

\usage{
	
	\S4method{fit}{depmix}(object, fixed=NULL, equal=NULL, conrows=NULL,
		conrows.upper=0, conrows.lower=0, method=NULL,em.control=list(tol=1e-8,
		crit=c("relative","absolute"),random.start=FALSE),verbose=TRUE,...)
	
	\S4method{summary}{depmix.fitted}(object,which="all")

	\S4method{fit}{mix}(object, fixed=NULL, equal=NULL, conrows=NULL,
		conrows.upper=0, conrows.lower=0, method=NULL,em.control=list(tol=1e-8,
		crit=c("relative","absolute"),random.start=FALSE),verbose=TRUE,...)
	
	\S4method{summary}{mix.fitted}(object,which="all")

	
}

\arguments{

	\item{object}{An object of class \code{(dep-)mix}.}
	
	\item{fixed}{Vector of mode logical indicating which parameters should 
		be fixed.}

	\item{equal}{Vector indicating equality constraints; see Details.}

	\item{conrows}{Rows of a general linear constraint matrix; see Details.}

	\item{conrows.upper, conrows.lower}{Upper and lower bounds for the
		linear constraints; see Details.}

	\item{method}{The optimization method; mostly determined by
		constraints.}
		
	\item{em.control}{Named list with control parameters for the EM
		algorithm (see Control parameters below).}
  
	\item{verbose}{Should optimization information be displayed on screen?}
	
	\item{which}{Should summaries be provided for "all" submodels? Options 
		are "prior", "response", and for fitted depmix models also "transition".}
  
	\item{...}{Further arguments passed on to the optimization methods.}

}

\details{ 

	Models are fitted by the EM algorithm if there are no constraints on
	the parameters. Aspects of the EM algorithm can be controlled through
	the \code{control.em} argument (see below). 
	Otherwise the general optimizers \code{solnp}, the
	default (from package \code{Rsolnp}) or \code{donlp2} (from package
	\code{Rdonlp2}) are used which handle general linear (in-)equality
	constraints.
	
	Three types of constraints can be specified on the parameters: fixed,
	equality, and general linear (in-)equality constraints.  Constraint
	vectors should be of length npar(object).  See help on
	\code{\link{getpars}} and \code{\link{setpars}} about the ordering of
	parameters. 
	
	The \code{equal} argument is used to specify equality constraints:
	parameters that get the same integer number in this vector are
	estimated to be equal. Any integers can be used in this way except 0
	and 1, which indicate fixed and free parameters, respectively. 

	Using \code{solnp} or \code{donlp2} , a Newton-Raphson scheme is employed
	to estimate parameters subject to linear constraints by imposing:
	
			bl <= A*x <= bu,
	
	where x is the parameter vector, bl is a vector of lower bounds, bu is
	a vector of upper bounds, and A is the constraint matrix.

	The \code{conrows} argument is used to specify rows of A directly, and
	the conrows.lower and conrows.upper arguments to specify the bounds on
	the constraints.  \code{conrows} is a matrix of npar(object) columns
	and one row for each constraint (a vector in the case of a single
	constraint).  Examples of these three ways of constraining parameters
	are provided below.
	
	Note that when specifying constraints that these should respect the
	fixed constraints inherent in e.g. the multinomial logit models for the
	initial and transition probabilities.  For example, the baseline
	category coefficient in a multinomial logit model is fixed on zero.
	
	\code{\link{llratio}} performs a log-likelihood ratio test on two
	\code{fit}'ted models; the first object should have the largest degrees
	of freedom (find out by using \code{\link{freepars}}).
	
}

\section{Control parameters for the EM algorithm}{
	
	Aspects of EM algorithm can be controlled by the \code{em.control}
	argument.  This named list currently observes the following parameters:
	
	\describe{
		\item{\code{tol}:}{sets the the tolerance level for convergence.}
	
		\item{\code{criterion}:}{sets the convergence criterion to either
		the relative change in the log-likelihood or the absolute change in
		the log-likelihood.  The relative likelihood criterion (the
		default) assumes convergence on iteration \eqn{i}{i} when
		\eqn{\frac{\log L(i) - \log L(i-1)}{\log L(i-1)} <
		tol}{\frac{\log L(i) - \log L(i-1)}{\log L(i-1)} <
		tol}.  The absolute likelihood criterion assumes convergence
		on iteration \eqn{i}{i} when \eqn{\log L(i) - \log L(i-1) <
		tol}{(logLik(i) - logLik(i-1)) < tol}.}
  
		\item{\code{random.start}:}{is used for a (limited) random
		initialization of the parameters.  In particular, if
		\code{random.start=TRUE}, the (posterior) state probabilities are
		randomized at iteration 0 (using a uniform distribution).  Random
		initialization is useful when no initial parameters can be given to
		distinguish between the states.  It is also useful for repeated
		estimation from different starting values.
		} 
	} 
}

\value{
	
	\code{fit} returns an object of class
	\code{\link{depmix.fitted}} which contains the
	original \code{depmix} object, and further has slots:
	
	\describe{
		\item{\code{message}:}{Convergence information.}
	
		\item{\code{conMa}t:}{The constraint matrix A, see Details.}
	
		\item{\code{posterior}:}{The posterior state sequence (computed
		with the viterbi algorithm), and the posterior probabilities (delta
		probabilities in Rabiner, 1989, notation).}
		
	}
	
	The print method shows the \code{message} along with the likelihood and
	AIC and BIC; the summary method prints the parameter estimates.

	Posterior densities and the viterbi state sequence can be accessed
	through \code{\link{posterior}}.
	
	As fitted models are depmixS4 models, they can be used as starting
	values for new fits, for example with constraints added. Note that 
	when refitting already fitted models, the constraints, if any, are not
	added automatically, they have to be added again. 

}

\references{
	
	Lawrence R. Rabiner (1989).  A tutorial on hidden Markov models and
	selected applications in speech recognition.  \emph{Proceedings of
	IEEE}, 77-2, p.  267-295.
	
}

\examples{

data(speed)

# 2-state model on rt and corr from speed data set 
# with Pacc as covariate on the transition matrix
# starting values for the transition pars (without 
# those EM does not get off the ground)
set.seed(1)
tr=runif(6)
trst=c(tr[1:2],0,tr[3:5],0,tr[6])
mod1 <- depmix(list(rt~1,corr~1),data=speed,transition=~Pacc,nstates=2,
	family=list(gaussian(),multinomial("identity")),trstart=trst)
# fit the model
fmod1 <- fit(mod1)
fmod1 # to see the logLik and optimization information
# to see the parameters
summary(fmod1)

# FIX SOME PARAMETERS

# get the starting values of this model to the optimized 
# values of the previously fitted model to speed optimization

pars <- c(unlist(getpars(fmod1)))

# constrain the initial state probs to be 0 and 1 
# also constrain the guessing probs to be 0.5 and 0.5 
# (ie the probabilities of corr in state 1)
# change the ones that we want to constrain
pars[1]=0
pars[2]=1 # this means the process will always start in state 2
pars[13]=0.5
pars[14]=0.5 # the corr parameters in state 1 are now both 0, corresponding the 0.5 prob
mod2 <- setpars(mod1,pars)

# fix the parameters by setting: 
free <- c(0,0,rep(c(0,1),4),1,1,0,0,1,1,1,1)
# fit the model
fmod2 <- fit(mod2,fixed=!free)

# likelihood ratio insignificant, hence fmod2 better than fmod1
llratio(fmod1,fmod2)

# NOW ADD SOME GENERAL LINEAR CONSTRAINTS

# set the starting values of this model to the optimized 
# values of the previously fitted model to speed optimization

pars <- c(unlist(getpars(fmod2)))
pars[4] <- pars[8] <- -4
pars[6] <- pars[10] <- 10
mod3 <- setpars(mod2,pars)

# start with fixed and free parameters
conpat <- c(0,0,rep(c(0,1),4),1,1,0,0,1,1,1,1)
# constrain the beta's on the transition parameters to be equal
conpat[4] <- conpat[8] <- 2
conpat[6] <- conpat[10] <- 3

fmod3 <- fit(mod3,equal=conpat)

llratio(fmod2,fmod3)

# above constraints can also be specified using the conrows argument as follows
conr <- matrix(0,2,18)
# parameters 4 and 8 have to be equal, otherwise stated, their diffence should be zero,
# and similarly for parameters 6 & 10
conr[1,4] <- 1
conr[1,8] <- -1
conr[2,6] <- 1
conr[2,10] <- -1

# note here that we use the fitted model fmod2 as that has appropriate 
# starting values
fmod3b <- fit(mod3,conrows=conr,fixed=!free) # using free defined above

data(balance)
# four binary items on the balance scale task

instart=c(0.5,0.5)
set.seed(1)
respstart=runif(16)
# note that ntimes argument is used to make this a mixture model
mod4 <- mix(list(d1~1,d2~1,d3~1,d4~1), data=balance, nstates=2,
	family=list(multinomial(),multinomial(),multinomial(),multinomial()),
	respstart=respstart,instart=instart)

fmod4 <- fit(mod4)

# add age as covariate on class membership by using the prior argument
instart=c(0.5,0.5,0,0) # we need the initial probs and the coefficients of age 
set.seed(2)
respstart=runif(16)
mod5 <- mix(list(d1~1,d2~1,d3~1,d4~1), data=balance, nstates=2,
	family=list(multinomial(),multinomial(),multinomial(),multinomial()),
	instart=instart, respstart=respstart, prior=~age, initdata=balance)

fmod5 <- fit(mod5)

# check the likelihood ratio; adding age significantly improves the goodness-of-fit
llratio(fmod5,fmod4)

}

\author{Ingmar Visser & Maarten Speekenbrink}

\keyword{methods}
