\name{fit}

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

	\S4method{fit}{mix}(object, fixed=NULL, equal=NULL, 
		conrows=NULL, conrows.upper=NULL, conrows.lower=NULL, 
		method=NULL, verbose=FALSE,
		emcontrol=em.control(),
		solnpcntrl=list(rho = 1, outer.iter = 400, inner.iter = 800, 
		delta = 1e-7, tol = 1e-8),
		donlpcntrl=donlp2Control(),
		...)
	
	\S4method{summary}{mix.fitted}(object,which="all")
	
	\S4method{summary}{depmix.fitted}(object,which="all")

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
	
	\item{verbose}{Should optimization information be displayed on screen?}
		
	\item{emcontrol}{Named list with control parameters for the EM
		algorithm (see \code{\link{em.control}}).}
  
	\item{solnpcntrl}{Control parameters passed to the 'rsolnp' optimizer;
	  see \code{\link{solnp}} for explanation and defaults used there.}
	
	\item{donlpcntrl}{Control parameters passed to 'donlp' optimizer; see
	\code{?donlp2Control} for explanation and defaults used there; this can 
		be used to tweak optimization but note that extra output is not 
		returned.}
	
	\item{which}{Should summaries be provided for "all" submodels? Options 
		are "prior", "response", and for fitted depmix models also "transition".}
  
	\item{...}{Further arguments passed on to the optimization methods.}

}

\details{ 

	Models are fitted by the EM algorithm if there are no constraints on the
	parameters.  Aspects of the EM algorithm can be controlled through the
	\code{emcontrol} argument; see details in \code{\link{em.control}}.
	Otherwise the general optimizers \code{solnp}, the default (from package
	\code{Rsolnp}) or \code{donlp2} (from package \code{Rdonlp2}) are used
	which handle general linear (in-)equality constraints.  These optimizers
	are selected by setting method='rsolnp' or method='donlp' respectively.
	
	Three types of constraints can be specified on the parameters: fixed,
	equality, and general linear (in-)equality constraints.  Constraint
	vectors should be of length \code{npar(object)}; note that this hence
	includes redundant parameters such as the base category parameter in
	multinomial logistic models which is always fixed at zero.  See help on
	\code{\link{getpars}} and \code{\link{setpars}} about the ordering of
	parameters.
	
	The \code{equal} argument is used to specify equality constraints:
	parameters that get the same integer number in this vector are
	estimated to be equal. Any integers can be used in this way except 0
	and 1, which indicate fixed and free parameters, respectively. 

	Using \code{solnp} (or \code{donlp2}), a Newton-Raphson scheme is employed
	to estimate parameters subject to linear constraints by imposing:
	
			bl <= A*x <= bu,
	
	where x is the parameter vector, bl is a vector of lower bounds, bu is
	a vector of upper bounds, and A is the constraint matrix.

	The \code{conrows} argument is used to specify rows of A directly, and
	the conrows.lower and conrows.upper arguments to specify the bounds on
	the constraints.  \code{conrows} must be a matrix of npar(object) columns
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

\value{
	
	\code{fit} returns an object of class
	\code{\link{depmix.fitted}} which contains the
	original \code{depmix} object, and further has slots:
	
	\describe{
		\item{\code{message}:}{Convergence information.}
	
		\item{\code{conMat}:}{The constraint matrix A, see Details.}
	
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
	
	Some of the below models for the \code{speed} data are reported in:
	
	Ingmar Visser, Maartje E. J. Raijmakers and Han L. J. van der Maas
	(2009).  Hidden Markov Models for Invdividual Time Series.  In: Jaan
	Valsiner, Peter C. M. Molenaar, M. C. D. P. Lyra, and N. Chaudhary
	(editors).  \emph{Dynamic Process Methodology in the Social and
	Developmental Sciences}, chapter 13, pages 269--289.  New York:
	Springer.
	
}

\examples{

data(speed)

# 2-state model on rt and corr from speed data set 
# with Pacc as covariate on the transition matrix
# ntimes is used to specify the lengths of 3 separate time-series
mod1 <- depmix(list(rt~1,corr~1),data=speed,transition=~Pacc,nstates=2,
	family=list(gaussian(),multinomial("identity")),ntimes=c(168,134,137))
# fit the model
set.seed(3)
fmod1 <- fit(mod1)
fmod1 # to see the logLik and optimization information
# to see the parameters
summary(fmod1)

# to obtain the posterior most likely state sequence, as computed by the 
# Viterbi algorithm
pst_global <- posterior(fmod1, type = "global")
# local decoding provides a different method for state classification:
pst_local <- posterior(fmod1,type="local")
identical(pst_global, pst_local)
# smoothing probabilities are used for local decoding, and may be used as 
# easily interpretable posterior state probabilities
pst_prob <- posterior(fmod1, type = "smoothing")

# testing viterbi states for new data
df <- data.frame(corr=c(1,0,1),rt=c(6.4,5.5,5.3),Pacc=c(0.6,0.1,0.1))
# define model with new data like above
modNew <- depmix(list(rt~1,corr~1),data=df,transition=~Pacc,nstates=2,
	family=list(gaussian(),multinomial("identity")))
# get parameters from estimated model
modNew <- setpars(modNew,getpars(fmod1))
# check the state sequence and probabilities
pst_new <- posterior(modNew, type="global")

# same model, now with missing data
\dontrun{
speed[2,1] <- NA
speed[3,2] <- NA

# 2-state model on rt and corr from speed data set 
# with Pacc as covariate on the transition matrix
# ntimes is used to specify the lengths of 3 separate series
mod1ms <- depmix(list(rt~1,corr~1),data=speed,transition=~Pacc,nstates=2,
	family=list(gaussian(),multinomial("identity")),ntimes=c(168,134,137))
# fit the model
set.seed(3)
fmod1ms <- fit(mod1ms)
}

# instead of the normal likelihood, we can also maximise the "classification" likelihood
# this uses the maximum a posteriori state sequence to assign observations to states
# and to compute initial and transition probabilities. 

fmod1b <- fit(mod1,emcontrol=em.control(classification="hard"))
fmod1b # to see the logLik and optimization information

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
pars[14]=0.5 # the corr parameters are now both 0.5
mod2 <- setpars(mod1,pars)

# fix the parameters by setting: 
free <- c(0,0,rep(c(0,1),4),1,1,0,0,1,1,1,1)
# fit the model
fmod2 <- fit(mod2,fixed=!free)

# likelihood ratio insignificant, hence fmod2 better than fmod1
llratio(fmod1,fmod2)


# ADDING SOME GENERAL LINEAR CONSTRAINTS

# set the starting values of this model to the optimized 
# values of the previously fitted model to speed optimization

\dontrun{

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

}

data(balance)
# four binary items on the balance scale task
mod4 <- mix(list(d1~1,d2~1,d3~1,d4~1), data=balance, nstates=2,
	family=list(multinomial("identity"),multinomial("identity"),
	multinomial("identity"),multinomial("identity")))

set.seed(1)
fmod4 <- fit(mod4)

\dontrun{

# add age as covariate on class membership by using the prior argument
mod5 <- mix(list(d1~1,d2~1,d3~1,d4~1), data=balance, nstates=2,
	family=list(multinomial("identity"),multinomial("identity"),
	multinomial("identity"),multinomial("identity")),
	prior=~age, initdata=balance)

set.seed(1)
fmod5 <- fit(mod5)

# check the likelihood ratio; adding age significantly improves the goodness-of-fit
llratio(fmod5,fmod4)

}

}

\author{Ingmar Visser & Maarten Speekenbrink}

\keyword{methods}
