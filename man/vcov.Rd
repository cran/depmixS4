\name{vcov}

\alias{vcov}
\alias{confint}
\alias{standardError}
\alias{hessian}

\alias{vcov,mix-method}
\alias{confint,mix-method}
\alias{standardError,mix-method}
\alias{hessian,mix-method}

\title{ Parameter standard errors }

\description{

	These functions provide standard errors for parameters of (dep-)mix models. 

}

\usage{
	\S4method{vcov}{mix}(object, fixed=NULL, equal=NULL, 
		conrows=NULL, conrows.upper=NULL, conrows.lower=NULL, tolerance=1e-6, 
		method="finiteDifferences", ...)	
	
	\S4method{standardError}{mix}(object, fixed=NULL, equal=NULL, 
		conrows=NULL, conrows.upper=NULL, conrows.lower=NULL, tolerance=1e-6, 
		method="finiteDifferences", ...)	
	
	\S4method{confint}{mix}(object, level=0.95, fixed=NULL, equal=NULL, 
		conrows=NULL, conrows.upper=NULL, conrows.lower=NULL, tolerance=1e-6, 
		method="finiteDifferences", ...)	

	\S4method{hessian}{mix}(object, tolerance=1e-6, 	
			method="finiteDifferences", ...)
}


\arguments{

  \item{object}{A (dep-)mix object; see \link{depmix} for details.}
  
  \item{fixed, equal}{These arguments are used to specify constraints on a model; see usage details here: \code{\link{fit}}.}
  
  \item{conrows}{These arguments are used to specify constraints on a model; see usage details here: \code{\link{fit}}.}
  
  \item{conrows.upper}{These arguments are used to specify constraints on a model; see usage details here: \code{\link{fit}}.}
  
  \item{conrows.lower}{These arguments are used to specify constraints on a model; see usage details here: \code{\link{fit}}.}
  
  \item{tolerance}{Threshold used for testing whether parameters are estimated on the boundary of the parameter space; if so, 
  they are ignored in these functions.}
  
  \item{method}{The method used for computing the Hessian matrix of the parameters; currently only a finite 
  differences method (using \code{fdHess} from package \code{\link[nlme]{nlme}}) is implemented and hence used by default.}
  
  \item{level}{The desired significance level for the confidence intervals.}
  
  \item{...}{Further arguments passed to other methods; currently not in use.}
  
}


\details{ 
	
	\code{vcov} computes the variance-covariance matrix of a (dep-)mix object, either fitted or not. 
	It does so by first constructing a Hessian matrix through the use of \code{hessian} and then 
	transforming this as described in Visser et al (2000), taking into account the linear constraints 
	that are part of the model. Currently, \code{hessian} has a single \code{method} using finite 
	differences to arrive at an approximation of the second order derivative matrix of the parameters. 
	
	\code{confint} and \code{standardError} use \code{vcov} to compute confidence intervals (the confidence
	level can be set through an argument) and standard errors respectively. The latter are computed first by 
	using \code{sqrt(diag(vcov))} and the confidence intervals are computed through the normal approximation. 
	
	If and when these methods are applied to \code{fit}'ted models, the linear constraint matrix is 
	obtained from the \code{mix.fitted} or \code{depmix.fitted} slot \code{lincon} (supplemented with 
	additional constraints if those are provided through the \code{equal} and other arguments to these
	functions). 
	
	All four functions exclude parameters that are estimated on or near (this can be controlled using 
	the \code{tolerance} argument) their boundary values. Setting this argument to zero can result in
	error as the \code{fdHess} function requires an environment around the parameter estimate that 
	provides proper log-likelihood values, which parameter on or over their boundary values are not
	guaranteed to provided. Fixed parameters are similarly ignored in these four functions. 
	
}

\value{
	
	\code{vcov} returns a named list with elements \code{vcov}, \code{elements}, and \code{lincon}. 
	\code{standardError} returns a \code{data.frame} with columns \code{par}, \code{elements}, 
	and \code{se}. \code{confint} returns a \code{data.frame} with columns \code{par}, 
	\code{elements}, and two columns for the lower and upper bounds of the confidence intervals 
	(with the column names indicating the \code{level} of the interval.)
	
	\item{vcov}{: The variance-covariance matrix of the parameters.}
	
	\item{elements}{: Vector of length \code{npar(object)} indicating which elements of 
	the parameter vector are included in computing the hessian, the variance-covariance matrix, 
	the standard errors and/or the confidence intervals.}
	
	
		\item{inc}{: 'inc'luded parameter.}
		\item{fix}{: 'fix'ed parameter.}
		\item{bnd}{: parameter estimated on the boundary.}
	
	\item{par}{: The values of the parameters.}
	
	\item{se}{: The values of the standard errors of the parameters.}
	
	\item{lower/upper}{: The lower and upper bounds of the confidence intervals; column names
	indicate the as in 0.5+/-level/2, using the \code{level} argument.}
	
}

\note{
	
	Note that the quality of the resulting standard errors is similar to those reported in 
	Visser et al (2000) for both bootstrap and the profile likelihood methods. In Visser et
	al (2000), the finite differences standard errors were somewhat less precise as they 
	relied on a very parsimonious but indeed less precise method for computing the finite
	differences approximation (computation time was a much scarcer resource at the time then
	it is now). 
	
}

\references{ 	

	Ingmar Visser, Maartje E. J. Raijmakers, and Peter C. M. Molenaar (2000).
	Confidence intervals for hidden Markov model parameters. \emph{British journal
	of mathematical and statistical psychology}, 53, p. 317-327. 

}


\examples{

	data(speed)

	# 2-state model on rt and corr from speed data set 
	# with Pacc as covariate on the transition matrix
	# ntimes is used to specify the lengths of 3 separate series
	mod1 <- depmix(list(rt~1,corr~1),data=speed,transition=~Pacc,nstates=2,
		family=list(gaussian(),multinomial("identity")),ntimes=c(168,134,137))
	
	# fit the model
	set.seed(3)
	fmod1 <- fit(mod1)
	
	vcov(fmod1)$vcov # $
	standardError(fmod1)
	confint(fmod1)
	
}

\author{ Ingmar Visser }

\keyword{ methods }

