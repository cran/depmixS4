\name{depmixS4-package}

\alias{depmixS4}
\alias{depmixS4-package}

\docType{package}

\title{

depmixS4 provides classes for specifying and fitting hidden Markov models

}

\description{

	\code{depmixS4} is a framework for specifying and fitting dependent
	mixture models, otherwise known as hidden or latent Markov models.
	Optimization is done with the EM algorithm or optionally with Rdonlp2
	when (general linear (in-)equality) constraints on the parameters need
	to be incorporated.  Models can be fitted on (multiple) sets of
	observations.  The response densities for each state may be chosen from
	the GLM family, or a multinomial.  User defined response densities are
	easy to add; for the latter an example is given for the ex-gauss distribution
	as well as the multivariate normal distribution. 
	
	Mixture or latent class (regression) models can also be fitted; these
	are the limit case in which the length of observed time series is 1 for
	all cases. 

}

\details{

	\tabular{ll}{
		Package: \tab depmixS4\cr
		Type: \tab Package\cr
		Version: \tab 0.3-0\cr
		Date: \tab 2010-01-30\cr
		License: \tab GPL\cr
	}

	Model fitting is done in two steps; first, models are specified through
	the \code{\link{depmix}} function (or the \code{\link{mix}} function for
	mixture and latent class models), which both use standard
	\code{\link{glm}} style arguments to specify the observed
	distributions; second, the model needs to be fitted by using the
	\code{\link{fit}} function; imposing constraints is done through the
	fit function.  Standard output includes the optimized parameters and
	the posterior densities for the states and the optimal state sequence.
	
	For full control and the possibility to add new response distributions,
	check the \code{\link{makeDepmix}} help page. 

}

\author{

Ingmar Visser & Maarten Speekenbrink

Maintainer: i.visser@uva.nl

}

\references{

	On hidden Markov models: Lawrence R. Rabiner (1989).  A tutorial on
	hidden Markov models and selected applications in speech recognition.
	\emph{Proceedings of IEEE}, 77-2, p.  267-295.

	On latent class models: A. L. McCutcheon (1987).  \emph{Latent class
	analysis}.  Sage Publications.

}

\keyword{ package }

\seealso{

	\code{\link{depmix}}, \code{\link{fit}}

}

\examples{
	# create a 2 state model with one continuous and one binary response
	data(speed)
	mod <- depmix(list(rt~1,corr~1),data=speed,nstates=2,family=list(gaussian(),multinomial()))
	# print the model, formulae and parameter values (ie the starting values)
	mod
}