\name{em.control}

\alias{em.control}

\title{Control parameters for the EM algorithm}

\description{Set control parameters for the EM algorithm.}

\usage{
	
	em.control(maxit = 500, tol = 1e-08, crit = "relative", random.start = TRUE)
	
}

\arguments{
	
	\item{maxit}{The maximum number of iterations.}
		
	\item{tol}{The tolerance level for convergence. See Details.}
	
	\item{crit}{Sets the convergence criterion to "relative" or "absolute" 
	change of the log-likelihood. See Details.}
  
	\item{random.start}{This is used for a (limited) random
	initialization of the parameters. See Details.}
	
}

\details{

The argument \code{crit} sets the convergence criterion to either the
relative change in the log-likelihood or the absolute change in the
log-likelihood.  The relative likelihood criterion (the default) assumes
convergence on iteration \eqn{i}{i} when 
\eqn{\frac{\log L(i) - \log L(i-1)}{\log L(i-1)} < tol}{ (log L(i) - log L(i-1))/(log L(i-1)) < tol}.  
The absolute likelihood criterion assumes convergence on iteration
\eqn{i}{i} when \eqn{\log L(i) - \log L(i-1) < tol}{(log L(i) - log L(i-1)) < tol}.  
Use \code{crit="absolute"} to invoke the latter
convergence criterion.  Note that in that case, optimal values of the 
tolerance parameter \code{tol} scale with the value of the log-likelihood (and these are not changed automagically). 

Argument \code{random.start} This is used for a (limited) random
initialization of the parameters.  In particular, if
\code{random.start=TRUE}, the (posterior) state probabilities are
randomized at iteration 0 (using a uniform distribution), i.e. the 
\eqn{\gamma} variables (Rabiner, 1989) are sampled from the Dirichlet
distribution with a (currently fixed) value of
\eqn{\alpha=0.1}; this results in values for each row of \eqn{\gamma}
that are quite close to zero and one; note that when these values are
chosen at zero and one, the initialization is similar to that used in
\code{kmeans}.  Random initialization is useful when no initial parameters can be
given to distinguish between the states.  It is also useful for repeated
estimation from different starting values.

}

\references{
	Lawrence R. Rabiner (1989).  A tutorial on hidden Markov models and
	selected applications in speech recognition.  \emph{Proceedings of
	IEEE}, 77-2, p.  267-295.
}

\value{
	
	\code{em.control} returns a list of the control parameters. 
	
}

\author{Ingmar Visser & Maarten Speekenbrink}

\keyword{methods}
