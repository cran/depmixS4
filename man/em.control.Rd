\name{em.control}

\alias{em.control}

\title{Control parameters for the EM algorithm}

\description{Set control parameters for the EM algorithm.}

\usage{
	
	em.control(maxit = 500, tol = 1e-08, crit = c("relative","absolute"), 
	random.start = TRUE, classification = c("soft","hard"))
	
}

\arguments{
	
	\item{maxit}{The maximum number of iterations.}
		
	\item{tol}{The tolerance level for convergence. See Details.}
	
	\item{crit}{Sets the convergence criterion to "relative" or "absolute" 
	change of the log-likelihood. See Details.}
  
	\item{random.start}{This is used for a (limited) random
	initialization of the parameters. See Details.}
	
	\item{classification}{Type of classification to states
	used. See Details.}
	
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
tolerance parameter \code{tol} scale with the value of the log-likelihood
(and these are not changed automagically). 

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

Argument \code{classification} is used to choose either soft (default) or
hard classification of observations to states. When using soft classification, observations
are assigned to states with a weight equal to the posterior probability of
the state. When using hard classification, observations are assigned to states
according to the maximum a posteriori (MAP) states (i.e., each observation
is assigned to one state, which is determined by the Viterbi algorithm in the
case of \code{depmix} models). As a result, the EM algorithm will find a local
maximum of the classification likelihood (Celeux & Govaert, 1992). 
Warning: hard classification is an experimental feature, 
especially for hidden Markov models, and its use is currently not advised.
}

\references{
	Lawrence R. Rabiner (1989).  A tutorial on hidden Markov models and
	selected applications in speech recognition.  \emph{Proceedings of
	IEEE}, 77-2, p.  267-295.
  
  Gilles Celeux and Gerard Govaert (1992). A classification EM algorithm 
  for clustering and two stochastic versions. \emph{Computational
Statistics and Data Analysis, 14}, p. 315-332.
}

\value{
	
	\code{em.control} returns a list of the control parameters. 
	
}

\author{Ingmar Visser & Maarten Speekenbrink}

\keyword{methods}

\examples{
# using "hard" assignment of observations to the states, we can maximise the
# classification likelihood instead of the usual marginal likelihood
data(speed)  
mod <- depmix(list(rt~1,corr~1),data=speed,nstates=2,
    family=list(gaussian(),multinomial("identity")),ntimes=c(168,134,137))
set.seed(1)
# fit the model by calling fit
fmod <- fit(mod,emcontrol=em.control(classification="hard"))
# can get rather different solutions with different starting values...
set.seed(3)
fmod2 <- fit(mod,emcontrol=em.control(classification="hard"))
}