\name{posterior}

\alias{posterior}

\alias{posterior,depmix.fitted-method}
\alias{posterior,mix.fitted-method}


\title{ Posterior states/classes }

\description{

	Return the posterior states for a fitted \code{(dep-)mix} object.  In
	the case of a latent class or mixture model these are the class
	probabilities.

}

\usage{
	\S4method{posterior}{depmix.fitted}(object)
	\S4method{posterior}{mix.fitted}(object)
}

\arguments{
  \item{object}{A (dep-)mix object.}
}

\value{

	\item{posterior}{: Returns a data.frame with nstates(object) + 1
	columns; the first column has the viterbi states, the other columns
	have the delta probabilities, see Rabiner (1989).}

}

\references{ 	

	Lawrence R. Rabiner (1989).  A tutorial on hidden Markov models and
	selected applications in speech recognition.  \emph{Proceedings of
	IEEE}, 77-2, p.  267-295. 

}

\author{ Ingmar Visser }

\details{
	See the \code{\link{fit}} help page for an example. 
}

\keyword{ methods }

