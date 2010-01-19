\name{speed}

\docType{data}

\alias{speed}

\title{Speed Accuracy Switching Data}

\description{

	This data set is a bivariate series of response times and accuracy
	scores of a single participant switching between slow/accurate
	responding and fast guessing on a lexical decision task.  The slow and
	accurate responding, and the fast guessing can be modelled using two
	states, with a switching regime between them.  The dataset further
	contains a third variable called Pacc, representing the relative
	pay-off for accurate responding, which is on a scale of zero to one.
	The value of Pacc was varied during the experiment to induce the
	switching.  This data set is a from participant A in experiment 
	1a from Dutilh et al (2010). 
	
}

\usage{data(speed)}

\format{
  A data frame with 439 observations on the following 3 variables.
  \describe{
	\item{\code{rt}}{a numeric vector of response times}
	\item{\code{corr}}{a numeric vector of accuracy scores (0/1)}
	\item{\code{Pacc}}{a numeric vector of the pay-off for accuracy}
  }
}

\source{
	
	Gilles Dutilh, Ingmar Visser, Eric-Jan Wagenmakers, & Han L. J. van der Maas.
	A Phase Transition Model for the Speed-Accuracy Trade-Off in 
	Response Time Experiments. \emph{2010, submitted for publication.}

	Please contact g.dutilh@uva.nl for preprints and more information. 
	
}

\examples{
data(speed)
## maybe str(speed) ; plot(speed) ...
}

\keyword{datasets}
