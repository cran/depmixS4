\name{speed}

\docType{data}

\alias{speed}

\title{Speed Accuracy Switching Data}

\description{

	This data set is a bivariate series of reaction times and accuracy
	scores of a single subject switching between slow and accurate
	responding and fast guessing on a lexical decision task.  The slow and
	accurate responding, and the fast guessing can be modelled using two
	states, with a switching regime between them.  The dataset further
	contains a third variable called Pacc, representing the relative
	pay-off for accurate responding, which is on a scale of zero to one.
	The value of Pacc was varied during the experiment to induce the
	switching.  This data set is a from a single subject from experiment 2
	in \cite{Van der Maas et al, 2005}.
	
}

\usage{data(speed)}

\format{
  A data frame with 439 observations on the following 3 variables.
  \describe{
	\item{\code{rt}}{a numeric vector}
	\item{\code{corr}}{a numeric vector}
	\item{\code{Pacc}}{a numeric vector}
  }
}

\source{

	Han L. J. Van der Maas, Conor V. Dolan and Peter C. M. Molenaar (2007),
	Phase Transitions in the Trade-Off between Speed and Accuracy in Choice
	Reaction Time Tasks.  \emph{Manuscript in preparation}. 

}

\examples{
data(speed)
## maybe str(speed) ; plot(speed) ...
}

\keyword{datasets}
