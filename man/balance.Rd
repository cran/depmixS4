\name{balance}

\docType{data}

\alias{balance}

\title{Balance Scale Data}

\description{Balance scale data of four distance items from 779 participants;
participants ages are included.}

\usage{data(balance)}

\format{
  
	A data frame with 779 observations on the following variables. The full 
	dataset is described and analyzed extensively in 
	Jansen & Van der Maas (2002). The trichotomous data are left, balance, right. 
	The dichotomous version of the data is scored correct, incorrect. 
  
  	\describe{
		\item{\code{sex}}{Participants sex.}
		\item{\code{agedays}}{Age in days.}
		\item{\code{age}}{Age in years.}
		\item{\code{t1}}{Trichotomously scored distance item.}
		\item{\code{t2}}{Trichotomously scored distance item.}
		\item{\code{t3}}{Trichotomously scored distance item.}
		\item{\code{t4}}{Trichotomously scored distance item.}
		\item{\code{d1}}{Dichotomously scored distance item.}
		\item{\code{d2}}{Dichotomously scored distance item.}
		\item{\code{d3}}{Dichotomously scored distance item.}
		\item{\code{d4}}{Dichotomously scored distance item.}
	}
	
}

\source{

	Brenda Jansen & Han van der Maas (2002). The development of
	children's rule use on the balance scale task. \emph{Journal 
	of experimental Child Psychology, 81}, p. 383-416.

}

\examples{
data(balance)
}

\keyword{datasets}
