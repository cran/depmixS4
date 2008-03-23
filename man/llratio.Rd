
\name{llratio}

\docType{methods}

\alias{llratio}

\alias{loglikelihoodratio}

\alias{show,llratio-method}

\title{Log likelihood ratio test on two fitted models}

\description{Performs a log likelihood ratio test on two fitted
\code{depmix} models.}

\usage{
	
	llratio(basemodel, constrainedmodel, ...)
	
}

\arguments{ 
	
	\item{basemodel}{Fitted model with a \code{logLik} method.}
	
	\item{constrainedmodel}{Fitted model with a \code{logLik} method.}
	
	\item{...}{Not currently used.}
	
}

\value{
	
	\code{llratio} returns an object of class \code{llratio} which has slots: 
	
	\item{value}{: Minus twice the loglikelihood difference.}
	
	\item{df}{: The degrees of freedom, ie the difference in number of freely
	estimated paraemters between the models.}
	
	The print method shows the value, the degrees of freedom and the
	corresponding p-value under the chisquared distribution. 
	
}

\author{Ingmar Visser}

\keyword{methods}