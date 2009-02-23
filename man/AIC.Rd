\name{AIC}

\alias{AIC}
\alias{AIC,depmix-method}
\alias{AIC,mix-method}

\alias{BIC}
\alias{BIC,depmix-method}
\alias{BIC,mix-method}

\title{Compute AIC and BIC for (dep-)mix objects}

\description{
	Compute AIC and BIC for \code{depmix} and \code{mix} objects.
}

\usage{
	AIC(object, ..., k = 2)
	BIC(object, ...)
}

\arguments{
  \item{object}{A \code{depmix} or \code{mix} model.}
  \item{\dots}{Not used currently.}
  \item{k}{The penalty factor which defaults to 2.}
}

\value{
	The value of the AIC or BIC respectively.
}

\details{
	The \code{n} that is used in the BIC is \code{sum(ntimes(object))}, where
	object is the model under consideration. 
}

\references{
	H. Akaike. Information theory and an extension of the maximum likelihood 
	principle. In B. N. Petrov and F. Csaki, editors, Second International
	Symposium on Information Theory, p. 267-281. Budapest: Academiai Kiado, 1973.

	G. Schwarz. Estimating the dimension of a model. Annals of Statistics, 
	volume 6, p. 4761-464, 1978. 
}

\author{Ingmar Visser}

\keyword{htest}
