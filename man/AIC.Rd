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

\references{}

\author{Ingmar Visser}

\keyword{htest}
