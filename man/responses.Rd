\name{responses}

\alias{responses}

\alias{BINOMresponse}
\alias{GAMMAresponse}
\alias{MULTINOMresponse}
\alias{MVNresponse}
\alias{NORMresponse}
\alias{POISSONresponse}

\alias{show,MVNresponse-method}

\title{Response models currently implemented in depmix.}

\description{Depmix contains a number of default response models. We provide a
brief description of these here.}

\section{BINOMresponse}{
	
	\code{BINOMresponse} is a binomial response model. It derives from the basic
	\code{\link{GLMresponse}} class.

  \describe{
    \item{y:}{The dependent variable can be either a
  binary vector, a factor, or a 2-column matrix, with successes and misses.}
    \item{x:}{The design matrix.}
    \item{Parameters:}{A named list with a single element ``coefficients'',
    which contains the GLM coefficients.}
  }
	
}

\section{GAMMAresponse}{

	\code{GAMMAresponse} is a model for a Gamma distributed response.
  It extends the basic \code{\link{GLMresponse}} class directly.

  \describe{
    \item{y:}{The dependent variable.}
    \item{x:}{The design matrix.}
    \item{Parameters:}{A named list with a single element ``coefficients'',
    which contains the GLM coefficients.}
  }

}

\section{MULTINOMresponse}{

	\code{MULTINOMresponse} is a model for a multinomial distributed response.
  It extends the basic \code{\link{GLMresponse}} class, although the
  functionality is somewhat different from other models that do so.

  \describe{
    \item{y:}{The dependent variable. This is a binary matrix with N rows and
    Y columns, where Y is the total number of categories.}
    \item{x:}{The design matrix.}
    \item{Parameters:}{A named list with a single element ``coefficients'',
    which is an \code{ncol(x)} by \code{ncol(y)} matrix which contains the GLM
    coefficients.}
  }

}

\section{MVNresponse}{

	\code{MVNresponse} is a model for a multivariate normal distributed response.

  \describe{
    \item{y:}{The dependent variable. This is a matrix.}
    \item{x:}{The design matrix.}
    \item{Parameters:}{A named list with a elements ``coefficients'',
    which contains the GLM coefficients, and ``Sigma'', which contains the
    covariance matrix.}
  }

}

\section{NORMresponse}{

	\code{NORMresponse} is a model for a normal (Gaussian) distributed response.
  It extends the basic \code{\link{GLMresponse}} class directly.

  \describe{
    \item{y:}{The dependent variable.}
    \item{x:}{The design matrix.}
    \item{Parameters:}{A named list with elements ``coefficients'',
    which contains the GLM coefficients, and ``sd'', which contains the
    standard deviation.}
  }

}

\section{POISSONresponse}{

	\code{POISSONresponse} is a model for a Poisson distributed response.
  It extends the basic \code{\link{GLMresponse}} class directly.

  \describe{
    \item{y:}{The dependent variable.}
    \item{x:}{The design matrix.}
    \item{Parameters:}{A named list with a single element ``coefficients'',
    which contains the GLM coefficients.}
  }

}

\examples{
	
	mod <- GLMresponse(rnorm(1000)~1)

	fit(mod)
	
	mod <- GLMresponse(sample(1:3,1000,rep=TRUE)~1,family=multinomial())
	
	fit(mod)
	
	colSums(mod@y)/1000
	
	x <- sample(0:1,1000,rep=TRUE)
	
	mod <- GLMresponse(sample(1:3,1000,rep=TRUE)~x,family=multinomial(),pstart=c(0.33,0.33,0.33,0,0,1))
	
	mod@y <- simulate(mod)
	
	fit(mod)
	
	colSums(mod@y[which(x==0),])/length(which(x==0))
	
	colSums(mod@y[which(x==1),])/length(which(x==1))

	x <- rnorm(1000)
	
	library(boot)
	
	p <- inv.logit(x)
	
	ss <- rbinom(1000,1,p)
	
	mod <- GLMresponse(cbind(ss,1-ss)~x,family=binomial())
	
	fit(mod)
	
	glm(cbind(ss,1-ss)~x, family=binomial)
	
	x <- abs(rnorm(1000,2))
	
	res <- rpois(1000,x)
	
	mod <- GLMresponse(res~x,family=poisson())
	
	fit(mod)
	
	glm(res~x, family=poisson)
	
	x=runif(1000,1,5)
	
	res <- rgamma(1000,x)
	
	## note that gamma needs proper starting values which are not
	## provided by depmixS4 (even with them, this may produce warnings)
	mod <- GLMresponse(res~x,family=Gamma(),pst=c(0.8,1/0.8))
	
	fit(mod)
	
	glm(res~x,family=Gamma)
	
	mn <- c(1,2,3)
	
	sig <- matrix(c(1,.5,0,.5,1,0,0,0,2),3,3)
	
	y <- mvrnorm(1000,mn,sig)
	
	mod <- MVNresponse(y~1)
	
	fit(mod)

	colMeans(y)
	
	var(y)
	
}

\author{Maarten Speekenbrink & Ingmar Visser}

\keyword{models}
