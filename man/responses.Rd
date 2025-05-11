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

	\code{MVNresponse} is a model for a multivariate normal distributed
	response.  See \code{\link{makeDepmix}} for an example of how to use
	this and other non-glm like distributions.

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
	
# binomial response model
x <- rnorm(1000)
p <- plogis(x)
ss <- rbinom(1000,1,p)
mod <- GLMresponse(cbind(ss,1-ss)~x,family=binomial())
fit(mod)
glm(cbind(ss,1-ss)~x, family=binomial)

# gamma response model
x=runif(1000,1,5)
res <- rgamma(1000,x)
## note that gamma needs proper starting values which are not
## provided by depmixS4 (even with them, this may produce warnings)
mod <- GLMresponse(res~x,family=Gamma(),pst=c(0.8,1/0.8))
fit(mod)
glm(res~x,family=Gamma)

# multinomial response model
x <- sample(0:1,1000,rep=TRUE)
mod <- GLMresponse(sample(1:3,1000,rep=TRUE)~x,family=multinomial(),pstart=c(0.33,0.33,0.33,0,0,1))
mod@y <- simulate(mod)
fit(mod)
colSums(mod@y[which(x==0),])/length(which(x==0))
colSums(mod@y[which(x==1),])/length(which(x==1))
# note that the response is treated as factor here, internal representation is in 
# dummy coded format:
head(mod@y)
# similar to the binomial model, data may also be entered in multi-column format 
# where the n for each row can be different
dt <- data.frame(y1=c(0,1,1,2,4,5),y2=c(1,0,1,0,1,0),y3=c(4,4,3,2,1,1))
m2 <- mix(cbind(y1,y2,y3)~1,data=dt,ns=2,family=multinomial("identity"))
fm2 <- fit(m2)
summary(fm2)

# multivariate normal response model
mn <- c(1,2,3)
sig <- matrix(c(1,.5,0,.5,1,0,0,0,2),3,3)
y <- mvrnorm(1000,mn,sig)
mod <- MVNresponse(y~1)
fit(mod)
colMeans(y)
var(y)

# normal (gaussian) response model
y <- rnorm(1000)
mod <- GLMresponse(y~1)
fm <- fit(mod)
cat("Test gaussian fit: ", all.equal(getpars(fm),c(mean(y),sd(y)),check.attributes=FALSE))

# poisson response model
x <- abs(rnorm(1000,2))
res <- rpois(1000,x)
mod <- GLMresponse(res~x,family=poisson())
fit(mod)
glm(res~x, family=poisson)

# this creates data with a single change point with Poisson distributed data
set.seed(3)
y1 <- rpois(50,1)
y2 <- rpois(50,2)
ydf <- data.frame(y=c(y1,y2))

# fit models with 1 to 3 states
m1 <- depmix(y~1,ns=1,family=poisson(),data=ydf)
fm1 <- fit(m1)
m2 <- depmix(y~1,ns=2,family=poisson(),data=ydf)
fm2 <- fit(m2)
m3 <- depmix(y~1,ns=3,family=poisson(),data=ydf)
fm3 <- fit(m3,em=em.control(maxit=500))

# plot the BICs to select the proper model
plot(1:3,c(BIC(fm1),BIC(fm2),BIC(fm3)),ty="b")

}
	
\author{Maarten Speekenbrink & Ingmar Visser}

\keyword{models}
