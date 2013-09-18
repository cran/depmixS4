\name{sp500}

\docType{data}

\alias{sp500}

\title{Standard & Poor's 500 index}

\description{

	This data set consists of (monthly) values of the S&P 500 stock exchange
	index. The variable of interest is the logarithm of the return values, i.e.,
	the logarithm of the ratio of indices, in this case the closing index is
	used. 
	
}

\usage{data(speed)}

\format{
  A data frame with 744 observations and 6 variables.
  \describe{
	\item{\code{Open}}{Index at the start of trading.}
	\item{\code{High}}{Highest index.}
	\item{\code{Low}}{Lowest index.}
	\item{\code{Close}}{Index at the close of trading.}
	\item{\code{Volume}}{The volume of trading.}
	\item{\code{logret}}{The log return of the closing index.}
  }
}

\source{
	
	Yahoo Data. 
}

\examples{

data(sp500)

# the data can be made with the following code (eg to include a longer or
# shorter time span)

\dontrun{

require(TTR)

# load SP500 returns
Sys.setenv(tz='UTC')

sp500 <- getYahooData('^GSPC',start=19500101,end=20120228,freq='daily')
ep <- endpoints(sp500, on="months", k=1)
sp500 <- sp500[ep[2:(length(ep)-1)]]
sp500$sp500_ret <- log(sp500$Close) - lag(log(sp500$Close))
sp500 <- na.exclude(sp500)

}

}

\keyword{datasets}
