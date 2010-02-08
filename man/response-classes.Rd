\name{response-classes}

\docType{class}

\alias{response-classes}

\alias{GLMresponse-class}
\alias{transInit-class}
\alias{MVNresponse-class}

\title{Class "GLMresponse" and class "transInit"}

\description{Specific instances of response models for \code{\link{depmix}}
models.}

\section{Slots}{
	
	Both \code{GLMresponse} and \code{transInit} contain the
	\code{response}-class. In addition to the slots of that class, these
	classes have the following slots: 
	
	\describe{
    	
		\item{\code{formula}:}{A formula that specifies the model.}

		\item{\code{family}:}{A family object specifying the link 
		function. See the \code{\link{GLMresponse}} help page for 
		possible options.}
	}
	
}

\details{
	
	The \code{GLMresponse}-class offers an interface to the
	\code{\link[stats]{glm}} functions that are subsequently used in fitting
	the \code{depmix} model of which the response is a part.
	
	The \code{transInit} is an extension of \code{response} that is used to
	model the transition matrix and the initial state probabilities by the
	use of a multinomial logistic model, the difference being that in fact
	the response is missing as the transitions between states are not
	observed.  This class has its own fit function which is an interface to
	the multinom function in \code{\link{nnet}}.
	
}

\section{Accessor Functions}{
	The following functions should be used for accessing the corresponding
	slots:
	\describe{
		\item{\code{npar}:}{The number of parameters of the model.}
		\item{\code{getdf}:}{The number of non-fixed parameters.}
	}
}

\author{Ingmar Visser & Maarten Speekenbrink}

\keyword{ methods }
