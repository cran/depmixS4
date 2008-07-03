\name{mix-class}

\docType{class}

\alias{mix-class}

\alias{nresp,mix-method}
\alias{nstates,mix-method}
\alias{ntimes,mix-method}

\title{Class "mix"}
\description{ A \code{\link{mix}} model. }

\section{Objects from the Class}{
Objects can be created by calls to \code{\link{mix}}.
}

\section{Slots}{

	\describe{ 

		\item{\code{response}:}{List of list of \code{response} objects.}

		\item{\code{prior}:}{\code{\link{transInit}} object; model for the
		prior probabilities, also unconditional probabilities}
		
		\item{\code{dens}:}{Array of dimension sum(ntimes)*nresp*nstates
		providing the densities of the observed responses for each state.}
		
		\item{\code{init}:}{Array of dimension
		\code{length(ntimes)}*nstates with the current predictions for the
		initial state probabilities.}
		
		\item{\code{nstates}:}{The number of states (classes) of the model.}
		
		\item{\code{nresp}:}{The number of independent responses.}
		
		\item{\code{ntimes}:}{A vector of 1's for each case; for internal
		use.}
		
		\item{\code{npars}:}{The total number of parameters of the model.
		This is not the degrees of freedom, ie there are redundancies in
		the parameters, in particular in the multinomial models for the
		transitions and prior.}

	}
}

\section{Accessor Functions}{
	The following functions should be used for accessing the corresponding
	slots:
	\describe{
		\item{\code{npar}:}{The number of parameters of the model.}
		\item{\code{nresp}:}{The number of responses.}
		\item{\code{nstates}:}{The number of states.}
		\item{\code{ntimes}:}{The vector of independent time series lengths.}
	}
}

\author{Ingmar Visser}

\examples{
	showClass("mix")
}

\keyword{classes}
