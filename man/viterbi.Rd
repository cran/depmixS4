\name{viterbi}

\alias{viterbi}
\alias{viterbi2}
\alias{viterbi.fb}

\title{Viterbi algorithm for decoding the most likely state sequence}

\description{

	Apply the Viterbi algorithm to compute the maximum a posteriori state sequence
	for a \code{mix} or \code{depmix} object.

}

\usage{

	viterbi(object, na.allow=TRUE)

}

\arguments{
  
	\item{object}{A \code{mix} or \code{depmix} object.}
  	
	\item{na.allow}{logical. If TRUE, the density of missing responses is set to 
	1, similar as in the \code{\link{forwardbackward}} algorithm. If FALSE, missing 
	values have NA as density values, and will result in an error.}

}

\value{
	
	\code{viterbi} returns a \code{data.frame} with in the first column 
	the maximum a posteriori state sequence. This is a vector with integers 
	corresponding to the index of the most likely hidden states. The remaining
	columns contain the normalized "delta" probabilities (see Details). 

}

\details{The Viterbi algorithm is used for global decoding of the hidden state 
sequence. Global decoding is based on the conditional probability
  \eqn{p(S_1, \ldots, S_T \mid Y_1, \ldots, Y_T)}, and consists of determining, 
  at each time point \eqn{t = 1, \ldots, T}:
  \deqn{s*_t = \arg \max_{i=1}^N p(S_1 = s*_1, \ldots, S_{t-1} = s*_{t-1}, S_t = i, S_{t+1} = s*_{t+1}, \ldots, S_T = s*_{T} \mid Y_1, \ldots, Y_T)}
  where \var{N} is the total number of states.
  
  The Viterbi algorithm is a dynamic programming algorithm that relies on 
  "delta" probabilities (see Rabiner, 1989), which are defined as the joint 
  probability of the most likely state sequence ending in state \var{i} at time \var{t},
  and all the observations up to time \var{t}. The implementation here normalizes 
  these probabilities on a time-point basis, dividing the delta probability
  by the sum of the delta probabilities for that time point for all possible states 
  \var{j} (including state \var{i})). The normalized delta probabilities for
  each state are returned in columns \code{2:(nstates(object) + 1)}, whilst
  column 1 contains the indices of the maximum a posteriori states.
  
}

\references{ 	

	Lawrence R. Rabiner (1989).  A tutorial on hidden Markov models and
	selected applications in speech recognition.  \emph{Proceedings of
	IEEE}, 77-2, p.  267-295. 

}

\author{Maarten Speekenbrink}

\examples{

data(speed)

# 2-state model on rt and corr from speed data set 
# with Pacc as covariate on the transition matrix
# ntimes is used to specify the lengths of 3 separate series
mod <- depmix(list(rt~1,corr~1),data=speed,transition=~Pacc,nstates=2,
	family=list(gaussian(),multinomial("identity")),ntimes=c(168,134,137))
fmod <- fit(mod)
# result of viterbi is stored in a depmix-fitted object in slot "posterior"
identical(viterbi(fmod),fmod@posterior)
}