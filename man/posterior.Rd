\name{posterior}

\alias{posterior}

\alias{posterior,depmix-method}
\alias{posterior,mix-method}
\alias{posterior,depmix.fitted-method}
\alias{posterior,mix.fitted-method}

\title{Posterior state/class probabilities and classification}

\description{

	Return posterior state classifications and/or
	probabilities for a fitted \code{(dep-)mix} object. In
	the case of a latent class or mixture model, states refer to the 
	classes/mixture components.
	
	There are different ways to define posterior state probabilities and the 
	resulting classifications. The 'type' argument can be used to specify the
	desired definition. The default is currently set to 'viterbi'.
	Other options are 'global' and 'local' for state classification, and 
	'filtering' and 'smoothing' for state probabilities. See Details for more 
	information.
}

\usage{
	\S4method{posterior}{depmix}(object, type = c("viterbi", "global", "local", "filtering", "smoothing"))
	\S4method{posterior}{depmix.fitted}(object, type = c("viterbi", "global", "local", "filtering", "smoothing"))
	\S4method{posterior}{mix}(object, type = c("viterbi", "global", "local", "filtering", "smoothing"))
	\S4method{posterior}{mix.fitted}(object, type = c("viterbi", "global", "local", "filtering", "smoothing"))
}

\arguments{
  \item{object}{A (fitted)(dep-)mix object.}
  \item{type}{character, partial matching allowed. The type of classification or posterior probability desired. }
}

\value{

  The return value of \code{posterior} depends on the value of the \code{type}
  argument:

	\item{type = 'viterbi'}{Returns a data.frame with \code{nstates(object) + 1}
	columns; the first column contains the states decoded through the Viterbi 
	algorithm, the remaining columns contain the (normalized) delta probabilities.}
	
	\item{type = 'global'}{Returns a vector which contains the states decoded 
	through the Viterbi algorithm.}
	
	\item{type = 'local'}{Returns a vector which contains the states decoded 
	as the maximum of the smoothing probabilities.}
	
	\item{type = 'filtering'}{Returns a matrix which contains the posterior 
	probabilities of each state, conditional upon the responses observed thus 
	far.}
	
	\item{type = 'smoothing'}{Returns a matrix which contains the posterior 
	probabilities of each state, conditional upon all the responses observed.}
	
	See Details for more information.
}

\references{ 	

	Lawrence R. Rabiner (1989). A tutorial on hidden Markov models and
	selected applications in speech recognition. \emph{Proceedings of
	IEEE}, 77-2, p.  267-295. 

}

\author{Maarten Speekenbrink & Ingmar Visser}

\details{

  After fitting a \code{mix} or \code{depmix} model, one is often interested 
  in determining the most probable mixture components or hidden states at each 
  time-point \var{t}. This is also called decoding the hidden states from the observed 
  data. There are at least two general ways to consider state classification:
  'global' decoding means determining the most likely state sequence, whilst 
  'local' decoding means determining the most likely state at each time point 
  whilst not explicitly considering the identity of the hidden states at other
  time points. For mixture models, both forms of decoding are identical.
  
  Global decoding is based on the conditional probability
  \eqn{p(S_1, \ldots, S_T \mid Y_1, \ldots, Y_T)}, and consists of determining, 
  at each time point \eqn{t = 1, \ldots, T}:
  \deqn{s*_t = \arg \max_{i=1}^N p(S_1 = s*_1, \ldots, S_{t-1} = s*_{t-1}, S_t = i, S_{t+1} = s*_{t+1}, \ldots, S_T = s*_{T} \mid Y_1, \ldots, Y_T)}
  where \var{N} is the total number of states. These probabilities and the 
  resulting classifications, are computed through the \code{\link{viterbi}} algorithm.
  Setting \code{type = 'viterbi'} returns a \code{data.frame} with the Viterbi-decoded 
  global state sequence in the first column, and the normalized "delta" probabilities
  in the remainining columns. These "delta" probabilities are defined as the joint 
  probability of the most likely state sequence ending in state \var{i} at time \var{t},
  and all the observations up to time \var{t}. The normalization of these joint
  probabilities is done on a time-point basis (i.e., dividing the delta probability
  by the sum of the delta probabilities for that time point for all possible states 
  \var{j} (including state \var{i})). These probabilities are not straightforward 
  to interpret. Setting \code{type = "global"} returns just a vector with the 
  Viterbi-decoded global state sequence.
  
  Local decoding is based on the smoothing probabilities 
  \eqn{p(S_t \mid Y_1, \ldots, Y_T)}, which are the "gamma" probabilities 
  computed with the \code{\link{forwardbackward}} algorithm. Local decoding then
  consists of determining, at each time point \eqn{t = 1, \ldots, T} 
  \deqn{s*_t = \arg \max_{i=1}^N p(S_t = i \mid Y_1, \ldots, Y_T)} 
  where \var{N} is the total number of states. Setting \code{type = "local"} returns
  a vector with the local decoded states. Setting \code{type = "smoothing"} returns
  the smoothing probabilities which underlie this classification. When considering
  the posterior probability of each state, the values returned by \code{type = "smoothing"}
  are most likely what is wanted by the user.
  
  The option \code{type = "filtering"} returns a matrix with the so-called filtering probabilities, 
  defined as \eqn{p(S_t \mid Y_1, \ldots, Y_t)}, i.e. the probability of a hidden 
  state at time \var{t} considering the observations up to and including time \var{t}.

	See the \code{\link{fit}} help page for an example. 
}

\note{The initial version of this function was a simple wrapper to return the value of the \code{posterior} slot in a \code{mix-fitted} or \code{depmix-fitted} object. The value of this slot is set by a call of the \code{viterbi} method. For backwards compatibility, the default value of the \code{type} argument is set to \code{"viterbi"}, which returns the same. As the "delta" probabilities returned as part of this may be misinterpreted, and may not be the desired posterior probabilities, the updated version of this method now allows for other return values, and the \code{type = "viterbi"} option should be considered depreciated.}

\examples{

data(speed)

# 2-state model on rt and corr from speed data set 
# with Pacc as covariate on the transition matrix
# ntimes is used to specify the lengths of 3 separate series
mod <- depmix(list(rt~1,corr~1),data=speed,transition=~Pacc,nstates=2,
	family=list(gaussian(),multinomial("identity")),ntimes=c(168,134,137))
fmod <- fit(mod)

# Global decoding:
pst_global <- posterior(fmod, type = "global")

# Local decoding:
pst_local <- posterior(fmod,type="local")

# Global and local decoding provide different results:
identical(pst_global, pst_local)

# smoothing probabilities are used for local decoding, and may be used as 
# easily interpretable posterior state probabilities
pst_prob <- posterior(fmod, type = "smoothing")

# "delta" probabilities from the Viterbi algorithm
pst_delta <- posterior(fmod, type="viterbi")[,-1]

# The smoothing and "delta" probabilities are different:
identical(pst_prob, pst_delta)

# Filtering probabilities are an alternative to smoothing probabilities:
pst_filt <- posterior(fmod, type = "filtering")

# The smoothing and filtering probabilities are different:
identical(pst_prob, pst_filt)
}

\keyword{ methods }

