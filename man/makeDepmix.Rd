
\name{depmix}

\docType{methods}

\alias{makeDepmix}

\title{Dependent Mixture Model Specifiction: the long way}

\description{
	
	\code{makeDepmix} creates an object of class \code{depmix}.  This
	function is meant for full control, e.g. specifying each response model
	and the transition and prior models 'by hand'.  For the default easier
	specification of models, please see \code{\link{depmix}}.  This
	function is meant for specifying one's own response models.
	
}

\usage{
	
	makeDepmix(response, transition, prior, ntimes = NULL, stationary = TRUE, 
    ...) 	
	
}

\arguments{
	
	\item{response}{A two-dimensional list of response models. See
	'Details'.}
	
	\item{transition}{A list of transition models, each created by a call
	to \code{\link{transInit}} on for possibilities of specifying such
	models.  The lenght of this list should be the nubmer of states of the
	model.}
	
	\item{prior}{The initial state probabilities model; created through a
	call to \code{\link{transInit}}.}

	\item{ntimes}{A vector specifying the lengths of individual, ie
		independent, time series. If not specified, the responses are
		assumed to form a single time series.}
	
	\item{stationary}{Logical indicating whether the transition models
		include time-varying covariates; used internally to determine the
		dimensions of certain arrays, notably \code{trDens}.}
		
	\item{...}{Not used currently.}
		
}

\details{ 

	The function \code{makeDepmix} creates an S4 object of class \code{depmix},
	which needs to be fitted using \code{\link{fit}} to optimize the
	parameters. This function is provided to have full control, eg by
	specifying one's own response models with distributions that are not
	provided. 
	
	The response model(s) should be created by call(s) to
	\code{\link{response}} or one's own created response models that should
	extend the response class and have the following methods: dens, predict
	and optionally fit.  The fit function should have an argument w,
	providing the weights.  If the fit function is not provided,
	optimization should be done by using Rdonlp (use method="donlp" in
	calling fit on the depmix model, this is currently not done
	automatically).  The first index of response models runs over the
	states of the model, and the seconde index over the responses to be
	modeled.
	
}

\value{
	
See the \code{\link{depmix}} help page for the return value, a
\code{depmix} object.

}

\author{Ingmar Visser & Maarten Speekenbrink}

\seealso{

	\code{\link{fit}}, \code{\link{transInit}}, \code{\link{response}},
	\code{\link{depmix-methods}} for accessor functions to \code{depmix}
	objects.

}

\examples{

# below example recreates the model from the depmix help page albeit in a
# roundabout way

data(speed)   

rModels <- list(
	list(
		GLMresponse(formula=rt~1,data=speed,family=gaussian(),pstart=c(5.52,.202)),
		GLMresponse(formula=corr~1,data=speed,family=multinomial(),pstart=c(0.5,0.5))
	),
	list(
		GLMresponse(formula=rt~1,data=speed,family=gaussian(),pstart=c(6.39,.24)),
		GLMresponse(formula=corr~1,data=speed,family=multinomial(),pstart=c(.1,.9))
	)
)

trstart=c(0.9,0.1,0.1,0.9)

transition <- list()
transition[[1]] <- transInit(~1,nstates=2,data=data.frame(1),pstart=c(trstart[1:2]))
transition[[2]] <- transInit(~1,nstates=2,data=data.frame(1),pstart=c(trstart[3:4]))

instart=c(0,1)
inMod <- transInit(~1,ns=2,ps=instart,data=data.frame(rep(1,3)))

mod <- makeDepmix(response=rModels,transition=transition,prior=inMod,ntimes=attr(speed,"ntimes"))

logLik(mod)

fm <- fit(mod)

fm 

summary(fm)

}

\keyword{methods}


