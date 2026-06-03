#
# Example regression tests for refactoring preparation.
#
# Suggested next unit-test targets:
# - parameter round trips: getpars(), setpars(), freepars(), fixed/equal/conrows
# - likelihood engines: fb(), forwardbackward(), viterbi(), posterior()
# - optimizers: fit() with EM, fixed parameters, equality constraints, and hard EM
# - response classes: GLMresponse and MVNresponse dens(), fit(), predict(), getpars()
# - data-shape handling: ntimes, multivariate responses, missing values, and weights
#
# The checks below intentionally exercise the examples and selected former
# \dontrun{} blocks against numeric baselines from the current implementation.

require(depmixS4)

check_equal <- function(value, expected, tolerance = 1e-6, label = NULL) {
	if(!isTRUE(all.equal(value, expected, tolerance = tolerance,
		check.attributes = FALSE))) {
		stop(paste("unexpected value", label), call. = FALSE)
	}
}

loglik <- function(object) {
	suppressWarnings(unname(as.numeric(logLik(object))))
}

quiet <- function(expr) {
	value <- NULL
	capture.output(value <- force(expr))
	value
}

run_example <- function(topic) {
	env <- new.env(parent = globalenv())
	quiet(suppressWarnings(utils::example(topic, package = "depmixS4", ask = FALSE,
		echo = FALSE, local = env, run.dontrun = FALSE,
		character.only = TRUE)))
	env
}

cat("Running Rd examples for regression coverage\n")

example_topics <- c(
	"balance", "depmix-methods", "fit", "depmix", "depmixS4-package",
	"em.control", "forwardbackward", "makeDepmix", "mix-class",
	"mix.fitted-class", "mix", "multistart", "posterior", "responses",
	"simulate", "sp500", "speed", "vcov", "viterbi"
)

example_envs <- lapply(example_topics, run_example)
names(example_envs) <- example_topics

check_equal(loglik(example_envs$depmix$fm), -296.107777497663,
	label = "depmix fm logLik")
check_equal(loglik(example_envs$depmix$fmsp), 1334.635143365299,
	label = "depmix fmsp logLik")
check_equal(loglik(example_envs$fit$fmod1), -248.972219690931,
	label = "fit fmod1 logLik")
check_equal(loglik(example_envs$fit$fmod2), -249.212895124180,
	label = "fit fmod2 logLik")
check_equal(loglik(example_envs$fit$fmod4), -1083.036278731958,
	label = "fit fmod4 logLik")
check_equal(example_envs$fit$pst_new, c(2, 1, 1),
	label = "fit new-data Viterbi states")
check_equal(dim(example_envs$fit$pst_prob), c(439, 2),
	label = "fit smoothing dimensions")
check_equal(loglik(example_envs$makeDepmix$fm1), -248.972219690931,
	label = "makeDepmix fm1 logLik")
check_equal(loglik(example_envs$makeDepmix$fm2), -297.952933433641,
	label = "makeDepmix fm2 logLik")
check_equal(loglik(example_envs$multistart$fmod2), -248.972216998551,
	label = "multistart fmod2 logLik")
check_equal(loglik(example_envs$posterior$fmod), -248.972216763550,
	label = "posterior fmod logLik")
check_equal(dim(example_envs$posterior$pst_prob), c(439, 2),
	label = "posterior smoothing dimensions")
check_equal(loglik(example_envs$vcov$fmod1), -248.972219690931,
	label = "vcov fmod1 logLik")
check_equal(loglik(example_envs$viterbi$fmod), -248.972218212895,
	label = "viterbi fmod logLik")

cat("Running fast local dontrun regression examples\n")

set.seed(3)
y1 <- rpois(50, 1)
y2 <- rpois(50, 2)
ydf <- data.frame(y = c(y1, y2))
m1 <- depmix(y ~ 1, ns = 1, family = poisson(), data = ydf)
set.seed(1)
fm1 <- quiet(fit(m1, verbose = FALSE))
m2 <- depmix(y ~ 1, ns = 2, family = poisson(), data = ydf)
set.seed(1)
fm2 <- quiet(fit(m2, verbose = FALSE))
m3 <- depmix(y ~ 1, ns = 3, family = poisson(), data = ydf)
set.seed(1)
fm3 <- quiet(fit(m3, em = em.control(maxit = 500), verbose = FALSE))
check_equal(c(loglik(fm1), loglik(fm2), loglik(fm3)),
	c(-159.305417842699, -145.527968416625, -145.414961863022),
	label = "Poisson change-point logLik")
check_equal(c(BIC(fm1), BIC(fm2), BIC(fm3)),
	c(323.216005871386, 314.081787763190, 341.486795771914),
	label = "Poisson change-point BIC")

dt <- data.frame(y1 = c(0, 1, 1, 2, 4, 5),
	y2 = c(1, 0, 1, 0, 1, 0),
	y3 = c(4, 4, 3, 2, 1, 1))
m2 <- mix(cbind(y1, y2, y3) ~ 1, data = dt, ns = 2,
	family = multinomial("identity"))
set.seed(1)
fm2 <- quiet(fit(m2, verbose = FALSE))
dm2 <- depmix(cbind(y1, y2, y3) ~ 1, data = dt, ns = 2,
	family = multinomial("identity"))
set.seed(1)
fdm2 <- quiet(fit(dm2, verbose = FALSE))
check_equal(c(loglik(fm2), loglik(fdm2)),
	c(-15.522894709188, -13.547386015719),
	label = "multicolumn multinomial logLik")

data(speed)
speed[2, 1] <- NA
speed[3, 2] <- NA
mod1ms <- depmix(list(rt ~ 1, corr ~ 1), data = speed, transition = ~Pacc,
	nstates = 2, family = list(gaussian(), multinomial("identity")),
	ntimes = c(168, 134, 137))
set.seed(3)
fmod1ms <- suppressWarnings(quiet(fit(mod1ms, verbose = FALSE)))
check_equal(loglik(fmod1ms), -247.430638392413,
	label = "missing-data fit logLik")

data(speed)
mod1 <- depmix(list(rt ~ 1, corr ~ 1), data = speed, transition = ~Pacc,
	nstates = 2, family = list(gaussian(), multinomial("identity")),
	ntimes = c(168, 134, 137))
set.seed(3)
fmod1 <- quiet(fit(mod1, verbose = FALSE))
pars <- c(unlist(getpars(fmod1)))
pars[1] <- 0
pars[2] <- 1
pars[13] <- 0.5
pars[14] <- 0.5
mod2 <- setpars(mod1, pars)
free <- c(0, 0, rep(c(0, 1), 4), 1, 1, 0, 0, 1, 1, 1, 1)
fmod2 <- quiet(fit(mod2, fixed = !free, verbose = FALSE))
pars <- c(unlist(getpars(fmod2)))
pars[4] <- pars[8] <- -4
pars[6] <- pars[10] <- 10
mod3 <- setpars(mod2, pars)
conpat <- c(0, 0, rep(c(0, 1), 4), 1, 1, 0, 0, 1, 1, 1, 1)
conpat[4] <- conpat[8] <- 2
conpat[6] <- conpat[10] <- 3
fmod3 <- suppressWarnings(quiet(fit(mod3, equal = conpat, verbose = FALSE)))
conr <- matrix(0, 2, 18)
conr[1, 4] <- 1
conr[1, 8] <- -1
conr[2, 6] <- 1
conr[2, 10] <- -1
fmod3b <- suppressWarnings(quiet(fit(mod3, conrows = conr,
	fixed = !free, verbose = FALSE)))
check_equal(c(loglik(fmod1), loglik(fmod2), loglik(fmod3), loglik(fmod3b)),
	c(-248.972219690931, -249.212895124180,
		-277.106880997915, -277.106880997915),
	label = "linear-constraint fit logLik")

data(balance)
mod4 <- mix(list(d1 ~ 1, d2 ~ 1, d3 ~ 1, d4 ~ 1), data = balance,
	nstates = 2, family = list(multinomial("identity"),
		multinomial("identity"), multinomial("identity"),
		multinomial("identity")))
set.seed(1)
fmod4 <- quiet(fit(mod4, verbose = FALSE))
mod5 <- mix(list(d1 ~ 1, d2 ~ 1, d3 ~ 1, d4 ~ 1), data = balance,
	nstates = 2, family = list(multinomial("identity"),
		multinomial("identity"), multinomial("identity"),
		multinomial("identity")), prior = ~age, initdata = balance)
set.seed(1)
fmod5 <- quiet(fit(mod5, verbose = FALSE))
check_equal(c(loglik(fmod4), loglik(fmod5)),
	c(-1083.036278731958, -951.291754995437),
	label = "balance prior logLik")

if(requireNamespace("gamlss", quietly = TRUE) &&
	requireNamespace("gamlss.dist", quietly = TRUE)) {
	if(!isClass("exgaus")) setClass("exgaus", contains = "response")
	if(!isGeneric("exgaus")) {
		setGeneric("exgaus", function(y, pstart = NULL, fixed = NULL, ...)
			standardGeneric("exgaus"))
	}
	setMethod("exgaus", signature(y = "ANY"),
		function(y, pstart = NULL, fixed = NULL, ...) {
			y <- matrix(y, length(y))
			x <- matrix(1)
			parameters <- list()
			npar <- 3
			if(is.null(fixed)) fixed <- as.logical(rep(0, npar))
			if(!is.null(pstart)) {
				if(length(pstart) != npar) {
					stop("length of 'pstart' must be ", npar)
				}
				parameters$mu <- pstart[1]
				parameters$sigma <- log(pstart[2])
				parameters$nu <- log(pstart[3])
			}
			new("exgaus", parameters = parameters, fixed = fixed, x = x,
				y = y, npar = npar)
		})
	setMethod("dens", "exgaus", function(object, log = FALSE) {
		gamlss.dist::dexGAUS(object@y, mu = predict(object),
			sigma = exp(object@parameters$sigma),
			nu = exp(object@parameters$nu), log = log)
	})
	if(!hasMethod("getpars", "response")) {
		setMethod("getpars", "response",
			function(object, which = "pars", ...) {
				switch(which, "pars" = unlist(object@parameters),
					"fixed" = object@fixed)
			})
	}
	setMethod("setpars", "exgaus",
		function(object, values, which = "pars", ...) {
			npar <- npar(object)
			if(length(values) != npar) {
				stop("length of 'values' must be", npar)
			}
			nms <- names(object@parameters)
			switch(which,
				"pars" = {
					object@parameters$mu <- values[1]
					object@parameters$sigma <- values[2]
					object@parameters$nu <- values[3]
				},
				"fixed" = {
					object@fixed <- as.logical(values)
				})
			names(object@parameters) <- nms
			object
		})
	setMethod("fit", "exgaus", function(object, w) {
		if(missing(w)) w <- NULL
		fitted <- gamlss::gamlss(object@y ~ 1, weights = w,
			family = gamlss.dist::exGAUS(),
			control = gamlss::gamlss.control(n.cyc = 100, trace = FALSE),
			mu.start = object@parameters$mu,
			sigma.start = exp(object@parameters$sigma),
			nu.start = exp(object@parameters$nu))
		pars <- c(fitted$mu.coefficients, fitted$sigma.coefficients,
			fitted$nu.coefficients)
		setpars(object, pars)
	})
	setMethod("predict", "exgaus", function(object) object@parameters$mu)
	data(speed)
	rt <- speed$rt
	rModels <- list(
		list(exgaus(rt, pstart = c(5, .1, .1)),
			GLMresponse(formula = corr ~ 1, data = speed,
				family = multinomial("identity"), pstart = c(0.5, 0.5))),
		list(exgaus(rt, pstart = c(6, .1, .1)),
			GLMresponse(formula = corr ~ 1, data = speed,
				family = multinomial("identity"), pstart = c(.1, .9)))
	)
	trstart <- c(0.9, 0.1, 0.1, 0.9)
	transition <- list()
	transition[[1]] <- transInit(~Pacc, nstates = 2, data = speed,
		pstart = c(trstart[1:2], 0, 0))
	transition[[2]] <- transInit(~Pacc, nstates = 2, data = speed,
		pstart = c(trstart[3:4], 0, 0))
	inMod <- transInit(~1, ns = 2, ps = c(0.5, 0.5),
		family = multinomial("identity"), data = data.frame(rep(1, 3)))
	mod <- makeDepmix(response = rModels, transition = transition,
		prior = inMod, ntimes = c(168, 134, 137), homogeneous = FALSE)
	fm3 <- quiet(fit(mod, emc = em.control(rand = FALSE), verbose = FALSE))
	check_equal(loglik(fm3), -224.4526, tolerance = 1e-4,
		label = "exgaus custom response logLik")
} else {
	cat("Skipping exgaus dontrun regression; gamlss packages unavailable\n")
}
