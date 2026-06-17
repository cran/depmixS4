#
# Source-driven example regression checks for refactoring preparation.
#
# Suggested next unit-test targets:
# - parameter round trips: getpars(), setpars(), freepars(), fixed/equal/conrows
# - likelihood engines: fb(), forwardbackward(), viterbi(), posterior()
# - optimizers: fit() with EM, fixed parameters, equality constraints, and hard EM
# - response classes: GLMresponse and MVNresponse dens(), fit(), predict(), getpars()
# - data-shape handling: ntimes, multivariate responses, missing values, and weights
#
# The examples are run through utils::example() so this file stays tied to the
# documented example code rather than carrying a second copy of it.

require(depmixS4)

loglik <- function(object) {
	suppressWarnings(unname(as.numeric(logLik(object))))
}

check_equal <- function(value, expected, tolerance = 1e-6, label = NULL) {
	comparison <- all.equal(value, expected, tolerance = tolerance,
		check.attributes = FALSE)
	if(!isTRUE(comparison)) {
		message <- "unexpected value"
		if(!is.null(label)) message <- paste(message, label)
		stop(paste(message, paste(comparison, collapse = "; "), sep = ": "),
			call. = FALSE)
	}
}

quiet <- function(expr) {
	value <- NULL
	capture.output(value <- suppressWarnings(force(expr)))
	value
}

run_example <- function(topic, run.dontrun = FALSE) {
	env <- new.env(parent = globalenv())
	quiet(utils::example(topic, package = "depmixS4", ask = FALSE,
		echo = FALSE, local = env, run.dontrun = run.dontrun,
		character.only = TRUE))
	env
}

check_metric <- function(env, name, expected, tolerance = 1e-6) {
	check_equal(loglik(get(name, envir = env)), expected, tolerance = tolerance,
		label = paste(name, "logLik"))
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

check_metric(example_envs$depmix, "fm", -296.107777497663)
check_metric(example_envs$depmix, "fmsp", 1334.635143365299)
check_metric(example_envs$fit, "fmod1", -248.972219690931)
check_metric(example_envs$fit, "fmod2", -249.212895124180)
check_metric(example_envs$fit, "fmod4", -1083.036278731958)
check_equal(example_envs$fit$pst_new, c(2, 1, 1),
	label = "fit new-data Viterbi states")
check_equal(dim(example_envs$fit$pst_prob), c(439, 2),
	label = "fit smoothing dimensions")
check_metric(example_envs$makeDepmix, "fm1", -248.972219690931)
check_metric(example_envs$makeDepmix, "fm2", -297.952933433641)
check_metric(example_envs$multistart, "fmod2", -248.972216998551)
check_metric(example_envs$posterior, "fmod", -248.972216763550)
check_equal(dim(example_envs$posterior$pst_prob), c(439, 2),
	label = "posterior smoothing dimensions")
check_metric(example_envs$vcov, "fmod1", -248.972219690931)
check_metric(example_envs$viterbi, "fmod", -248.972218212895)

cat("Running source-driven dontrun example topics\n")

dontrun_fit <- run_example("fit", run.dontrun = TRUE)
check_metric(dontrun_fit, "fmod1", -248.972219690931)
check_metric(dontrun_fit, "fmod1ms", -247.430638392413)
check_metric(dontrun_fit, "fmod4", -1083.036278731958)
check_metric(dontrun_fit, "fmod5", -951.291754995437)

fmod2_pars <- getpars(dontrun_fit$fmod2)
check_equal(fmod2_pars[c(1, 2, 13, 14)], c(0, 1, 0.5, 0.5),
	label = "fixed-parameter fit constraints")
check_metric(dontrun_fit, "fmod2", -249.212895124180,
	tolerance = 1e-4)

constrained_ll <- c(fmod3 = loglik(dontrun_fit$fmod3),
	fmod3b = loglik(dontrun_fit$fmod3b))
check_equal(constrained_ll, c(fmod3 = -277.106880997915,
	fmod3b = -277.106880997915), tolerance = 1e-4,
	label = "linear-constraint fit logLik")
check_equal(constrained_ll["fmod3"], constrained_ll["fmod3b"],
	tolerance = 1e-4, label = "equivalent linear-constraint specifications")

fmod3_pars <- getpars(dontrun_fit$fmod3)
fmod3b_pars <- getpars(dontrun_fit$fmod3b)
check_equal(fmod3_pars[c(4, 6)], fmod3_pars[c(8, 10)],
	label = "equal-pattern fit constraints")
check_equal(fmod3b_pars[c(4, 6)], fmod3b_pars[c(8, 10)],
	label = "constraint-matrix fit constraints")

if(requireNamespace("gamlss", quietly = TRUE) &&
	requireNamespace("gamlss.dist", quietly = TRUE)) {
	dontrun_makeDepmix <- run_example("makeDepmix", run.dontrun = TRUE)
	check_metric(dontrun_makeDepmix, "fm1", -248.972219690931)
	check_metric(dontrun_makeDepmix, "fm2", -297.952933433641)
	check_metric(dontrun_makeDepmix, "fm3", -224.452639279067,
		tolerance = 1e-4)
} else {
	cat("Skipping makeDepmix dontrun regression; gamlss packages unavailable\n")
}
