# Classes for simulated mix and depmix models

setClass("mix.sim",
	contains="mix",
	representation(
		states="matrix"
	)
)

setClass("depmix.sim",
	contains="depmix",
	representation(
		states="matrix"
	)
)

