# Class for simulated depmix model

setClass("depmix.sim",
  contains="depmix",
  representation(
    states="matrix"
  )
)