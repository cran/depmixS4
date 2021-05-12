# M. Speekenbrink, 4-5-2021
#
# Using type argument to allow different forms of posterior probabilities
# and state sequences.
#
# Default of type is currently set to 'viterbi' for backwards compatibility. 
# This return value may not be what is expected from the user, and hence we
# provide a warning when type is not explicitly set and may change the default
# in later releases.

setMethod("posterior","mix",
          function(object, type = c("viterbi", "global", "local", "filtering", "smoothing")) {
            if(identical(type,c("viterbi", "global", "local", "filtering", "smoothing"))) {
              warning("Argument 'type' not specified and will default to 'viterbi'. This default may change in future releases of depmixS4. Please see ?posterior for alternative options.")
            }
            type <- match.arg(type)
            switch(type,
                   viterbi = viterbi(object),
                   global = viterbi(object)[,1],
                   local = apply(forwardbackward(object)$gamma, 1, which.max),
                   filtering = forwardbackward(object)$alpha,
                   smoothing = forwardbackward(object)$gamma)
          }
)


setMethod("posterior","mix.fitted",
          function(object, type = c("viterbi", "global", "local", "filtering", "smoothing")) {
            if(identical(type,c("viterbi", "global", "local", "filtering", "smoothing"))) {
              warning("Argument 'type' not specified and will default to 'viterbi'. This default may change in future releases of depmixS4. Please see ?posterior for alternative options.")
            }
            type <- match.arg(type)
            switch(type,
                   viterbi = object@posterior,
                   global = object@posterior[,1],
                   local = apply(forwardbackward(object)$gamma, 1, which.max),
                   filtering = forwardbackward(object)$alpha,
                   smoothing = forwardbackward(object)$gamma)
          }
)

setMethod("posterior","depmix",
          function(object, type = c("viterbi", "global", "local", "filtering", "smoothing")) {
            if(identical(type,c("viterbi", "global", "local", "filtering", "smoothing"))) {
              warning("Argument 'type' not specified and will default to 'viterbi'. This default may change in future releases of depmixS4. Please see ?posterior for alternative options.")
            }
            type <- match.arg(type)
            switch(type,
                   viterbi = viterbi(object),
                   global = viterbi(object)[,1],
                   local = apply(forwardbackward(object)$gamma, 1, which.max),
                   filtering = forwardbackward(object)$alpha,
                   smoothing = forwardbackward(object)$gamma)
          }
)

setMethod("posterior","depmix.fitted",
          function(object, type = c("viterbi", "global", "local", "filtering", "smoothing")) {
            if(identical(type,c("viterbi", "global", "local", "filtering", "smoothing"))) {
              warning("Argument 'type' not specified and will default to 'viterbi'. This default may change in future releases of depmixS4. Please see ?posterior for alternative options.")
            }
            type <- match.arg(type)
            switch(type,
                   viterbi = object@posterior,
                   global = object@posterior[,1],
                   local = apply(forwardbackward(object)$gamma, 1, which.max),
                   filtering = forwardbackward(object)$alpha,
                   smoothing = forwardbackward(object)$gamma)
          }
)