
# 
# Started by Ingmar Visser 26-2-2008
# 
# Usage: go to trunk directory and source("depmixNew-test4.R")
# 

# 
# BALANCE SCALE data example with age as covariate on class membership
# 

# library(depmixS4) 

# setwd("/Users/ivisser/Documents/projects/depmixProject/depmixNew/rforge/depmix/trunk/")

# 
# optimization speed profile: case 1: latent class data
# 

require(depmixS4)

data(balance)

# now fit some latent class models
trstart=c(1,0,0,1) # as this is a latent class model, the transition are not optimized
instart=c(0.5,0.5)
set.seed(1)
respstart=runif(16)
# note that ntimes argument is used to make this a mixture model
mod1 <- depmix(list(d1~1,d2~1,d3~1,d4~1), data=balance, nstates=2,
	family=list(multinomial(),multinomial(),multinomial(),multinomial()),
	respstart=respstart,trstart=trstart,instart=instart,
	ntimes=rep(1,nrow(balance)))

logLik(mod1)

# mod1 <- fit(mod)

gc()
Rprof(file="lca1")
mod1 <- fit(mod1)
Rprof(NULL)
summaryRprof("lca1")


# 
# optimization speed profile: case 1: latent class data with cov on prior
# 

data(balance)

instart=c(0.5,0.5,0,0)
respstart=c(rep(c(0.1,0.9),4),rep(c(0.9,0.1),4))
trstart=c(1,0,0,1)
mod2 <- depmix(list(d1~1,d2~1,d3~1,d4~1), data=balance, nstates=2,
	family=list(multinomial(),multinomial(),multinomial(),multinomial()),
	trstart=trstart, instart=instart, respstart=respstart,
	ntimes=rep(1,nrow(balance)), prior=~age, initdata=balance)

gc()
Rprof(file="lca2")
mod2 <- fit(mod2)
Rprof(NULL)
summaryRprof("lca2")


# 
# speed data, no covariates
# 

data(speed)

set.seed(1)
trstart=runif(4)

mod <- depmix(list(rt~1,corr~1),data=speed,nstates=2,family=list(gaussian(),multinomial()), trstart=trstart)

logLik(mod)

gc()
Rprof("speed1")
mod1 <- fit(mod)
Rprof(NULL)
summaryRprof("speed1")


trstart=c(0.899,0.101,0,0.01,0.084,0.916,0,0)
instart=c(0.5,0.5)
resp <- c(5.52,0.202,0.472,0.528,6.39,0.24,0.098,0.902)

mod <- depmix(list(rt~1,corr~1),data=speed,transition=~Pacc,nstates=2,family=list(gaussian(),multinomial()),
	respstart=resp,trstart=trstart,instart=instart)

logLik(mod)

gc()
Rprof("speed2")
mod1 <- fit(mod)
Rprof(NULL)
summaryRprof("speed2")

summaryRprof("lca1")
summaryRprof("lca2")
summaryRprof("speed1")
summaryRprof("speed2")

