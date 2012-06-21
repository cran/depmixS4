em.control <- 
function(maxit=500,tol=1e-8,crit="relative",random.start=TRUE) {
	return(list(maxit=maxit,tol=tol,crit=crit,random.start=random.start))
}
