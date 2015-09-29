# convert a parameter vector to alphas
# v = a paramter vector contain phi and p
mdmAlphas <- function(v,total=FALSE) {
	if(is.vector(v)) {
		v <- t(v)
	}
	phi <- v[,1]
	p <- v[,-1,drop=FALSE]
	at <- ((1.0-phi)/phi)
	a <- p * at
	colnames(a) <- paste("a", seq_len(ncol(a)),sep="")
	if(total) {
		a <- cbind(a,aa=at)
	}
	rownames(a) <- NULL
	structure(a, class="mdmAlphas")
}

#Convert proportions or dispersion params into a paramater vector
mdmParams <- function(phi, p=NULL) {
	if(inherits(phi, "mdmParams")) {
		return(phi)
	}
	if(!is.null(p)) {
		if(is.vector(p)) {
			p <- t(p)
		}
		p <- p/rowSums(p)
	} else {
		if(is.vector(phi)) {
			a <- t(phi)
		} else {
			a <- phi
		}
		A <- rowSums(a)
		phi <- 1.0/(A+1.0)
		p <- a/A		
	}

	colnames(p) <- paste("p", seq_len(ncol(p)),sep="")
	v <- cbind(phi,p)
	class(v) <- "mdmParams"
	v
}

`[.mdmParams` <- function(x, i, j, ...) {
  y <- NextMethod(.Generic)
  class(y) <- .Class
  y
}
