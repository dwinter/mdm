#' Generate random samples of Direchelet Multinomial observations
#' @param n  number of observations 
#' @param m  observation sizes vector possibly of length one.
#' @param p  matrix or vector of proportions
#' @param phi Overdispersion paramater. Value must be in [0,1], where phi=0 
#' is 'pure' Multinomial distribuion and increasing values of phi lead to
#' increasingly over-dispersed distributions. 
#' @param scale Scale parameters, vector of matrix.

#' @return A matrix of n-rows, each row representing a random draw from 
#' @export
#' @examples
#' #simulate reads with 1% sequencing error
#' genotypes <- matrix( 0.01/3, nrow=4, ncol=4)
#' diag(genotypes) <- 0.99
#' round(genotyes,3)
#' rdm(4, 40, phi=0, p=genotypes)
#' rdm(4, 40, phi=0.01, p=genotypes)



rdm <- function(n, m, p, phi=NULL, scale=NULL) {
    if ((is.null(scale) + is.null(p)) != 1){
        stop("Must specify either (but not both) 'scale' (alpha) or 'p' vector")
    }
    if(is.null(scale)){
	    params <- mdmParams(phi,p)
    
    	if(any(params[,1] < 0.0 | 1.0 < params[,1])) {
		    stop("phi must be in [0,1].")
	    }
        if(any(params[,-1] <= 0.0 | 1.0 <= params[,-1])) {
	    	stop("p must be in (0,1).")
    	}
	    params[params[,1] < .Machine$double.eps/2,1] <- .Machine$double.eps/2
    	p <- params[,-1,drop=FALSE]
    
    }
    else{
       params <- mdmParams(scale)
    }
    alphas <- mdmAlphas(params)
	# choose initial 
	y <- mc2d::rmultinomial(n,1,p)
	# update params, conditional on what has occurred
	ny <- nrow(y)
	na <- nrow(alphas)
	if(na != ny) {
		n1 <- ny %/% na
		n2 <- ny %% na
		u <- rep(seq_len(na),n1)
		if(n2 > 0) {
			u <- c(u,seq_len(n2))
		}
		a <- y + alphas[u,]
	} else {
		a <- y + alphas
	}
	# choose following
	y+ mc2d::rmultinomial(n,m-1,mc2d::rdirichlet(n,a))	
}

#' Generate random samples from a mixture of Direchelet Multinomials
#' @param n  number of observations 
#' @param m  observation sizes vector possibly of length one.
#' @param p  matrix or vector of proportions
#' @param f mixture proportions
#' @param phi Overdispersion paramater. Value must be in [0,1], where phi=0 
#' is 'pure' Multinomial distribuion and increasing values of phi lead to
#' increasingly over-dispersed distributions. 
#' @param scale Scale parameters, vector of matrix.

#' @return A matrix of n-rows, each row representing a random draw from 
#' @export
#' @examples
#' #simulate reads with 1% sequencing error
#' #from AT-biased genome
#' genotypes <- matrix( 0.01/3, nrow=4, ncol=4)
#' diag(genotypes) <- 0.99
#' nfreq <- c(0.4, 0.1, 0.1, 0.4)
#' rmdm(10, 40, phi=0, f=nfreq p=genotypes)
#' rmdm(10, 40, phi=0.01,f=nfreq,  p=genotypes)


rmdm <- function(n, m, p, f, phi=NULL, scale=NULL ) {
    if ((is.null(scale) + is.null(p)) != 1){
        stop("Must specify either (but not both) 'scale' (alpha) or 'p' proportions")
    }
    if(is.null(scale)){
	    params <- mdmParams(phi,p)
    
    	if(any(params[,1] < 0.0 | 1.0 < params[,1])) {
		    stop("phi must be in [0,1].")
	    }
        if(any(params[,-1] <= 0.0 | 1.0 <= params[,-1])) {
	    	stop("p must be in (0,1).")
    	}
	    params[params[,1] < .Machine$double.eps/2,1] <- .Machine$double.eps/2
    	p <- params[,-1,drop=FALSE]
    
    }
    else{
        params <- mdmParams(scale)
    }
    alphas <- mdmAlphas(params)
	k <- nrow(p)
	if(length(f) != k) {
		stop("The length of 'f' and number of rows in params must be equal.")
	}
	
	# generate the mixture
	q <-  mc2d::rmultinomial(1, n, f)
	mix <- rep.int(seq_len(k), q)
	mix <- sample(mix)
	# generate the samples
	x <- rdm(n,m,params[mix,])
	# use the rownames to store the mixture categories
	rownames(x) <- mix
	x
}