make_params <- function(phi = NULL, p=NULL, scale=NULL){
    if ((is.null(scale) + is.null(phi)) != 1){
        stop("Must specify either (but not both) 'scale' (alpha) or 'phi'")
    }
    if(is.null(scale)){
        if(is.null(p)){
            stop("Must specify 'p' probability vector/matrix when using phi parameterization")
        }
	    params <- mdmParams(phi,p)
    } else   params <- mdmParams(scale)
        if(any(params[,1] < 0.0 | 1.0 < params[,1])) {
        stop("phi must be in [0,1].")
    }
    if(any(params[,-1] <= 0.0 | 1.0 <= params[,-1])) {
        stop("p must be in (0,1).")
    }

    params[params[,1] < .Machine$double.eps/2,1] <- .Machine$double.eps/2
    params
}



#' Generate random samples of Direchelet Multinomial observations
#' @param n  number of observations 
#' @param m  observation sizes vector possibly of length one.
#' @param p  matrix or vector of proportions
#' @param phi overdispersion paramater. Value must be in [0,1], where phi=0 
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



rdm <- function(n, m, phi=NULL, p=NULL, scale=NULL) {

    params <- make_params(phi, p, scale) 
    
    p <- params[,-1,drop=FALSE]
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

#' Calculate the likelihood of a given Direchelet Multiomial model
#' @param p  matrix or vector of proportions
#' @param phi overdispersion paramater. Value must be in [0,1], where phi=0 
#' is 'pure' Multinomial distribuion and increasing values of phi lead to
#' increasingly over-dispersed distributions. 
#' @param scale Scale parameters, vector of matrix.
#' @export

ddm <- function(x, phi=NULL, scale=NULL, p=NULL, log=TRUE){
   params <- make_params(phi, p, scale)
   summ <- mdmSumStats(x)
   params <- params[1,c(TRUE, summ$mask)]
   res <- mdmSingleLogLikeCore(s, params)
   if(!log){
      return(exp(res))
   }
  res 
}




#' Generate random samples from a mixture of Dirchelet Multinomials
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
#' rmdm(10, 40, phi=0,    f=nfreq, p=genotypes)
#' rmdm(10, 40, phi=0.01, f=nfreq, p=genotypes)


rmdm <- function(n, m,  f, phi=NULL, p=NULL, scale=NULL ) {

    params <- make_params(phi, p, scale) 
	k <- nrow(params)
	if(length(f) != k) {
		stop("The length of 'f' and number of rows in params must be equal.")
	}
	
	# generate the mixture
	q <-  mc2d::rmultinomial(1, n, f)
	mix <- rep.int(seq_len(k), q)
	mix <- sample(mix)
	# generate the samples
    alphas <- mdmAlphas(params)
	x <- rdm(n,m,scale=alphas[mix,])
	# use the rownames to store the mixture categories
	rownames(x) <- mix
	x
}



#' Calculate the likelihood of a given mixture of Direchelet Multiomials model
#' @param p  matrix or vector of proportions
#' @param phi overdispersion paramater. Value must be in [0,1], where phi=0 
#' is 'pure' Multinomial distribuion and increasing values of phi lead to
#' increasingly over-dispersed distributions. 
#' @param scale Scale parameters, vector of matrix.
#' @param f mixture proportions
#' @export

dmdm <- function(x, phi=NULL, scale=NULL, p=NULL, f, log=TRUE){
   params <- make_params(phi, p, scale)
   r <- mdmAugmentData(x)
   mask <- c(TRUE,r$mask)
   params <- params[,mask]
   res <- mdmLogLikeCore(r$y,r$w,f, params)
   if(!log){
      return(exp(res))
   }
  res 
}


