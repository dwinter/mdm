hash <- function(x) apply(x, 1, paste, collapse="")

#' Calcualte the probabiliy that each of a set of observation belongs each
#' component in mixture of Direchelet Multinomials
#'@export
#'@param x matrix of counts, with on row per obseravation
#'@param fit fitted model from which component probabilities will be calculated
get_component_probs <- function(x, fit){
    LL <- dmdm(x=x, phi=fit$params[,1], p=fit$params[,-1], f=fit$f)
    rownames(LL$w) <- hash(unique(x))
    LL$w[hash(x),]
}

#'@export
AIC.mdm_model <- function(x){
    k <- length(x$freeParams)
    2*k - 2 * x$ll
}

#'@export
BIC.mdm_model <- function(x){
    k <- length(x$freeParams)
    -2 * x$ll + k * log(x$nobs)
}

#'@export
logLik.mdm_model<- function(x) x$ll

#'@export
coef.mdm_model <- function(x) x$params

#'@export
print.mdmd_model <- function(x, ...){
    mtype <- if( !is.null(x[["f"]]) ) "mixture" else ""
    cat(paste("Direchelet multinomial", mtype, "model\n\n"))
    cat("Log likelihood:\t", logLik(x), "\n")
    if(mtype == "mixture"){
        cat("Mixture proportions:\t [", paste(round(x$f,2), collapse=", "), "]\n", sep="")
    }
    cat("BIC:\t", BIC(x), "\n\n")
    cat("Paramater estimates:\n")
    print(x$params)
    cat("\n\n")
}
