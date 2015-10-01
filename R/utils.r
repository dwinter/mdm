match_rows <- function(big_m, u_m){
 apply(big_m ,1, function(r) which(apply(u_m, 1, identical, x=r)))
}
#' Calcualte the probabiliy that each of a set of observation belongs each
#' component in mixture of Dirichelet Multinomials
#'@export
#'@param x matrix of counts, with on row per obseravation
#'@param fit fitted model from which component probabilities will be calculated
get_component_probs <- function(x, fit){
    if(is.null(fit[["f"]])){
        stop("Model has only one component")
    }
    LL <- dmdm(x=x, phi=fit$params[,1], p=fit$params[,-1], f=fit$f)
    LL$w[match_rows(x, unique(x)),]
}

#'@export
AIC.mdm_model <- function(object, ..., k=NULL){
    if(is.null(object)){
        k <- length(object$freeParams)
    }
    2*k - 2 * object$ll
}

#'@importFrom stats4 BIC
#'@export
BIC.mdm_model <- function(object, ...){
    k <- length(object$freeParams)
    -2 * object$ll + k * log(object$nobs)
}

#'@export
logLik.mdm_model<- function(object, ...) object$ll

#'@export
coef.mdm_model <- function(object, ...) object$params

#'@export
print.mdm_model <- function(x, ...){
    mtype <- if( !is.null(x[["f"]]) ) "mixture" else ""
    cat(paste("Dirichelet multinomial", mtype, "model\n\n"))
    cat("Log likelihood:\t", logLik(x), "\n")
    if(mtype == "mixture"){
        cat("Mixture proportions:\t [", paste(round(x$f,2), collapse=", "), "]\n", sep="")
    }
    cat("Paramater estimates:\n")
    cat("AIC:\t", AIC(x), "\n")
    cat("BIC:\t", BIC(x), "\n\n")
    print(x$params)
    cat("\n\n")
}
