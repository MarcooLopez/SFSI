
#====================================================================
# Add the value 'lambda' to the diagonal of 'A' matrix
#====================================================================
penalize_cov <- function(A, lambda, inplace = FALSE,
                         n = NULL, uplo = NULL, byrow = NULL)
{
  if(missing(lambda)){
    stop("Argument 'lambda' is missing, with no default")
  }

  dm <- dim(A)
  if(length(dm) == 2L){
    if(dm[1] != dm[2]){
      stop("Input 'A' must be a squared matrix")
    }
    n <- dm[1]
    type <- "full"
  }else{
    if(is.null(n)){
      if('n' %in% names(attributes(A))){
        n <- attr(A, "n")
      }else{
        stop("Dimension 'n' of the matrix must be provided")
      }
    }
    if(length(A) != (n*(n+1)/2)){
      stop("Input 'A' must contain n(n+1)/2 entries when is not a matrix")
    }

    if('uplo' %in% names(attributes(A))){
      type <- attr(A, "uplo")
    }else{
      if(is.null(uplo)){
        stop("Parameter 'uplo' must be specified: either 'upper' or 'lower'")
      }
      type <- uplo
    }

    if('byrow' %in% names(attributes(A))){
      byrow <- attr(A, "byrow")
    }else{
      if(is.null(byrow)){
        stop("Parameter 'byrow' whether data is stored by row must be must be specified")
      }
    }
  }

  stopifnot(length(lambda)==1 | length(lambda)==n)

  type <- ifelse(type=="full",0L,ifelse(type=="upper",1L,2L))
  #isBigMatrix <- bigmemory::is.big.matrix(A)
  isBigMatrix <- FALSE

  if(isBigMatrix){
    message(" Routine 'penalize_cov' is not implemented yet for 'bigmatrix'")
    #stopifnot(bigmemory::typeof(A) == "double")

  }else{
    #if(storage.mode(A) != "double") storage.mode(A) <- "double"

    #dyn.load("c_penalize_cov.so")
    if(inplace){
      out <- .Call('R_penalize_cov',n, A, lambda, type, byrow)
    }else{
      out <- A[]
      tmp <- .Call('R_penalize_cov',n, out, lambda, type, byrow)
    }
    #dyn.unload("c_penalize_cov.so")
  }

  invisible(out)
}
