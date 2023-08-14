
#====================================================================
# Add the value 'lambda' to the diagonal of 'A' matrix
#====================================================================
# a = 1; lambda = 0; inplace = FALSE; rows <- cols <- NULL; verbose = TRUE
penalize_cov <- function(A, a = 1, lambda, inplace = FALSE,
                         rows = NULL, cols = NULL, drop = TRUE,
                         unpack = TRUE, verbose = TRUE)
{
  eps <- .Machine$double.eps

  tmp <- info_matrix(A, check=TRUE)
  if(tmp$nrows != tmp$ncols){
    stop("Input 'A' must be squared matrix")
  }
  n <- tmp$nrows
  type <- tmp$type
  diag <- tmp$diag
  byrow <- tmp$byrow

  if(verbose & type!="full"){
    message(" Penalizing a ",n,"x",n," matrix packed ",ifelse(byrow,"row","column"),
            "-wise with the '",type,"' triangular data")
  }

  if(length(a) == 1){
    a <- as.vector(a)
  }
  if(length(lambda) == 1){
    lambda <- as.vector(lambda)
  }

  dm1 <- dim(a)
  dm2 <- dim(lambda)
  scale <- TRUE
  if((length(dm1)==2L) & (length(dm2)==2L)){
    if(!all(dm1==dm2)){
      stop("Inputs 'a' and 'lambda' must be matrices of the same dimension")
    }

  }else{
    if(length(dm1) == 2L){ # a is matrix, lambda is vector
      if(dm1[1] != dm1[2]){
        stop("Input 'a' must be a squared matrix")
      }
      if(length(lambda) == 1L){
        lambda <- matrix(lambda, nrow=dm1[1], ncol=dm1[2])
      }else{
        if(length(lambda) == dm1[1]){
          lambda <- diag(lambda)
          if(verbose){ message(" A diagonal matrix of dimension ",dm1[1]," for 'lambda' is used") }
        }else{
          stop("Incompatible dimensions between 'a' and 'lambda'")
        }
      }

    }else{
      if(length(dm2) == 2L){ # a is vector, lambda is matrix
        if(dm2[1] != dm2[2]){
          stop("Input 'lambda' must be a squared matrix")
        }
        scale <- any(abs(a-1) > eps)
        if(length(a) == 1L){
          a <- matrix(a, nrow=dm2[1], ncol=dm2[2])
        }else{
          stop("Input 'a' must be an scalar when 'lambda' is a matrix")
        }

      }else{ # a and lambda are vectors
        na <- 1
        if(length(a) != 1){
          stop("Input 'a' must be an scalar")
        }
        if(!(length(lambda)==1 | length(lambda)==n)){
          stop("Input 'lambda' must be an scalar or a vector of length equal to dim(A)=",n)
        }
        scale <- abs(a-1) > eps
      }
    }
  }

  isKronecker <- (length(dim(a))==2L & length(dim(lambda))==2L)
  na <- ifelse(isKronecker,nrow(a),1)
  if(length(dim(a))==2L){
    stopifnot(isSymmetric(a))
  }

  if(!is.null(rows)){
    stopifnot(all(!is.na(rows)))
    if(any(rows<1) | any(rows>(n*na))){
      stop("Input 'rows' must be integers between 1 and dim(A)",ifelse(na==1,"",paste0("x",na)))
    }
  }
  if(!is.null(cols)){
    stopifnot(all(!is.na(cols)))
    if(any(cols<1) | any(cols>(n*na))){
      stop("Input 'cols' must be integers between 1 and dim(A)",ifelse(na==1,"",paste0("x",na)))
    }
  }

  #isBigMatrix <- bigmemory::is.big.matrix(A)
  isBigMatrix <- FALSE

  if(isBigMatrix){
    message(" Routine 'penalize_cov' is not implemented yet for 'bigmatrix'")
    #stopifnot(bigmemory::typeof(A) == "double")

  }else{
    #dyn.load("c_penalize_cov.so")
    if(isKronecker)
    {
      if(inplace){
        stop("Routine cannot be implemented 'inplace' for a kronecker product")
      }
      return(.Call('R_penalize_cov_kronecker',na,a,n,A,lambda,rows-1,cols-1,
                                              drop,type,diag,byrow))

    }else{
      return(.Call('R_penalize_cov',a,n,A,lambda,rows-1,cols-1,drop,type,diag,
                                    byrow,scale,inplace,unpack,verbose))
    }
    #dyn.unload("c_penalize_cov.so")
  }
}
