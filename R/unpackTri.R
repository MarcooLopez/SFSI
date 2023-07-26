
#====================================================================
# Convert a packed upper/lower triangle matrix into a symmetric matrix
#====================================================================
# n = uplo = diag = byrow = NULL
unpackTri <- function(A, n = NULL, uplo = NULL, byrow = NULL)
{
   if(is.null(n)){
     if('n' %in% names(attributes(A))){
       n <- attr(A, "n")
     }else{
       stop("Dimension 'n' of the matrix must be provided")
     }
   }
   if(length(A)==(n*(n+1)/2) | length(A)==(n*(n-1)/2)){
     diag <- ifelse(length(A)==n*(n+1)/2, TRUE, FALSE)
   }else{
     stop("Input 'A' must contain either n(n+1)/2 or n(n-1)/2 entries")
   }

   if('uplo' %in% names(attributes(A))){
     uplo <- attr(A, "uplo")
   }else{
     if(is.null(uplo)){
       stop("Parameter 'uplo' must be specified: either 'upper' or 'lower'")
     }
   }
   uplo <- ifelse(uplo=="upper",1L,ifelse(uplo=="lower",2L,NA))

   if('byrow' %in% names(attributes(A))){
     byrow <- attr(A, "byrow")
   }else{
     if(is.null(byrow)){
       stop("Parameter 'byrow' whether data is stored by row must be must be specified")
     }
   }

   #dyn.load("c_triangular.so")
   return(.Call('R_unpacktri', n, A, uplo, diag, byrow))
   #dyn.unload("c_triangular.so")
}
