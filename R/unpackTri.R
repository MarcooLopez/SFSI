
#====================================================================
# Convert a packed upper/lower triangle matrix into a symmetric matrix
#====================================================================
# n = uplo = diag = byrow = NULL
unpackTri <- function(A, rows = NULL, cols = NULL,
                      drop = TRUE, verbose = TRUE)
{
   if(length(dim(A)) == 2L){
     stop("Input 'A' must be a vector with attributes")
   }

   tmp <- info_matrix(A, check=TRUE)
   n <- tmp$nrows
   uplo <- tmp$type
   diag <- tmp$diag
   byrow <- tmp$byrow

   if(verbose){
     message(" Unpacking a ",n,"x",n," matrix packed ",ifelse(byrow,"row","column"),
             "-wise with the '",uplo,"' triangular data")

     if(!diag){
       message(" Data in 'A' does not contain a diagonal, it is assumed to be 1")
     }
   }

   if(!is.null(rows)){
     stopifnot(all(!is.na(rows)))
     if(any(rows<1) | any(rows>n)){
       stop("Input 'rows' must be integers between 1 and dim(A)=",n)
     }
   }
   if(!is.null(cols)){
     stopifnot(all(!is.na(cols)))
     if(any(cols<1) | any(cols>n)){
       stop("Input 'cols' must be integers between 1 and dim(A)=",n)
     }
   }

   #dyn.load("c_triangular.so")
   return(.Call('R_unpack_tri', n, A, rows-1, cols-1, drop, uplo, diag, byrow))
   #dyn.unload("c_triangular.so")
}
