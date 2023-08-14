
#====================================================================
#
#====================================================================
# n = uplo = diag = byrow = NULL
subsetTri <- function(A, rowcols, verbose = TRUE)
{
   if(length(dim(A)) == 2L){
     stop("Input 'A' must be a vector")
   }

   tmp <- info_matrix(A, check=TRUE)
   n <- tmp$nrows
   uplo <- tmp$type
   diag <- tmp$diag
   byrow <- tmp$byrow

   if(verbose){
     message(" Subsetting a ",n,"x",n," matrix packed ",ifelse(byrow,"row","column"),
             "-wise with the '",uplo,"' triangular data")
   }

   stopifnot(all(!is.na(rowcols)))
   if(any(rowcols<1) | any(rowcols>n)){
    stop("Input 'rowcols' must be integers between 1 and dim(A)=",n)
   }

   #dyn.load("c_triangular.so")
   return(.Call('R_subset_tri',n, A, rowcols-1, uplo, diag, byrow))
   #dyn.unload("c_triangular.so")
}
