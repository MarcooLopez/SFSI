
# Covariance matrix to distance matrix

cov2dist <- function(A, a = 1, inplace = FALSE)
{
   tmp <- info_matrix(A, check=TRUE)
   if(tmp$nrows != tmp$ncols){
     stop("Input 'A' must be squared matrix")
   }
   n <- tmp$nrows
   type <- tmp$type
   diag <- tmp$diag
   byrow <- tmp$byrow

   type <- ifelse(type=="full",0L,ifelse(type=="upper",1L,2L))

   #isBigMatrix <- bigmemory::is.big.matrix(A)
   isBigMatrix <- FALSE

   if(isBigMatrix){
     message(" Routine 'cov2dist' is not implemented yet for 'bigmatrix'")
     #stopifnot(bigmemory::typeof(A) == "double")

   }else{
     #dyn.load("c_cov2dist.so")
     if(inplace){
       tmp <- .Call('R_cov2dist', n, a, A, type, byrow)
       #return(tmp)
     }else{
       out <- A[]
       tmp <- .Call('R_cov2dist', n, a, out, type, byrow)
       return(out)
     }
     #dyn.unload("c_cov2dist.so")
  }
}
