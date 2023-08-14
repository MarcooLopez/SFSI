
#====================================================================
# Extract the lower triangle section from a squared matrix
#====================================================================
# diag = TRUE; byrow = FALSE
lowerTri <- function(A, diag = TRUE, byrow = FALSE)
{
    dm <- dim(A)
    if((sum(dm)/2)^2 != length(A)){
       stop("Input 'A' must be a squared matrix")
    }
    n <- dm[1]

    #isBigMatrix <- bigmemory::is.big.matrix(A)
    isBigMatrix <- FALSE
    #dyn.load("c_triangular.so")
    if(isBigMatrix){
      message(" Routine 'lowerTri' is not implemented yet for 'bigmatrix'")
      return(NULL)

    }else{
      return(.Call('R_lower_tri', n, A, diag, byrow))
    }
    #dyn.unload("c_triangular.so")

}
