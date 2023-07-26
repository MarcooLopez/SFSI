
#====================================================================
# Extract the lower triangle section from a squared matrix
#====================================================================
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

    }else{
      #if(storage.mode(A) != "double") storage.mode(A) <- "double"
      out <- .Call('R_lowertri', n, A, diag, byrow)
    }
    #dyn.unload("c_triangular.so")

    attr(out, "uplo") <- "lower"
    attr(out, "n") <- n
    attr(out, "include.diag") <- diag
    attr(out, "byrow") <- byrow

    return(out)
}
