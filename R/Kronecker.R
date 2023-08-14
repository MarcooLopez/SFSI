
#====================================================================
# Kronecker product between matrices A and B
#====================================================================
# rows <- cols <- NULL; verbose=TRUE
Kronecker <- function(A, B, rows = NULL, cols = NULL,
                      drop = TRUE, unpack = TRUE, verbose = TRUE)
{
  eps <- .Machine$double.eps

  if((length(dim(A))!=2) & !isTri(A)){
    A <- matrix(A, ncol=1L)
  }
  if((length(dim(B))!=2) & !isTri(B)){
    B <- matrix(B, ncol=1L)
  }

  info <- info_matrix(A, check=FALSE)
  if(any(unlist(lapply(info,is.null)))){
    stop("'A' must be either a full matrix or vector containing the triangular part")
  }
  nrow1 <- info$nrows
  ncol1 <- info$ncols
  type1 <- info$type
  diag1 <- info$diag
  byrow1 <- info$byrow

  info <- info_matrix(B, check=FALSE)
  if(any(unlist(lapply(info,is.null)))){
    stop("'B' must be either a full matrix or vector containing the triangular part")
  }
  nrow2 <- info$nrows
  ncol2 <- info$ncols
  type2 <- info$type
  diag2 <- info$diag
  byrow2 <- info$byrow

  if(verbose & type1!="full"){
    message(" 'A' is a ",nrow1,"x",nrow1," matrix packed ",ifelse(byrow1,"row","column"),
            "-wise with the '",type1,"' triangular data")
    if(!diag1){
      message(" Data in 'A' does not contain a diagonal, it is assumed to be 1")
    }
  }
  if(verbose & type2!="full"){
    message(" 'B' is a ",nrow2,"x",nrow2," matrix packed ",ifelse(byrow2,"row","column"),
            "-wise with the '",type2,"' triangular data")
    if(!diag2){
      message(" Data in 'B' does not contain a diagonal, it is assumed to be 1")
    }
  }

  if(!is.null(rows)){
    stopifnot(all(!is.na(rows)))
    if(any(rows<1) | any(rows>(nrow1*nrow2))){
      stop("Input 'rows' must be integers between 1 and nrow(A)*nrow(B)=",nrow1*nrow2)
    }
  }
  if(!is.null(cols)){
    stopifnot(all(!is.na(cols)))
    if(any(cols<1) | any(cols>(ncol1*ncol2))){
      stop("Input 'cols' must be integers between 1 and ncol(A)*ncol(B)=",ncol1*ncol2)
    }
  }

  # Check if triangular matrices are from the same type
  if((type1!="full") & (type2!="full")){
    if((type1!=type2) | (byrow1!=byrow2)){
      stop("Both 'A' and 'B' must be of the same triangular type")
    }

    if(diag1 != diag2){
      stop("Both 'A' and 'B' must either have or not have diagonal data")
    }
  }

  #isBigMatrix <- bigmemory::is.big.matrix(A)
  isBigMatrix <- FALSE

  if(isBigMatrix){
    message(" Routine 'Kronecker' is not implemented yet for 'bigmatrix'")
    #stopifnot(bigmemory::typeof(A) == "double")

  }else{
    #dyn.load("c_kronecker.so")
    return(.Call('R_kronecker',nrow1, ncol1, A, nrow2, ncol2, B, rows-1, cols-1,
                               drop, type1, byrow1, diag1, type2, byrow2, diag2,
                               unpack, verbose))
    #dyn.unload("c_kronecker.so")
  }

}
