
# index.row = NULL; index.col = NULL; verbose = TRUE
# Read a file as binary
readBinary <- function(file = paste0(tempdir(),"/file.bin"),
                       index.row = NULL, index.col = NULL,
                       verbose = TRUE)
{
  if(!file.exists(file)){
    stop("File '",file,"' does not exist")
  }

  # Read lines
  #dyn.load("c_functions.so")
  X <- .Call("R_readBinFile", file,
                          as.integer(index.row),
                          as.integer(index.col))
  #dyn.unload("c_functions.so")

  n <- X[[1]]
  p <- X[[2]]
  vartype <- X[[3]]
  size <- X[[4]]
  nError <- X[[5]]
  #type <- ifelse(vartype==1L,"integer",ifelse(vartype==2L,"logical","double"))
  if(nError == 0L){
    X <- X[[6]]
    type <- storage.mode(X)
    if(verbose){
      tmp <- c(Gb=1E9,Mb=1E6,Kb=1E3,b=1E0)
      sz <- object.size(X)/tmp[min(which(object.size(X)/tmp>1))]
      message("Loaded file '",file,"'")
      message("  nrow=",n,", ncol=",p,", type=",type,", size=",size," bytes, object.size=",round(sz,2),names(sz))
    }

  }else{
    X <- NULL
  }
  return(X)
}
