
# Save a file as binary

saveBinary <- function(X, file = paste0(tempdir(),"/file.bin"),
                       precision.format = c("double","single"),
                       verbose = TRUE)
{
  precision.format <- match.arg(precision.format)

  dm <- dim(X)
  if(length(dm) == 2L){
    X <- as.matrix(X)
  }else{
    dm <- c(length(X),0)
    #stop("Object 'X' must be a matrix")
  }

  type <- storage.mode(X)
  if(!type %in% c("logical","double","integer")){
     stop("'storage.mode(X)' must be either 'integer', 'logical', or 'double'")
  }

  unlink(file)

  doubleprecision <- as.logical(precision.format=="double")

  #dyn.load("c_functions.so")
  out <- .Call('R_writeBinFile', file, dm[1], dm[2], X, doubleprecision)
  #dyn.unload("c_functions.so")

  if(verbose){
    size <- out[[4]]
    tmp <- c(Gb=1E9,Mb=1E6,Kb=1E3,b=1E0)
    sz <- file.size(file)/tmp[min(which(file.size(file)/tmp>1))]
    message(" Saved file '",file,"'")
    message("  nrow=",out[[1]],", ncol=",out[[2]],", type=",type,", size=",size," bytes, file.size=",round(sz,2),names(sz))
  }
}
