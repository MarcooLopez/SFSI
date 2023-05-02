#====================================================================
# Prunne elements based on pairwise similarity R (correlation, R^2, etc)
#====================================================================

Prune <- function(R, threshold=0.95, D=NULL, d.max=NULL, MAF=NULL, verbose=FALSE)
{
    nc <- nchar(ncol(R))
    PrFn <- function(x, nc=3)sprintf(paste0('%',nc,'d'),x)

    pruneIn <- c()
    pruneOut <- c()

    if(!is.null(d.max) & is.null(D)){
      message("Maximum distance 'd.max' is ignored when matrix 'D=NULL'")
    }
    if(is.null(d.max) & !is.null(D)){
      stop("A 'd.max' value must be provided along with the distance matrix 'D'")
    }

    if(is.null(D)){
      A <- (R > threshold)
    }else{
      A <- (R > threshold) & (D <= d.max)
    }
    diag(A) <- FALSE
    colnames(A) <- rownames(A) <- NULL

    if(verbose){
      tmp <- ifelse(is.null(d.max),"",paste0(" within ",d.max," bp"))
      cat("Pruning ", ncol(R)," subjects",tmp," ...\n",sep="")
    }

    if(is.null(colnames(R))){
      ID <- 1:ncol(R)
    }else{
      ID <- colnames(R)
    }

    if(!is.null(MAF)){
      stopifnot(length(MAF)==nrow(R))
      names(MAF) <- ID
    }

    CON <- rowSums(A)
    nConn <- CON[]
    remain <- 1:ncol(R)
    mc <- nchar(max(CON))

    cont <- 0
    if(any(CON==0)){
        tmp <- which(CON==0)
        pruneIn <- ID[tmp]
        remain <- remain[-tmp]
        CON <- CON[-tmp]
        ID <- ID[-tmp]
        if(verbose){
          cont <- cont + 1
          cat("--------------------------------------------------------\n")
          cat(" S",PrFn(cont),". nConn=",PrFn(0,mc), ". In: n=",PrFn(length(tmp),nc),
              ". Out: n=",PrFn(0,mc),". Remain: n=",PrFn(length(remain),nc),"\n",sep="")
        }
    }

    CON0 <- rep(0,length(CON))
    cond <- any(CON>0)
    while(cond){
      #CON <- CON-CON0
      if(any(CON > 0)){
        tmp <- which(CON==max(CON))
        if(is.null(MAF)){
          keep <- tmp[1]
        }else{
          keep <- tmp[which.max(MAF[ID[tmp]])]
        }
        pruneIn <- c(pruneIn, ID[keep])
        remove <- which(A[remain[keep],remain])
        # range(R[ID[c(keep)],ID[remove]])
        pruneOut <- c(pruneOut, ID[remove])

        if(verbose){
          cont <- cont + 1
          cat("--------------------------------------------------------\n")
          cat(" S",PrFn(cont),". nConn=",PrFn(CON[keep],mc),". In: i=",PrFn(keep,nc),
                ". Out: n=",PrFn(length(remove),mc),". Remain: n=",
                PrFn(ncol(R)-length(pruneOut)-length(pruneIn),nc),"\n",sep="")
        }

        CON0 <- rowSums(A[remain[-remove], remain[remove], drop=FALSE])
        remain <- remain[-remove]
        CON <- CON[-remove]
        ID <- ID[-remove]
        CON <- CON-CON0
        cond=(length(remain) > 1)
      }else{
         cond <- FALSE
         tmp <- which(!ID %in% pruneIn)
         if(length(tmp)>0){
          pruneIn <- c(pruneIn, ID[tmp])
          if(verbose){
            cont <- cont + 1
            cat("--------------------------------------------------------\n")
            cat(" S",PrFn(cont),". nConn=",PrFn(max(CON),mc),
                ". In: n=",PrFn(length(tmp),nc),". Out: n=",PrFn(0,mc),". Remain: n=",
                PrFn(ncol(R)-length(pruneOut)-length(pruneIn),nc),"\n",sep="")
          }
        }
      }
    }
    return(list(nConn=nConn, prune.in=pruneIn, prune.out=pruneOut))
}
