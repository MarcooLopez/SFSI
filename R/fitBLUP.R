
fitBLUP <- function(y, X = NULL, Z = NULL, K = NULL, U = NULL,
                    d = NULL, varU = NULL, varE = NULL,
                    intercept = TRUE, BLUP = TRUE,
                    method = c("REML","ML"), interval = c(1E-9,1E9),
                    tol = 1E-8, maxiter = 1000, n.regions = 10,
                    verbose = TRUE)
{
  method <- match.arg(method)
  dmin <- .Machine$double.eps
  ratio <- NULL
  # y=yNA[,]; X=X0; K=G0; method="REML"; BLUP=FALSE

  if(length(dim(y)) == 2L){
    y <- as.matrix(y)
  }else{
    y <- matrix(y, ncol=1L)
  }

  n <- nrow(y) # Number of total observations
  q <- ncol(y)
  trn_list <- get_common_trn(y)

  # Track if the training set is the same across all response variables.
  # If not, SVD is performed for each group with common trn set
  commonTRN <- (length(trn_list)==1L)

  if((q > 1L) & !commonTRN){
    if(verbose){
      message(" ",length(trn_list)," different training sets were found for the response variable.")
      message(" Eigenvalue decomposition is applied to each common training set")
    }
    #message(" This process can be time-consuming for large sample size")
  }

  BLUE <- TRUE
  if(is.null(X)){ # Design matrix for fixed effects (including an intercept)
    X <- stats::model.matrix(~1,data=data.frame(rep(1,n)))
    if(!intercept){
      if(verbose){
        message(" No intercept is estimated. Response is assumed to have mean zero")
      }
      BLUE <- FALSE
    }
  }else{
    if(length(dim(X))==2L){
      X <- as.matrix(X)
    }else{
      X <- stats::model.matrix(~X)
      if(ncol(X)>2)  colnames(X)[-1] <- substr(colnames(X)[-1],2,nchar(colnames(X)[-1]))
    }
  }
  stopifnot(n.regions > 0)
  stopifnot(nrow(X) == n)
  p <- ncol(X)
  isREML <- (method=="REML")

  isEigen <- FALSE
  if(is.null(U) | is.null(d))
  {
    G <- U <- d <- NULL
    if(is.null(Z) & is.null(K)){
      ratio <- 0   # Only GLS (only fixed effects)

    }else{
      if(is.null(Z)){  # Case G=K
        if(length(dim(K)) != 2L) stop("Object 'K' must be a matrix")

      }else{
        if(length(dim(Z)) != 2L) stop("Object 'Z' must be a matrix")
        stopifnot(nrow(Z) == n)

        if(is.null(K)){
          G <- tcrossprod(Z)  # G = ZKZ' = ZZ' with K=I
        }else{
          G <- tcrossprod(Z,tcrossprod(Z,K))  # G = ZKZ'
        }
      }
    }

  }else{
    isEigen <- TRUE
    G <- K <- Z <- NULL
    tmp <- unlist(lapply(trn_list,function(x)length(x$trn)))
    if(any(tmp != n)){
       stop("No 'NA' values are allowed when parameters 'U' and 'd' are provided")
    }
  }

  # If varE and varU are provided
  if(is.null(ratio) & !is.null(varU) & !is.null(varE)){
    stopifnot(length(varU) == length(varE))
    ratio <- varU/varE
    if(length(ratio)==1L){
       ratio <- rep(ratio, q)
    }else{
      if(length(ratio) != q){
        stop("Length of 'varU' and 'varE' must be equal to the number of columns of 'y'")
      }
    }
  }

  bounds <- exp(seq(log(interval[1]), log(interval[2]), length=n.regions+1))

  out <- vector("list", q)
  conty <- 0

  # Perform the analysis for all traits
  if(verbose & q>1L){
    pb = utils::txtProgressBar(style=3)
  }
  for(tr in 1:length(trn_list))
  {
    trn <- trn_list[[tr]]$trn
    nTRN <- length(trn)

    if(!isEigen){
      if(is.null(Z) & is.null(K)){
        # SVD of a diagonal matrix
        d <- rep(1, nTRN)
        U <- matrix(0, ncol=nTRN, nrow=nTRN)
        for(i in 1:nTRN){
          U[i,nTRN-i+1] <- 1
        }

      }else{
        if(is.null(Z)){
          tmp <- eigen(K[trn,trn], symmetric=TRUE)
        }else{
          tmp <- eigen(G[trn,trn], symmetric=TRUE)
        }
        d <- tmp$values
        U <- tmp$vectors
        rm(tmp)
      }
    }
    stopifnot(ncol(U) == length(d))

    iy <- trn_list[[tr]]$iy  # columns of y with a common trn set

    for(k in seq_along(iy))
    {
      ind <- iy[k]
      conty <- conty + 1
      ytrn <- as.vector(y[trn, ind])
      stopifnot(all(!is.na(ytrn)))
      ratio0 <- NULL
      if(!is.null(ratio)) ratio0 <- ratio[ind]

      #dyn.load("c_blup.so")
      res <- .Call('R_solve_mixed', n, ratio0, trn-1, ytrn, X, Z, K,
                    U, d, bounds, tol, maxiter, dmin, isREML,
                    isEigen, BLUE, BLUP)
      #dyn.unload("c_blup.so")

      bHat0 <- res[[10]]
      if(BLUE){
        names(bHat0) <- colnames(X)
      }
      if(verbose & q>1L){
        utils::setTxtProgressBar(pb, conty/q)
      }

      out[[ind]] <- list(ind=ind, solution=res[[1]], status=res[[4]],
                         convergence=(res[[5]]>0), nDsmall=res[[6]],
                         varU=res[[7]], varE=res[[8]], h2=res[[9]],
                         bHat=bHat0, yHat=res[[11]], uHat=res[[12]])
    }
  }

  if(verbose & q>1L) {
    close(pb)
  }

  # Checkpoint
  if(any(seq(q) != unlist(lapply(out,function(x) x$ind)) )){
      stop("Some sub-processes failed. Something went wrong during the analysis")
  }

  if(verbose){
    nDsmall <- unlist(lapply(out,function(x)x$nDsmall))
    tmp <- ifelse(any(nDsmall>0),ifelse(sum(nDsmall>0)==1,nDsmall[nDsmall>0],
                  paste(range(nDsmall[nDsmall>0]),collapse="-")),NA)
    if(any(nDsmall>0) | !is.na(tmp)){
      message(ifelse(max(commonTRN),nDsmall,tmp)," eigenvalue(s) are very small.",
              " The corresponding eigenvector(s) were ignored")
    }

    status <- unlist(lapply(out, function(x)x$status))
    status <- status[!is.na(status)]
    if(any(status > 0)){
      for(k in 1:4){
        index <- which(status==k)
        if(length(index)>0){
          msg <- switch(k,
            '1'=paste0("The log Likelihood function is horizontal. No search for ratio varU/varE\n",
                     " was performed and was set to varU/varE=",interval[1]),
            '2'=paste0("Algorithm to find ratio varU/varE did not converge after ",
                     maxiter," iterations.\n Results are doubtful"),
            '3'=paste0("Ratio varU/varE is around the lower bound ",interval[1],
                     "\n Results might be doubtful"),
            '4'=paste0("Ratio varU/varE is around the upper bound ",interval[2],
                     "\n Results might be doubtful"))
          if(q > 1L){
            tmp <- ifelse(length(index)<=15, paste(index,collapse=","),
                paste0(paste(index[1:3],collapse=","),",...,",paste(index[length(index)-(3:1)],collapse=",")))
            msg <- paste(msg, "for",length(index),"trait(s):",tmp)
          }
          message(msg)
        }
      }
    }
  }

  out <- list(varE=unlist(lapply(out,function(x)x$varE)),
              varU=unlist(lapply(out,function(x)x$varU)),
              h2=unlist(lapply(out,function(x)x$h2)),
              b=do.call(cbind,lapply(out,function(x)x$bHat)),
              yHat=do.call(cbind,lapply(out,function(x)x$yHat)),
              u=do.call(cbind,lapply(out,function(x)x$uHat)),
              convergence=unlist(lapply(out,function(x)x$convergence)),
              method=method)

  if(q == 1L){
    out$b <- drop(out$b)
    out$u <- as.vector(out$u)
    out$yHat <- as.vector(out$yHat)
  }

  return(out)
}
