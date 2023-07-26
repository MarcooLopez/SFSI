
SSI <- function(y, X = NULL, b = NULL, Z = NULL, K,
                trn_tst = NULL, varU = NULL, varE = NULL,
                intercept = TRUE, alpha = 1, lambda = NULL,
                nlambda = 100, lambda.min = .Machine$double.eps^0.5,
                common.lambda = TRUE, subset = NULL, tol = 1E-4,
                maxiter = 500, method = c("REML","ML"),
                name = NULL, save.at = NULL, mc.cores = 1L,
                precision.format = c("single","double"),
                verbose = TRUE)
{
  method <- match.arg(method)
  precision.format <- match.arg(precision.format)
  # y=Y0[,k]; K=G0; method="REML"; precision.format='double'
  # varU=varU0; varE=varE0; subset=c(2,5); save.at=prefix
  if(length(dim(y)) == 2L){
    y <- as.matrix(y)
  }else{
    y <- matrix(y, ncol=1L)
  }
  dimnames(y) <- NULL
  n <- nrow(y)
  ntraits <- ncol(y)   # Number of traits

  if(is.null(trn_tst)){
    trn_tst <- matrix(TRUE, nrow=n, ncol=ntraits)
  }else{
    if(length(dim(trn_tst)) == 2L){
      trn_tst <- as.matrix(trn_tst)
    }else{
      trn_tst <- matrix(trn_tst, ncol=1L)
    }
    stopifnot(length(trn_tst) == (n*ntraits))
    if(storage.mode(trn_tst) %in% c("double","integer")){
      if(any( !unique(as.vector(trn_tst)) %in% c(0,1) )){
        stop("Input 'trn_tst' must contain only 0, 1, or NA")
      }
      trn_tst <- (trn_tst == 1)
    }
    stopifnot(storage.mode(trn_tst) == "logical")
  }
  dimnames(trn_tst) <- NULL

  MAP <- data.frame(index=seq(ntraits*n),reshape2::melt(trn_tst))
  MAP$value <- c("tst","trn")[as.numeric(MAP$value)+1]
  colnames(MAP) <- c("index","i","j","set")
  MAP$index_set <- NA

  tst.index <- data.frame(i=rep(seq(n),ntraits), j=rep(seq(ntraits),each=n))
  trn <- which(as.vector(trn_tst))
  MAP[trn,"index_set"] <- seq_along(trn)

  if(any(!trn_tst)){  # If any testing set
    tst <- which(as.vector(!trn_tst))
    MAP[tst,"index_set"] <- seq_along(tst)
    tst.index <- data.frame(tst.index[tst,], index=1:length(tst))
  }else{
    if(verbose){
      message(" No testing set was found. The SSI will be fitted to the entire data set")
    }
    tst.index <- data.frame(tst.index, index=seq(nrow(tst.index)))
    MAP[trn,"set"] <- "trn_tst"
    tst <- trn[]
  }
  nTST <- length(tst)
  rownames(tst.index) <- NULL

  storage.mode(K) <- "double"

  # Design matrix for fixed effects
  BLUE <- TRUE
  if(is.null(X)){   # Only an intercept
    X <- model.matrix(~1, data=data.frame(rep(1,n)))
    if(!intercept){
      if(verbose){
        message(" No intercept is estimated. Response is assumed to have mean zero")
      }
      BLUE <- FALSE
    }

  }else{
    if(length(dim(X)) == 2L){
      X <- as.matrix(X)
    }else{
      X <- stats::model.matrix(~X)
      if(ncol(X) > 2L)  colnames(X)[-1] <- substr(colnames(X)[-1],2,nchar(colnames(X)[-1]))
    }
  }

  if(!is.null(Z)){
    if(length(dim(Z)) != 2L){
       stop("Object 'Z' must be a matrix with nrow(Z)=n and ncol(Z)=nrow(K)")
    }
    K <- tcrossprod(Z, tcrossprod(Z,K))   # Z%*%K%*%t(Z)
  }

  if(length(dim(K)) != 2L | (length(K) != n^2)){
    stop("Product Z K Z' must be a squared matrix with number of rows\n",
         "      (and columns) equal to the number of elements in 'y'")
  }

  labels <- NULL
  if(has_names(K)){
    labels <- rownames(K)
  }

  if(ntraits == 1L)  # Single-trait case
  {
    if(is.null(varU) | is.null(varE))
    {
      # Get variance components and estimate fixed effects
      res <- fitBLUP(y[trn,], X=X[trn, ,drop=FALSE], K=K[trn,trn],
                     intercept=intercept, method=method,
                     BLUP=FALSE, verbose=FALSE)
      if(res$convergence){
        varU <- res$varU
        varE <- res$varE
        b <- res$b

        if(verbose){
          message(" 'varU' and 'varE' were estimated using a BLUP model")
        }
      }else{
        stop("Convergence was not reached in the 'GEMMA' algorithm.\n",
             "       Please provide variance components' estimates")
      }
    }else{   # Only estimate fixed effects as GLS
      if(is.null(b) & BLUE){
        b <- fitBLUP(y[trn,], X=X[trn, ,drop=FALSE], K=K[trn,trn],
                     BLUP=FALSE, varU=varU, varE=varE, verbose=FALSE)$b
      }
    }
    if(BLUE){
      b <- t(t(b))
    }

  }else{   # Multitrait case
    trn0 <- which(apply(trn_tst,1,all))  # Common TRN set
    if(is.null(varU) | is.null(varE)){
      if(length(trn0) == 0){
         stop("No common training to all response variables was found.\n",
              "       Please provide 'varU' and 'varE' estimates")
      }
      res <- getGenCov(y[trn0,], X=X[trn0, ,drop=FALSE], K=K[trn0,trn0],
                       pairwise=TRUE, intercept=intercept, verbose=FALSE)
      varU <- res$varU
      varE <- res$varE
      b <- res$b

      if(verbose){
        message(" 'varU' and 'varE' were pairwise estimated using BLUP",
                " with n=",length(trn0)," common records")
      }

    }else{
      if(is.null(b) & BLUE){
        if(length(trn0) == 0){
           stop("No common training for all response variables was found.\n",
                "       Please provide an estimate for 'b'")
        }
        b <- fitBLUP(y[trn0,], X=X[trn0, ,drop=FALSE], K=K[trn0,trn0], BLUP=FALSE,
                     varU=diag(varU), varE=diag(varE), verbose=FALSE)$b
      }
    }

    if(length(dim(varU)) != 2L | length(dim(varE)) != 2L){
      stop("'varU' and 'varE' must be matrices of dimension equal to ncol(y)=",ncol(y))
    }
  }

  h2 <- varU/(varU + varE)
  if(any(diag(t(h2)) < 0.001) & verbose){
     message(" The 'heritability' is too small. Results may be affected")
  }

  # Adjusted training phenotypes
  Xb <- c()
  if(BLUE){
    if(nrow(b) != ncol(X)){
       stop("The length of 'b' must be the same as the number of columns of 'X'")
    }
    for(k in 1:ntraits){
      Xb <- c(Xb,as.vector(X%*%b[,k]))
    }
    yTRN <- matrix(y[trn]-Xb[trn], nrow=1)
  }else{
    yTRN <- matrix(y[trn], nrow=1)
  }

  # Getting K <- varU*G and H <- varU*G + varE*I
  K <- kronecker(varU, K)
  if(ntraits == 1L){
    H <- penalize_cov(K, lambda=varE, inplace=FALSE)
  }else{
    H <- K + kronecker(varE, diag(n))
  }
  H <- H[trn, trn]
  K <- K[trn, tst, drop=FALSE]

  if(is.null(lambda)){
    if(common.lambda){
      sdx <-  sqrt(diag(H))
      Cmax <- max(abs(sweep(K, 1L, sdx, FUN="/"))/alpha)
      Cmax <- ifelse(alpha > .Machine$double.eps, Cmax, 5)
      lambda <- exp(seq(log(Cmax),log(lambda.min),length=nlambda))
    }
  }else{
    if(length(dim(lambda))==2L){
      if(ncol(lambda)>1L & ncol(lambda)<nTST){
         stop("Number of columns of 'lambda' must be equal to 'length(tst)'")
      }
    }else{
      lambda <- matrix(lambda, ncol=1)
    }
    if(any(apply(lambda, 2L, function(x) any(diff(x)>0)))){
        stop("Object 'lambda' must be a matrix of decreasing numbers")
    }
  }

  name <- ifelse(is.null(name),"SSI",name)

  # Split the testing set into subsets. Only the subset provided will be fitted
  if(is.null(subset)){
    subset_index <- fileID <- NULL
    tmp <- ""
  }else{
     if(!is.numeric(subset) & length(subset) != 2L){
       stop("Object 'subset' must be a 2-elements vector")
     }
     sets <- sort(rep(1:subset[2],ceiling(nTST/subset[2]))[1:nTST])
     subset_index <- do.call(rbind,lapply(1:subset[2], function(k){
       tmp <- which(sets == k)
       data.frame(subset=k,n=length(tmp),from=min(tmp),to=max(tmp))
     }))
     index <- which(sets == subset[1])
     fileID <- index[]
     tmp <- paste0(" of ",length(tst))
     tst <- tst[index]
     K <- K[,index,drop=FALSE]
  }

  if(verbose){
    message(" Fitting ",ifelse(ntraits>1,paste("Multi-trait SSI for",ntraits,"traits"),"SSI"),
            " with nTST=",length(tst),tmp," and nTRN=",length(trn))
  }

  # If 'save.at' is not NULL
  if(!is.null(save.at)){
    prefix <- basename(tempfile(pattern=""))
    if(is.null(subset)){
      outfile <- paste0(save.at,"output.RData")
    }else{
      prefix <- paste0("subset_",subset[1],"_of_",subset[2],"_")
      outfile <- paste0(save.at,prefix,"output.RData")
    }
    unlink(outfile)
  }

  out <- solveEN(Sigma=H, Gamma=K, alpha=alpha, lambda=lambda,
                 nlambda=nlambda, lambda.min=lambda.min,
                 common.lambda=common.lambda, tol=tol,
                 maxiter=maxiter, save.at=save.at, fileID=fileID,
                 verbose=verbose, mc.cores=mc.cores,
                 precision.format=precision.format)

  if(length(tst) == 1L){
    nsup0 <- matrix(out$nsup, nrow=1)
    lambda0 <- matrix(out$lambda, nrow=1)
  }else{
    nsup0 <- do.call(rbind, out$nsup)
    lambda0 <- do.call(rbind, out$lambda)
  }

  u <- fitted.LASSO(out, yTRN)
  dimnames(u) <- list(tst, paste0("SSI.",1:ncol(u)))
  if(BLUE){
    yHat <- sweep(u, 1L, Xb[tst], FUN="+")
  }else{
    yHat <- u[]
  }

  out <- list(n=n, p=out$p, q=out$q, ntraits=ntraits,
              labels=labels, name=name, nlambda=nlambda,
              y=y, Xb=Xb, u=u, yHat=yHat,
              b=b, varU=varU, varE=varE, h2=h2,
              trn=trn, tst=tst, MAP=MAP, tst.index=tst.index, alpha=alpha,
              nsup = nsup0, lambda = lambda0,
              beta = out$beta,
              file_beta=out$file_beta,
              fileID=out$fileID,
              precision.format=precision.format
            )

  if(!is.null(subset)){
    out$subset <- subset_index
  }

  class(out) <- c("SSI")

  # Save outputs if 'save.at' is not NULL
  if(!is.null(save.at)){
    save(out, file=outfile)
    #out <- NULL
  }

  return(out)
}
