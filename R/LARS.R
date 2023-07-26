
LARS <- function(Sigma, Gamma, method=c("LAR","LASSO"),
                 nsup.max = NULL, eps = .Machine$double.eps*100,
                 scale = TRUE, mc.cores = 1L, save.at = NULL,
                 precision.format = c("double","single"),
                 fileID = NULL, verbose = FALSE)
{
  precision.format <- match.arg(precision.format)
  method <- match.arg(method)

  if(length(dim(Gamma)) == 2L){
    Gamma <- as.matrix(Gamma)
  }else{
    Gamma <- matrix(Gamma, ncol=1L)
  }
  dimnames(Sigma) <- NULL
  p <- nrow(Gamma)
  q <- ncol(Gamma)

  if((sum(dim(Sigma))/2)^2 != p^2){
    stop("Object 'Sigma' must be a p*p squared matrix where p=nrow(Gamma)")
  }

  storage.mode(Sigma) <- "double"
  storage.mode(Gamma) <- "double"

  if(scale){
    sdx <-  sqrt(diag(Sigma))
    cov2cor2(Sigma, inplace=TRUE)     # Equal to Sigma=cov2cor(Sigma) but faster
    Gamma <- sweep(Gamma, 1L, sdx, FUN="/")
  }else{
    sdx <- rep(1, p)
  }

  nsup.max <- ifelse(is.null(nsup.max), p, nsup.max)
  flagsave <- as.logical(!is.null(save.at))
  isLASSO <- as.logical(method=="LASSO")
  verbose2 <- as.logical(verbose & q==1L)
  mc.cores <- ifelse(q==1L & mc.cores>1L, 1L, mc.cores)
  doubleprecision <- as.logical(precision.format=="double")

  compApply <- function(ind)
  {
    rhs <- as.vector(Gamma[,ind])

<<<<<<< HEAD
    if(flagsave){
      filename <- paste0(file_beta,fileID[ind],".bin")
=======
    im <- inactive <- seq(p)
    ignores <- NULL
    if(is.null(dfmax))  dfmax <- p
    beta <- Matrix::Matrix(0, nrow=p, ncol=dfmax*8)
    lambda <- double(dfmax*8)
    active <- NULL
    drops <- FALSE
    Sign <- NULL
    R <- NULL
    k <- 0
    time <- proc.time()[3]
    while((length(active) < dfmax) & (length(active) < (p-length(ignores))))
    {
      covar <- rhs[inactive]
      Cmax <- max(abs(covar))
      if(Cmax < eps*100){
        if(verbose & q==1L) cat(" Max |corr| = 0; exiting...\n")
        break
      }
      k <- k+1
      lambda[k] <- ifelse(isfloat, float::dbl(Cmax), Cmax)
      if(!any(drops))
      {
        # new <- abs(covar) >= Cmax-eps
        new <- abs(abs(covar) - Cmax) <= eps
        covar <- covar[!new]
        new <- inactive[new]
        for(inew in new)
        {
          R <- upDateR(Sigma[inew,inew],R,as.vector(Sigma[inew,active]),eps=eps)
          if(attr(R,"rank")==length(active))
          {
            nR <- seq(length(active))
            R <- R[nR,nR,drop=FALSE]
            attr(R,"rank") <- length(active)
            ignores <- c(ignores,inew)
            if(verbose & q==1L){
              cat(" LARS Step ",k,":\t Variable ", inew,"\tcollinear; dropped for good\n",sep="")
            }
          }else{
            active <- c(active,inew)
            Sign <- c(Sign,sign(rhs[inew]))
            if(verbose & q==1L){
              cat("--------------------------------------------------------------\n")
              tmp <- proc.time()[3]
              cat(paste(textPrint,"=",c(sprintf("%*d",Sp,k),sprintf('%.*f',4,tmp-time),
                  sprintf("%*d",Sp,inew)))," added\n",sep="")
              time <- tmp
            }
          }
        }
      }
      Gi1 <- float::backsolve(R,backsolvet(R,Sign))
      A <- 1/sqrt(sum(Gi1*Sign))
      w <- A*Gi1

      if((length(active) >= (p-length(ignores)))){
          gamhat <- Cmax/A
      }else{
        # a <- drop(w %*% Sigma[active, -c(active,ignores),drop=FALSE])
        a <- float::crossprod(Sigma[active, -c(active,ignores),drop=FALSE],w)[,1]
        gam <- c((Cmax-covar)/(A-a),(Cmax+covar)/(A+a))
        gamhat <- min(gam[gam > eps],Cmax/A)
      }
      if(method == "LAR-LASSO")
      {
        dropid <- NULL
        b1 <- beta[active,k]
        z1 <- -b1/w
        zmin <- min(z1[z1 > eps],gamhat)
        if(zmin < gamhat){
            gamhat <- zmin
            drops <- z1 == zmin
        }else drops <- FALSE
      }
      beta[,k+1] <- beta[,k]
      beta[active,k+1] <- beta[active,k+1] + gamhat*w
      # rhs <- rhs - gamhat*Sigma[,active,drop=FALSE]%*%w
      rhs <- rhs - gamhat*float::crossprod(Sigma[active, ,drop=FALSE],w)[,1]
      if(method == "LAR-LASSO" && any(drops))
      {
        dropid <- seq(drops)[drops]
        for(id in rev(dropid))
        {
          if(verbose & q==1L){
            cat("--------------------------------------------------------------\n")
            tmp <- proc.time()[3]
            cat(paste(textPrint,"=",c(sprintf("%*d",Sp,k+1),sprintf('%.*f',4,tmp-time),
              sprintf("%*d",Sp,active[id])))," dropped\n",sep="")
            time <- tmp
          }
          R <- downDateR(R,id)
        }
        dropid <- active[drops]
        beta[dropid,k+1] <- 0
        active <- active[!drops]
        Sign <- Sign[!drops]
      }
      inactive <- im[-c(active, ignores)]
    }
    beta <- beta[,seq(k+1) ,drop = FALSE]
    lambda  <-  c(lambda[seq(k)],0)
    if(scale){
       beta <- sweep(beta, 1L, sdx, FUN="/")
    }
    df <- do.call(c,lapply(1:ncol(beta),function(i) sum(abs(beta[,i])>0)))

    if(!is.null(X)){
      yHat <- X%*%as.matrix(beta)
>>>>>>> f172c692471fce78aa4eef4bd2994fe6869e1be8
    }else{
      filename <- NULL
    }

    #dyn.load("c_lasso.so")
    res <- .Call("R_lars", Sigma, rhs, eps, nsup.max,
                 scale, sdx, isLASSO, filename,
                 doubleprecision, verbose2)
    #dyn.unload("c_lasso.so")

    if(verbose & q>1L){
      cat(1,file=con,append=TRUE)
      utils::setTxtProgressBar(pb, nchar(scan(con,what="character",quiet=TRUE))/q)
    }

    list(ind=ind, beta=res[[1]], lambda=res[[2]], nsup=res[[3]])
  }

  tmpdir0 <- tempdir()

  file_beta <- NULL
  if(flagsave){
    stopifnot(is.character(save.at))
    file_beta <- paste0(save.at,"beta_i_")

    if(!file.exists(dirname(file_beta))){
      dir.create(dirname(file_beta),recursive=TRUE)
    }

    if(is.null(fileID)){
      fileID <- seq(q)
    }else{
      stopifnot(length(fileID)==q)
    }
  }

  # Run the analysis for 1:ncol(Gamma)
  if(verbose & q>1L){
     pb = utils::txtProgressBar(style=3)
     con <- tempfile(tmpdir=tmpdir0)
  }
  if(mc.cores == 1L){
    out = lapply(X=seq(q),FUN=compApply)
  }else{
    out = parallel::mclapply(X=seq(q),FUN=compApply,mc.cores=mc.cores)
  }
  if(verbose & q>1L) {
    close(pb); unlink(con)
  }

  # Checkpoint
  if(any(seq(q) != unlist(lapply(out,function(x) x$ind)) )){
    stop("Some sub-processes failed. Something went wrong during the analysis.")
  }

  out <- list(p=p, q=q, method=method,
              nlambda=unlist(lapply(out, function(x)length(x$lambda))),
              lambda = lapply(out, function(x)x$lambda),
              nsup = lapply(out, function(x)x$nsup),
              beta = lapply(out, function(x)x$beta)
            )

  if(q == 1L){
    out$nsup <- out$nsup[[1]]
    out$lambda <- out$lambda[[1]]
    out$beta <- as.matrix(out$beta[[1]])
  }

  if(flagsave){
    out$file_beta <- gsub("i_[0-9]+.bin$","i_\\*.bin",
                            normalizePath(paste0(file_beta,fileID[1],".bin")))
    out$fileID <- fileID
    out$beta <- NULL
  }

  class(out) <- "LASSO"
  return(out)
}
