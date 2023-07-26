
solveEN <- function(Sigma, Gamma, alpha = 1, lambda = NULL,
                    nlambda = 100, lambda.min = .Machine$double.eps^0.5,
                    lambda.max = NULL, common.lambda = TRUE,
                    nsup.max = NULL, scale = TRUE, tol = 1E-5,
                    maxiter = 1000, mc.cores = 1L, save.at = NULL,
                    precision.format = c("double","single"),
                    fileID = NULL, verbose = FALSE)
{
    precision.format <- match.arg(precision.format)
    alpha <- as.numeric(alpha)
    scale <- as.logical(scale)
    tol <- as.numeric(tol)
    maxiter <- as.integer(maxiter)

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

    if(alpha<0 | alpha>1){ stop("Parameter 'alpha' must be a number between 0 and 1")}
    stopifnot(maxiter>0)
    if(tol<.Machine$double.eps){
      stop("Parameter 'tol' must be > 0")
    }

    storage.mode(Sigma) <- "double"
    storage.mode(Gamma) <- "double"
    nsup.max <- ifelse(is.null(nsup.max), p, as.integer(nsup.max))

    if(scale){
      sdx <-  sqrt(diag(Sigma))
      cov2cor2(Sigma, inplace=TRUE)     # Equal to Sigma=cov2cor(Sigma) but faster
      Gamma <- sweep(Gamma, 1L, sdx, FUN="/")
    }else{
      sdx <- rep(1, p)
    }

    if(is.null(lambda)){
      if(!is.null(lambda.max)){
        if(length(lambda.max)>1L & length(lambda.max)!=q){
            stop("Length of 'lambda.max' must be the same as ncol(Gamma)=",q)
        }
      }
      if(common.lambda){
        if(is.null(lambda.max)){
         lambda.max <- ifelse(alpha > .Machine$double.eps,max(abs(Gamma)/alpha),5)
        }
        lambda <- matrix(exp(seq(log(lambda.max),log(lambda.min),length=nlambda)), ncol=1)
      }else{
        lambda <- do.call(cbind,lapply(1:q,function(k){
          if(is.null(lambda.max)){
            lambda.maxi <- ifelse(alpha > .Machine$double.eps,max(abs(Gamma[,k])/alpha),5)
          }else{
            lambda.maxi <- lambda.max[k]
          }
          exp(seq(log(lambda.maxi),log(lambda.min),length=nlambda))
        }))
      }
    }else{
      if(length(dim(lambda))==2L){
        if(ncol(lambda) > 1 & ncol(lambda) < q){
          stop("Number of columns of 'lambda' must be equal to the number of columns in 'Gamma'")
        }
      }else{
        lambda <- matrix(lambda, ncol=1)
      }
      if(any(apply(lambda, 2L, function(x) any(diff(x) > 0)))){
        stop("Object 'lambda' must be a matrix of decreasing numbers")
      }
    }
    storage.mode(lambda) <- "double"
    nlambda <- nrow(lambda)

    flagsave <- as.logical(!is.null(save.at))
    verbose2 <- as.logical(verbose & q==1L)
    mc.cores <- ifelse(q==1L & mc.cores>1L, 1L, mc.cores)
    doubleprecision <- as.logical(precision.format=="double")

    compApply <- function(ind)
    {
      if(ncol(lambda)>1L){
        lambda0 <- lambda[,ind]
      }else{
        lambda0 <- lambda[,1]
      }

      if(flagsave){
        filename <- paste0(file_beta,fileID[ind],".bin")
      }else{
        filename <- NULL
      }

      #dyn.load("c_lasso.so")
      res <- .Call('R_updatebeta', Sigma, Gamma[,ind],
                  lambda0, alpha, tol, maxiter, nsup.max,
                  scale, sdx, filename,
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
      out = lapply(X=seq(q), FUN=compApply)
    }else{
      out = parallel::mclapply(X=seq(q), FUN=compApply, mc.cores=mc.cores)
    }
    if(verbose & q>1L) {
      close(pb); unlink(con)
    }

    # Checkpoint
    if(any(seq(q) != unlist(lapply(out,function(x) x$ind)) )){
      stop("Some sub-processes failed. Something went wrong during the analysis.")
    }

    out <- list(p=p, q=q, nlambda=nlambda,
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
