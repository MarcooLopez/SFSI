
#====================================================================
# Obtain a list with common trn elements for all columns in 'y'
#====================================================================
get_common_trn <- function(y){
  trn_list0 <- lapply(1:ncol(y), function(j)which(!is.na(y[,j])))
  ID <- 1:ncol(y)
  trn_list <- list()
  cont <- 0
  while(length(ID)>0){
    trn0 <- trn_list0[[ID[1]]]
    index <- c()
    if(length(ID)>1){
      for(k in 2:length(ID)){
        tmp <- trn_list0[[ID[k]]]
        if(length(trn0)==length(tmp)){
          if(all(trn0==tmp)){
            index <- c(index,k)
          }
        }
      }
    }
    if(length(index)>0){
      common <- ID[c(1,index)]
      ID <- ID[-c(1,index)]
    }else{
      common <- ID[1]
      ID <- ID[-1]
    }
    cont <- cont + 1
    trn_list[[cont]] <- list(iy=common, trn=trn0)
  }

  tmp <- unlist(lapply(trn_list,function(x)x$iy))
  if(any(seq(ncol(y)) != tmp[order(tmp)])){
    stop("There was an error in obtaining common training sets")
  }
  trn_list
}

#====================================================================
# Labels and breaks for the nsup axis
#    x = -log(lambda)
#    y = nsup
#====================================================================
get_breaks <- function(x, y, nbreaks.x=6, ymin=1)
{
  fm <- stats::smooth.spline(x,y)

  tmp <- cbind(x,stats::fitted(fm))
  tmp <- tmp[tmp[,2]>=ymin,]
  tmp <- tmp[abs(tmp[,2]-min(tmp[,2]))<1E-8,,drop=F]
  xmin <- mean(tmp[,1], na.rm=T)

  tmp <- seq(min(x), ifelse(xmin<=0,0,2*xmin), length=1000)
  tt <- stats::predict(fm,tmp)
  xmin <- tt$x[min(which(tt$y >= ymin))]

  breaks.x <- seq(xmin, max(x), length=nbreaks.x)
  breaks.y <- stats::predict(fm, breaks.x)$y

  return(list(breaks.x=breaks.x,breaks.y=breaks.y))
}

#====================================================================
#====================================================================
map_set <- function(n, q, x=NULL, y=NULL, labels=NULL,
                    xlab="x", ylab="y")
{

  MAP <- data.frame(index=seq(n*q), expand.grid(i=seq(n), j=seq(q)))
  MAP$i_j <- paste0(MAP$i,"_",MAP$j)
  MAP$set <- NA
  MAP$index_set <- NA

  if(!is.null(y)){
    y <- as.integer(y)
    stopifnot(all(y %in% MAP$index))
    MAP[y,"set"] <- ylab
    MAP[y,"index_set"] <- seq_along(y)
  }

  if(!is.null(x)){
    x <- as.integer(x)
    stopifnot(all(x %in% MAP$index))
    MAP[x,"set"] <- xlab
    MAP[x,"index_set"] <- seq_along(x)
  }

  if(is.null(labels)){
    MAP$label <- MAP$i
  }else{
    stopifnot(length(labels) == n)
    MAP$label <- labels[MAP$i]
  }
  MAP$label_j <- paste0(MAP$label,"_",MAP$j)
  MAP
}

#====================================================================
#====================================================================
get_summary_nsup <- function(object, map=NULL, eps = .Machine$double.eps){

  ntraits <- object$ntraits
  nlambda <- object$nlambda
  names_nsup <- paste0("nsup_",1:ntraits)
  if(is.null(map)){
    map <- map_set(object$n, ntraits,
                   x=object$trn, y=object$tst,
                   xlab="trn", ylab="tst")
  }
  map_trn <- map[object$trn,]
  map_tst <- map[object$tst,]

  out <- expand.grid(SSI=paste0("SSI.",seq(nlambda)), trait=seq(ntraits))
  out <- data.frame(out, matrix(NA, nrow=nrow(out), ncol=ntraits))
  colnames(out) <- c("SSI","trait",names_nsup)
  ID <- paste(out$SSI,out$trait)

  nsup0 <- matrix(NA, nrow=ntraits, ncol=ntraits)

  for(ilambda in 1:nlambda){
    SSIname <- paste0("SSI.",ilambda)
    b <- coef.SSI(object, ilambda=ilambda)
    for(i in 1:ntraits){
      b0 <- b[map_tst$j==i,,drop=F]
      for(j in 1:ntraits){
        nsup0[i,j] <- mean(apply(b0[,map_trn$j==j,drop=F],1,function(x)sum(abs(x)>eps)))
      }
    }
    out[match(paste(SSIname,seq(ntraits)), ID),names_nsup] <- nsup0
  }
  return(out)
}

#====================================================================
#====================================================================

has_names <- function(A){
  if(length(dim(A)) == 2L){
    out <- length(unlist(dimnames(A))) == sum(dim(A))
  }else{
     out <- FALSE
  }

  out
}

#====================================================================
#====================================================================

capitalize <- function(string){
  substr(string,1,1) <- toupper(substr(string,1,1))
  string
}

#====================================================================
#====================================================================

circleFun <- function(center=c(0,0), radius=1, n=200){
    tt <- seq(0, 2*pi, length.out = n)
    xx <- center[1] + radius * cos(tt)
    yy <- center[2] + radius * sin(tt)
    return(data.frame(x = xx, y = yy))
}

#====================================================================
# Obtain layout to for the net.plot function
#====================================================================
get_net <- function(X, MAP, symmetric,
                    K = NULL, p.radius = 1.75,
                    delta = .Machine$double.eps)
{

  axis_labels <- namesK <- NULL

  n <- length(unique(MAP$i))
  q <- length(unique(MAP$j))

  if(!is.null(K)){
    if(has_names(K)){
      stopifnot(all(rownames(K)==colnames(K)))
      namesK <- rownames(K)
    }
  }

  xxx <- which(MAP$set=="x")
  yyy <- which(MAP$set=="y")
  if(symmetric | length(yyy)==0){
    yyy <- xxx[]
  }

  # List of modules: Group of nodes (primary nodes) connected to other nodes (secondary nodes)
  modules <- vector("list",length(yyy))
  for(k in 1:nrow(X)){
    modules[[k]] <- NA
    if(symmetric){
      modules[[k]] <- xxx[which(abs(X[k,k:ncol(X)]) > delta) + k -1]
    }else{
      modules[[k]]  <- xxx[which(abs(X[k, ]) > delta)]
    }
  }

  if(is.null(K)){
    if(symmetric){
      eee <- c()
      gr <- igraph::make_empty_graph(n = n)
      for(j in 1:nrow(X)){
        index <- modules[[j]][-1]
        if(length(index[!is.na(index)]) > 0) eee <- c(eee, as.vector(rbind(j,index)))
      }
      gr <- igraph::add_edges(gr, eee)
      xy <- list(igraph::layout_with_fr(gr, dim=2, niter=2000))
    }else{
      gr <- igraph::make_empty_graph(n = n)
      eee <- c()
      MAP0 <- MAP[yyy,]
      uniquei <- unique(MAP0$i)
      for(j in 1:length(uniquei)){
        i0 <- uniquei[j]
        map0 <- MAP0[MAP0$i==i0, ]
        con0 <- c()
        for(k in 1:nrow(map0)){
          yi <- which(yyy == map0$index[k])
          tmp <- unique(MAP[modules[[yi]],'i'])
          con0 <- c(con0, tmp)
        }
        con0 <- con0[!is.na(con0)]
        if(length(con0) > 0){
          eee <- c(eee, as.vector(rbind(i0,con0)))
        }
      }
      gr <- igraph::add_edges(gr, eee)
      tmp <- igraph::layout_with_fr(gr, dim=2, niter=2000)
      xy <- vector('list',q)
      for(j in 1:q){
        xy[[j]] <- tmp
      }
    }

  }else{
    tmp <- svd(K, nu=2, nv=0)
    d <- tmp$d
    xy <- tmp$u[,1:2]
    xy <- lapply(1:q,function(k) xy)

    expvarPC <- paste0(" (",sprintf('%.1f',100*d/sum(d)),"%)")
    if(length(expvarPC) < 2) expvarPC <- NULL
    axis_labels <- paste0("PC ",1:2,expvarPC[1:2])
  }

  # Center and  scale to -1, 1
  xy <- lapply(xy,function(tt){
    mm0 <- apply(tt,2, function(x) (max(x)+min(x))/2)
    tt <- sweep(tt, 2L, mm0, FUN="-")
    apply(tt,2,function(x)x/max(abs(x)))
  })

  radius <- unlist(lapply(xy,function(tt){
    max(apply(tt,2,function(x)diff(range(x))/2)) # Radius
  }))

  if(q == 1L){
    xy <- do.call(rbind,xy)
    mid_point <- matrix(c(x=sum(range(xy[,1]))/2, y=sum(range(xy[,1]))/2),nrow=1)

  }else{
    radius0 <- p.radius*max(radius)
    deg <- seq(0,q-1)*(360/q)
    rad <- deg*pi/180   #  Convert to radians
    mid_point <- cbind(x=radius0*cos(rad), y=radius0*sin(rad))
    xy <- do.call(rbind,lapply(1:q,function(k) sweep(xy[[k]],2L,mid_point[k,],FUN="+")))
  }

  colnames(xy) <- c("x","y")

  return(list(xy=xy, axis_labels=axis_labels,
              mid_point=mid_point, radius=radius,
              modules=modules, index_module=yyy))
}

#====================================================================
#====================================================================
.onAttach <- function(libname, pkgname) {
  packageStartupMessage("
  |=======================================================================|
  |    ._______. ._______. ._______. ._______.                            |
  |    | ._____| | ._____| | ._____| |__. .__|                            |
  |    | |_____. | |___.   | |_____.    | |                               |
  |    |_____. | | .___|   |_____. |    | |      Authors:                 |
  |    ._____| | | |       ._____| | .__| |__.     Marco Lopez-Cruz       |
  |    |_______| |_|       |_______| |_______|     Gustavo de los Campos  |
  |                                                                       |
  |    Sparse Family and Selection Index. Version 1.3.0 (Jul 25, 2023)    |
  |    Type 'citation('SFSI')' to know how to cite SFSI                   |
  |    Type 'help(package='SFSI',help_type='html')' to see help           |
  |    Type 'browseVignettes('SFSI')' to see documentation                |
  |    Type 'demo(package='SFSI')' to see demos                           |
  |                                                                       |
  |=======================================================================|
  ")
}
