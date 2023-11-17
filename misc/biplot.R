


# scores.group <- loadings.group <- NULL;
# scores.color <- "gray70"; loadings.color <- "orange3";
# method="princomp"; title <- "Biplot"; linewidth = 0.5
# text.size=2.5; point.size=1.3; x.disp = 0.05; y.disp = 0
# scores.label <- scores.point <- loadings.arrow <- TRUE; text.position="top"

# =================================================================
# Make a biplot from a SVD showing scores and loadings
# =================================================================
biplot2 <- function(X, center = TRUE, scale = FALSE,
                    title = "PCA - biplot", linewidth = 0.5,
                    scores.group = NULL, loadings.group = NULL,
                    scores.color = "darkorange4", loadings.color = "royalblue4",
                    text.size = 2.5, point.size = 1.5, arrow.text.size = 3.2,
                    text.position = c("top","right","bottom","left","center","topleft","topright","bottomleft","bottomright"),
                    x.disp = 0.6, y.disp = 0.6, arrow.size = 3.1,
                    scores.point = TRUE, scores.label = TRUE,
                    show.origin = TRUE, origin.color = "gray65",
                    expand.x = 0.06, expand.y = 0.06,
                    loadings.arrow = TRUE)
{
	text.position <- match.arg(text.position)

	tmp <- svd(scale(X, center=center, scale=scale))
	U <- tmp$u[,1:2]
	V <- tmp$v[,1:2]
	SD <- tmp$d
	datU <- data.frame(names=rownames(X),PC1=U[,1],PC2=U[,2])
	datV <- data.frame(names=colnames(X),PC1=V[,1],PC2=V[,2])

	VAR <- SD^2
	tmp <- (100*VAR/sum(VAR))
	xlab0 <- paste0("PC1 (",sprintf('%.1f',tmp[1]),"%)")
	ylab0 <- paste0("PC2 (",sprintf('%.1f',tmp[2]),"%)")

  	# Center and  scale to -1, 1
  	scalePC <- function(x){
    	mm <-  (max(x)+min(x))/2
    	x2 <- x - mm
    	x2/max(abs(x2))
  	}
  	datU$PC1 <- scalePC(datU$PC1)
  	datU$PC2 <- scalePC(datU$PC2)

  	datV$PC1 <- scalePC(datV$PC1)
  	datV$PC2 <- scalePC(datV$PC2)

	theme0 <- theme(plot.title = element_text(hjust=0.5),
					panel.grid.major=element_blank(),
					panel.grid.minor=element_blank())

  	# Scores setting
  	gpname1 <- ""; guide_label <- guide_point <- 'none'
  	if(is.null(scores.group)){
    	datU <- data.frame(datU, group=factor("none"))
		stopifnot(length(scores.color)==1)
		names(scores.color) <- "none"
  	}else{
		stopifnot(all(rownames(X) %in% rownames(scores.group)))
    	datU <- data.frame(datU, group=as.character(scores.group[rownames(X),1]))
		datU$group <- factor(datU$group)
		if(length(scores.color) != nlevels(datU$group)){
			scores.color <- RColorBrewer::brewer.pal(8,"Dark2")[1:nlevels(datU$group)]
			names(scores.color) <- levels(datU$group)
		}else{
			if(is.null(names(scores.color))){
				names(scores.color) <- levels(datU$group)
			}
		}
		gpname1 <- colnames(scores.group)[1]
		guide_point <- guide_legend(override.aes=aes(label="", size=3.2, alpha=1))
  	}

  	# Loadings setting
	gpname2 <- ""; guide_arrow <- 'none'
	if(is.null(loadings.group)){
		datV <- data.frame(datV, group=factor("none"))
		stopifnot(length(loadings.color)==1)
		names(loadings.color) <- "none"
	}else{
		stopifnot(all(colnames(X) %in% rownames(loadings.group)))
		datV <- data.frame(datV, group=as.character(loadings.group[colnames(X),1]))
		datV$group <- factor(datV$group)
		if(length(loadings.color) != nlevels(datV$group)){
			loadings.color <- gg_color_hue(nlevels(datV$group))
			#loadings.color <- RColorBrewer::brewer.pal(8,"Set2")[1:nlevels(datV$group)]
			names(loadings.color) <- levels(datV$group)
		}else{
			if(is.null(names(loadings.color))){
				names(loadings.color) <- levels(datV$group)
			}
		}
		gpname2 <- colnames(loadings.group)[1]
		guide_arrow <- guide_legend(override.aes=aes(label="", linewidth=1))
	}

  	# Get text position for arrows
  	datV <- data.frame(datV,hjust=NA,vjust=NA)
  	for(i in 1:nrow(datV)){
		pos <- text_arrow_pos(datV$PC1[i], datV$PC2[i], min.angle=20)
    	datV[i,c("hjust","vjust")] <- get_hvjust(pos, x.disp=0.3, y.disp=0.3)
  	}

  	# Making the plot
  	pp <- ggplot(datU, aes(PC1,PC2))
	if(show.origin){
		pp <- pp + geom_hline(yintercept=0, linetype="dashed", color=origin.color, linewidth=0.5) +
		           geom_vline(xintercept=0, linetype="dashed", color=origin.color, linewidth=0.5)
	}
  	if(scores.label){
		if(scores.point){
			text.position <- ifelse(text.position=="center","top",text.position)
		}else{
			text.position <- "center"
		}
		tmp <- get_hvjust(text.position, x.disp=x.disp, y.disp=y.disp)
    	pp <- pp + geom_text(aes(label=names, color=group), size=text.size,
												 hjust=tmp[1], vjust=tmp[2])
  	}
  	if(scores.point){
    	pp <- pp + geom_point(aes(fill=group),shape=21, size=point.size)
  	}else{
		# Plot point with alpha=0 if labels are plotted for a better legend
		if(scores.label){
			pp <- pp + geom_point(aes(fill=group), shape=21, alpha=0)
		}
	}
	pp <- pp +
	      scale_color_manual(gpname1, values=scores.color, guide='none') +
				scale_fill_manual(gpname1, values=scores.color, guide=guide_point) +
  			ggnewscale::new_scale_fill() + ggnewscale::new_scale_color()

	if(loadings.arrow){
		pp <- pp +
          geom_segment(data=datV, aes(x=0,y=0,xend=PC1,yend=PC2,color=group),
                       linewidth=linewidth,
                       arrow = arrow(length=unit(arrow.size, "pt"))) +
          geom_text(data=datV, aes(label=names, vjust=vjust, hjust=hjust, color=group),
                    size=arrow.text.size) +
					scale_color_manual(gpname2,values=loadings.color, guide=guide_arrow)
	}

  	pp <- pp +
			scale_x_continuous(expand=c(expand.x,expand.x)) +
			scale_y_continuous(expand=c(expand.x,expand.x)) +
			labs(x=xlab0, y=ylab0, title=title) +
			theme_bw() + theme0

 	pp
}

# =================================================================

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# =================================================================

get_hvjust <- function(position, x.disp=0, y.disp=0){
  switch(position,
		'top'         = c(0.5,0.0-y.disp),
		'right'       = c(0.0-x.disp,0.5),
		'bottom'      = c(0.5,1.0+y.disp),
		'left'        = c(1.0+x.disp,0.5),
		'center'      = c(0.5,0.5),
		'topright'    = c(0.0-x.disp,0.0-y.disp),
		'topleft'     = c(1.0+x.disp,0.0-y.disp),
		'bottomleft'  = c(1.0+x.disp,1.0+y.disp),
		'bottomright' = c(0.0-x.disp,1.0+y.disp)
  )
}

# =================================================================

get_quadrant <- function(x, y){
	 angle <- abs(atan(x/y)*180/pi)
	 if(x>=0 & y>=0) return(1)
	 if(x<0 & y>=0) return(2)
	 if(x<0 & y<0) return(3)
	 if(x>=0 & y<0) return(4)
}

# =================================================================

text_arrow_pos <- function(x, y, min.angle=10)
{
	#angle <- atan(y/x)*180/pi
	qq <- get_quadrant(x,y)
	angle <- atan2(y,x)*180/pi

	if(qq == 1){ # [x=1;y=0.0001]  [x=0.0001;y=1]
		pos <- ifelse(angle<min.angle,"right",ifelse(abs(angle-90)<min.angle,"top","topright"))
	}
	if(qq == 2){ # [x=-0.0001;y=1]  [x=-1;y=0.0001]
		angle <- angle - 90
		pos <- ifelse(angle<min.angle,"top",ifelse(abs(angle-90)<min.angle,"left","topleft"))
	}
	if(qq == 3){ # [x=-1;y=-0.0001]  [x=-0.0001;y=-1]
		angle <- 360+angle
		angle <- angle - 180
		pos <- ifelse(angle<min.angle,"left",ifelse(abs(angle-90)<min.angle,"bottom","bottomleft"))
	}
	if(qq == 4){ # [x=0.0001;y=-1]  [x=1;y=-0.0001]
		angle <- 360+angle
		angle <- angle - 270
		pos <- ifelse(angle<min.angle,"bottom",ifelse(abs(angle-90)<min.angle,"right","bottomright"))
	}
	return(pos)
}
