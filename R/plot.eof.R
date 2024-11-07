#' Function to plot the EOF pattern (spatial/temporal)
#'
#' This functions is a rather easy plotting function to see the spatial & temporal
#' representation of the first n calculated EOF modes.
#' Typically a maximum of 15 modes should not be exceeded.
#'
#' @param eof.out An output object from the \code{spatial.PCA} function
#' @param n The number of corresponding EOF-pattern/PC-time series to show.

#' @export
plot.eof = function(eof.out,n){
  opar = graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(opar))
  # get PC time series
  PCs_temporal = eof.out$sp.PCA$x
  # get spatial pattern
  EOFs_spatial = eof.out$raster
  if(n>15){
    warning("'n' should be chosen to be below/equal to 15 for a reasonable display!")
  }

  if(dim(EOFs_spatial)[3] < n){
    warning("'n' is chosen smaller than the number of available EOF fields.\n
             Please make sure to increase the 'var.theshold' in 'spatial.PCA' to
             allow displaying all 'n' EOF-patterns.")
    warning(paste0("Setting 'n' to: ",dim(EOFs_spatial)[3],"!"))
    n = dim(EOFs_spatial)[3]
  }
  graphics::par(mfrow = c(1,1),mar = c(2,4,1,3))
  # create plot layout
  tmp = matrix(1:(2*n),nrow = 2)
  if(n > 3){
    frst.row = tmp[,1:ceiling(n/2)]
    scnd.row = tmp[,(1+ceiling(n/2)):n]
    if(ncol(frst.row) > ncol(scnd.row)){
      scnd.row = cbind(scnd.row,c(max(scnd.row)+1,max(scnd.row)+2))
    }
    mat = rbind(frst.row,scnd.row)
  } else if(n >= 9){
    indx.split = split(1:n,f = rep(1:3,each = ceiling(n/3))[1:n])
    frst.row = tmp[,indx.split[[1]]]
    scnd.row = tmp[,indx.split[[2]]]
    thrd.row = tmp[,indx.split[[3]]]
    if(ncol(scnd.row) > ncol(thrd.row)){
      thrd.row = cbind(thrd.row,c(max(thrd.row)+1,max(thrd.row)+2))
    }
    mat = rbind(frst.row,scnd.row,thrd.row)
  } else {
    mat = tmp
  }
  graphics::layout(mat,
         height = rep(c(1.5,1),nrow(mat)/2))

  for (i in 1:n) {
    terra::image(EOFs_spatial[[i]],
                 col = grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(10,
                                                                                "RdYlBu")))(100),
                 las = 1,xlab = "",ylab = "lat")
    terra::contour(x = EOFs_spatial[[i]], add = TRUE,
                   col = "gray30",labcex = 0.8,vfont = c("sans serif","bold"))
    maps::map(add = TRUE, fill = TRUE, col = "gray90")
    graphics::box()
    graphics::title(paste0("EOF", i),adj = 0)
    graphics::mtext("lon",side = 1,line = 1.5,cex = 0.8)
    # add time series
    plot(PCs_temporal[,i],type = "l",xlab = "year",
         ylab = "",las = 1)
    graphics::title(paste0("PC",i),adj = 0,cex = 0.9)
  }
}
