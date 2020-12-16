# ================================================================ #
# spatial.PCA()
# ---------------------------------------------------------------- #
# Function to perform spatial PCA on a raster dataset (S-mode PCA)
# plot relevant PCs and write results to matrix and raster format
# ---------------------------------------------------------------- #
# author: Bernhard Kuehn
# last change: 2018-01-30
# ================================================================ #

#' EOF analysis (S-mode PCA) for a spatio-temporal field
#'
#' @description Wrapper function for \code{prcomp(...)} on a spatio-temporal field organised as raster
#' @param x Either a raster \code{brick} for scalar valued spatio-temporal fields (e.g. SST) with each layer corresponding to a time step or a two-element \code{list} with each element being a \code{raster} for vector valued fields (e.g. u & v of Currents).
#'  If layer names contain a reference to a date it will be used for writing it out.
#' @param center A logical value indicating whether the variables should be centered around the mean.
#' @param scale. A logical value indicating whether to scale each gridcell by its variance, if \code{FALSE} the PCA is calculated on the covariance matrix, otherwise (\code{TRUE}) on the correlation matrix.
#' @param spatial.extent Region to which the spatial field should be cropped, either \code{NULL} for no cropping, an object of class \code{extent} or the character string \code{"North Sea"} to cut it to the extent(c(xmin = -7.01,xmax = 12.9,ymin = 48.9,ymax = 61)) corresponding to the North Sea region.
#' @param plot A logical value if the spatial EOF fields should be returned in a plot after calculation.
#' @param var.threshold A value, indicating up to which percentage of explained variance the spatial fields should be returned (default = 0.95 (aka 95% variance explained)).
#' @return A \code{list} with the following elements:
#' \itemize{
#'   \item raster - A \code{raster} for scalar valued input fields or a \code{list} of \code{rasters} for vector valued fields. Contains the spatial-EOF fields.
#'   \item sp.PCA - A list returned by \link[stats]{prcomp} containing the output of the PCA with
#'    \cr \code{x} - containing the PC timeseries and
#'    \cr \code{rotation} - the eigenvectors (aka EOF-weights per grid-point).
#'   \item time - the time associated to the PCs
#'   }
#' @examples
#' # load AHOI SST data
#' library(raster)
#' data(sst.ahoi)
#' # perform EOF-analysis
#' SST.EOF = spatial.PCA(sst.ahoi,center = TRUE,scale. = FALSE,spatial.extent = NULL,plot = FALSE)
#' # plot first EOF with PC timeseries
#' layout(matrix(rep(c(1,2),each = 2),nrow = 2,ncol = 2,byrow = TRUE),heights = c(0.65,0.35))
#' par(mar = c(4,4,3,3))
#' image(SST.EOF$raster[[1]],las = 1,xlab = "lon",ylab = "lat")
#' maps::map(add = TRUE,fill = TRUE,col = "gray80");box()
#' title("EOF1",adj = 0)
#' plot(SST.EOF$time,SST.EOF$sp.PCA$x[,1],type = "l",las = 1,xlab = "time",ylab = "")
#' title("PC1",adj = 0)

#' @export
spatial.PCA = function(x,center = TRUE , scale. = FALSE, spatial.extent = "North Sea",plot = FALSE,var.threshold = 0.95)  {

  # ----------------------------------------------------------------------------------- #
  # 1. step: prepare raster data as input for PCA
  #   - crop to spatial extent
  #   - convert raster to format: matrix(time,space) and save position of raster cells
  #   - x = raster brick with time dimension as additional layers
  # ----------------------------------------------------------------------------------- #

  # # x can be a raster layer or a list of raster-layers if two variables are calculated simultaniously (e.g. u and v)
  if(class(x)=="list") {
    if(all(lapply(x,class)%in%c("raster","RasterBrick","RasterStack"))) {

      # define spatial extent or write codeword
      if(is.character(spatial.extent) == TRUE) {
        if(spatial.extent == "North Sea") {
          # define spatial extent
          spatial.extent = raster::extent(c(xmin = -7.01,
                                            xmax = 12.9,
                                            ymin = 48.9,
                                            ymax = 61))
        } else {
          stop("Only 'North Sea' as string allowed! Otherwise use NULL for no cropping or obj. of class(extent)")
        }
      } else if(is.null(spatial.extent)){
        spatial.extent = raster::extent(x[[1]])
      }
      # crop to spatial extent
      x.crop = vector("list",length(x))
      for(i in 1:length(x)) {
        x.crop[[i]] = raster::crop(x[[i]],y = spatial.extent)
      }

      # combine files to one matrix
      #x.add = do.call(addLayer,x)
      # get geo.position
      geo.pos = vector("list",length(x))
      for(i in 1:length(geo.pos)) {
        geo.pos[[i]] = raster::xyFromCell(x.crop[[i]],1:raster::ncell(x.crop[[i]]))
      }
      if(do.call(all.equal.numeric,geo.pos)) {
        geo.pos = geo.pos[[1]]
      } else {
        stop("Rasters in list have not the same spatial resolution!")
      }

      # convert to matrix
      print("1. convert raster to matrix...")
      tmp.mat = vector("list",length(x))
      for(i in 1:length(x)) {
        tmp.mat[[i]] = t(raster::as.matrix(x.crop[[i]]))
      }

      mat = do.call(cbind,tmp.mat); rm(tmp.mat)
      # bring files together in correct form

    }else{
      stop("Not all elements of list(x) are raster-files!")
    }
  } else {
     # define spatial extent or write codeword
  if(is.character(spatial.extent) == TRUE) {
    if(spatial.extent == "North Sea") {
      # define spatial extent
    spatial.extent = raster::extent(c(xmin = -7.01,
                                      xmax = 12.9,
                                      ymin = 48.9,
                                      ymax = 61))
    } else {
      stop("Only 'North Sea' as string allowed! Otherwise use NULL for no cropping or obj. of class(extent)")
    }
  } else if(is.null(spatial.extent)){
    spatial.extent = raster::extent(x)
  }
  # crop to spatial extent
  x.crop = raster::crop(x,y = spatial.extent)

  # convert to matrix
  print("1. convert raster to matrix...")
  mat = t(raster::as.matrix(x.crop))

  # save geo.pos of raster grid cells
  geo.pos = raster::xyFromCell(x.crop,1:raster::ncell(x.crop))
  }

  # check for NA column
  NA.cols = which(apply(mat, 2, function(x) all(is.na(x)))) # | var(x)== 0 maybe also if var == 0, but not applicable for currents

  # remove NA data, but only if there are NA.cols, otherwise not
  if(length(NA.cols >0)){
  mat = mat[,-NA.cols]
  geo.pos = geo.pos[-NA.cols,]
  }
  # -------------- #
  # 2. perform PCA
  # -------------- #
  print("2. perform PCA...")
  sp.PCA = stats::prcomp(mat, center = center, scale. = scale.)

  # explained variance
  expl.var <- sp.PCA$sdev^2 / sum(sp.PCA$sdev^2) # explained variance
  cum.expl.var <- cumsum(expl.var) # cumulative explained variance

  # get all PCs with var > var.threshold
  indx.PCs = which(cum.expl.var>=var.threshold)[1]

  # ------------------------------ #
  # 3. convert back to raster format
  # ------------------------------ #
  print("3. convert back to raster-format...")
  # create long table format with lat, long and loadings of first axes that explain 95% of the variance
  sp.PCA.df = data.frame(geo.pos,sp.PCA$rotation[,1:indx.PCs])
  # rasterize
  if(class(x)=="list") {
    split.df = split(sp.PCA.df,f = rep(1:length(x),each = nrow(geo.pos)))
    sp.PCA.raster = vector("list",length = length(x))
    for(i in 1:length(x)) {
      sp.PCA.raster[[i]] = raster::rasterFromXYZ(split.df[[i]])
      }
    } else{
      sp.PCA.raster = raster::rasterFromXYZ(sp.PCA.df)
  }

  # ------------ #
  # 4. plot data
  # ------------ #

  if(plot == TRUE) {
    opar = graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(opar))
    print("4. plot PCA...")
    # devide plot window in optimal way
    plot.x.y = indx.PCs/(15:1)[which(((indx.PCs/(15:1))%%1 == 0))[1:2]]
    plot.x.y[is.na(plot.x.y)] <- 1
    graphics::par(mfrow = c(plot.x.y[1],plot.x.y[2]),mar = c(5,4,3,3))
    for(i in 1:indx.PCs) {
       raster::image(sp.PCA.raster[[i]], col = grDevices::colorRampPalette(RColorBrewer::brewer.pal(10,"RdYlBu"))(100),las = 1)
       sp::plot(raster::rasterToContour(x = sp.PCA.raster[[i]]), add=TRUE, col="gray30")
       maps::map(add = TRUE,fill = TRUE,col = "gray90")
       graphics::box()
       graphics::title(paste0("PC",i), sub = paste("% Var expl:",100 *summary(sp.PCA)$importance["Proportion of Variance",i]))
    }
  }

  # ---- #
  # time
  # ---- #
  if(class(x) == "list"){
    x = x[[1]]
  }
  year = as.numeric(substr(x@data@names,start = 2,stop = 5))
  month = as.numeric(substr(x@data@names,start = 7,stop = 8))
  day = as.numeric(substr(x@data@names,start = 10,stop = 11))

  if(all(is.na(day)) & all(is.na(month))) {
     time = year
  } else if(all(is.na(day)) ){
    time = as.Date(paste0(year,"-",month,"-",15),format = "%Y-%m-%d")
    } else {
    time = as.Date(paste0(year,"-",month,"-",day),format = "%Y-%m-%d")
  }

  # -------- #
  # read out
  # -------- #

    sp.PCA.all = list(raster = sp.PCA.raster,sp.PCA = sp.PCA,time = time)

    return(sp.PCA.all)

}
