# =========================================================== #
# wrapper for the field anomaly function in the sinkr package
# works with raster (brick) data sets
# changes: - incorporate terra-functionalities
#          - incorporate way to remove the mean based on some additional anomaly raster
# author: Bernhard Kuehn
# first version: 2019.09.19
# last change: 2023.08.18
# =========================================================== #

#' Field anomaly for a spatio-temporal field
#'
#' @description Wrapper function for \code{fieldAnomaly(...)} from the package sinkr to calculate the daily/monthly field anomaly for a raster brick.
#' @param rst1 A raster either a \code{RasterBrick} from the raster package or \code{SpatRaster} from terra with each layer corresponding to a time step.
#' @param rst2.anomalies An optional second raster, where anomalies already has been calculated on.
#'                       Used for correcting the first raster with the same mean anomalies e.g. as a historical reference period to base the correction on.
#' @param time A vector of the class "Date", "POSIXct", or "POSIXlt" corresponding to each layer in the raster.
#' @param level A character string. The temporal resolution over which to calculate the anomalies. Either \code{julian} / \code{day} for daily or \code{month} for monthly anomalies.
#'              Will consider "Day of the year" from as.POSIXlt(x)$yday for \code{julian} and  unique format(x,"%m%d") strings for \code{day}.
#' @return A raster \code{brick} containing the field anomalies.
#' @references Taylor, M. (2017). sinkr: Collection of functions with emphasis in multivariate data analysis.
#'             R package version 0.6. \href{https://github.com/marchtaylor/sinkr}{https://github.com/marchtaylor/sinkr}
#' @details The field anomalies are calculated by subtracting the daily/monthly means from each spatial location of a field.
#' @examples
#' # load AHOI SST data
#' library(raster)
#' data(sst.ahoi)
#' # calculate anomalies
#' sst.anom = fieldAnomaly.raster(sst.ahoi,
#'                                time = as.Date(sub("X","",names(sst.ahoi)),format = "%Y.%m.%d"),
#'                                level = "month")
#' # plot for comparison
#' par(mfrow = c(1,2))
#' plot(sst.anom[[1]],las = 1)
#' title("SST-anomaly",adj = 0)
#' plot(sst.ahoi[[1]],las = 1)
#' title("SST-raw",adj = 0)


#' @export
fieldAnomaly.raster = function(rst1,rst2.anomalies = NULL,time,level = "month"){

  # rst1 needs to have the dim(lat,long,times)
  if((class(rst1) != "RasterBrick"| class(rst1) != "SpatRaster") & length(dim(rst1)) != 3) {
    stop("Input data 'rst1' is not of class 'Rasterbrick' or 'SpatRaster'!")
  }
  if(!is.null(rst2.anomalies)){
    if((class(rst2.anomalies) != "RasterBrick"| class(rst2.anomalies) != "SpatRaster") & length(dim(rst2.anomalies)) != 3) {
      stop("Input data 'rst2' is not of class 'Rasterbrick' or 'SpatRaster'!")
    }
  }


  # calculate anomalies
  if(class(rst1) == "RasterBrick"){
    mat= t(raster::as.matrix(rst1)) # rows = temporal dimension # cols = spatial dimension
  } else if(class(rst1) == "SpatRaster"){
    mat= t(terra::as.matrix(rst1))
  }
  if(length(time) != nrow(mat)) {
    stop(paste("Supplied time vector needs to be the same length as 3rd dim. in the raster rst1!",
               "\n 3rd dim raster:",nrow(mat),"\n length time:", length(time)))
  }

  #based on function from sinkr-pkg
  if (level == "month") {
    levs <- as.POSIXlt(time)$mon
    ulevs <- sort(unique(levs))
    levs_lookup <- lapply(ulevs, function(x, y) which(y == x), levs)
    names(levs_lookup) <- ulevs
  }
  if (level == "julian") {
    levs <- as.POSIXlt(time)$yday
    ulevs <- sort(unique(levs))
    levs_lookup <- lapply(ulevs, function(x, y) which(y == x), levs)
    names(levs_lookup) <- ulevs
  }
  if (level == "day") {
    levs <- format(time, format = "%m%d")
    ulevs <- sort(unique(levs))
    levs_lookup <- lapply(ulevs, function(x, y) which(y == x), levs)
    names(levs_lookup) <- ulevs
  }
  mat.anom <- mat * NaN
  if(!is.null(rst2.anomalies)){
    means = attributes(rst2.anomalies)$grid.means
  } else{
    means = matrix(NA,nrow = ncol(mat),ncol = length(levs_lookup))
  }
  for (j in seq(levs_lookup)) {
    if(is.null(rst2.anomalies)){
      means[,j] = apply(as.matrix(mat[levs_lookup[[j]], ]),
                        2, mean, na.rm = TRUE)
    }
    mat.anom[levs_lookup[[j]], ] <- t(t(as.matrix(mat[levs_lookup[[j]],
                                                      ])) - means[,j])
  }
  if(class(rst1) == "RasterBrick"){
    monthly.anom = raster::rasterFromXYZ(data.frame(raster::xyFromCell(rst1,cell = 1:raster::ncell(rst1)),t(mat.anom)))
    # crs
    raster::crs(monthly.anom) = raster::crs(rst1)
    attributes(monthly.anom)$grid.means = means
  } else if(class(rst1) == "SpatRaster"){
    monthly.anom = terra::rast(data.frame(terra::xyFromCell(rst1,cell = 1:raster::ncell(rst1)),t(mat.anom)),
                               type = "xyz",crs = terra::crs(rst1))
    attributes(monthly.anom)$grid.means = means
  }

  return(monthly.anom)
}
