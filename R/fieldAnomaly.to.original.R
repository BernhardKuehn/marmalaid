# =========================================================== #
# a function that back-calculates the original field from
# field anomalies (produced with the marmalaid package)
# author: Bernhard Kuehn
# last change: 2023-08-28
# =========================================================== #

#' Add back field anomaly to a spatio-temporal field
#'
#' @description The counterpart of the \code{fieldAnomaly.raster(...)} function to add back the anomalies and get back the original field.
#' @param rst1 A raster either a \code{RasterBrick} from the raster package or \code{SpatRaster} from terra with each layer corresponding to a time step.
#' @param time A vector of the class "Date", "POSIXct", or "POSIXlt" corresponding to each layer in the raster.
#' @param level A character string. The temporal resolution over which to calculate the anomalies. Either \code{julian} / \code{day} for daily or \code{month} for monthly anomalies.
#'              Will consider "Day of the year" from as.POSIXlt(x)$yday for \code{julian} and  unique format(x,"%m%d") strings for \code{day}.
#' @return A raster either of form \code{RasterBrick} or \code{SpatRaster} containing the original field.
#' @references Taylor, M. (2017). sinkr: Collection of functions with emphasis in multivariate data analysis.
#'             R package version 0.6. \href{https://github.com/marchtaylor/sinkr}{https://github.com/marchtaylor/sinkr}
#' @details The field anomalies are calculated by subtracting the daily/monthly means from each spatial location of a field.
#'          Therefore in order to get back to the original field the mean on each location is added back to the anomalies.
#' @examples
#' library(raster)
#' data("sst.ahoi")
#'
#' # calculate anomalies
#' time.ahoi = as.Date(sub("X","",names(sst.ahoi)),format = "%Y.%m.%d")
#' sst.anom = fieldAnomaly.raster(sst.ahoi,
#'                               time = time.ahoi,
#'                               level = "month")
#' # back calculate the original field
#' original.field = fieldAnomaly.to.original(sst.anom,time = time.ahoi,level = "month")
#'
#' # check if they are the same
#' all(cellStats(original.field,mean) -cellStats(sst.ahoi,mean)== 0)
#' # correct

#' @export
fieldAnomaly.to.original = function(rst1,time,level = "month"){

  # rst1 needs to have the dim(lat,long,times)
  if((class(rst1) != "RasterBrick"| class(rst1) != "SpatRaster") & length(dim(rst1)) != 3) {
    stop("Input data 'rst1' is not of class 'Rasterbrick' or 'SpatRaster'!")
  }

  # get the means
  means = attributes(rst1)$grid.means

  # back-calculate
  if(class(rst1) == "RasterBrick"){
    mat.anom = t(raster::as.matrix(rst1)) # rows = temporal dimension # cols = spatial dimension
  } else if(class(rst1) == "SpatRaster"){
    mat.anom = t(terra::as.matrix(rst1))
  }
  if(length(time) != nrow(mat.anom)) {
    stop(paste("Supplied time vector needs to be the same length as 3rd dim. in the raster rst1!",
               "\n 3rd dim raster:",nrow(mat),"\n length time:", length(time)))
  }

  # additional check if the spatial dimensions of the means and the anomaly match
  if(dim(means)[1] != dim(mat.anom)[2]){
    stop("Spatial dimension of the supplied raster and the attached grid-means does not match!")
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
  mat <- mat.anom * NaN
  for (j in seq(levs_lookup)) {
    mat[levs_lookup[[j]], ] <- t(t(as.matrix(mat.anom[levs_lookup[[j]],
                                                      ])) + means[,j])
  }
  if(class(rst1) == "RasterBrick"){
    original = raster::rasterFromXYZ(data.frame(raster::xyFromCell(rst1,cell = 1:raster::ncell(rst1)),t(mat)))
    # crs
    raster::crs(original) = raster::crs(rst1)
    } else if(class(rst1) == "SpatRaster"){
    original = terra::rast(data.frame(terra::xyFromCell(rst1,cell = 1:raster::ncell(rst1)),t(mat)),
                               type = "xyz",crs = terra::crs(rst1))
    }
  return(original)
}
