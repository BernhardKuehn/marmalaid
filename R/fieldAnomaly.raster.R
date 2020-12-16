# =========================================================== #
# wrapper for the field anomaly function in the sinkr package
# works with raster (brick) data sets
# author: Bernhard Kuehn
# changed on 2019.09.19
# =========================================================== #

#' Field anomaly for a spatio-temporal field
#'
#' @description Wrapper function for \code{fieldAnomaly(...)} from the package sinkr to calculate the monthly field anomaly for a raster brick.
#' @param brick A raster \code{brick} with each layer corresponding to a time step.
#' @param time A vector of the class "Date", "POSIXct", or "POSIXlt" corresponding to each layer in the raster \code{brick}.
#' @param level = "month" The temporal resolution over which to calculate the anomalies. See \code{fieldAnomaly} for details.
#' @return A raster \code{brick} containing the field anomalies.
#' @references Taylor, M. (2017). sinkr: Collection of functions with emphasis in multivariate data analysis.
#'             R package version 0.6. \href{https://github.com/marchtaylor/sinkr}{https://github.com/marchtaylor/sinkr}
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
fieldAnomaly.raster = function(brick,time,level = "month"){

  # brick needs to have the dim(lat,long,times)
  if(class(brick) != "RasterBrick" & length(dim(brick)) != 3) {
    stop("Input data is not of class 'RasterBrick'!")
  }

  # calculate anomalies
  mat= t(raster::as.matrix(brick)) # rows = temporal dimension # cols = spatial dimension
  if(length(time) != nrow(mat)) {
    stop(paste("Supplied time vector needs to be the same length as 3rd dim. in the raster brick!",
               "\n 3rd dim raster:",nrow(mat),"\n length time:", length(time)))
  }

  #based on function from sinkr-pkg
  if (level == "month") {
      levs <- as.POSIXlt(time)$mon
      ulevs <- unique(levs)
      levs_lookup <- lapply(ulevs, function(x, y) which(y == x), levs)
      names(levs_lookup) <- ulevs
    }
    if (level == "julian") {
      levs <- as.POSIXlt(time)$yday
      ulevs <- unique(levs)
      levs_lookup <- lapply(ulevs, function(x, y) which(y == x), levs)
      names(levs_lookup) <- ulevs
    }
    if (level == "day") {
      levs <- format(time, format = "%m%d")
      ulevs <- unique(levs)
      levs_lookup <- lapply(ulevs, function(x, y) which(y == x), levs)
      names(levs_lookup) <- ulevs
    }
  mat.anom <- mat * NaN
    for (j in seq(levs_lookup)) {
      mat.anom[levs_lookup[[j]], ] <- t(t(as.matrix(mat[levs_lookup[[j]],
                                                    ])) - apply(as.matrix(mat[levs_lookup[[j]], ]), 2,
                                                                mean, na.rm = TRUE))
    }
  monthly.anom = raster::rasterFromXYZ(data.frame(raster::xyFromCell(brick,cell = 1:raster::ncell(brick)),t(mat.anom)))

  # crs
  raster::crs(monthly.anom) = raster::crs(brick)
  return(monthly.anom)
}
