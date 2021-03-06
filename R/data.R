#' SST data of the North Sea
#'
#' @source Thuenen Institute for Sea fisheries, "AHOI : A physical-statistical model of hydrography for fishery and ecology studies"
#'  \url{https://www.thuenen.de/en/sf/projects/ahoi-a-physical-statistical-model-of-hydrography-for-fishery-and-ecology-studies/}
#' @description A raster dataset containing the averaged temperature over the first 10m in the North Sea for the years 1990 - 2014.
#' @format A raster brick:
#' \describe{
#'  \item{dimensions}{latitude, longitude, number cells, time}
#'  \item{resolution}{spatial resolution of the data set}
#'  \item{extent}{spatial extension of the data set}
#'  \item{crs}{Coordinate reference system}
#'  \item{names}{Names of each layer in the time dimension - contains the time step, here denoted with XYYYY.MM.DD}
#' }
#' @references Núñez-Riboni I, Akimova A (2015) "Monthly maps of optimally interpolated in situ hydrography in the North Sea from 1948 to 2013." J Mar Syst 151:15–34
#' @details This dataset contains the averaged temperature over the first 10m (in the following refered to as sea surface temperature aka "SST") in the North Sea between -5.1 to 10.1 longitude and 50.5 to 62.1 latitude.
#' It is based on an extraction of the AHOI dataset (v19.02), which spans temperature & salinity data from 1948 - 2014 at 54 vertical depth layers at a resolution of 0.2°x 0.2°.
#' @usage data(sst.ahoi)
"sst.ahoi"

#' Salinity data of the North Sea
#'
#' @source Thuenen Institute for Sea fisheries, "AHOI : A physical-statistical model of hydrography for fishery and ecology studies"
#'  \url{https://www.thuenen.de/en/sf/projects/ahoi-a-physical-statistical-model-of-hydrography-for-fishery-and-ecology-studies/}
#' @description A raster dataset containing the averaged salinity over the first 10m in the North Sea for the years 1990 - 2014.
#' @format A raster brick:
#' \describe{
#'  \item{dimensions}{latitude, longitude, number cells, time}
#'  \item{resolution}{spatial resolution of the data set}
#'  \item{extent}{spatial extension of the data set}
#'  \item{crs}{Coordinate reference system}
#'  \item{names}{Names of each layer in the time dimension - contains the time step, here denoted with XYYYY.MM.DD}
#' }
#' @references Núñez-Riboni I, Akimova A (2015) "Monthly maps of optimally interpolated in situ hydrography in the North Sea from 1948 to 2013." J Mar Syst 151:15–34
#' @details This dataset contains the averaged salinity over the first 10m  in the North Sea between -5.1 to 10.1 longitude and 50.5 to 62.1 latitude.
#' It is based on an extraction of the AHOI dataset (v19.02), which spans temperature & salinity data from 1948 - 2014 at 54 vertical depth layers at a resolution of 0.2°x 0.2°.
#' @usage data(salt.ahoi)
"salt.ahoi"
