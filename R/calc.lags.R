# ---------------------------------------- #
# calculate lagged versions of time-series
# ---------------------------------------- #

#' Calculate (multiple) Lags of a time series
#'
#' @description Function to calculate lagged versions of time-series.
#' @param ts A vector of a time series.
#' @param lags A vector of lags (positive or negative).
#' @return A \code{data.frame} containing the ts at different lags of the same length as the original ts. NAs are used to fill in the missing values resulting from the shift.
#' @export
calc.lags =  function(ts,lags) {
  ts.lag = data.frame(matrix(NA,nrow = length(ts),ncol  = length(lags)))

  for(i in 1:length(lags)) {
    lag = lags[i]
    if(lag >0) {
      ts.lag[,i] = dplyr::lag(ts,abs(lag))

    } else if(lag <0) {
      ts.lag[,i] = dplyr::lead(ts,abs(lag))
    }else{
      ts.lag[,i] = ts
    }
  }
  colnames(ts.lag) = paste0("lag",lags)
  return(ts.lag)
}
