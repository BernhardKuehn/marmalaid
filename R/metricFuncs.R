# =============== #
# various Metrics
# =============== #

#' @name metricfuns
#' @aliases RRSE
#' @aliases RMSE
#' @aliases RMSPE
#' @aliases MAE
#' @aliases PseudoR2
#'
#' @title Various metric functions for performance evaluation
#' @description Various metric functions to assess the fit of a statistical model in a regression context.
#' @param data A \code{dataframe} containing a column with numeric observations (\code{obs}) and the prediction (\code{pred}).
#' @details {\code{RMSE} calculates the "Root mean squared error". See \href{https://en.wikipedia.org/wiki/Root-mean-square_deviation}{Wikipedia} for details. \cr
#'          \code{MAE} calculates the "Mean absolute error". See \href{https://en.wikipedia.org/wiki/Mean_absolute_error}{Wikipedia} for details. \cr
#'          \code{PseudoR2} calculates a "R2"-like metric as cor(obs,pred)^2. \cr
#'          \code{RMSPE} calculates the "Root mean squared percentage error". See \href{https://stats.stackexchange.com/questions/413249/what-is-the-correct-definition-of-the-root-mean-square-percentage-error-rmspe}{Stackexchange} for details. \cr
#'          \code{RRSE} calculates the "Root relative squared error", which is the RMSE relative to the RMSE of a naive prediction by using the deviations to the mean.
#'}
#' @note \code{metricfuns} is a generic name for the functions documented.


# -------------------------------------------------- #
# root relative squared error (RRSE)  for regression
# -------------------------------------------------- #
#' @rdname metricfuns
#' @export
RRSE = function(data){
  sq.residuals = (data$obs - data$pred)^2
  sq.naive.prediction = (data$obs-mean(data$obs))^2
  RRSE = sqrt(sum(sq.residuals)/sum(sq.naive.prediction))
  return(RRSE)
}

# -------------------------------------------------- #
# Root mean square percentage (fraction error) error
# -------------------------------------------------- #
#' @rdname metricfuns
#' @export
RMSPE = function(data){
  x.rel = (data$pred/data$obs) - 1
  RMSPE = sqrt(mean(x.rel^2))#*100
  return(RMSPE)
}

# ---- #
# RMSE
# ---- #
#' @rdname metricfuns
#' @export
RMSE = function(data){
  sqrt(mean((data$pred - data$obs)^2))
}

# --- #
# MAE
# --- #
#' @rdname metricfuns
#' @export
MAE = function(data){
  mean(abs(data$pred - data$obs))
}

# --------- #
# Pseudo-R2
# --------- #
#' @rdname metricfuns
#' @export
PseudoR2 = function(data){
  stats::cor(data$pred,data$obs)^2
}


