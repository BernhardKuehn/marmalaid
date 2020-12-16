#' Calculates Model metrics (RMSE,MAE,Pseudo-R2) for output of feature selection NSGA-II
#'
#' @description Function to calculate Model metrics (RMSE,MAE,Pseudo-R2) for output of feature selection NSGA-II
#' @param NSGA.II.obj The output object from the NSGA-II feature selection.
#' @param fitness.func A function used as fitness function in the feature selection process.
#' @return A \code{data.frame} specifying the solutions within the pareto front, with the following columns:
#' \itemize{
#'   \item ID - the index in 'NSGA.II.obj$pareto.varnames' to which each solution corresponds to
#'   \item metric - the metric used to evaluate the fitness of the solution.
#'   \item params - a number between 0 and 1 representing the scaled number of parameters used for feature selection
#'   \item nr.of.params - An integer. Specifies the number of covariates used in the model.
#'   \item dist.to.CoordOrigin - distance of the solution within the pareto front to the coordinate origin (0,0)
#'   \item RMSE - the RMSE (with the same level of Cross validation, which was used during the feature selection)
#'   \item PseudoR2 - the pseudo R2 (with the same level of Cross validation, which was used during the feature selection)
#'   \item MAE - the MAE (with the same level of Cross validation, which was used during the feature selection)
#'   }

#' @export
Model.metrics.regression_NSGAII = function(NSGA.II.obj,fitness.func){

  # get best solutions and fitness function
  feat.pareto = NSGA.II.obj$pareto.individuals
  pareto.solution = NSGA.II.obj$pareto.solution
  x = NSGA.II.obj$fitted.dataset$x
  y = NSGA.II.obj$fitted.dataset$y
  seed = NSGA.II.obj$algorithm.params$seed
  folds = NSGA.II.obj$algorithm.params$CV["folds"]
  reps = NSGA.II.obj$algorithm.params$CV["reps"]

  # ------------------------------------------------------------------------ #
  # various regression metrics that can be plugged into the fitness function
  # ------------------------------------------------------------------------ #

  # ---- #
  # RMSE
  # ---- #

  RMSE = function(data){
    sqrt(mean((data$pred - data$obs)^2))
  }

  # --------- #
  # Pseudo-R2
  # --------- #

  PseudoR2 = function(data){
    stats::cor(data$pred,data$obs)^2
  }

  # --- #
  # MAE
  # --- #

  MAE = function(data){
    mean(abs(data$pred - data$obs))
  }

  metric.funcs = list(RMSE = RMSE, PseudoR2 = PseudoR2,MAE = MAE)

  # --------------------- #
  # calculate metrics
  # --------------------- #

 all.scores = list()
  for(i in 1:nrow(feat.pareto)){
    if(sum(feat.pareto[i,]) == 0){
      all.scores[[i]]  = NA
    } else{
      data = fitness.func(feat = feat.pareto[i,],
                          x = x,
                          y = y,
                          folds = folds,
                          reps = reps,
                          seed = seed,
                          metric.func = NULL,
                          only.obs.vs.pred = TRUE)
      # metrics
      all.scores[[i]] = sapply(metric.funcs,function(x) x(data = data))

    }
  }

  all.scores = do.call(rbind,all.scores)

  # bring together with pareto.solutions
  out = cbind(pareto.solution,all.scores)
  rownames(out) = rownames(pareto.solution)

  # return
  return(data.frame(out))
}
