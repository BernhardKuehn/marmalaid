# ============================ #
# Multi-obj. fitness functions
# ============================ #

#' @name fitnessfuns
#' @aliases MO.fitness.func.RF.regression
#' @aliases MO.fitness.func.extraTrees.ranger
#'
#' @title Random Forest & extreme randomized Trees multi-obj. fitness functions for use with NSGA-II feature selection
#' @description Various multi-obj. fitness functions for use with NSGA-II feature selection
#' @param feat The binary coded feature vector that denotes which gene (feature) is on (1) or off (0).
#' @param x The dataframe containing all the features to run feature selection on.
#' @param y Numeric vector of the response variable.
#' @param seed The setting for the seed.
#' @param metric.func Slot for the metric function to be included. The format is metric.func(data) see \link[marmalaid]{metricfuns} for details.
#' @param mtry Number of variables to possibly split at in each node. Default is the square root of the Number of Variables.
#' @param nTree Number of Trees to build in Forest for Random Forest/ Extreme Randomized Trees algorithm.
#' @param folds Number of folds (k) for n-repeated k-fold crossvalidation.
#' @param reps  Number of repetitions (n) for n-repeated k-fold crossvalidation.
#' @param only.obs.vs.pred A logical value. Setting to return only observations vs. predictions without any metric-function being considered.
#' @param return.scaled.metric  A logical value. Setting if the metric should be scaled by a naive forecast (mean of all observations).
#' @details \code{MO.fitness.func.RF.ranger}  - Fitness function for the "randomForest" - algorithm from the \link[ranger]{ranger} package.
#' \cr \code{MO.fitness.func.extraTrees.ranger}  - Fitness function for the "extraTrees" - algorithm from the \link[ranger]{ranger} package.
#' @importFrom foreach %do%
#' @note \code{fitnessfuns} is a generic name for the functions documented.
#'
# # --------------------------------------------------------------- #
# # Random forest for a regression problem using the ranger package
# # --------------------------------------------------------------- #
#' @rdname fitnessfuns
#' @export
MO.fitness.func.RF.ranger = function(feat,
                                      x,y,
                                      seed,
                                      metric.func,
                                      mtry = NULL,
                                      nTree = 200,
                                      folds = 3,
                                      reps = 20,
                                      only.obs.vs.pred = FALSE,
                                      return.scaled.metric = TRUE){

  if(exists(".Random.seed")){
    old.seed = .Random.seed
  } else{
    set.seed(NULL)
    old.seed = .Random.seed
  }
  on.exit(assign(".Random.seed",value = old.seed,envir = .GlobalEnv))

  # select subset
  x.subset = x[,which(feat == 1),drop = FALSE]

  set.seed(seed)
  cc.folds = caret::createMultiFolds(y = y, k = folds,times = reps)

  # calculate fraction of features
  frac.feat = sum(feat)/length(feat)
  # control for empty set
  if(frac.feat == 0){
    return(c(metric = Inf,params =  frac.feat))
  }

  if(is.null(mtry)){
    mtry= floor(sqrt(sum(feat)))
  }

  # perform cv
  i = 1:length(cc.folds)
  out = foreach::foreach(i = i,.packages = c("caret","ranger"),.multicombine = TRUE) %do% {

    # create train & test
    y.train = y[cc.folds[[i]]]
    y.test = y[-cc.folds[[i]]]
    x.train = x.subset[cc.folds[[i]],,drop = FALSE]
    x.test = x.subset[-cc.folds[[i]],,drop = FALSE]

    xy.train = data.frame(y.train = y.train,x.train)
    xy.test = data.frame(y.test = y.test,x.test)

    # train classifier (RandomForest)
    set.seed(seed)
    rf.out = ranger::ranger(y.train~.,data = xy.train,num.trees =  nTree,probability = FALSE,mtry = mtry)
    out = data.frame(pred = stats::predict(rf.out,x.test)$predictions,obs = y.test)
  }

  # combine all obs and pred to one
  data = data.table::rbindlist(out)

  # option to only return the predicted vs. observed CV results
  if(only.obs.vs.pred == TRUE | is.null(metric.func)){
    return(data)
  }

  # calculate metric over all folds and repetitions (to minimize the metric)
  metric = metric.func(data)

  # mean
  if(return.scaled.metric == TRUE){
    naive.metric = metric.func(data = data.frame(pred = mean(y,na.rm = TRUE),obs = y))
    rel.metric = metric/naive.metric
    # control for predictions, which are worse than the mean
    rel.metric = ifelse(rel.metric>1,1,rel.metric)
  } else{
    rel.metric = metric
  }

  return.out = c(metric = rel.metric,params =  frac.feat)
  return(return.out)
}

# ------------------------ #
# extreme randomized trees
# ------------------------ #
#' @rdname fitnessfuns
#' @export
MO.fitness.func.extraTrees.ranger = function(feat,x,y,
                                             seed,
                                             metric.func,
                                             mtry = NULL,
                                             nTree = 200,
                                             folds = 3,
                                             reps = 20,
                                             only.obs.vs.pred = FALSE,
                                             return.scaled.metric = TRUE){

  if(exists(".Random.seed")){
    old.seed = .Random.seed
  } else{
    set.seed(NULL)
    old.seed = .Random.seed
  }
  on.exit(assign(".Random.seed",value = old.seed,envir = .GlobalEnv))

  # select subset
  x.subset = x[,which(feat == 1),drop = FALSE]

  set.seed(seed)
  cc.folds = caret::createMultiFolds(y = y, k = folds,times = reps)

  # calculate fraction of features
  frac.feat = sum(feat)/length(feat)
  # control for empty set
  if(frac.feat == 0){
    return(c(metric = Inf,params =  frac.feat))
  }

  if(is.null(mtry)){
    mtry= floor(sqrt(sum(feat)))
  }

  # perform cv
  i = 1:length(cc.folds)
  out = foreach::foreach(i = i,.packages = c("caret","ranger"),.multicombine = TRUE) %do% {

    # create train & test
    y.train = y[cc.folds[[i]]]
    y.test = y[-cc.folds[[i]]]
    x.train = x.subset[cc.folds[[i]],,drop = FALSE]
    x.test = x.subset[-cc.folds[[i]],,drop = FALSE]

    xy.train = data.frame(y.train = y.train,x.train)
    xy.test = data.frame(y.test = y.test,x.test)

    # extreme randomized Trees regression with ranger
    extraT.out = ranger::ranger(y.train~.,data = xy.train,num.trees = nTree,mtry = mtry,splitrule = "extratrees",num.random.splits = 1)

    out = data.frame(pred = stats::predict(extraT.out,x.test)$predictions,obs = y.test)
  }

  # combine all obs and pred to one
    data = data.table::rbindlist(out)

  # option to only return the predicted vs. observed CV results
  if(only.obs.vs.pred == TRUE | is.null(metric.func)){
    return(data)
  }

  # calculate metric over all folds and repetitions (to minimize the metric)
  metric = metric.func(data)

  # mean
  if(return.scaled.metric == TRUE){
    # return a metric scaled between 0 and 1
    naive.metric = metric.func(data = data.frame(pred = mean(y,na.rm = TRUE),obs = y))
    rel.metric = metric/naive.metric
    # control for predictions, which are worse than the mean
    rel.metric = ifelse(rel.metric>1,1,rel.metric)
  } else{
    rel.metric = metric
  }
  return.out = c(metric = rel.metric,params =  frac.feat)
  return(return.out)
}
