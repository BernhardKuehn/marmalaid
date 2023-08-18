# =============================================== #
# function to calculate spatial-temporal SOMs
# to extract pattern from spatial-temporal fields
# author: Bernhard Kuehn
# last change: 2019-10-11
# =============================================== #

#' Self Organising Map (SOM) for a spatio-temporal field
#'
#' @description Wrapper function for \code{\link[kohonen:supersom]{som}} to perform a S-mode SOM on a spatio-temporal field organised as raster
#' @param x A raster \code{brick} with each layer corresponding to a time step specified in the argument \code{time}.
#' @param time the time step associated with each layer of the raster
#' @param seed Seed settings for reproducibility, either \code{numeric} or \code{NULL}.
#' @param parallel List, containing two elements: parallel - logical, if analysis needs to be run in parallel mode, cores - number of CPU-cores to run analysis on
#' @param reps Number of repetitions to calculate for the SOM
#' @param return.all Logical. If \code{TRUE} all SOMs (specified in 'reps' besides the 'best' one) are returned.
#' @param plot A logical value if a simple plot with the convergence of the SOM algorithm and the BMU (best matching unit) time series should be plotted.
#' @param ... additional arguments specified in \code{\link[kohonen:supersom]{som}}.
#' @return A named \code{list} with the following elements:
#' \itemize{
#'   \item SOM.out - The 'best' SOM mapping found in regards to both QE and TE. Object of \code{class(kohonen)}. See \code{\link[kohonen:supersom]{som}} for details.
#'   \item SOM.quality - The quality of the mapping assessed with both QE and TE.
#'   \item SOM.raster - \code{raster} of the spatial-SOM fields.
#'   \item Freq.pattern - the frequency of occurence of the spatial pattern within the BMU time series
#'   \item BMU.ts - a \code{data.frame} showing the time and the associated BMU
#'   \item return.all.SOMs - if \code{return.all.SOMs} is \code{TRUE}, a list with all stored SOM objects, otherwise \code{NULL}.
#'   }
#' @examples
#' library(raster)
#' data(sst.ahoi)
#' # calculate field anomaly
#' sst.anom = fieldAnomaly.raster(rst1 = sst.ahoi,
#'                             time = as.Date(sub("X","",names(sst.ahoi)),format = "%Y.%m.%d"),
#'                             level = "month")
#' # calc. mean over summer season
#'sst.JJA = calc.mean.over.Month(sst.anom,
#'                              time = as.Date(sub("X","",names(sst.anom)),
#'                                              format = "%Y.%m.%d"),
#'                               month = c(6,7,8),shiftYear = FALSE)
#' # calculate spatial SOM
#' SST.SOM = spatial.SOM(x = sst.JJA,time = as.numeric(sub("X","",names(sst.JJA))),
#'                     plot = FALSE,seed = 1234,
#'                     parallel = c(parallel = TRUE,cores = 2),
#'                     reps = 50,return.all = FALSE,
#'                     grid = kohonen::somgrid(xdim = 3,
#'                                             ydim = 3,
#'                                             topo = "rectangular",
#'                                             neighbourhood.fct = "bubble",
#'                                             toroidal = FALSE),
#'                     mode = "online",rlen = 1000)
#' plot(SST.SOM$SOM.raster)
#'
#' @importFrom foreach %dopar%

#' @export
spatial.SOM = function(x,time,plot = TRUE,seed = NULL,parallel = list(parallel = FALSE,cores = NULL),reps = 100,return.all = FALSE,...) {

  # ---------------------------------------- #
  # seed management to control stochasticity
  # ---------------------------------------- #

  if(is.null(seed)){
    warning("No internal seed set! Solution might be different for each run...")
  } else{
    if(exists(".Random.seed")){
      old.seed = .Random.seed
    } else{
      set.seed(NULL)
      old.seed = .Random.seed
    }
    on.exit(assign(".Random.seed",value = old.seed,envir = .GlobalEnv))
    # generate various seeds for the repetitions
    set.seed(seed)
    sim.seeds = stats::runif(reps,min = 1,max = 1E9)
  }

  # -------------------------------------------------------------- #
  # prepare spatio-temporal fields (raster files) for SOM-analysis
  # -------------------------------------------------------------- #

  # # x can be a raster layer or a list of raster-layers if two variables are calculated simultaniously (e.g. u and v)
  if(class(x)=="list") {
    if(all(lapply(x,class)%in%c("raster","RasterBrick","RasterStack"))) {

      # combine files to one matrix

      # get geo.position
      geo.pos = vector("list",length(x))
      for(i in 1:length(geo.pos)) {
        geo.pos[[i]] = raster::xyFromCell(x[[i]],1:raster::ncell(x[[i]]))
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
        tmp.mat[[i]] = t(raster::as.matrix(x[[i]]))
      }

      mat = do.call(cbind,tmp.mat); rm(tmp.mat)
      # bring files together in correct form

    }else{
      stop("Not all elements of list(x) are raster-files!")
    }
  } else {
    # convert to matrix
    print("1. convert raster to matrix...")
    mat = t(raster::as.matrix(x))

    # remove "land"-datapoints
    geo.pos = raster::xyFromCell(x,1:raster::ncell(x))
  }

  # check for NA column
  NA.cols = which(apply(mat, 2, function(x) all(is.na(x)))) # | var(x)== 0 maybe also if var == 0, but not applicable for currents

  # remove NA data, but only if there are NA.cols, otherwise not
  if(length(NA.cols >0)){
    mat = mat[,-NA.cols]
    geo.pos = geo.pos[-NA.cols,]
  }

  # --------------- #
  # 2. calculate SOM
  # --------------- #

  print("2. perform SOM...")

  if(parallel[[1]] == TRUE){

    if(is.null(parallel[[2]])){
      n.cores = parallel::detectCores()-1
    } else{
      n.cores = parallel[[2]]
    }
    cl = parallel::makeCluster(n.cores)
    doParallel::registerDoParallel(cl)

    # run SOM calculation in parallel mode
    cat("Running in parallel mode...","\n")
    tmp.store = foreach::foreach(i = 1:reps,.packages = "kohonen") %dopar% {
      if(!is.null(seed)){
        set.seed(sim.seeds[i])
      }
      SOM = kohonen::som(X = mat,...) #som(X = mat,somgrid(3,3,"rectangular"),
                #rlen = 10,mode = "online",dist.fcts = "euclidean") #
    }
    parallel::stopCluster(cl)
  } else{
    cat("Running in single core mode...","\n")

    tmp.store = rep(list(NA),reps)
    for(i in 1:reps){
      if(!is.null(seed)){
        set.seed(sim.seeds[i])
      }
      tmp.store[[i]] = kohonen::som(X = mat,...)
    }
  }

  # get best SOM-representation
  if(reps == 1){
    best.SOM = tmp.store[[1]]
    SOM.quality = NULL
    warning("Robustness of SOM was not evaluated!")
  } else{
    cat("Checking quality of the solutions...","\n")

    # calculate quality indices for the mapping
    QE = lapply(tmp.store,marmalaid::SOM.quant.error)
    TE = lapply(tmp.store,marmalaid::SOM.topo.error)

    # evaluate quality
    all.TEs = sapply(TE,function(x) x$TE)
    all.QEs = do.call(c,QE)

    # best SOM
    # get the one with the best QE where TE is 0
    indx.TE = which(all.TEs == min(all.TEs))
    indx.QE = which(all.QEs == min(all.QEs[indx.TE]))

    # best SOM is just some random SOM
    best.SOM = tmp.store[[indx.QE]]
    SOM.quality = data.frame(QE = all.QEs,TE = all.TEs)
  }

  # -------------- #
  # prepare output
  # -------------- #

  print("3. convert back to raster-format...")

  # pattern into raster
  SOM.space = best.SOM$codes[[1]]

  # create long table format with lat, long and loadings of first axes that explain 95% of the variance
  SOM.space.df = data.frame(geo.pos,t(best.SOM$codes[[1]]))
  # rasterize
  if(class(x)=="list") {
    split.df = split(SOM.space.df,f = rep(1:length(x),each = nrow(geo.pos)))
    SOM.raster = vector("list",length = length(x))
    for(i in 1:length(x)) {
      SOM.raster[[i]] = raster::rasterFromXYZ(split.df[[i]])
    }
  } else{
    SOM.raster = raster::rasterFromXYZ(SOM.space.df)
  }

  #frequency of occurence
  freq = prop.table(table(best.SOM$unit.classif))[order(prop.table(table(best.SOM$unit.classif)),decreasing = TRUE)]

  if(is.vector(freq)){
    freq.pattern = data.frame(Pattern = names(freq),Freq.percent = freq*100)
  } else{
    freq.pattern = data.frame(freq*100)
    names(freq.pattern) = c("Pattern","Freq.percent")
  }


  if(plot == TRUE) {
    opar = graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(opar))
    graphics::par(mfrow = c(1,2))
    # check convergence
    graphics::plot(best.SOM, type="changes")

    # # plot spatial pattern
    # plot(SOM.raster)

    # time series of pattern occurence
    graphics::plot(time,best.SOM$unit.classif,type = "l",xlab = "time",ylab = "BMU",main = "Pattern time series",las = 1)
    graphics::points(time,best.SOM$unit.classif,col = "gray50",pch = 16,cex = 0.5)
  }

  if(return.all == FALSE){
    return.all.SOMs = NULL
  } else{
    return.all.SOMs = tmp.store
  }

  # organise output in a list-structure
  output = list(SOM.out = best.SOM,
                SOM.quality = SOM.quality,
                SOM.raster = SOM.raster,
                Freq.Pattern = freq.pattern,
                BMU.ts = data.frame(time = time,ts = best.SOM$unit.classif),
                return.all.SOMs = return.all.SOMs)

  return(output)
}
