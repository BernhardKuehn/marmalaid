#' @name SOM.qualityMetrics
#' @aliases SOM.quant.error
#' @aliases SOM.topo.error
#'
#' @title Quality metrics for the map representation of a SOM
#' @description Various quality metrics to assess the quality of the SOM-mapping returned by \code{\link[kohonen:supersom]{som}}.
#' @param SOM.obj An object of \code{class(kohonen)} produced via \code{\link[kohonen:supersom]{som}}.
#' @details \code{SOM.quant.error}  -  Calculates the Quantisation error (QE) for the SOM mapping
#' \cr \code{SOM.topo.error}  -  Calculates the Topographic error (TE) (Kiviluoto 1996) and a combined index of QE & TE based on the work of Kaski & Lagus 1996
#'
#' @note \code{SOM.qualityMetrics} is a generic name for the functions documented.
#' @return
#' ## \code{SOM.quant.error}
#' A \code{numeric} value for the QE
#' ## \code{SOM.topo.error}
#' A named \code{list} with following elements: \itemize{
#'                        \item TE - a numeric value for the Topographic error
#'                        \item C.index - a numeric value for the combined metric of QE & TE error from Kaski & Lagus 1996
#'                         \item Mapping.of.data - A \code{data.frame} showing the best matching unit (BMU), the second best matching unit (SNDMU), the C.index and an integer value denoting if BMU and SNDMU are neighbours (1) or not (0) for every data point.
#' }
#' @examples
#' data(sst.ahoi)
#' # calculate field anomaly
#' sst.anom = fieldAnomaly.raster(brick = sst.ahoi,
#'                             time = as.Date(sub("X","",names(sst.ahoi)),format = "%Y.%m.%d"),
#'                             level = "month")
#' # calc. mean over summer season
#' sst.JJA = calc.mean.over.Month(sst.anom,
#'                              time = as.Date(sub("X","",names(sst.anom)),
#'                                              format = "%Y.%m.%d"),
#'                               month = c(6,7,8),shiftYear = FALSE)
#' # calculate spatial SOM
#' SST.SOM = spatial.SOM(x = sst.JJA,time = as.numeric(sub("X","",names(sst.JJA))),
#'                     plot = FALSE,seed = 1234,
#'                     parallel = c(parallel = TRUE,cores = 2),
#'                     reps = 50,return.all = FALSE,
#'                     grid = kohonen::somgrid(xdim = 2,ydim = 3,
#'                                             topo = "rectangular",
#'                                             neighbourhood.fct = "bubble",
#'                                             toroidal = FALSE),
#'                     mode = "online",rlen = 1000)
#' # calculate QE
#' SOM.quant.error(SST.SOM$SOM.out)
#' # calculate TE & C-index
#' SOM.topo.error(SST.SOM$SOM.out)

#' @references  Kiviluoto K (1996) "Topology preservation in self-organizing maps." In: Proceedings of International Conference on Neural Networks (ICNN’96). IEEE, p 294–299 \cr
#'              Kaski S, Lagus K (1996) "Comparing self-organizing maps." In: von der Malsburg C , von Seelen W , Vorbrüggen JC, Sendhoff B (eds) Artificial Neural Networks — ICANN 96. ICANN 1996. Lecture Notes in Computer Science, vol 1112. Springer, Berlin, Heidelberg.p 809–814 \cr
#'              Pötzlbauer G (2004) "Survey and comparison of quality measures for self-organising maps." In: Proc. 5th Workshop on Data Analysis (WDA 2004).p 67–82
#'
#'
#' @rdname SOM.qualityMetrics
#' @export
SOM.quant.error = function(SOM.obj){
  # mean of the distances of the fitted SOM-obj.
  QE = sum(SOM.obj$distances)/length(SOM.obj$distances)
  return(QE)
}

#' @rdname SOM.qualityMetrics
#' @export
SOM.topo.error = function(SOM.obj){

  # data
  data = SOM.obj$data[[1]]

  # pattern the data was mapped to
  pattern = SOM.obj$codes[[1]]

  # distance of data to the BMUs
  data.to.Nodes = data.frame(BMU = rep(NA,nrow(data)),
                             SNDMU = rep(NA,nrow(data)),
                             C.index = rep(NA,nrow(data)))

  for(i in 1:nrow(data)){
    dp = data[i,]
    dist.dp = rep(NA,nrow(pattern))
    # calculate distance of each data point to each Node
    for(j in 1:nrow(pattern)){
      mat = rbind(dp,pattern[j,])
      dist.dp[j] = stats::dist(mat)
    }
    # calculate best and second best matching pattern and distance to best
    best.mu = which(dist.dp == min(dist.dp))
    d.bmu.x = dist.dp[best.mu]

    dist.dp[best.mu] = NA
    snd.mu = which(dist.dp == min(dist.dp,na.rm = TRUE))

    # if distance to two or more nodes is equal take the first one
    if(length(snd.mu)>1){
      snd.mu = snd.mu[1]
    }

    # calculate distance between best and second best matching unit
    dist.mat.BMU = as.matrix(stats::dist(pattern,upper = TRUE,diag = TRUE)) # object.distances(SOM.obj,type = "codes")
    d.bmu.sndmu = dist.mat.BMU[best.mu,snd.mu]

    # calculate combined index
    C.index = d.bmu.x + d.bmu.sndmu

    # write in data.frame
    data.to.Nodes[i,"BMU"] = best.mu
    data.to.Nodes[i,"SNDMU"] = snd.mu
    data.to.Nodes[i,"C.index"]  = C.index
  }

  # distance between BMUs on the grid
  dists.2D = kohonen::unit.distances(SOM.obj$grid, toroidal = SOM.obj$grid$toroidal)
  # check if BMU & SNDMU are neighbours
  data.to.Nodes$neighbours = apply(data.to.Nodes,1,function(x) ifelse(dists.2D[x[1],x[2]] <= 1,0,1))

  # calculate topographic error
  TE = sum(data.to.Nodes$neighbours)/nrow(data.to.Nodes)
  # Calculate mean C-index over all nodes
  C.index.mean = mean(data.to.Nodes$C.index)

  return(list(TE = TE,
              C.index = C.index.mean,
              Mapping.of.data =  data.to.Nodes))
}
