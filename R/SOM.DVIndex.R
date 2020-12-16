# ============================================================= #
# function "SOM.DVIndex"
# Implementation of the DVIndex of Shen et al. 2005
# to identify the correct number of clusters in a dataset
# for SOM algorithm (in the kohonen package)
# takes as input a list of SOM objects from the kohonen package

# author: Bernhard Kuehn
# last change: 2020-01-22
# ============================================================== #

#' Dynamic Validity Index (DVI) of Shen et al. (2005) for SOMs
#'
#' @description Function to calculate the DVI for SOMs with various grid-sizes.
#' @param list.SOM.obj A \code{list} of SOM objects of \code{class(kohonen)} with different map sizes to be evaluated with the DVI.
#' @param c.param Control parameter to weight the influence of within cluster density and in-between cluster seperability.
#' @param control.for.empty.nodes Logical value, denoting if empty clusters should be discarded in the calculation of the in-between distance (\code{TRUE}) or not (\code{FALSE}).
#' @return A \code{data.frame} with the following columns:
#' \itemize{
#'   \item xdim - the y-dimension of the SOM mapping.
#'   \item ydim - the x-dimension of the SOM mapping
#'   \item size - the total number of clusters (product of xdim and ydim)
#'   \item intra.ratio - A measure of the inner cluster density.
#'   \item inter.ratio - A measure of the seperability between clusters.
#'   \item DVI - The dynamic validity index, being a combination of intra - and inter ratio.
#'   }
#' @details This function implements the DVI of Shen et al. (2005) for SOMs, which was originally proposed for kmeans clustering. It identifies the optimal number of clusters to be a compromise between intra
#' compactness of clusters and inter-seperatedness. Clusters should have a small distance to members of their own cluster (intra-ratio) as well as having a larger distance to members of other clusters (inter-ratio).
#' The DVI is now the sum of both the normalised intra- and inter ratios (optional scaling with a correction parameter is possible) and has a minimum at the prefered number of clusters. Since larger map sizes of SOMs can lead to empty clusters
#' (nodes that have no datapoint attached to them), a correction with the option \code{control.for.empty.nodes = TRUE} is possible, that only considers clusters with data attached to them.
#' @references Shen J, Chang SI, Lee ES, Deng Y, Brown SJ (2005) "Determination of cluster number in clustering microarray data." Appl Math Comput 169:1172–1185
#' @examples
#' library(raster)
#' data(sst.ahoi)
#  # calculate field anomaly
#' sst.anom = fieldAnomaly.raster(brick = sst.ahoi,
#'                             time = as.Date(sub("X","",names(sst.ahoi)),format = "%Y.%m.%d"),
#'                              level = "month")
#' # calc. mean over summer season
#' sst.JJA = calc.mean.over.Month(sst.anom,
#'                              time = as.Date(sub("X","",names(sst.anom)),
#'                                             format = "%Y.%m.%d"),
#'                              month = c(6,7,8),shiftYear = FALSE)
#'
#' # grid sizes
#' grid.size = data.frame(expand.grid(1:5,1:5))
#' names(grid.size) = c("xdim","ydim")
#' SOMs.store = list()
#' for(i in 1:nrow(grid.size)){
#'  cat(i,"\n")
#'  # calculate spatial SOM
#'  SOMs.store[[i]] = spatial.SOM(x = sst.JJA,time = as.numeric(sub("X","",names(sst.JJA))),
#'                               plot = FALSE,seed = NULL,
#'                               parallel = c(parallel = FALSE,cores = NA),
#'                               reps = 1,return.all = FALSE,
#'                               grid = kohonen::somgrid(xdim = grid.size$xdim[i],
#'                                                       ydim = grid.size$ydim[i],
#'                                                       topo = "rectangular",
#'                                                       neighbourhood.fct = "bubble",
#'                                                       toroidal = FALSE),
#'                               mode = "online",rlen = 1000)
#' }
#' # calculate DVI
#' DVI.SST.JJA = SOM.DVIndex(list.SOM.obj = lapply(SOMs.store,function(x) x$SOM.out),
#'                          c.param = 1,
#'                          control.for.empty.nodes = TRUE)
#' DVI.SST.JJA = DVI.SST.JJA[order(DVI.SST.JJA$size),]
#' # plot DVI
#' par(mfrow = c(1,1))
#' plot(1:nrow(DVI.SST.JJA),DVI.SST.JJA$DVI,type = "l",xaxt = "n",
#'      xlab = " grid size",ylab  = "DVI",las = 1)
#' points(1:nrow(DVI.SST.JJA),DVI.SST.JJA$DVI)
#' abline(v = which(DVI.SST.JJA$DVI == min(DVI.SST.JJA$DVI,na.rm = TRUE)),
#'        col = "red",lty = 2)
#' axis(side = 1,at = 1:nrow(DVI.SST.JJA),labels = DVI.SST.JJA$size)
#' title("DVI",adj = 0)


#' @export
SOM.DVIndex = function(list.SOM.obj,c.param = 1,control.for.empty.nodes = FALSE){

  # check if all elements in the list are of class "kohonen"
  if(!all(sapply(list.SOM.obj,class) %in% "kohonen")){
    stop("Check if all elements in 'list.SOM.obj' are SOM results from the kohonen package!")
  }

  eucl.dist = function(x1,x2){
    sqrt(sum((x1 - x2) ^ 2))
  }

  # get SOM architecture
  SOM.dims = lapply(list.SOM.obj,
                    function(x) data.frame(xdim = x$grid$xdim,ydim = x$grid$ydim,size = x$grid$xdim*x$grid$ydim))

  # get data
  data = list.SOM.obj[[1]]$data[[1]]
  # get codes
  codes = lapply(list.SOM.obj,function(x) x$codes[[1]])
  # get BMUs
  BMUs = lapply(list.SOM.obj,function(x) x$unit.classif)

  intra.distances = Inter = list()
  for(i in 1:length(list.SOM.obj)){
    code = codes[[i]]
    dims = SOM.dims[[i]]
    BMU = BMUs[[i]]
    intra.distances[[i]] = rep(NA,dims$size)
    sq.dist.to.all.others = list()
    for(j in 1:dims$size){
      # take the values corresponding to the respective cluster
      indx.cluster = which(BMU == j)
      partition = data[indx.cluster,]
      if(length(partition) == 0){
        intra.distances[[i]][j] = NA
        if(control.for.empty.nodes == TRUE){
           # control for empty nodes
          code[j,] = NA
        }
      } else if(is.vector(partition)){
        # calculate the the sum of the intra-cluster distances (squared distances to the node)
        intra.distances[[i]][j] = sum(eucl.dist(x1 = partition,x2 = code[j,])^2)
      } else{
        intra.distances[[i]][j] = sum(apply(partition,1,function(x) eucl.dist(x1 = x,x2 = code[j,]))^2)
      }
      # calculate inter-cluster distances (distance between nodes)
      if(length(code[-j,]) == 0){
        sq.dist.to.all.others = NA
      } else if(is.vector(code[-j,])){
        sq.dist.to.all.others[[j]] = stats::dist(code)^2
      } else{
        sq.dist.to.all.others[[j]] = apply(code[-j,],1,function(x) eucl.dist(x1 = code[j,],x2 = x))^2
      }
    }
    sum.dist.between.codes = sapply(sq.dist.to.all.others,sum,na.rm = TRUE)
    if(control.for.empty.nodes == TRUE){
        sum.dist.between.codes = sum.dist.between.codes[sum.dist.between.codes>0]
    }
    Inter[[i]] = (max(unlist(sq.dist.to.all.others),na.rm = TRUE)/min(unlist(sq.dist.to.all.others),na.rm = TRUE))*sum(1/sum.dist.between.codes)
  }

  # calculate the mean
  Intra = sapply(intra.distances,mean,na.rm = TRUE) # to take care of the empty clusters (that can happen during SOM formation)
  Inter = do.call(c,Inter)

  max.intra = max(Intra,na.rm = TRUE)
  max.inter =  max(Inter,na.rm = TRUE)
  # calc ratios
  intra.ratio = Intra/max.intra
  inter.ratio = Inter/max.inter
  out = data.frame(do.call(rbind,SOM.dims),
                   intra.ratio = intra.ratio,
                   inter.ratio = inter.ratio,
                   DVI = intra.ratio+c.param*inter.ratio)
  return(out)
}
