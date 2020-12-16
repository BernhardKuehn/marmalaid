# =============================================================== #
# function to calculate seasonal means over defined period (month)
# =============================================================== #

#' Seasonal mean for spatio-temporal field
#'
#' @description Function to calculate the mean over subsequent month in a raster brick.
#' @param raster A raster \code{brick} with each layer corresponding to a time step.
#' @param time A vector of dates corresponding to each layer in the raster \code{brick}.
#' @param month A vector of numbers between 1 and 12 denoting the month over which to calculate the average (e.g. 3,4,5). Only subsequent month are allowed (e.g \code{c(12,1,2)} works; \code{c(1,3,5)} does not!).
#' @param shiftYear A logical value. How to handle averages over subsequent month but different years (e.g. Dec.-Jan.-Feb.). If \code{TRUE} the time returned corresponds to the older year (e.g. 1989), otherwise to the younger one (1990).
#' @return A raster \code{brick} containing the avg. spatial fields per specified monthly period per year.
#' @examples
#'## For an artificial toy dataset
#' library(raster)
#' r = raster::raster(extent(c(-5,5,40,60)),
#'            res = c(0.1,0.1))
#' date = seq(as.Date("1990-01-01"),
#'            as.Date("2000-12-31"),by = "month")
#' br = raster::brick(lapply(seq_along(date), function(i) setValues(r,runif(ncell(r)))))
#' calc.mean.over.Month(raster = br, time = date,month = c(12,1,2),shiftYear = TRUE)
#' calc.mean.over.Month(raster = br, time = date,month = c(12,1,2),shiftYear = FALSE)
#' ## For SST of the North Sea
#' data(sst.ahoi)
#' SST.winter = calc.mean.over.Month(raster = sst.ahoi,
#'                                   time = as.Date(sub("X","",names(sst.ahoi)),format = "%Y.%m.%d"),
#'                                   month = c(12,1,2),shiftYear = TRUE)
#' plot(SST.winter,1:6)
#'
#' @export
calc.mean.over.Month = function(raster,time,month,shiftYear = TRUE) {

  # check if times and raster are the same size
  if(dim(raster)[3] != length(time)) {
    stop("Raster and time-object are not the same length!")
  }

  # order month (based on shortest distance)
  shortest.dist.month = function(a,b){
    dist1 = (12-max(a,b))+min(a,b)
    dist2 = (b-a)
    # get shortest distance
    indx = which(abs(c(dist1,dist2)) == min(abs(c(dist1,dist2))))[1]
    # direction
    if(indx == 1 & dist2 <0){
      return(c(dist1,dist2)[indx])
    } else if(indx == 1 & dist2 >0){
      return(-1*c(dist1,dist2)[indx])
    }
    return(c(dist1,dist2)[indx])
  }
  month.order = month[order(month)]
  if(any(abs(diff(month.order)) !=1)){
    a = TRUE
    while(a == TRUE){

      nr = month.order[1]
      diff1 = rep(NA,length(month.order)-1)
      for(i in 1:length(month.order)){
        diff1[i] = shortest.dist.month(nr,month.order[i])
      }
      new.order = order(diff1)
      month.order.new = month.order[new.order]

      diff.sum = diff.sum.new = rep(NA,length(month.order)-1)
      for(i in 1:(length(month.order)-1)){
        diff.sum[i] = shortest.dist.month(month.order[i],month.order[i+1])
        diff.sum.new[i] = shortest.dist.month(month.order.new[i],month.order.new[i+1])
      }
      if(sum(abs(diff.sum))< sum(abs(diff.sum.new)) | all(month.order == month.order.new)){
        a = FALSE
      } else{
        month.order = month.order.new
      }
    }
    # check if ordered month are subsequent
    if(any(abs(diff.sum) != 1)){
      stop("Only subsequent months allowed!")
    }
  }
  month = month.order
  cat("Calc mean over:",month.name[month],"...","\n")

  # drop !month
  drop.month = c(1:12)[!(1:12 %in% month)]
  drop.indx = which(as.numeric(format(time,"%m")) %in% drop.month)
  raster.drop = raster::dropLayer(raster,drop.indx)
  time.drop = time[-drop.indx]

  if(any(abs(diff(month))>1)) {
    # apply mean over respective month
    reps = rle(as.numeric(format(time,"%m")) %in% month)$lengths * rle(as.numeric(format(time,"%m")) %in% month)$values
    reps = reps[reps>0]
    indx = rep(1:length(reps),reps)
    e = 0
    for(i in unique(indx)) {
      s = ifelse(i == 1,1,e+1)
      e = cumsum(reps)[i]
      cat("year:", format(time.drop[s:e],"%Y"),"|","month:",substr(format(time.drop[s:e],"%b"),1,1),"\n")
    }
    month.mean = raster::stackApply(raster.drop,indx,mean) # calculate mean over respective month
    if(reps[1] < length(month)) {
      if(shiftYear == TRUE) { # added on the 2019-10-11 to control the time-stamp if month span over different years
        names(month.mean) = c(as.numeric(unique(format(time.drop,"%Y")))-1,max(unique(format(time.drop,"%Y"))))
      } else{
        if(reps[length(reps)]< length(month)) {
          month.mean = month.mean[[-length(reps)]]
        }
        names(month.mean) = as.numeric(unique(format(time.drop,"%Y")))
      }

    } else {
      names(month.mean) =  unique(format(time.drop,"%Y"))
    }
  } else {
    # calculate mean over respective month
    month.indx = which(as.numeric(format(time,"%m")) %in% month)
    indx = as.numeric(as.factor(as.numeric(format(time,"%Y"))[month.indx]))
    for(jj in unique(indx)){
      ii = which(indx == jj)
      cat("year:", format(time[month.indx[ii]],"%Y"),"|","month:",
          substr(format(time[month.indx[ii]],"%b"),1,1),"\n")
    }
    month.mean = raster::stackApply(raster.drop,indx,mean)
    names(month.mean) = unique(format(time,"%Y")[month.indx])
  }
  return(month.mean)
}

