# ===================================================== #
# function to merge multiple data.frames/lists together
# ===================================================== #

#' Merge multiple tables
#'
#' @description Function to merge various data.frames/list of data.frames to one.
#' @param fileA A string denoting the name of a data.frame
#' @param listFileB A character vector denoting data.frames or lists containg data.frames, which are merged with fileA.
#' @param by Variable name by which to merge.
#' @param by.firstCol A logical value. If merging should be done by the first column in each data.frame regardless of the name.
#' @return A merged \code{data.frame}.
#' @export
mult.merge = function(fileA,listFileB,by,by.firstCol = FALSE) {
  # function to merge multiple lists and dataframes
  a = get(fileA)

  if(missing(by) & by.firstCol == TRUE) {
    by.x = colnames(a)[1]
  } else {
    by.x = by.y = by
  }
  for(i in listFileB) {
    print(i)
    b = get(i)
    if(class(b) == "data.frame") {
      if(by.firstCol == TRUE) {
        by.y = colnames(b)[1]
      }
      a = merge(a,b,by.x = by.x,by.y = by.y,all.x = TRUE)
    } else if(class(b) == "list") {
      for(j in 1:length(b)) {
        bb = b[[j]]
        if(!is.null(names(b))){
          names(bb) = c(names(bb)[1],paste(names(b)[j],names(bb)[-1],sep = "_"))
          #by.y = paste(names(b)[j],by,sep = "_")
        }
        if(class(bb) != "data.frame") {
          stop("List elements are not Data.frames!")
        }
        if(by.firstCol == TRUE) {
          by.y = colnames(bb)[1]
        }
        a = merge(a,bb,by.x = by.x,by.y = by.y,all.x = TRUE)
      }
    } else {
      stop("File not in appropriate format (List or data.frame)!")
    }
  }
  return(a)
}
