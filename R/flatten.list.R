# =============================================================== #
# function to flatten a deeply nested list with various levels
# to just one shallow one
# =============================================================== #

#' Function to flatten a deeply nested list to a shallow one
#'
#' @description Function to flatten a deeply nested list to a shallow one.
#' @param x A deeply nested \code{list}.
#' @param verbose Toggle printing of different layer nestings to console on or off if function steps through them.
#' @return A shallow \code{list}
#' @details A function that flattens a deeply nested list and returns a shallow flattened one with only one level of nesting.
#'          Be aware that the single arguments are returned somewhat out of order, since the degree of nesting can vary between elements.
#' @examples
#'# For an artificial toy dataset
#' lst = list(a = list(a = 1:10,b = list(1,2,3),c = list(a = list(1:100,NA,list(rep(NA,100))))),
#'            b = list(data.frame(id = 1:4,val = c(2,NA,3,5)),
#'                     list(letters,data.frame(m = month.abb,no = 1:length(month.abb)))))
#' lst_unnested = flatten.list(lst)
#' print(lst_unnested)
#'
#' @export
flatten.list = function(x,verbose = TRUE){
  tmp = x
  # function to flatten differently nested lists to one shallow one
  go.deeper = TRUE
  unnested.list = list()
  n = 1
  while(go.deeper == TRUE){
    if(verbose == TRUE){
      cat(n," | ")
    }
    # check if all are lists
    if(all(sapply(tmp,class) == "list")){
      tmp = unlist(tmp,recursive = FALSE)
    } else{
      # remove these from list & store in unnested list
      indx.no.list = which(sapply(tmp,class) != "list")
      unnested.list[[n]] = tmp[indx.no.list]
      tmp[indx.no.list] = NULL
      n = n+1
    }
    # check length of tmp
    if(length(tmp)==0){
      go.deeper = FALSE
    }
  }
  return(unlist(unnested.list,recursive = FALSE))
}
