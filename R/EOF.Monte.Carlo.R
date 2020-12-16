# ----------------------------------------------------------------------- #
# permute columns of pca-input to get a statistically valid number of PCs
# for further analysis
# last change: 2020-06-30
# ----------------------------------------------------------------------- #

#' Finds Nr. of PCs to retain by means of Monte-Carlo analysis
#'
#' @description A form of Horns parallel analysis to compare the eigenvalues of a PCA with those of a null-model (randomly permuted colums).
#' @param raster A raster brick from which to calculate the spatial-PCA.
#' @param center If the mean should be removed prior to analysis.
#' @param scale If each grid point should be scaled by its variance.
#' @param size Nr. of permutations to generate for the Null-model.
#' @param CI.size Nr. of bootstrap samples for calculating confidence-intervals.
#' @param plot A logical value, indicating if the eigenvalues of the Null-model, PCA-analysis and CIs are plotted.
#' @importFrom foreach %dopar%
#' @return A \code{list} with the following elements:
#' \itemize{
#'   \item eigen.data - Eigenvalues of the original dataset
#'   \item eigen.data.bootstrap - Eigenvalues of the bootstrapped dataset
#'   \item eigen0 - Eigenvalues of the Null-model
#'   \item Nr.PCs - the number of PCs to retain based on comparing the eigenvalues of the Null-model with the bootstrapped eigenvalues.
#'   }
#' @details ...
#' @references ...

#' @export
EOF.Monte.Carlo = function(raster,center = TRUE,scale = FALSE,size = 100,CI.size = 1000,plot = TRUE){

  # S-mode EOF-analysis
  if(class(raster) == "list"){
    tmp.mat = vector("list",length(raster))
    for(i in 1:length(raster)) {
      tmp.mat[[i]] = raster::as.matrix(raster[[i]])
      complete = stats::complete.cases(tmp.mat[[i]])
      tmp.mat[[i]] = tmp.mat[[i]][complete,]
    }
    data = do.call(rbind,tmp.mat)
    coords = as.data.frame(raster::coordinates(raster[[1]])[complete,])
    rm(tmp.mat)
  } else{
    data = raster::as.matrix(raster)
    complete = stats::complete.cases(data)
    coords = as.data.frame(raster::coordinates(raster)[complete,])
    data = data[complete,]
  }

data = t(data)

# scale or not scale
if(scale == TRUE){
  data = apply(data,2,scale)
}

# calculate EOF
eof.out = stats::prcomp(data,center = center,scale. = scale)
boot.out = boot_prcomp_parallel(valores = data,size = CI.size,alpha = 0.05,center = center,scale = scale,parallel = TRUE)

# calculate Monte Carlo simulation and compare to EOF or original data
eof0 = EOF.null.Monte.Carlo(mat = data,center = center, scale = scale,reps = size)
eigen.data = eof.out$sdev^2
eigen0 = eof0$eigen0

# return number of EOFs selected
Nr.PCs = which(apply(eigen0,2,mean)>boot.out$eigen_quantiles[,1])[1]-1

if(plot == TRUE){
  opar = graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(opar))
  graphics::par(mfrow = c(1,1),mar = c(4,4,3,3),oma = c(0,0,0,0))
  graphics::plot(eigen.data,type = "b",xlab = "eigenvalues",pch = 16,ylab = "",ylim = c(min(boot.out$eigen_quantiles[,1]),max(boot.out$eigen_quantiles[,2])))
  graphics::polygon(c(rev(1:nrow(boot.out$eigen_quantiles)),1:nrow(boot.out$eigen_quantiles)),
          c(rev(boot.out$eigen_quantiles[,1]),boot.out$eigen_quantiles[,2]),col=scales::alpha("gray50",0.1),border=NA)
  graphics::lines(boot.out$eigen_quantiles[,1],col = "gray50",lty = 1)
  graphics::lines(boot.out$eigen_quantiles[,2],col = "gray50",lty = 1)
  graphics::lines(apply(eigen0,2,mean),type = "b",col  ="red",pch = 16)
  graphics::abline(v = Nr.PCs,lty = 2,col = "gray50")
  graphics::text(x = Nr.PCs+(0.025*length(eigen.data)),y = max(boot.out$eigen_quantiles[,2])*0.8,Nr.PCs,col = "gray50")
  graphics::title("PCs selected",adj = 0)
  graphics::legend("topright",legend = c("data set","Null model","Nr. PCs\n selected"),title = "Eigenvalues",cex = 0.9,
         pch = c(16,16,NA),lty = c(NA,NA,2),col = c("black","red","gray50"),border = NA,bty = "n")
}

return(list(eigen.data = eigen.data,
            eigen.data.bootstrap = boot.out,
            eigen0 = eigen0,
            Nr.PCs = Nr.PCs))
}


# function to calculate CI of Eigenvalues
# parallel version

boot_prcomp_parallel = function (valores, size = 1000, alpha = 0.05,center = TRUE,scale = TRUE,parallel = FALSE) {
  if (is.matrix(valores) == FALSE) {
    return("Object is not a matrix")
  }
  else {
    lengs = length(valores[, 1])
    hlengs = length(valores[1, ])
    store_bootstrap = matrix(0, size, lengs)
    elements_aux = NULL
    indexes = seq(1, lengs, 1)
    work_matrix = matrix(0, lengs, hlengs)
    empirical_quantiles1 = matrix(0, lengs - 1, 2)
    store_eigen = matrix(0, size, lengs)
    empirical_quantiles2 = matrix(0, lengs, 2)
    if(parallel == FALSE){
      for (sim in 1:size) {

        resample_indexes = sample(indexes, length(indexes),
                                  replace = TRUE)
        for (ips in 1:length(resample_indexes)) {
          work_matrix[ips, ] = valores[resample_indexes[ips],
                                       ]
        }
        #cova = (lengs - 1) * cov(work_matrix)/(lengs)
        #auto_va = unlist(eigen(cova)$values)
        # faster way to do this: svd
        #auto_va = svd(cova,nu = hlengs)$d[1:hlengs]
        auto_va = stats::prcomp(work_matrix,center = TRUE,scale = scale)$sdev^2
        cum_eigen = cumsum(auto_va)
        total_var = sum(auto_va)
        percent = (cum_eigen/total_var)
        store_bootstrap[sim, ] = percent
        store_eigen[sim, ] = auto_va
      }
    } else{
      nc = parallel::detectCores()-1
      cl = parallel::makeCluster(nc)
      doParallel::registerDoParallel(cl)

      auto_va = foreach::foreach(i = 1:size,.combine = rbind) %dopar% {
        resample_indexes = sample(indexes, length(indexes),
                                  replace = TRUE)
        for (ips in 1:length(resample_indexes)) {
          work_matrix[ips, ] = valores[resample_indexes[ips],]
        }
        auto_va = stats::prcomp(work_matrix,center = center,scale = scale)$sdev^2
      }
      parallel::stopCluster(cl)

      cum_eigen = apply(auto_va,1,cumsum)
      total_var = apply(auto_va,1,sum)
      percent = (cum_eigen/total_var)
      store_bootstrap = percent
      store_eigen = auto_va
    }
    store_bootstrap = store_bootstrap[, -lengs]
    for (iter in 1:(lengs - 1)) {
      empirical_quantiles1[iter, 2] = stats::quantile(store_bootstrap[,
                                                               iter], 1 - alpha/2)
      empirical_quantiles1[iter, 1] = stats::quantile(store_bootstrap[,
                                                               iter], alpha/2)
    }
    for (iter in 1:(lengs)) {
      empirical_quantiles2[iter, 2] = stats::quantile(store_eigen[,
                                                           iter], 1 - alpha/2)
      empirical_quantiles2[iter, 1] = stats::quantile(store_eigen[,
                                                           iter], alpha/2)
    }
    colnames(empirical_quantiles1) = c(paste(alpha/2 * 100,
                                             "%"), paste((1 - alpha/2) * 100, "%"))
    colnames(empirical_quantiles2) = c(paste(alpha/2 * 100,
                                             "%"), paste((1 - alpha/2) * 100, "%"))
    lista = list(proportions_quantiles = empirical_quantiles1,
                 proportions_used = store_bootstrap, eigen_quantiles = empirical_quantiles2,
                 eigen_used = store_eigen)
    return(lista)
  }
}

# function to reshuffle columns
EOF.null.Monte.Carlo = function(mat,center = TRUE,scale = FALSE,reps = 1000){
  eigen0 = matrix(NA,reps,nrow(mat))
  EOF0 = list()
  for(j in 1:reps){
    reshuffle.mat = matrix(NA,nrow(mat),ncol(mat))
    for(i in 1:ncol(mat)){
      shuffle = sample(mat[,i],nrow(mat),replace = FALSE)
      reshuffle.mat[,i] = shuffle
    }
    out = stats::prcomp(reshuffle.mat,center = center,scale. = scale)
    eigen0[j,] = out$sdev^2
    EOF0[[j]] = out$rotation
  }
  return(list(eigen0 = eigen0,EOF0 = EOF0))
}

