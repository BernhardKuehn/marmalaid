# ====================================================================== #
# NSGA.II_wrapper(), plot_NSGAII()
# flexible wrapper-function to perform multi-objective feature selection
# using NSGA-II with user-defined fitness functions
# author: Bernhard Kuehn
# last change: 2019-12-13
# ====================================================================== #

#' NSGA-II Wrapper function to perform multi-objective (2) feature selection
#'
#' @description Performs Bi-objective Feature selection via a specified fitness function
#' @param x A \code{data.frame} containing all features that should be tried in the feature selection process.
#' @param y The response variable.
#' @param seed Defines the seed used for the internal loop to make sure to get the same fit with the same features used.
#' @param ga.input A named \code{list} list(fitness.func, metric.func) containing the fitness function & the metric used to evaluate the fit (e.g. RMSE).
#' @param crossover.rate Fraction of crossover between parent solutions.
#' @param mutation.rate Probability of a mutation (bitflip) within the feature vector.
#' @param pop.size Parent population size for the genetic algorithm.
#' @param offspring.size Offspring population size for the genetic algorithm.
#' @param max.iter Maximum number of iterations the algorithm should run.
#' @param stop.criterion Maximum number of iterations without improvement of the solution.
#' @param remove.overlap A logical value. If overlapping solutions should be removed (see Nojima et al. 2005)?
#' @param initialize.equal A logical value. If initial population should be drawn at random (\code{FALSE}) or if it should contain individuals drawn uniformly along the number of features dimension (\code{TRUE}).
#' @param init.limit A number. Maximum number of genes to activate in the initial population. This is rather useful if the model fitting is slow for models with lots of features (so you want to limit it for computational purposes).
#' @param allowZeroGenes A logical value. If individuals containing only zeros (no features) should be included in the search (\code{TRUE}) or not (\code{FALSE}).
#' @param ref.point Vector of length 2, defining the Reference point for the Hypervolume Indicator.
#' @param n.cores Number of CPU-cores to parallelize the external loop (parallelize number of individuals in the population). If \code{NULL} the function will try all possible number of cores and chooses the fastest (the evaluation will however take some time, based on the machine).
#' @param memoisation A logical value. If memoisation should be used (\code{TRUE}) or not (\code{FALSE}). If memoisation is used, function calls with the same input parameters (aka same features) are stored in a Cache and returned if solution is visited again.
#'
#' @return A \code{list} with the following elements:
#' \itemize{
#'   \item generation.fitness - \code{data.frame} of tracked fitness per iteration
#'   \item generation.populations - A list of all individuals within the population per iteration
#'   \item generation.offspring - A list of all offspring individuals per iteration
#'   \item pareto.individuals - individuals of the pareto-optimal solution
#'   \item pareto.varnames - Names of features included in the pareto-optimal solution
#'   \item pareto.solution - Summary table of the final pareto-solution
#'   \item algorithm.params - Input parameters and other summarized in a list to make output more reproducible
#'   \item fitted.dataset - \code{list} of \code{x} and \code{y} used for model fitting
#' }
#' @details This function implements the non-dominated sorting genetic algorithm (NSGA-II) feature selection for the bi-objective case. It simultanously minimizes the Nr. of features and
#'          a measure of performance (e.g. RMSE) for a machine learning model (e.g. Random Forest).
#' @importFrom ecr evaluateFitness

#' @export
NSGA.II_wrapper = function(x,y,seed,
                           ga.input,
                           memoisation = F,
                           mutation.rate = 0.1,
                           crossover.rate = 0.8,
                           remove.overlap = T,
                           pop.size = 50,
                           offspring.size = NULL,
                           max.iter = 150,
                           initialize.equal = T,
                           init.limit = NULL,
                           stop.criterion = NULL,
                           n.cores = NULL,
                           allowZeroGenes = F,
                           ref.point = c(1,1)){

  if(exists(".Random.seed")){
    old.seed = .Random.seed
  } else{
    set.seed(NULL)
    old.seed = .Random.seed
  }
  on.exit(assign(".Random.seed",value = old.seed,envir = .GlobalEnv))

  # .................................. #
  # record the time taken for analysis
  start.time = Sys.time()
  # .................................. #

  # .............. #
  # parameters:
  n.feat = ncol(x)
  # If offspring.size == NULL, set it to population size
  if(is.null(offspring.size)){
    offspring.size = pop.size
  }
  # get stuff from the ga.input
  metric.func = ga.input$metric.func

  if(memoisation == T){
    # add memoisation (cache) that stores already visited solutions, thereby increasing the computational speed
    fitness.func = R.cache::addMemoization(f = ga.input$fitness.func)
    on.exit(R.cache::clearCache(prompt = F),add = T)
  } else{
    fitness.func = ga.input$fitness.func
  }

  # .............. #

  # ......................... #
  # define the control object

  ctrl = ecr::initECRControl(fitness.fun = function(feat,seed) fitness.func(feat,x = x,seed = seed,y = y,metric.func = metric.func),
                             n.objectives = 2L, minimize = c(TRUE,TRUE))

  # define mutation function
  ctrl = ecr::registerECROperator(ctrl, "mutate",ecr::mutBitflip)
  # define crossover function
  ctrl = ecr::registerECROperator(ctrl,"recombine",ecr::recUnifCrossover)
  # mating rules
  ctrl = ecr::registerECROperator(ctrl , "selectForMating",ecr::selSimple) # tournament selection or simple selection?
  # survival (unique to NSGA-II algorithm: non dominated)
  ctrl = ecr::registerECROperator(ctrl ,"selectForSurvival",ecr::selNondom)

  #define the population (binomial)
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  if(initialize.equal == T){
    if(!is.null(init.limit)){
      nr.params = sample(1:init.limit,pop.size,replace = T)
    } else{
        nr.params = sample(1:n.feat,pop.size,replace = T)
      }
    population = rep(list(rep(0,n.feat)),pop.size)
    for(i in 1:length(nr.params)){
      indx = sample(1:n.feat,nr.params[i],replace = F)
      population[[i]][indx] = 1
    }
  } else{
    population = ecr::genBin(pop.size, n.feat)
    if(!is.null(init.limit)){
      feat.limit = function(n.dim,init.limit){
        indx.feat = sample(1:n.dim, size = sample(1:init.limit,1), replace = F)
        feat.vec = rep(0,n.dim)
        feat.vec[indx.feat] = 1
        return(feat.vec)}

      population = ecr::gen(expr = feat.limit(n.dim = n.feat, init.limit = init.limit),n = pop.size)
    }

  }

  # check if 0 gene-solution exists
  if(allowZeroGenes == F){
    if(any(sapply(population,sum)== 0)){
      indx = which(sapply(population,sum)== 0)
      # create random individual
      population[[indx]] = sample(c(0,1),n.feat,replace = T)
    }
  }

  #evaluate fitness for the first population

  # (parallelize fitness evaluation on n cores)
  on.exit(parallelMap::parallelStop(),add = T)
  if(is.null(n.cores)){
    cat("Checking number of cores for optimal parallelisation...","\n")

    # evaluate best nr. of cores yourself... (might take some time)
    system.cores = parallel::detectCores()
    # if number of cores is large (e.g. on a server) > 30 it will throw an error
    if(system.cores > 30){
      stop(paste("Automatically checking number of cores not advised if your machine has more than 30!\n", "Please specify 'n.cores'!\n"))
    }
    cat(paste("Number of cores detected:",system.cores),"\n")

    time.taken = rep(NA,system.cores-1)
    cat("Checking...","\n")
    for(ii in 2:system.cores){
      start = Sys.time()
      parallelMap::parallelStartSocket(cpus = ii, level = "ecr.evaluateFitness",show.info = F)
      RNGkind("L'Ecuyer-CMRG")
      parallelMap::parallelExport("x","y","fitness.func","metric.func","seed",
                     level = "ecr.evaluateFitness",show.info = F)
      fitness = ecr::evaluateFitness(inds = population,control =  ctrl,seed = seed)
      parallelMap::parallelStop()
      time.taken[ii-1] = difftime(time1 = Sys.time(),time2 = start,units = "secs")
      cat(paste("# cores:",ii,"|| time taken:",round(time.taken[ii-1],digits = 2)),"\n")
    }
    n.cores = (2:system.cores)[which(time.taken == min(time.taken))]
    cat(paste("Select number of cores to be:",n.cores),"\n")
  } else{
    parallelMap::parallelStartSocket(cpus = n.cores, level = "ecr.evaluateFitness")
    RNGkind("L'Ecuyer-CMRG")
    parallelMap::parallelExport("x","y","fitness.func","metric.func","seed",
                   level = "ecr.evaluateFitness")
    fitness = ecr::evaluateFitness(inds = population,control =  ctrl,seed = seed)
    parallelMap::parallelStop()
  }
  # ........................... #

  # ............................................... #
  # create a logger (to keep track of the evolution)
  # logs the population, time and the fitness

  # logger for multi-objective optimisation using the hypervolume indicator to calculate fitness of a solution set

  ref.point = ref.point # maximum of each parameter minimised

  logger = ecr::initLogger(control = ctrl,
                      log.stats = list(fitness = list("HV" = list(
                        fun = ecr::computeHV,
                        pars = list(ref.point = ref.point)))),
                      log.pop = T, # to store the populations in each run
                      init.size = max.iter +1L)

  ecr::updateLogger(log = logger,
               population = population,
               fitness = fitness,
               n.evals = offspring.size)

  # create list to store the offspring
  store.offspring = rep(list(NA),max.iter)
  # ............................................... #

  # ---------------------- #
  # simulate the evolution
  # ---------------------- #

  for (i in seq_len(max.iter)) {
    # create offspring
    offspring = ecr::recombinate(control = ctrl ,
                                 inds = population,
                                 fitness = fitness ,
                                 lambda = offspring.size,
                                 p.recomb = crossover.rate)

    # mutate offspring
    offspring = ecr::mutate(control = ctrl ,
                            inds = offspring ,
                            p.mut = mutation.rate)

    # control for 0-genes individuals
    if(allowZeroGenes == F){
      if(any(sapply(offspring,sum)== 0)){
        indx = which(sapply(offspring,sum)== 0)
        # create random individual
        for(zeroGenes in indx){
          # take an offspring and activate one gene
          select.offspring = sample(1:length(offspring),1,replace  = T)
          replace.offspring = offspring[[select.offspring]]
          rnd.gene = sample(1:n.feat,1,replace  = T)
          replace.offspring[rnd.gene] = 1
          offspring[[zeroGenes]] =  replace.offspring
        }
      }
    }

    # remove overlapping individuals
    if(remove.overlap == T){
      pop.tmp = do.call(rbind,population)
      offspr.tmp = do.call(rbind,offspring)
      # unique offspring
      offspr.tmp = unique(offspr.tmp)
      # check if there are offspring solutions that are already in the population
      dubs = rep(NA,nrow(offspr.tmp))
      for(mm in 1: nrow(offspr.tmp)){
        dubs[mm] = any(apply(pop.tmp,1,function(x) identical(x,offspr.tmp[mm,])))
      }
      # remove dublicates
      index.rm = which(dubs == T)
      if(!pracma::isempty(index.rm)){
        offspr.tmp = offspr.tmp[-index.rm,]
      }
      # create list again
      offspring = split(offspr.tmp, seq(nrow(offspr.tmp)))
    }

    # print to console
    cat(paste("NSGA-II | iter =", i, "| pop.size =",
              length(population), "| offspr.size =", length(offspring)))
    cat("\n")
    utils::flush.console()

    # parallelize fitness evaluation on n cores
    parallelMap::parallelStartSocket(cpus = n.cores, level = "ecr.evaluateFitness")
    RNGkind("L'Ecuyer-CMRG")
    parallelMap::parallelExport("x","y","fitness.func","metric.func","seed",
                   level = "ecr.evaluateFitness") # need to export own functions, so that parallel tasks can find them

    # evaluate fitness of offspring
    fitness.o = ecr::evaluateFitness(inds = offspring ,control =  ctrl,seed = seed)
    parallelMap::parallelStop()

    # Survival (non-dominated): select a number of mu individuals from the parent and offspring generation that matches the initial population size
    sel = ecr::replaceMuPlusLambda(control = ctrl ,population =  population ,
                              offspring = offspring ,
                              fitness = fitness ,
                              fitness.offspring = fitness.o)
    # select the new population
    population = sel$population
    # select the fitness of the new population
    fitness = sel$fitness

    ecr::updateLogger(logger,population,
                 fitness = fitness, n.evals = length(offspring))

    # store offspring
    store.offspring[[i]] = list(offspring = offspring,fitness.o = fitness.o)


    #check if stopping criterion is reached
    if(!is.null(stop.criterion)){
      if(i == 1){
        best.fitness = ecr::getStatistics(logger)$fitness.HV[1]
        n = 1
      }else{
        new.fitness = ecr::getStatistics(logger)$fitness.HV[i]
        if(best.fitness == new.fitness){
          n = n+1
        } else {
          n = 1
        }
        best.fitness = new.fitness
        # stoping criterion
        if(n >= stop.criterion){
          break()
        }
      }
    }
  }

  # ----------- #
  # wrap output
  # ----------- #

  # get the statistics from the logger
  stats.fitness = ecr::getStatistics(logger)
  stats.populations = ecr::getPopulations(logger)

  # get the pareto front
  Pareto.set = ecr::initParetoArchive(control = ctrl)
  ecr::updateParetoArchive(Pareto.set,inds = population,fitness = fitness)

  pareto.front = t(ecr::getFront(Pareto.set))
  fit.individuals = do.call(rbind,ecr::getIndividuals(Pareto.set))

  # get unique solutions
  indx.uniq = which(duplicated(pareto.front) == F)
  pareto.front.uniq = data.frame(ID = indx.uniq,pareto.front[indx.uniq,])

  # get individuals
  pareto.individuals = fit.individuals[indx.uniq,]

  # nr. of parameters
  if(is.null(dim(pareto.individuals))){ # control if pareto-set contains only one solution
    pareto.front.uniq$nr.of.params = sum(pareto.individuals)
    pareto.varnames = list(names(x)[which(pareto.individuals == 1)])
  } else{
    pareto.front.uniq$nr.of.params = apply(pareto.individuals,1,sum)
    # names of variables in pareto.front
    pareto.varnames = list()
    for(j in 1:nrow(pareto.individuals)){
      pareto.varnames[[j]] = names(x)[which(pareto.individuals[j,] == 1)]
    }
  }

  # calculate euclidean distance to c(0,0)
  pareto.front.uniq$dist.to.CoordOrigin = apply(pareto.front.uniq[,c(2,3)],1,function(x) sqrt(sum((x - c(0,0)) ^ 2)))

  # add specifications of the model fitting and stuff
  algorithm.params = list(seed = seed,
                          fitness.func = ga.input$fitness.func,
                          CV = c(folds = formals(ga.input$fitness.func)$folds,reps = formals(ga.input$fitness.func)$reps),
                          metric = ga.input$metric.func,
                          ga.parameters = data.frame(mutation.rate = mutation.rate,
                                                     crossover.rate = crossover.rate,
                                                     pop.size = pop.size,
                                                     offspring.size = offspring.size,
                                                     only.unique.offspring = remove.overlap,
                                                     initialize.equal = initialize.equal,
                                                     max.iter = max.iter,
                                                     iter.realised = i,
                                                     ref.point.x = ref.point[1],ref.point.y = ref.point[2]),
                          n.cores = n.cores,
                          time.taken = difftime(time1 = Sys.time(),time2 = start.time,units = "hours"))



  # return list with individuals, pareto-optimal solutions,
  # fitness, the tracked populations and names of the variables in pareto solution
  return(
    list(generation.fitness = stats.fitness, # tracked fitness levels
         generation.populations = stats.populations, # tracked populations
         generation.offspring = store.offspring, # tracked offspring (to get an idea of the search space visited)
         pareto.individuals = pareto.individuals, # individuals of the pareto-optimal solution
         pareto.varnames = pareto.varnames, # names of the variables selected in the pareto-solutions
         pareto.solution = pareto.front.uniq, # pareto-front
         algorithm.params = algorithm.params, # input parameters and others summarized in a list to make output more reproducible
         fitted.dataset = list(x = x,y = y)
    )
  )
}

# ----------------- #
# plotting function
# ----------------- #

#' plot function for output of \code{NSGA.II_wrapper}
#'
#' @param output.NSGAII An obj. (list) containing results of the NSGA-II feature selection.
#' @description A summary plot for the NSGA-II feature selection.
#' @return A plot showing the evolution with the Hyper-volume Indicator as well as the Pareto-Front.
#' @export
plot_NSGAII = function(output.NSGAII){
  opar = graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(opar))
  graphics::par(mfrow = c(1,2),mar = c(4,4,3,2))
  # first plot
  graphics::plot(output.NSGAII$generation.fitness$gen,output.NSGAII$generation.fitness$fitness.HV,
       las = 1,xlab = "generations",ylab = "fitness (Hypervolume Ind.)",lwd = 0.7,
       col = "darkgoldenrod4",type = "l")
  graphics::points(output.NSGAII$generation.fitness$gen,output.NSGAII$generation.fitness$fitness.HV,cex = 1,pch = 16,
         col = "darkgoldenrod3")
  graphics::title("NSGA-optimisation",adj = 0)
  # second plot
  tmp = output.NSGAII$pareto.solution[order(output.NSGAII$pareto.solution$params,decreasing = T),] # order solutions

  graphics::plot(tmp$metric,tmp$params,col = "deeppink2",lwd = 0.7,las = 1,
       xlab = "metric",type = "l",
       ylab = "fraction parameters",xlim = c(0,max(output.NSGAII$pareto.solution$metric,output.NSGAII$pareto.solution$params)+0.1),
       ylim = c(0,max(output.NSGAII$pareto.solution$metric,output.NSGAII$pareto.solution$params)+0.1))
  graphics::points(output.NSGAII$pareto.solution$metric,output.NSGAII$pareto.solution$params,
         pch = 16,col = "darkmagenta")
  graphics::title("Pareto-Front",adj = 0)
}
