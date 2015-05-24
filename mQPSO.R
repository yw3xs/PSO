# My implementation of multi-swarm quantum PSO 
# Reference: Multi-swarm Optimization in Dynamic Environments -- Tim Blackwell and Jurgen Branke2

################################################################################
##                          InitializeX                                       ##
################################################################################
# Purpose  : To initialize X for the main algorithm, generate a matrix with dim(particle, nparam)
#            given the lower and upper bound of the search space.

# 'particle'   : number of particles in the swarm
# 'nparam'     : number of parameters to be estimated in the model
# 'lower', 'upper' : lower and upper bound of the parameter in each dimension ('nparam'-D vector)

InitializeX <- function(particle, lower, upper){
  nparam = length(lower)
  lowerM <- matrix(rep(lower, particle), nrow=particle, byrow=TRUE)
  upperM <- matrix(rep(upper, particle), nrow=particle, byrow=TRUE)
  X <- matrix(runif(nparam*particle, min=0, max=1), nrow=particle, ncol=nparam)
  X <- lowerM + (upperM-lowerM)*X 
} # 'InitializeX' end

################################################################################
##                          clamp                                             ##
################################################################################
# Purpose  : To constrain the position of each particle within upper/lower bound

# 'X'  : nparam-D velocity/position vector/matrix
# 'upper/lower' : upper and lower bound of searching space

clamp <- function(X, upper, lower){
  if(is.matrix(X)){
    n <- ncol(X)        # number of parameters (columns in X matrix)
    for (i in 1:n){
      pos <- which(X[,i]>upper[i])
      X[pos,i] <- upper[i]
      pos <- which(X[,i]<lower[i])
      X[pos,i] <- lower[i]
    }  
    return(X)
  }
  else if(is.vector(X)){
    pos <- which(X>upper)
    X[pos] <- upper[pos]
    pos <- which(X<lower)
    X[pos] <- lower[pos] 
    return(X)
  }
} # 'clamp' end

################################################################################
##                          clamp.sphere                                      ##
################################################################################
# Purpose  : To constrain the position of each particle within a sphere centered around a point

# 'X'  : matrix of points
# 'radius' : radius of the sphere
# 'center' : the center of the sphere

clamp.sphere <- function(X, radius, center){
    n <- nrow(X)        # number of parameters (columns in X matrix)
    for (i in 1:n){
      if(distance(X[i,],center) > radius){
        X[i,] = (X[i,] - center) * radius / distance(X[i,],center) + center
      }
    }  
    return(X)
} # 'clamp.sphere' end 



################################################################################
##                          Initialize quantum particles(box)                 ##
################################################################################
# purpose : initialize/update quantum particles in a box with length rcloud*(upper-lower) centered around swarm attractor

# center : position of swarm attractor,
# qpar : number of quantum particles
# rcloud : percentage of span of searching space
# upper/lower : range of searching space

quantum.box <- function(center, qpar, rcloud, upper, lower){
  nparam = length(center)
  X <- matrix(rep(NA, qpar*nparam), nrow = qpar, ncol = nparam)
  for (i in 1:nparam){
    X[,i]<-runif(qpar, min = -(upper[i]-lower[i])*rcloud, max=(upper[i]-lower[i])*rcloud) + center[i]        
  }
  return (X)  
} # 'quantum.box' end

################################################################################
##                          Initialize quantum particles(sphere)              ##
################################################################################
# purpose : initialize/update quantum particles in a sphere with radius rcloud*(upper-lower) centered around swarm attractor.

# center : position of swarm attractor,
# qpar : number of quantum particles
# rcloud : percentage of span of searching space
# upper/lower : range of searching space
quantum.sphere <- function(center, qpar, rcloud, upper, lower) {
  
  # dimension of 'center' (number of parameters)
  nparam = length(center)
  radius = rcloud * max(upper - lower)
  
  # Step 1. Direction 
  x <- matrix(rnorm(nparam * qpar, mean=0, sd=1), nrow = qpar)
  l <- sqrt( rowSums(x*x) )
  
  # Step 2. Random Radius
  r <- runif(qpar)
  
  x <- r * radius * x / l
  
  # Centering the random point at 'G'
  centermatrix <- matrix(rep(center, qpar), nrow = qpar, byrow = TRUE)
  return( x +  centermatrix)  
} # 'quantum.sphere' end


################################################################################
##                          Distance calculation                              ##
################################################################################
# purpose : calculate the euclidean distance of two points in the searching space

# point1, point2
distance <- function(point1, point2){
  return (sqrt(sum((point1 - point2)^2)))
} # 'distance' end


################################################################################
##                       Swarm diameter calculation                           ##
################################################################################
# purpose : calculate the diameter of neutral particles in a swarm

# point1, point2
diameter <- function(swarm){
  size <- nrow(swarm)
  DISMAT <- matrix(NA, nrow = size, ncol = size)
  for(i in 1:size){
    for (j in 1:size){
      DISMAT[i,j] = distance(swarm[i,], swarm[j,])
    }
  }
  return (max(DISMAT))
} # 'diameter' end


################################################################################
##                          Steepest descent                                  ##
################################################################################
# Purpose  : Given a function and a starting point, calculate the steepest descent direction and randomly choose a step size to move to the next point. repeat for LMAX iterations

SD <- function(fun,                            
               startP,
               LMAX, 
               lower = -Inf,
               upper = Inf,              
               stepsize = NULL,              
               grad_tol = 1e-7
){
  result <- list()       # for returning results (best position and # of iterations)
  iter <- 0              # communication iterations (gradient and line search)
  gradconv <- FALSE
  P <- startP
  while(LMAX > 0){        
    LMAX = LMAX -1    
    require(numDeriv)    
    grad <- grad(fun,P)  # one communication iteration     
    iter <- iter + 1
    if (sqrt(sum(grad^2)) < grad_tol) {
      message("gradient achieves 0!")
      result[[1]] <- P
      result[[2]] <- fun(P)
      result[[3]] <- iter
      gradconv <- TRUE
      result[[4]] <- gradconv
      return (result)
    }
    ref <- fun(P)
    require(ppls)     # for the function normalize.vector
    direction <- normalize.vector(grad)    
    stepsize <- seq(from = 5, to = 1, by=-1)
    stepsize <- c(stepsize, 2^seq(from = 0, to = -40, by = -1))
    newP <- list()
    newFit <- NULL
    for (i in 1:length(stepsize)){
      newP[[i]] <- P - stepsize[i] * direction
      newP[[i]] <- clamp(newP[[i]], upper, lower)
      newFit[i] <- fun(newP[[i]])
    }
    iter <- iter + 1
    bestfit <- newFit[which.min(newFit)]
    bestP <- newP[[which.min(newFit)]]
    #message("best stepsize : ", stepsize[which.min(newFit)])
    P <- bestP
  } 
  result[[1]] <- P
  result[[2]] <- bestfit
  result[[3]] <- iter
  result[[4]] <- gradconv
  return (result)
} # 'SD' end


################################################################################
##                          Method: mQPSO                                     ##
################################################################################
mQPSO <- function(fun,                            
                  nparam,                                      
                  lower = -Inf,
                  upper = Inf,
                  minmax = c("min","max"),
                  control = list()
                 ){
  ############################################################################  
  # 0)                             Parameter Setup                           #
  ###################################################$$$$##################### 
  minmax <- match.arg(minmax)
  con <- list(
    swarm = 5,
    npar=5, 
    qpar=10,
    maxit=100, 
    sdmaxit = 10,
    maxfn=Inf,
    c1_cf= 2.05,
    c2_cf= 2.05,
    cf= 0.729843788,
    rexecl= 1.5,
    rcloud= .15,   
    reltol = 1e-8,
    parallel = TRUE,    
    verbose = TRUE,
    dmax = 2
  )
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  if (length(noNms <- namc[!namc %in% nmsC])) 
    warning("[Unknown names in control: ", paste(noNms, collapse = ", "), " (not used) !]")         
  swarm <- con[["swarm"]] 
  npar <- con[["npar"]]
  qpar <- con[["qpar"]]  
  maxit <- con[["maxit"]] 
  sdmaxit <- con[["sdmaxit"]] 
  maxfn <- con[["maxfn"]]
  c1_cf <- con[["c1_cf"]] 
  c2_cf <- con[["c2_cf"]]   
  cf <- con[["cf"]]
  rexecl <- con[["rexecl"]]
  rcloud <- con[["rcloud"]]
  reltol <- con[["reltol"]] 
  parallel <- as.logical(con[["parallel"]])
  verbose <- con[["verbose"]]
  dmax <- con[["dmax"]]
  ############################################################################  
  # 1)            Initialisation of X and V in each swarm                    #
  ############################################################################
  if (parallel){    
    require(parallel)
    if (verbose) message("                               ")
    if (verbose) message("[ Parallel initialization ... ]")
    nnodes.pc <- parallel::detectCores()
    if (verbose) message("[ Number of cores/nodes detected: ", nnodes.pc, " ]")
    cl <- makeForkCluster(nnodes = nnodes.pc)        
  }
  allswarm <- list()
  allv <- list()
  for (i in 1:swarm){    
    swarmname <- paste("swarm",i,sep="")
    xtmp <- InitializeX(npar, lower, upper) 
    allswarm[[i]] <- xtmp
    names(allswarm)[i] <- swarmname    
    vname <- paste("V",i,sep="")
    vtmp <- (InitializeX(npar, lower, upper)-xtmp)/2 
    allv[[i]] <- vtmp
    names(allv)[i] <- vname    
  }
  
  ############################################################################  
  # 2)                              Main algorithm                           #
  ############################################################################
  
  if(minmax == "max") fun = -fun  
  f.x <- list()
  Pbest <- list()
  Pfit <- list()
  Gfit <- numeric()
  Gbest <- list()
  relconv <- rep(FALSE, swarm)
  iter <- rep(0, swarm)
  
  for (i in 1:swarm){
    ## first evaluation
    if(parallel){
      f.x[[i]] <- parallel::parRapply(cl= cl, x=allswarm[[i]], FUN=fun)  
    }
    else { f.x[[i]] <- apply(allswarm[[i]], MARGIN = 1, fun) }
    Pbest[[i]] <- allswarm[[i]][1:npar,]   # personal best positions of neutral particles
    Pfit[[i]] <- f.x[[i]][1:npar]         # personal best function fit of neutral particles
    Gfit[i] <- f.x[[i]][which.min(f.x[[i]])]      # global best function fit
    Gbest[[i]] <- allswarm[[i]][which.min(f.x[[i]]),]  # global best position
    swarmD <- 9999
    while (swarmD > dmax && iter[i] < maxit){  
    #while (!relconv[i] && iter[i] < maxit){
      #Pfit.prior <- Pfit[[i]]
      iter[i] <- iter[i] + 1    
      ### update neutral particles
      c1M <- matrix(runif(nparam*npar, 0, c1_cf), nrow = npar, ncol = nparam)
      c2M <- matrix(runif(nparam*npar, 0, c2_cf), nrow = npar, ncol = nparam)
      allv[[i]] <- cf*(allv[[i]] + c1M * (Pbest[[i]] - allswarm[[i]][1:npar,]) + c2M * (matrix(rep(Gbest[[i]],npar), nrow = npar, byrow = TRUE) - allswarm[[i]][1:npar,]))     
      allswarm[[i]][1:npar,] <- clamp(allswarm[[i]][1:npar,] + allv[[i]], upper,lower)
           
      ### evaluate the function fit at new positions
      if(parallel){
        f.x[[i]] <- parallel::parRapply(cl= cl, x=allswarm[[i]], FUN=fun)  
      }
      else { f.x[[i]] <- apply(allswarm[[i]], MARGIN = 1, fun) }
      
      tmp <- ( which( f.x[[i]][1:npar] < Pfit[[i]] ))
      if(length(tmp) > 0){
        Pfit[[i]][tmp] <- f.x[[i]][tmp]
        Pbest[[i]][tmp,] <- allswarm[[i]][tmp,]
      }
      tmp <- which.min(f.x[[i]])
      if(f.x[[i]][tmp] < Gfit[i]){
        Gbest[[i]] <- allswarm[[i]][tmp,]
        Gfit[i] <- f.x[[i]][tmp]      
      }
      
      swarmD <- diameter(allswarm[[i]])  
      #print(allswarm[[i]])
      #print(swarmD)
#       relerror <- max(abs(Pfit[[i]] - Pfit.prior)/Pfit.prior)
#       if(verbose) cat("Relative error at iteration ", iter, ", swarm", i, " : ", relerror, "\n")
#       print(relerror)
#       relconv[i] <- (relerror < reltol)
    } # end while    
    cat("\n##################################\n")
    cat("Swarm attractor: ", Gbest[[i]], " best fit : ", Gfit[[i]], "diameter: ", swarmD)  
  }  # end for

  Center <- list()
  Cfit <- rep(9999,swarm)
  sol <- numeric()   # best solution
  fit <- 9999   # best fit
  stagit = 0 # number of iteration when best solution doesn't change
  while(max(iter) < maxit && stagit < 15){
    ### SD iterations  
    SD.result <- list()    
    gradconv <- rep(FALSE, swarm)
    
    for (i in 1:swarm){
      SD.result[[i]] <- SD(fun, Gbest[[i]], sdmaxit, lower,upper) 
      #Gbest[[i]] <- SD.result[[i]][[1]]
      #Gfit[i] <- SD.result[[i]][[2]]
      if(SD.result[[i]][[2]] < Cfit[i]){
        Center[[i]] <- SD.result[[i]][[1]]
        Cfit[i] <- SD.result[[i]][[2]]
        iter[i] <- iter[i] + SD.result[[i]][[3]]
        gradconv[i] <- SD.result[[i]][[4]]                  
      }
    }
    # find the best center
    if(fit != Cfit[which.min(Cfit)]){
      sol <- Center[which.min(Cfit)]
      fit <- Cfit[which.min(Cfit)]
    }
    else {stagit = stagit+1}
    
    
    ### find the worst swarms and reinitialize the particles
    resetworst <- which.max(Cfit)
    xtmp <- InitializeX(npar, lower, upper) 
    allswarm[[resetworst]] <- xtmp
    vtmp <- (InitializeX(npar, lower, upper)-xtmp)/2 
    allv[[resetworst]] <- vtmp        
    Cfit[resetworst] <- 9999
    ## check if two swarms are optimizing the same peak
    resetcompeting <- numeric()
    for(j in (1:swarm)[-resetworst]){
      for(k in (1:swarm)[-resetworst]){
        if (j < k){          
          if (distance(Center[[j]], Center[[k]]) < rexecl){
            if (Cfit[[j]] < Cfit[[k]]){
              # record the index of the swarm to be reset
              resetcompeting <- c(resetcompeting, k)
            }
            else resetcompeting <- c(resetcompeting, j)          
          }                
        }            
      }
    }
    resetcompeting <- unique(resetcompeting)
    for (i in resetcompeting) {    
      xtmp <- InitializeX(npar, lower, upper)
      allswarm[[i]] <- xtmp
      vtmp <- (InitializeX(npar, lower, upper)-xtmp[1:npar,])/2 
      allv[[i]] <- vtmp
      Cfit[i] <- 9999
    }
    resetindex <- c(resetworst, resetcompeting)
    
    ### construct a sphere for quantum particles around every gbest
    for (i in (1:swarm)[-resetindex]){
      #allswarm[[i]] <- rbind(allswarm[[i]], clamp(quantum.sphere(Gbest[[i]],qpar,rcloud,upper,lower),upper,lower))      
      allswarm[[i]] <- clamp(quantum.sphere(Center[[i]],npar+qpar,rcloud,upper,lower),upper,lower)
      allv[[i]] <- (InitializeX(npar, lower, upper)-InitializeX(npar, lower, upper))/2
    }
    ### run PSO iterations on all swarms
    for (i in 1:swarm){
      if(i %in% resetindex){
        ## first evaluation
        if(parallel){
          f.x[[i]] <- parallel::parRapply(cl= cl, x=allswarm[[i]], FUN=fun)  
        }
        else { f.x[[i]] <- apply(allswarm[[i]], MARGIN = 1, fun) }
        Pbest[[i]] <- allswarm[[i]][1:npar,]   # personal best positions of neutral particles
        Pfit[[i]] <- f.x[[i]][1:npar]         # personal best function fit of neutral particles
        Gfit[i] <- f.x[[i]][which.min(f.x[[i]])]      # global best function fit
        Gbest[[i]] <- allswarm[[i]][which.min(f.x[[i]]),]  # global best position
        relconv[i] = FALSE
        #while (!relconv[i] && iter[i] < maxit){  
        swarmD <- 99999
        while (swarmD > dmax && iter[i] < maxit){  
        #  Pfit.prior <- Pfit[[i]]
          iter[i] <- iter[i] + 1    
          ### update neutral particles
          c1M <- matrix(runif(nparam*npar, 0, c1_cf), nrow = npar, ncol = nparam)
          c2M <- matrix(runif(nparam*npar, 0, c2_cf), nrow = npar, ncol = nparam)
          allv[[i]] <- cf*(allv[[i]] + c1M * (Pbest[[i]] - allswarm[[i]][1:npar,]) + c2M * (matrix(rep(Gbest[[i]],npar), nrow = npar, byrow = TRUE) - allswarm[[i]][1:npar,]))      
          allswarm[[i]][1:npar,] <- clamp(allswarm[[i]][1:npar,] + allv[[i]], upper, lower)    
          
          ### evaluate the function fit at new positions
          if(parallel){
            f.x[[i]] <- parallel::parRapply(cl= cl, x=allswarm[[i]], FUN=fun)  
          }
          else { f.x[[i]] <- apply(allswarm[[i]], MARGIN = 1, fun) }
          
          tmp <- ( which( f.x[[i]][1:npar] < Pfit[[i]] ))
          if(length(tmp) > 0){
            Pfit[[i]][tmp] <- f.x[[i]][tmp]
            Pbest[[i]][tmp,] <- allswarm[[i]][tmp,]
          }
          tmp <- which.min(f.x[[i]])
          if(f.x[[i]][tmp] < Gfit[[i]]){
            Gbest[[i]] <- allswarm[[i]][tmp,]
            Gfit[i] <- f.x[[i]][tmp]      
          }
          swarmD <- diameter(allswarm[[i]])      
#            relerror <- max(abs(Pfit[[i]] - Pfit.prior)/Pfit.prior)
#            if(verbose) cat("Relative error at iteration ", iter, ", swarm", i, " : ", relerror, "\n")
#            relconv[i] <- (relerror < reltol)
        } # end while      
      } # end if   
      else{
        ## first evaluation
        if(parallel){
          f.x[[i]] <- parallel::parRapply(cl= cl, x=allswarm[[i]], FUN=fun)  
        }
        else { f.x[[i]] <- apply(allswarm[[i]], MARGIN = 1, fun) }
        
        Pbest[[i]] <- allswarm[[i]][1:npar,]   # personal best positions of neutral particles
        Pfit[[i]] <- f.x[[i]][1:npar]         # personal best function fit of neutral particles                
        Gfit[i] <- f.x[[i]][which.min(f.x[[i]])]      # global best function fit
        Gbest[[i]] <- allswarm[[i]][which.min(f.x[[i]]),]  # global best position
        
#         if(f.x[[i]][which.min(f.x[[i]])] < Gfit[i]){
#           Gfit[i] <- f.x[[i]][which.min(f.x[[i]])]      # global best function fit
#           Gbest[[i]] <- allswarm[[i]][which.min(f.x[[i]]),]  # global best position
#         }
       # relconv[i] <- FALSE
        swarmD <- 99999
        while (swarmD > dmax && iter[i] < maxit){ 
      #  while (!relconv[i] && iter[i] < maxit){  
       #   Pfit.prior <- Pfit[[i]]
          iter[i] <- iter[i] + 1    
          ### update neutral particles
          c1M <- matrix(runif(nparam*npar, 0, c1_cf), nrow = npar, ncol = nparam)
          c2M <- matrix(runif(nparam*npar, 0, c2_cf), nrow = npar, ncol = nparam)
          allv[[i]] <- cf*(allv[[i]] + c1M * (Pbest[[i]] - allswarm[[i]][1:npar,]) + c2M * (matrix(rep(Gbest[[i]],npar), nrow = npar, byrow = TRUE) - allswarm[[i]][1:npar,]))      
          allswarm[[i]][1:npar,] <- clamp(clamp.sphere(allswarm[[i]][1:npar,] + allv[[i]],rcloud * max(upper - lower), Center[[i]]), upper, lower)
      
          allswarm[[i]] <- rbind(allswarm[[i]][1:npar,], clamp(quantum.sphere(Center[[i]],qpar,rcloud,upper,lower), upper, lower))                      
          ### evaluate the function fit at new positions
          if(parallel){
            f.x[[i]] <- parallel::parRapply(cl= cl, x=allswarm[[i]], FUN=fun)  
          }
          else { f.x[[i]] <- apply(allswarm[[i]], MARGIN = 1, fun) }
          
          tmp <- ( which( f.x[[i]][1:npar] < Pfit[[i]] ))
          if(length(tmp) > 0){
            Pfit[[i]][tmp] <- f.x[[i]][tmp]
            Pbest[[i]][tmp,] <- allswarm[[i]][tmp,]
          }
          tmp <- which.min(f.x[[i]])
          if(f.x[[i]][tmp] < Gfit[[i]]){
            Gbest[[i]] <- allswarm[[i]][tmp,]
            Gfit[i] <- f.x[[i]][tmp]      
          }      
          swarmD <- diameter(allswarm[[i]][1:npar,])
#            relerror <- max(abs(Pfit[[i]] - Pfit.prior)/Pfit.prior)          
#            relconv[i] <- (relerror < reltol)
        } # end while            
      } # end else      
    }  # end for
    
  }  # end while
cat("\n center \n")
print(sol)

cat("\n Cfit \n")
print(fit)
  print(iter)

}