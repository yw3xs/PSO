# My implementation of Gradient PSO (reference: A new gradient based particle swarm optimization algorithm for accurate computation of global minimum)
################################################################################
##                          InitializeX                                       ##
################################################################################
# Purpose  : To initialize X for the main algorithm, generate a matrix with dim(particle, nparam)
#            given the lower and upper bound of the search space.

# 'particle'   : number of particles in the swarm
# 'nparam'     : number of parameters to be estimated in the model
# 'lower', 'upper' : lower and upper bound of the parameter in each dimension ('nparam'-D vector)

InitializeX <- function(nparam, particle, lower, upper){
  lowerM <- matrix(rep(lower, particle), nrow=particle, byrow=TRUE)
  upperM <- matrix(rep(upper, particle), nrow=particle, byrow=TRUE)
  X <- matrix(runif(nparam*particle, min=0, max=1), nrow=particle, ncol=nparam)
  X <- lowerM + (upperM-lowerM)*X 
} # 'InitializeX' end

################################################################################
##                          clamp                                             ##
################################################################################
# Purpose  : To constrain the position of each particle within upper/lower bound

# 'X'  : nparam-D position vector/matrix
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

# ################################################################################
# ##                          Gradient calculation                              ##
# ################################################################################
# # Purpose  : To calculate the gradient of a function at a point X using finite difference method
# 
# # 'fun'   : objective function
# # 'X' : gradient at the point X
# # epsilon : finite difference
# gradient <- function(fun, X, epsilon){
#   grad <- NULL  
#   tmp <- X
#   for (i in 1:length(X)){
#     tmp[i] <- X[i] + epsilon
#     diff <- fun(tmp) - fun(X)
#     grad[i] <- diff / epsilon        
#     tmp <- X
#   }
#   return (grad)
# } # 'gradient' end

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
               grad_tol = 1e-5
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
      result[[2]] <- iter
      gradconv <- TRUE
      result[[3]] <- gradconv
      return (result)
    }
    ref <- fun(P)
    require(ppls)       # for the function normalize.vector
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
    P <- bestP
  } 
  result[[1]] <- P
  result[[2]] <- iter
  result[[3]] <- gradconv
  return (result)
} # 'SD' end


# ################################################################################
# ##                          Newton Raphson's method (NR)                      ##
# ################################################################################
# # Purpose  : Given a function and a starting point, calculate the steepest descent direction and randomly choose a step size to move to the next point. repeat for LMAX iterations
# 
# NR <- function(fun,                            
#                startP,
#                LMAX, 
#                lower = -Inf,
#                upper = Inf,              
#                stepsize = NULL,              
#                grad_tol = 1e-5
# ){
#   result <- list()       # for returning results (best position and # of iterations)
#  
#   sditer = 0
#   P <- startP
#   while(sditer < LMAX){
#     sditer = sditer + 1    
#     require(numDeriv)    
#     grad <- grad(fun,P)  # one communication iteration    
#     cat("\n grad: ", grad)
#     hessian <- hessian(fun,P)
#     cat("\n hessian: ", hessian)
#     direction <- solve(hessian)%*%grad
#     ref <- fun(P)
#     
#     stepsize <- seq(from = 5, to = 1, by=-1)
#     stepsize <- c(stepsize, 2^seq(from = 0, to = -40, by = -1))
#     newP <- list()
#     newFit <- NULL
#     for (i in 1:length(stepsize)){
#       newP[[i]] <- P - stepsize[i] * direction
#       newP[[i]] <- clamp(newP[[i]], upper, lower)
#       newFit[i] <- fun(newP[[i]])
#     }
# 
#     bestfit <- newFit[which.min(newFit)]
#     bestP <- newP[[which.min(newFit)]]
#     #message("best stepsize : ", stepsize[which.min(newFit)])
#     P <- bestP
#   } 
#   result[[1]] <- P
#   result[[2]] <- sditer
#   return (result)
# } # 'SD' end



################################################################################
##                          Method : Gradient PSO                             ##
################################################################################

GPSO <- function(fun,                            
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
    particle=40,
    maxit = 1000,
    psomaxit=20,       # number of pso iteration before going to steepest descent
    sdmaxit=5,         # number of steepest descent iteration before going to pso  
    maxfn=Inf,
    c1_w= 0.5+log(2), 
    c2_w= 0.5+log(2),  
    parallel = TRUE,
    verbose = TRUE
  )
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  if (length(noNms <- namc[!namc %in% nmsC])) 
    warning("[Unknown names in control: ", paste(noNms, collapse = ", "), " (not used) !]")         
  particle <- con[["particle"]] 
  maxit <- con[["maxit"]] 
  psomaxit <- con[["psomaxit"]]
  sdmaxit <- con[["sdmaxit"]]
  maxfn <- con[["maxfn"]]
  c1_w <- con[["c1_w"]] 
  c2_w <- con[["c2_w"]]
  parallel <- as.logical(con[["parallel"]])
  verbose <- con[["verbose"]]
  ############################################################################  
  # 1)                      Initialisation of X and V                        #
  ############################################################################ 
  if (parallel){    
    require(parallel)
    if (verbose) message("                               ")
    if (verbose) message("[ Parallel initialization ... ]")
    nnodes.pc <- parallel::detectCores()
    if (verbose) message("[ Number of cores/nodes detected: ", nnodes.pc, " ]")
    cl <- makeForkCluster(nnodes = nnodes.pc)        
  }
  
  if(length(upper) != nparam || length(lower) != nparam) stop("The dimension of lower bound or upper bound is not equal to number of parameters in the model!")
  X <- InitializeX(nparam,particle,lower,upper)
  V <- (InitializeX(nparam,particle,lower,upper)-X)/2
  ############################################################################  
  # 2)                              Main algorithm                           #
  ############################################################################
  
  ## first evaluation
  if(minmax == "max") fun = -fun  
  if(parallel){
    f.x <- parallel::parRapply(cl= cl, x=X, FUN=fun)  
  }
  else { f.x <- apply(X,MARGIN = 1,fun) }
  Pbest <- X                              # personal best positions
  Pfit <- f.x                             # personal best function fit  
  Gfit <- Pfit[which.min(Pfit)]           # global best function fit
  Gbest <- Pbest[which.min(Pfit),]        # global best position
  
  iter <- 0  
  Best <- NULL
  gradconv = FALSE 
  
  cf.upper <- .9
  cf.lower <- .7
  cf <- seq(from = cf.upper, to = cf.lower, by = -(cf.upper-cf.lower) / maxit)
  
  ## start the main loop
  
  #while (iter < maxit && !gradconv){
  while (iter < maxit - psomaxit){
    for (psoiter in 1:psomaxit){    
      iter <- iter+1 
      #Pfit.prior <- Pfit       
      ### vectorize the computation      
      c1M <- matrix(runif(nparam*particle, 0, c1_w), nrow = particle, ncol = nparam)
      c2M <- matrix(runif(nparam*particle, 0, c2_w), nrow = particle, ncol = nparam)
      
      V <- cf[iter] * (V + c1M * (Pbest - X) + c2M * (matrix(rep(Gbest,particle), nrow = particle, byrow = TRUE) - X))
      X <- X + V      
      X <- clamp(X, upper, lower) 
      
      if(parallel){
        f.x <- parallel::parRapply(cl= cl, x=X, FUN=fun)  
      }
      else { f.x <- apply(X,MARGIN = 1,fun) }
      
      tmp <- ( which( f.x < Pfit ))
      if(length(tmp) >0){
        Pfit[tmp] <- f.x[tmp]
        Pbest[tmp,] <- X[tmp,]
      }
      tmp <- which.min(f.x)
      if(f.x[tmp] < Gfit){
        Gbest <- X[tmp,]
        Gfit <- f.x[tmp]   
      }
      Best[iter] <- Gfit         
    }  # end of psomaxit loop
    
    SD.result <- SD(fun, Gbest, sdmaxit, lower,upper)
    Gbest.sd <- SD.result[[1]]
    iter <- iter + SD.result[[2]]
    gradconv <- SD.result[[3]]
    message("One set of SD finished")
    if (fun(Gbest.sd) < Gfit) {Gbest <- Gbest.sd;Gfit <-fun(Gbest.sd); message("improved!")}
    Best[iter] <- Gfit 
  }  # end of while loop
  parallel::stopCluster(cl) 
  
  param <- Gfit
  fit <- Gbest
  cat("\n##################################\n")
  cat("Optimal parameters: ", param, " best fit : ", fit)  
  cat("\n best value :", Best)
  test <- apply(X,MARGIN = 1,fun)
  cat("\n best value of the particles in the swarm: ",min(test))
}


