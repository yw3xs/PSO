# My implementation of dynamic PSO (part of the code is copied from the R package hydroPSO and PSO)

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
# Purpose  : To constrain the velocity/position of each particle within VMAX/LMAX (if larger than VMAX/LMAX, then set to VMAX/LMAX)

# 'X'  : nparam-D velocity/position vector/matrix
# 'XMAX' : nparam-D velocity/position vector/matrix, limit of velocity/position in each dimension

clamp <- function(X, XMAX){
  pos <- which( abs(X) > XMAX )
  if ( length(pos) > 0 ) 
    X[pos] <- sign(X[pos])*abs(XMAX[pos])    
  return(X)    
} # 'V.clamp' end

################################################################################
##                          Method 0: standard PSO                            ##
################################################################################

PSO <- function(fun,                            
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
    maxit=1000, 
    maxfn=Inf,
    c1_w= 0.5+log(2), 
    c2_w= 0.5+log(2),     
    use.IW= TRUE, IW.w = 1/(2*log(2)),    
    lambda= 1,       
    reltol=sqrt(.Machine$double.eps), 
    parallel = TRUE,
    verbose = TRUE
  )
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  if (length(noNms <- namc[!namc %in% nmsC])) 
    warning("[Unknown names in control: ", paste(noNms, collapse = ", "), " (not used) !]")         
  particle <- con[["particle"]] 
  maxit <- con[["maxit"]] 
  maxfn <- con[["maxfn"]]
  c1_w <- con[["c1_w"]] 
  c2_w <- con[["c2_w"]]   
  use.IW <- as.logical(con[["use.IW"]])
  IW.w <- con[["IW.w"]]
  lambda <- con[["lambda"]]
  reltol <- con[["reltol"]] 
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
  LMAX <- matrix(rep(upper - lower, particle), nrow = particle, byrow = TRUE) 
  VMAX <- lambda*LMAX   
  X <- InitializeX(nparam,particle,lower,upper)
  V <- (InitializeX(nparam,particle,lower,upper)-X)/2
  V <- clamp(V, VMAX)  
  
  ############################################################################  
  # 2)                              Main algorithm                           #
  ############################################################################
  
  ## first evaluation
  funeval <- 0                    # number of function evaluation
  if(minmax == "max") fun = -fun  
  if(parallel){
    f.x <- parallel::parRapply(cl= cl, x=X, FUN=fun)  
  }
  else { f.x <- apply(X,MARGIN = 1,fun) }
  funeval <- funeval + particle
  Pbest <- X                      # personal best positions
  Pfit <- f.x                     # personal best function fit
  G.index <- which.min(Pfit)      # global best index
  Gfit <- Pfit[G.index]           # global best function fit
  Gbest <- Pbest[G.index,]        # global best position
  
  iter <- 0  
  relconv = FALSE 
  
  ## start the main loop
  
  while (!relconv && iter < maxit){  
    Pfit.prior <- Pfit
    iter <- iter+1    
    ### vectorize the computation
    V <- IW.w*V
    c1M <- matrix(runif(nparam*particle, 0, c1_w), nrow = particle, ncol = nparam)
    c2M <- matrix(runif(nparam*particle, 0, c2_w), nrow = particle, ncol = nparam)
    V <- V + c1M * (Pbest - X) + c2M * (matrix(rep(Gbest,particle), nrow = particle, byrow = TRUE) - X)
    V <- clamp(V, VMAX)
    X <- X + V
    X <- clamp(X, LMAX)
    
    if(parallel){
      f.x <- parallel::parRapply(cl= cl, x=X, FUN=fun)  
    }
    else { f.x <- apply(X,MARGIN = 1,fun) }
    funeval <- funeval + particle
    
    tmp <- ( which( f.x < Pfit ))
    if(length(tmp) >0){
      Pfit[tmp] <- f.x[tmp]
      Pbest[tmp,] <- X[tmp,]
    }
    tmp <- which.min(f.x)
    if(f.x[tmp] < Gfit){
      Gbest <- X[tmp,]
      Gfit <- f.x[tmp]
      G.index <- tmp      
    }
    
    relerror <- max(abs(Pfit - Pfit.prior)/Pfit.prior)
    if(verbose) cat("Relative error at iteration ", iter, " : ", relerror, "\n")
    relconv <- (relerror < reltol)
  }
  parallel::stopCluster(cl) 
  if (relconv) cat("reltol converged after ", iter, "iterations and ", funeval," function evaluations!")
  param <- Gfit
  fit <- Gbest
  cat("\n##################################\n")
  cat("Optimal parameters: ", param, " best fit : ", fit)  
}







################################################################################
##                          Method 1: resetting PSO                           ##
################################################################################

PSO.resetting <- function(fun1,
                          fun2,                            
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
    maxit=1000, 
    maxfn=Inf,
    c1_w= 0.5+log(2), 
    c2_w= 0.5+log(2), 
    c1_cf= 2.05,
    c2_cf = 2.05,
    use.IW= FALSE, IW.w = 1/(2*log(2)),
    use.CF= TRUE,  
    lambda= 1,
    switch_iter = sample(10:40,1),
    reset_iter = 20,        
    reltol=sqrt(.Machine$double.eps), 
    verbose = TRUE
    )
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  if (length(noNms <- namc[!namc %in% nmsC])) 
    warning("[Unknown names in control: ", paste(noNms, collapse = ", "), " (not used) !]")         
  particle <- con[["particle"]] 
  maxit <- con[["maxit"]] 
  maxfn <- con[["maxfn"]]
  c1_w <- con[["c1_w"]] 
  c2_w <- con[["c2_w"]] 
  c1_cf <- con[["c1_cf"]]
  c2_cf <- con[["c2_cf"]]
  use.IW <- as.logical(con[["use.IW"]])
  IW.w <- con[["IW.w"]]
  use.CF <- as.logical(con[["use.CF"]])
  lambda <- con[["lambda"]]
  reltol <- con[["reltol"]]                   
  switch_iter <- con[["switch_iter"]]
  reset_iter <- con[["reset_iter"]]
  verbose <- con[["verbose"]]
  ############################################################################  
  # 1)                      Initialisation of X and V                   s     #
  ############################################################################ 
  
  if(length(upper) != nparam || length(lower) != nparam) stop("The dimension of lower bound or upper bound is not equal to number of parameters in the model!")
  LMAX <- upper - lower
  VMAX <- lambda*LMAX      
  X <- InitializeX(nparam,particle,lower,upper)
  V <- (InitializeX(nparam,particle,lower,upper)-X)/2
  V <- t(apply(X = V, MARGIN = 1, FUN = V.clamp, VMAX))    
  
  ############################################################################  
  # 2)                              Main algorithm                           #
  ############################################################################
  
  ## first evaluation
  funeval <- 0                    # number of function evaluation
  if(minmax == "max") {fun1 = -fun1; fun2 = -fun2}  
  FUN = fun1
  f.x <- apply(X,MARGIN = 1,FUN)  
  funeval <- funeval + particle
  Pbest <- X                      # personal best positions
  Pfit <- f.x                     # personal best function fit
  G.index <- which.min(Pfit)      # global best index
  Gfit <- f.x[G.index]            # global best function fit
  Gbest <- X[G.index,]            # global best position
  
  iter <- 0  
  relconv = FALSE 
  
  ## start the main loop
  
  while (!relconv && iter < maxit){  
    Pfit.prior <- Pfit
    iter <- iter+1    
    index <- 1:particle
    for (i in index){    
      V[i,] <- IW.w*V[i,] 
      V[i,] <- V[i,] + runif(nparam, 0, c1_w)*(Pbest[i,]-X[i,])
      if (i != G.index) V[i,] <- V[i,] + runif(nparam, 0, c2_w)*(Gbest-X[i,])
      V[i,] <- V.clamp(V[i,], VMAX) 
      X[i,] <- X[i,]+V[i,]
      ## Check bounds
      temp <- (X[i,]<lower)
      if (any(temp)) {
        X[i,temp] <- lower[temp]
        V[i,temp] <- 0
      }
      temp <- X[i,]>upper
      if (any(temp)) {
        X[i,temp] <- upper[temp]
        V[i,temp] <- 0
      }
      ## evaluate the fit at new position
      f.x[i] <- FUN(X[i,])
      funeval <- funeval + 1
      if (verbose) cat("At iteration ", iter, " Particle ", i, ":\n (", X[i,], ") fit: ", f.x[i],"\n")
      ## update best values
      if (f.x[i] < Pfit[i]){
        Pbest[i,] <- X[i,]
        Pfit[i] <- f.x[i]
        if (Pfit[i] < Gfit) {
          G.index <- i
          Gfit <- Pfit[i]
          Gbest <- X[i,]        
        }
      }
    }
    
    relerror <- max(abs(Pfit - Pfit.prior)/Pfit.prior)
    print(Pfit)
    print(Pfit.prior)
    if(verbose) cat("Relative error at iteration ", iter, " : ", relerror, "\n")
  #  if(iter != switch_iter + 1){  
      relconv <- (relerror < reltol)
   # }
        
    if(iter %% reset_iter == 0){
      cat("######reset particles#########\n")
      Pbest <- X                      
      Pfit <- f.x              
      G.index <- which.min(Pfit)     
      Gfit <- Pfit[G.index]          
      Gbest <- X[G.index,]    
      print(X)
    }
    
    if(iter == switch_iter){
      FUN = fun2    
      cat("######switch particles#########\n")
      Pfit <- apply(Pbest, MARGIN = 1, FUN)
#       funeval <- funeval + particle
#       G.index <- which.min(Pfit)     
#       Gfit <- Pfit[G.index]          
#       Gbest <- Pbest[G.index,] 
      print(X)
    }    
  }
  if (relconv) cat("reltol converged after ", iter, "iterations and ", funeval," function evaluations!")
  param <- Gfit
  fit <- Gbest
  cat("##################################\n")
  cat("Optimal parameters: ", param, " best fit : ", fit)  
}




################################################################################
##                          Method 2:  PSO with sentry                        ##
################################################################################
PSO.sentry <- function(fun1,
                       fun2,                            
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
    maxit=1000, 
    maxfn=Inf,
    c1_w= 0.5+log(2), 
    c2_w= 0.5+log(2), 
    c1_cf= 2.05,
    c2_cf = 2.05,
    use.IW= FALSE, IW.w = 1/(2*log(2)),
    use.CF= TRUE, 
    lambda= 1,
    switch_iter = sample(10:40,1),
    reltol=sqrt(.Machine$double.eps),
    verbose = TRUE
  )
  nmsC <- names(con)
  con[(namc <- names(control))] <- control
  if (length(noNms <- namc[!namc %in% nmsC])) 
    warning("[Unknown names in control: ", paste(noNms, collapse = ", "), " (not used) !]")         
  particle <- con[["particle"]] 
  maxit <- con[["maxit"]] 
  maxfn <- con[["maxfn"]]
  c1_w <- con[["c1_w"]] 
  c2_w <- con[["c2_w"]] 
  c1_cf <- con[["c1_cf"]]
  c2_cf <- con[["c2_cf"]]
  use.IW <- as.logical(con[["use.IW"]])
  IW.w <- con[["IW.w"]]
  use.CF <- as.logical(con[["use.CF"]])
  lambda <- con[["lambda"]]
  reltol <- con[["reltol"]]                   
  switch_iter <- con[["switch_iter"]]
  verbose <- con[["verbose"]]
  ############################################################################  
  # 1)                      Initialisation of X and V                   s     #
  ############################################################################ 
  
  if(length(upper) != nparam || length(lower) != nparam) stop("The dimension of lower bound or upper bound is not equal to number of parameters in the model!")
  LMAX <- upper - lower
  VMAX <- lambda*LMAX      
  X <- InitializeX(nparam,particle,lower,upper)
  V <- (InitializeX(nparam,particle,lower,upper)-X)/2
  V <- t(apply(X = V, MARGIN = 1, FUN = V.clamp, VMAX)) 
  sentry <- InitializeX(nparam,1,lower,upper)
  
  ############################################################################  
  # 2)                              Main algorithm                           #
  ############################################################################
  
  ## first evaluation
  funeval <- 0                    # number of function evaluation
  if(minmax == "max") {fun1 = -fun1; fun2 = -fun2}  
  FUN = fun1
  f.x <- apply(X,MARGIN = 1,FUN)  
  funeval <- funeval + particle
  Pbest <- X                      # personal best positions
  Pfit <- f.x                     # personal best function fit
  G.index <- which.min(Pfit)      # global best index
  Gfit <- f.x[G.index]            # global best function fit
  Gbest <- X[G.index,]            # global best position
  f.sen <- FUN(sentry)
  funeval <- funeval + 1
  
  iter <- 0  
  relconv = FALSE 
  
  ## start the main loop
  
  while (!relconv && iter < maxit){  
    Pfit.prior <- Pfit
    f.sen.prior<- f.sen
    iter <- iter+1    
    index <- 1:particle
    for (i in index){    
      if(use.IW){
        V[i,] <- IW.w*V[i,] 
        V[i,] <- V[i,] + runif(nparam, 0, c1_w)*(Pbest[i,]-X[i,])
        if (i != G.index) V[i,] <- V[i,] + runif(nparam, 0, c2_w)*(Gbest-X[i,])
        V[i,] <- V.clamp(V[i,], VMAX)
        X[i,] <- X[i,]+V[i,]
      }
      if(use.CF){       
        V[i,] <- V[i,] + runif(nparam, 0, c1_cf)*(Pbest[i,]-X[i,])
        if (i != G.index) V[i,] <- V[i,] + runif(nparam, 0, c2_cf)*(Gbest-X[i,])
        V[i,] <- V.clamp(V[i,], VMAX)
        CF = 2/(c1_cf+c2_cf-2+sqrt((c1_cf+c2_cf)^2-4*(c1_cf+c2_cf)))
        X[i,] <- X[i,] + CF * V[i,]
      }
      ## Check bounds
      temp <- (X[i,]<lower)
      if (any(temp)) {
        X[i,temp] <- lower[temp]
        V[i,temp] <- 0
      }
      temp <- X[i,]>upper
      if (any(temp)) {
        X[i,temp] <- upper[temp]
        V[i,temp] <- 0
      }
      ## evaluate the fit at new position
      f.x[i] <- FUN(X[i,])
      funeval <- funeval + 1
      if (verbose) cat("At iteration ", iter, " Particle ", i, ":\n", X[i,], " fit: ", f.x[i],"\n")
      ## update best values
      if (f.x[i] < Pfit[i]){
        Pbest[i,] <- X[i,]
        Pfit[i] <- f.x[i]
        if (Pfit[i] < Gfit) {
          G.index <- i
          Gfit <- Pfit[i]
          Gbest <- X[i,]        
        }
      }
    }
    
    relerror <- max(abs(Pfit - Pfit.prior)/Pfit.prior)
    if(verbose) cat("Relative error at iteration ", iter, " : ", relerror, "\n")
   # if(iter != switch_iter + 1){  
      relconv <- (relerror < reltol)
  #  }
    
    f.sen <- FUN(sentry)
    funeval <- funeval + 1
    if(f.sen != f.sen.prior){
      Pbest <- X                      # personal best positions
      Pfit <- f.x                     # personal best function fit
      G.index <- which.min(Pfit)      # global best index
      Gfit <- Pfit[G.index]           # global best function fit
      Gbest <- Pbest[G.index,]        # global best position  
      cat("######reset particles#########\n")
      print(X)
    }
    
    if(iter == switch_iter){
      FUN = fun2    
      cat("######switch particles#########\n")
      Pfit <- apply(Pbest, MARGIN = 1, FUN)
      funeval <- funeval + particle
      G.index <- which.min(Pfit)     
      Gfit <- Pfit[G.index]          
      Gbest <- Pbest[G.index,] 
      print(X)
    }    
  }
  if (relconv) cat("reltol converged!")
  param <- Gfit
  fit <- Gbest
  cat("##################################\n")
  cat("Optimal parameters: ", param, " best fit : ", fit)  
}
  