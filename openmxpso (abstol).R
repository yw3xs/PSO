# ---------------------------------------------------------------------
# Program: UnivariateStd-OpenMx100214.R
#  Author: Steven M. Boker
#   Date: Sun Feb 14 12:13:20 EST 2010
#
# This program fits a univariate model to the multiData simulated data.

# ----------------------------------
# Read libraries and set options.
require(OpenMx)
require(hydroPSO)
options(digits = 15)
# ----------------------------------
# Read the data and print descriptive statistics.
data(multiData1)
# ----------------------------------
# Build an OpenMx univariate regression models using y and x1, x2, x3, x4 separately 

manifests <- c("x", "y")
multiData1Cov1 <- cov(multiData1[,c(1,5)])
multiData1Cov2 <- cov(multiData1[,c(2,5)])
multiData1Cov3 <- cov(multiData1[,c(3,5)])
multiData1Cov4 <- cov(multiData1[,c(4,5)])

multiData1Cov5 <- matrix(c(13.1387075, 11.7619982, 11.7619982, 12.1608437),2,2)

dimnames(multiData1Cov1) <- list(c("x", "y"), c("x", "y"))
dimnames(multiData1Cov2) <- list(c("x", "y"), c("x", "y"))
dimnames(multiData1Cov3) <- list(c("x", "y"), c("x", "y"))
dimnames(multiData1Cov4) <- list(c("x", "y"), c("x", "y"))
dimnames(multiData1Cov5) <- list(c("x", "y"), c("x", "y"))

                                 
# To change dimnames of each matrix, the following lapply is not working, why?
# because the whole paradigm of R is to "pass by value". Environments are typically the only objects that can be passed by reference in R. The normal way to handle it - just return the modified value. 
# See http://stackoverflow.com/questions/8419877/modify-variable-within-r-function
# lapply(list(multiData1Cov1, multiData1Cov2, multiData1Cov3, multiData1Cov4), function(data){dimnames(data) <<- list(c("x","y"), c("x","y"))})

uniRegModel <- function(data, optimize = FALSE){
  model <- mxModel("Univariate Regression of y on x",
                       type="RAM",
                       manifestVars=manifests,
                       mxPath(from="x", to="y", arrows=1, 
                              free=TRUE, values=.2, labels="b"),
                       mxPath(from=manifests, arrows=2, 
                              free=TRUE, values=.8, labels=c("VarX", "VarE")),
                       mxData(observed=data, type="cov", numObs=500)
  )
  if(!optimize)  model <- mxModel(model = model, mxComputeOnce('fitfunction', 'fit'))
  return (model)
}
model1 <- uniRegModel(multiData1Cov1, TRUE)
model2 <- uniRegModel(multiData1Cov2, TRUE)
model3 <- uniRegModel(multiData1Cov3, TRUE)
model4 <- uniRegModel(multiData1Cov4, TRUE)
model5 <- uniRegModel(multiData1Cov5, TRUE)

# find the optimal value
ModelOut1 <- mxRun(model1, suppressWarnings=TRUE)
model1fit <- ModelOut1$output$fit     
# [1] 1312.9853580061
# b      VarX      VarE 
# 0.6691782 1.1387074 1.6509311 
ModelOut2 <- mxRun(model2, suppressWarnings=TRUE)
model2fit <- ModelOut2$output$fit
# [1] 1453.908
# b      VarX      VarE 
# 0.669178184696338 1.138707360984092 1.650931103515756 
ModelOut3 <- mxRun(model3, suppressWarnings=TRUE)
model3fit <- ModelOut3$output$fit
# [1] 1499.465
# b      VarX      VarE 
# 0.6404838 2.1069835 1.2965179 
ModelOut4 <- mxRun(model4, suppressWarnings=TRUE)
model4fit <- ModelOut4$output$fit
# [1] 1526.066
# b      VarX      VarE 
# 0.6359235 2.5581003 1.1263511 
ModelOut5 <- mxRun(model5, suppressWarnings=TRUE)
model5fit <- ModelOut5$output$fit 
# [1] 2527.40484377555
#  b     VarX     VarE 
#  0.89521726817729 13.13870136144152  1.63129984228987 






model1 <- uniRegModel(multiData1Cov1)
model2 <- uniRegModel(multiData1Cov2)
model3 <- uniRegModel(multiData1Cov3)
model4 <- uniRegModel(multiData1Cov4)
model5 <- uniRegModel(multiData1Cov5)

modeleval <- function(par, model){
  model.labels <- names(omxGetParameters(model))
  model <- omxSetParameters(model = model,labels = model.labels, values = par)
  modelout <- mxRun(model)
  return (modelout$output$fit)  
}

# My implementation (from R package PSO)
#----------------------------------------------------------------
#sink('analysis-output.txt')
npar <- length(omxGetParameters(model1))       # number of parameters in the model
swarm <- 10   # swarm size
lower = 0.0001   # lower bound
upper = 15    # upper bound
abstol = 1e-3  # absolute tolerance, optimal value known in advance.
reltol = 1e-5  # relative tolerance, see hydroPSO, no requirement for optimal value

switch_iter <- 20 # After how many iterations the model will be changed?
reset_iter <- 15  # After how many iterations the particles will be reset?

#initialization
mrunif <- function(npar,swarm,lower,upper){
  return(matrix(runif(npar*swarm,lower,upper),nrow = npar, ncol = swarm))
}

X <- mrunif(npar,swarm,lower,upper)
V <- (mrunif(npar,swarm,lower,upper)-X)/2

# experiment: starting from model 1, after 1000 (?) iterations, change to model 2 and so on
# Method 1: after 500 (?) itertion, reset the swarm (set pbest to the current position)

# first evaluations

f.x <- apply(X,2,modeleval, model=model1)
P <- X   # personal best positions
f.p <- f.x  # personal best function fit
P.improved <- rep(FALSE,swarm)
i.best <- which.min(f.p)  # global best index
f.gbest <- f.x[i.best]   # global best function fit
G <- X[, i.best]  # global best position
error <- abs(f.gbest - model1fit)


## Iterations
iter <- 1


percent <- .15 # percentage of informants for each particle
w <- 1/(2*log(2)) # inertia
c.p <- .5+log(2)
c.g <- .5+log(2)
model = model1
modelfit = model1fit
while (error > abstol){

  iter <- iter+1
  #links <- matrix(runif(swarm*swarm,0,1)<=percent,swarm,swarm)  # 
  #diag(links) <- TRUE
  index <- 1:swarm
  for (i in index){
#    j <- which(links[,i])[which.min(f.p[links[,i]])]
    V[,i] <- w*V[,i] 
    V[,i] <- V[,i]+runif(npar,0,c.p)*(P[,i]-X[,i])
    if (i!=i.best) V[,i] <- V[,i]+runif(npar,0,c.g)*(G-X[,i])
    X[,i] <- X[,i]+V[,i]
    ## Check bounds
    temp <- X[,i]<lower
    if (any(temp)) {
      X[temp,i] <- lower
      V[temp,i] <- 0
    }
    temp <- X[,i]>upper
    if (any(temp)) {
      X[temp,i] <- upper
      V[temp,i] <- 0
    }
    f.x[i] <- modeleval(X[,i], model)
    cat("parameter:", X[,i], "fit: ", f.x[i],"\n")
    if (f.x[i]<f.p[i]){
      P[,i] <- X[,i]
      f.p[i] <- f.x[i]
      if (f.p[i]<f.gbest) {
        i.best <- i
        f.gbest <- f.p[i]
        G <- X[, i.best]        
      }
    }
  }
  

#   d <- X-G
#   d <- sqrt(max(colSums(d*d)))
#   if (d<reltol) {
#     X <- mrunif(npar,swarm,lower,upper)
#     V <- (mrunif(npar,swarm,lower,upper)-X)/2
#   }
  
  
 # init.links <- f.p[i.best]==error
  error <- abs(f.gbest - modelfit)
  cat("error at iteration ", iter, ":", error, "\n")
  if(iter %% reset_iter == 0){
    P <- X   # personal best positions
    f.p <- f.x  # personal best function fit
    i.best <- which.min(f.p)  # global best index
    f.gbest <- f.x[i.best]   # global best function fit
    G <- X[, i.best]  # global best position  
    cat("######reset particles#########\n")
    print(X)
  }
  
#   if((iter %/% switch_iter) %% 2 == 0 ){
#     model = model1
#   }
#   else{
#     model = model5
#   }
  if(iter == switch_iter){
    model = model5
    modelfit = model5fit
    cat("######switch particles#########\n")
    print(X)}
}
  






if (error<=abstol) cat("converged!")
par <- P[,i.best]
modeleval(par, model)

#sink()
#---------------------

# check values: uniRegModelOut

# expectVal <- c(0.669178, 1.138707, 1.650931)
# 
# expectSE <-c(0.053902, 0.07209, 0.104518)
# 
# expectMin <- 1312.985
# 
# omxCheckCloseEnough(expectVal, uniRegModelOut$output$estimate, 0.001)
# 
# omxCheckCloseEnough(expectSE, 
#                     as.vector(uniRegModelOut$output$standardError), 0.001)
# 
# omxCheckCloseEnough(expectMin, uniRegModelOut$output$minimum, 0.001)
# 
# 



