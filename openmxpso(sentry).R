# pso on dynamic environments (sentry particle)
# 1. UnivariateRaw-OpenMx (good)
require(OpenMx)
data(multiData1)

manifests <- c("x", "y")
Data1 <- multiData1[,c(1,5)]
Data2 <- multiData1[,c(4,5)]
names(Data1) <- c("x","y")
names(Data2) <- c("x","y")

uniRegModelRaw <- function(data, optimize = FALSE){
  model <- mxModel("FIML Univariate Regression of y on x",
                   type="RAM",
                   manifestVars=manifests,
                   mxPath(from="x", to="y", arrows=1, 
                          free=TRUE, values=.2, labels="b1"),
                   mxPath(from=manifests, 
                          arrows=2, free=TRUE, values=.8, 
                          labels=c("VarX", "VarE")),
                   mxPath(from="one", to=manifests, 
                          arrows=1, free=TRUE, values=.1, 
                          labels=c("MeanX", "MeanY")),
                   mxData(observed=data, type="raw")
  )
  if(!optimize)  model <- mxModel(model = model, mxComputeOnce('fitfunction', 'fit'))
  return (model)
}
model1 <- uniRegModelRaw(Data1, TRUE)
model2 <- uniRegModelRaw(Data2, TRUE)
# find the optimal value
ModelOut1 <- mxRun(model1, suppressWarnings=TRUE)
model1fit <- ModelOut1$output$fit     
model1estimate <- ModelOut1$output$estimate
# fit: 3151.4916549273
# b1              VarX              VarE             MeanX             MeanY 
# 0.669178103010170 1.136430307702897 1.647628906697996 0.984894604002143 3.189369536976615 
ModelOut2 <- mxRun(model2, suppressWarnings=TRUE)
model2fit <- ModelOut2$output$fit     
model2estimate <- ModelOut2$output$estimate
# fit: 3364.99912689123
# b1              VarX              VarE             MeanX             MeanY 
# 0.635923340142127 2.552983287073350 1.124098309147817 4.055540415024234 1.269426021674500 

model1 <- uniRegModelRaw(Data1)
model2 <- uniRegModelRaw(Data2)

modeleval <- function(par, model){
  model.labels <- names(omxGetParameters(model))
  model <- omxSetParameters(model = model,labels = model.labels, values = par)
  modelout <- mxRun(model)
  return (modelout$output$fit)  
}

#####################################################################
npar <- length(omxGetParameters(model1))       # number of parameters in the model
swarm <- 20   # swarm size
lower = 0.0001   # lower bound
upper = 8   # upper bound
reltol = 1e-8  # relative tolerance, see hydroPSO, no requirement for optimal value
relconv = FALSE 


switch_iter <- sample(10:40,1) # After how many iterations the model will be changed?

#initialization
mrunif <- function(npar,swarm,lower,upper){
  return(matrix(runif(npar*swarm,lower,upper),nrow = npar, ncol = swarm))
}

X <- mrunif(npar,swarm,lower,upper)
V <- (mrunif(npar,swarm,lower,upper)-X)/2
sentry <- mrunif(npar,1,lower,upper)

# first evaluations

f.x <- apply(X,2,modeleval, model=model1)
f.sen <- modeleval(sentry, model1)
P <- X   # personal best positions
f.p <- f.x  # personal best function fit
P.improved <- rep(FALSE,swarm)
i.best <- which.min(f.p)  # global best index
f.gbest <- f.x[i.best]   # global best function fit
G <- X[, i.best]  # global best position
#error <- abs(f.gbest - model1fit)


## Iterations
iter <- 1


percent <- .15 # percentage of informants for each particle
w <- 1/(2*log(2)) # inertia
c.p <- .5+log(2)
c.g <- .5+log(2)
testmodel = model1

#modelfit = model1fit
while (!relconv){  
  f.p.prior <- f.p
  f.sen.prior<- f.sen
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
    f.x[i] <- modeleval(X[,i], testmodel)
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
  #error <- abs(f.gbest - modelfit)
  
  relerror <- max(abs(f.p - f.p.prior)/f.p.prior)#/f.p.prior[which.max(abs(f.p - f.p.prior))]
  print(f.p)
  print(f.p.prior)
  cat("relative error at iteration ", iter, ":", relerror, "\n")
  if(iter != switch_iter + 1){  
    if(relerror < reltol) relconv <- TRUE
  }
  
  f.sen <- modeleval(sentry, testmodel)
  if(f.sen != f.sen.prior){
    P <- X   # personal best positions
    f.p <- f.x  # personal best function fit
    i.best <- which.min(f.p)  # global best index
    f.gbest <- f.x[i.best]   # global best function fit
    G <- X[, i.best]  # global best position  
    cat("######reset particles#########\n")
    print(X)
  }
  
  
  
  if(iter == switch_iter){
    testmodel = model2    
    cat("######switch particles#########\n")
    f.p <- apply(P,2,modeleval, model=testmodel)
    print(X)}
   
}


if (relconv) cat("reltol converged!")
par <- P[,i.best]
modeleval(par, testmodel)





