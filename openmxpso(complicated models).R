# pso on dynamic environments
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


# 2. ACEDuplicateMatrices.R (good)

require(OpenMx)

varNames <- c('x','y')

dataMZ <- mxData(matrix(c(1,.8,.8,1), nrow = 2, ncol=2, 
                        dimnames = list(varNames,varNames)), type="cov",
                 numObs=100)
dataDZ <- mxData(matrix(c(1,.5,.5,1), nrow = 2, ncol=2,
                        dimnames = list(varNames,varNames)), type="cov",
                 numObs=100)

X <- mxMatrix("Full",.6,free=TRUE, labels='param1', nrow=1, ncol=1, name="X")
Y <- mxMatrix("Full",.6,free=TRUE, labels='param2', nrow=1, ncol=1, name="Y")
Z <- mxMatrix("Full",.6,free=TRUE, labels='param3', nrow=1, ncol=1, name="Z")
h <- mxMatrix("Full",.5,free=FALSE, nrow=1, ncol=1, name="h")
A <- mxAlgebra(X * t(X), name="A")
C <- mxAlgebra(Y * t(Y), name="C")
E <- mxAlgebra(Z * t(Z), name="E")
cMZ <- mxAlgebra(rbind(cbind(A+C+E,A+C),cbind(A+C,A+C+E)), 
                 name="cMZ")
cDZ <- mxAlgebra(rbind(cbind(A+C+E,h%x%A+C),cbind(h%x%A+C,A+C+E)),
                 name="cDZ")

objMZ <- mxExpectationNormal("cMZ", dimnames = varNames)
objDZ <- mxExpectationNormal("cDZ", dimnames = varNames)

modelMZ <- mxModel("modelMZ", dataMZ, X,Y,Z,A,C,E,cMZ, objMZ, mxFitFunctionML())
modelDZ <- mxModel("modelDZ", dataDZ, X,Y,Z,h,A,C,E,cDZ, objDZ, mxFitFunctionML())

twin <- mxAlgebra(modelMZ.objective + modelDZ.objective, name="twin")
obj <- mxFitFunctionAlgebra("twin")

model <- mxModel("both", twin, obj, modelMZ, modelDZ)

modelOut <- mxRun(model)

expectedACE <- c(.6, .2, .2)
observedACE <- c(modelOut$modelMZ.A$result, 
                 modelOut$modelMZ.C$result, modelOut$modelMZ.E$result)

omxCheckCloseEnough(expectedACE, observedACE, epsilon = 10 ^ -4)
 
################ ACE model2
varNames <- c('x2','y2')
dataMZ2 <- mxData(matrix(c(10,3.8,3.8,5), nrow = 2, ncol=2, 
                        dimnames = list(varNames,varNames)), type="cov",
                 numObs=100)
dataDZ2 <- mxData(matrix(c(6,2.5,2.5,3), nrow = 2, ncol=2,
                        dimnames = list(varNames,varNames)), type="cov",
                 numObs=100)

X2 <- mxMatrix("Full",1.6,free=TRUE, labels='param1', nrow=1, ncol=1, name="X2")
Y2 <- mxMatrix("Full",1.6,free=TRUE, labels='param2', nrow=1, ncol=1, name="Y2")
Z2 <- mxMatrix("Full",1.6,free=TRUE, labels='param3', nrow=1, ncol=1, name="Z2")
h2 <- mxMatrix("Full",.5,free=FALSE, nrow=1, ncol=1, name="h2")
A2 <- mxAlgebra(X2 * t(X2), name="A2")
C2 <- mxAlgebra(Y2 * t(Y2), name="C2")
E2 <- mxAlgebra(Z2 * t(Z2), name="E2")
cMZ2 <- mxAlgebra(rbind(cbind(A2+C2+E2,A2+C2),cbind(A2+C2,A2+C2+E2)), 
                 name="cMZ2")
cDZ2 <- mxAlgebra(rbind(cbind(A2+C2+E2,h2%x%A2+C2),cbind(h2%x%A2+C2,A2+C2+E2)),
                 name="cDZ2")

objMZ2 <- mxExpectationNormal("cMZ2", dimnames = varNames)
objDZ2 <- mxExpectationNormal("cDZ2", dimnames = varNames)

modelMZ2 <- mxModel("modelMZ2", dataMZ2, X2,Y2,Z2,A2,C2,E2,cMZ2, objMZ2, mxFitFunctionML())
modelDZ2 <- mxModel("modelDZ2", dataDZ2, X2,Y2,Z2,h2,A2,C2,E2,cDZ2, objDZ2, mxFitFunctionML())

twin2 <- mxAlgebra(modelMZ2.objective + modelDZ2.objective, name="twin2")
obj2 <- mxFitFunctionAlgebra("twin2")

model2 <- mxModel("both2", twin2, obj2, modelMZ2, modelDZ2)

modelOut2 <- mxRun(model2)

expectedACE2 <- c(0, 3.15, 2.85)
observedACE2 <- c(modelOut2$modelMZ2.A2$result, 
                 modelOut2$modelMZ2.C2$result, modelOut2$modelMZ2.E2$result)

################expected results##############
modelOut2$output$estimate
param1       param2       param3 
1.489984e-07 1.774824e+00 1.688194e+00 
modelOut2$output$fit
[1] 1041.692

modelOut$output$fit
[1] 266.376
param1    param2    param3 
0.7745967 0.4472135 0.4472135 

c(modelOut2$modelMZ2.A2$result, 
     modelOut2$modelMZ2.C2$result, modelOut2$modelMZ2.E2$result)
[1] 2.220054e-14 3.150000e+00 2.850000e+00
c(modelOut$modelMZ.A$result, 
  +   modelOut$modelMZ.C$result, modelOut$modelMZ.E$result)
[1] 0.6000000 0.1999999 0.1999999
##########################
model <- mxModel(model, mxComputeOnce('fitfunction', 'fit'))
model2 <- mxModel(model2, mxComputeOnce('fitfunction', 'fit'))

modeleval <- function(par, model){
  model.labels <- names(omxGetParameters(model))
  model <- omxSetParameters(model = model,labels = model.labels, values = par)
  modelout <- mxRun(model)
  return (modelout$output$fit)  
}













#########################################################
# pso implementation
#########################################################
npar <- length(omxGetParameters(model))       # number of parameters in the model
swarm <- 20   # swarm size
lower = 0.0001   # lower bound
upper = 3   # upper bound
reltol = 1e-8  # relative tolerance, see hydroPSO, no requirement for optimal value
relconv = FALSE 

switch_iter <- 20 # After how many iterations the model will be changed?
reset_iter <- 15  # After how many iterations the particles will be reset?

#initialization
mrunif <- function(npar,swarm,lower,upper){
  return(matrix(runif(npar*swarm,lower,upper),nrow = npar, ncol = swarm))
}

X <- mrunif(npar,swarm,lower,upper)
V <- (mrunif(npar,swarm,lower,upper)-X)/2


# first evaluations

f.x <- apply(X,2,modeleval, model=model)
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
testmodel = model

#modelfit = model1fit
while (!relconv){  
  f.p.prior <- f.p
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
  
  relerror <- max(abs(f.p - f.p.prior))/f.p.prior[which.max(abs(f.p - f.p.prior))]
  cat("relative error at iteration ", iter, ":", relerror, "\n")
  if(iter != switch_iter + 1){  
    if(relerror < reltol) relconv <- TRUE
  }
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
    testmodel = model2    
    cat("######switch particles#########\n")
    f.p <- apply(P,2,modeleval, model=testmodel)
    print(X)}
}


if (relconv) cat("reltol converged!")
par <- P[,i.best]
modeleval(par, testmodel)





