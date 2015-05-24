# Test functions for dynamic PSO

################################################################################
##                          IntroSEM_UnivariateRaw                            ##
################################################################################
require(OpenMx)
options(digits = 15)
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

### find the optimal value

ModelOut1 <- mxRun(model1, suppressWarnings=TRUE)
model1fit <- ModelOut1$output$fit     
model1estimate <- ModelOut1$output$estimate
model1fit
model1estimate
# fit: 3151.4916549273
# b1              VarX              VarE             MeanX             MeanY 
# 0.669178103010170 1.136430307702897 1.647628906697996 0.984894604002143 3.189369536976615 
ModelOut2 <- mxRun(model2, suppressWarnings=TRUE)
model2fit <- ModelOut2$output$fit     
model2estimate <- ModelOut2$output$estimate
model2fit
model2estimate
# fit: 3364.99912689123
# b1              VarX              VarE             MeanX             MeanY 
# 0.635923340142127 2.552983287073350 1.124098309147817 4.055540415024234 1.269426021674500 

### set the model to MxComputeonce
model1 <- uniRegModelRaw(Data1)
model2 <- uniRegModelRaw(Data2)


fun1 <- function(par){
  model.labels <- names(omxGetParameters(model1))
  model <- omxSetParameters(model = model1,labels = model.labels, values = par)
  modelout <- mxRun(model)
  return (modelout$output$fit)  
}

fun2 <- function(par){
  model.labels <- names(omxGetParameters(model2))
  model <- omxSetParameters(model = model2,labels = model.labels, values = par)
  modelout <- mxRun(model)
  return (modelout$output$fit)  
}

PSO.resetting(fun1 = fun1, fun2 = fun2, nparam = 5, lower = c(-2,.00001,.00001,-3,-3), upper = c(2, 4, 4, 5, 5), minmax = "min", control = list(particle = 40, verbose = FALSE,parallel=TRUE))
PSO(fun1,nparam = 5,lower = c(-2,.00001,.00001,-3,-3), upper = c(2, 4, 4, 5, 5), minmax = "min", control = list(particle = 40, verbose = FALSE, parallel=TRUE))
GPSO(fun1,nparam = 5,lower = c(-2,.00001,.00001,-3,-3), upper = c(2, 4, 4, 5, 5), minmax = "min", control = list(particle = 40, verbose = FALSE, parallel=TRUE))
##################################################
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

model3 <- mxModel("both", twin, obj, modelMZ, modelDZ)

model3Out <- mxRun(model3)

expectedACE <- c(.6, .2, .2)
observedACE <- c(model3Out$modelMZ.A$result, 
                 model3Out$modelMZ.C$result, model3Out$modelMZ.E$result)

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

model4 <- mxModel("both2", twin2, obj2, modelMZ2, modelDZ2)

model4Out <- mxRun(model4)

expectedACE2 <- c(0, 3.15, 2.85)
observedACE2 <- c(model4Out$modelMZ2.A2$result, 
                  model4Out$modelMZ2.C2$result, model4Out$modelMZ2.E2$result)

################expected results##############
model4Out$output$estimate
# param1       param2       param3 
# 1.489984e-07 1.774824e+00 1.688194e+00 
model4Out$output$fit
# [1] 1041.692

model3Out$output$fit
# [1] 266.376
model3Out$output$estimate
# param1    param2    param3 
# 0.7745967 0.4472135 0.4472135 

c(model4Out$modelMZ2.A2$result, 
  model4Out$modelMZ2.C2$result, model4Out$modelMZ2.E2$result)
# [1] 2.220054e-14 3.150000e+00 2.850000e+00
c(model3Out$modelMZ.A$result, 
     model3Out$modelMZ.C$result, model3Out$modelMZ.E$result)
# [1] 0.6000000 0.1999999 0.1999999
##########################
model3 <- mxModel(model3, mxComputeOnce('fitfunction', 'fit'))
model4 <- mxModel(model4, mxComputeOnce('fitfunction', 'fit'))

fun3 <- function(par){
  model.labels <- names(omxGetParameters(model3))
  model <- omxSetParameters(model = model3,labels = model.labels, values = par)
  modelout <- mxRun(model)
  return (modelout$output$fit)  
}

fun4 <- function(par){
  model.labels <- names(omxGetParameters(model4))
  model <- omxSetParameters(model = model4,labels = model.labels, values = par)
  modelout <- mxRun(model)
  return (modelout$output$fit)  
}

PSO(fun3,nparam = 3,lower = rep(.00001,3), upper = rep(3, 3), minmax = "min", control = list(particle = 40, verbose = FALSE, parallel=TRUE))
PSO.resetting(fun3,fun4,nparam = 3,lower = rep(.00001,3), upper = rep(3, 3), minmax = "min", control = list(particle = 40, verbose = FALSE, parallel=TRUE))



funtest <- function(par){
  model.labels <- names(omxGetParameters(grpModel))
  model <- omxSetParameters(model = grpModel,labels = model.labels, values = par)
  modelout <- mxRun(model)
  return (modelout$output$fit)  
}



