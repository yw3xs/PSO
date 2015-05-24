# hydroPSO and openmx
require(OpenMx)
require(hydroPSO)
#############################################################
# 1. IntroSEM-UnivariateStd.R
data(multiData1)
manifests <- c("x1", "y")
multiData1Cov <- cov(multiData1[,c(1,5)])

uniRegModel <- mxModel("Univariate Regression of y on x1",
                       type="RAM",
                       manifestVars=manifests,
                       mxPath(from="x1", to="y", arrows=1, 
                              free=TRUE, values=.2, labels="b1"),
                       mxPath(from=manifests, arrows=2, 
                              free=TRUE, values=.8, labels=c("VarX1", "VarE")),
                       mxData(observed=multiData1Cov, type="cov", numObs=500),
                       mxComputeOnce('fitfunction', 'fit')
)

uniRegModelOut <- mxRun(uniRegModel, suppressWarnings=TRUE)
uniRegModelOut$output$fit

modeleval <- function(par){
  model.labels <- names(omxGetParameters(uniRegModel))
  uniRegModel <- omxSetParameters(model = uniRegModel,labels = model.labels, values = par)
  modelout <- mxRun(uniRegModel)
  return (modelout$output$fit)  # 1312.985 is the known optimal value
}

# use package hydroPSO

hydroPSO(fn = 'modeleval', lower = rep(0.0001,3), upper = rep(2,3))

# expected results

expectVal <- c(0.669178, 1.138707, 1.650931)

expectSE <-c(0.053902, 0.07209, 0.104518)

expectMin <- 1312.985

omxCheckCloseEnough(expectVal, uniRegModelOut$output$estimate, 0.001)

omxCheckCloseEnough(expectSE, 
                    as.vector(uniRegModelOut$output$standardError), 0.001)

omxCheckCloseEnough(expectMin, uniRegModelOut$output$minimum, 0.001)



#################################################################
# 2. IntroSEM-UnivariateRaw.R

data(multiData1)

manifests <- c("x1", "y")

uniRegModelRaw <- mxModel("FIML Univariate Regression of y on x1",
                          type="RAM",
                          manifestVars=manifests,
                          mxPath(from="x1", to="y", arrows=1, 
                                 free=TRUE, values=.2, labels="b1"),
                          mxPath(from=manifests, 
                                 arrows=2, free=TRUE, values=.8, 
                                 labels=c("VarX1", "VarE")),
                          mxPath(from="one", to=manifests, 
                                 arrows=1, free=TRUE, values=.1, 
                                 labels=c("MeanX1", "MeanY")),
                          mxData(observed=multiData1, type="raw"),
                          mxComputeOnce('fitfunction', 'fit')
)
uniRegModelRawOut <- mxRun(uniRegModelRaw, suppressWarnings=TRUE)
uniRegModelRawOut$output$fit

modeleval <- function(par){
  model.labels <- names(omxGetParameters(uniRegModelRaw))
  uniRegModelRaw <- omxSetParameters(model = uniRegModelRaw,labels = model.labels, values = par)
  modelout <- mxRun(uniRegModelRaw)
  return (modelout$output$fit)  
}

hydroPSO(fn = 'modeleval', lower = c(-10,0.0001,0.0001,-10,-10), upper = c(10,2,2,10,10))

# check values: uniRegModelRawOut

expectVal <- c(0.669179, 1.13643, 1.647629, 0.984894, 3.189368)

expectSE <-c(0.053849, 0.071873, 0.104204, 0.047674, 0.078154)

expectMin <- 3151.492

#######################################################
# 3. IntroSEM-BivariateStd.R

data(multiData1) 

manifests <- c("x1", "x2", "y")
multiData1Cov <- cov(multiData1[,c(1,2,5)])

biRegModel <- mxModel("Bivariate Regression of y on x1 and x2",
                      type="RAM",
                      manifestVars=manifests,
                      mxPath(from=c("x1","x2"), to="y", 
                             arrows=1, 
                             free=TRUE, values=.2, labels=c("b1", "b2")),
                      mxPath(from=manifests, 
                             arrows=2, 
                             free=TRUE, values=.8, 
                             labels=c("VarX1", "VarX2", "VarE")),
                      mxPath(from="x1", to="x2",
                             arrows=2, 
                             free=TRUE, values=.2, 
                             labels=c("CovX1X2")),
                      mxData(observed=multiData1Cov, type="cov", numObs=500),
                      mxComputeOnce('fitfunction', 'fit')
)
biRegModelOut <- mxRun(biRegModel, suppressWarnings=TRUE)

modeleval <- function(par){
  model.labels <- names(omxGetParameters(biRegModel))
  biRegModel <- omxSetParameters(model = biRegModel,labels = model.labels, values = par)
  modelout <- mxRun(biRegModel)
  return (modelout$output$fit)  
}

hydroPSO(fn = 'modeleval', lower = c(-10,-10,0.0001,-10,0.0001,0.0001), upper = c(10,10,2,10,2,2))

# 4. simpleconstraint.R (not working)

constraint <- mxConstraint(A > B, name = 'AdominatesB')

# Constrain matrix 'K' to be equal to matrix 'limit'
model <- mxModel(model="con_test",
                 mxMatrix(type="Full", nrow=2, ncol=2, free=TRUE, name="K"),
                 mxMatrix(type="Full", nrow=2, ncol=2, free=FALSE, name="limit", values=1:4),
                 mxConstraint(K == limit, name = "Klimit_equality"),
                 mxAlgebra(min(K), name="minK"),
                 mxFitFunctionAlgebra("minK"),
                 mxComputeOnce('fitfunction', 'fit')
)
ModelOut <- mxRun(model, suppressWarnings=TRUE)

modeleval <- function(par){
  model.labels <- names(omxGetParameters(model))
  model <- omxSetParameters(model = model,labels = model.labels, values = par)
  modelout <- mxRun(model)
  return (modelout$output$fit)  
}

hydroPSO(fn = 'modeleval', lower = c(-10,-10,-10,-10), upper = c(10,10,10,10))

