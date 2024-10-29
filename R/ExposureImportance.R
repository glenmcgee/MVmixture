require(randomForest)

ExposureImportance = function(obj, exposures,
                              nMC = 250,
                              nSamp=NULL) {

  n = dim(obj$u)[2]
  if(is.null(nSamp)){
    nSamp = n
  }
  p <- obj$const$L[[1]]
  K <- obj$const$K
  X <- obj$const$X[[1]]

  ## Create indices for group of interest and for remaining exposures
  groupJ = exposures
  groupNegJ = (1:p)[-groupJ]

  ## For each exposure in groupJ, we need predictions given the other exposures
  predJ = matrix(NA, n, length(groupJ))
  diffJ = matrix(NA, n, length(groupJ))

  for (jj in 1 : length(groupJ)) {
    ## Need a model for the this exposure given the others
    modJJ = randomForest(x = X[,groupNegJ], y = X[,groupJ[jj]])

    ## Extract predictions and residuals
    predJ[,jj] = modJJ$predicted
    diffJ[,jj] = predJ[,jj] - X[,groupJ[jj]]
  }

  ## Get conditional variance matrix
  varJ = var(diffJ)
  cholVarJ = chol(varJ)

  ## nSamp is typically going to be n unless n is very large
  samp = sample(1:n, nSamp, replace=FALSE)

  ## Store predictions once X_j is integrated out
  predStore = matrix(NA, nSamp, K)

  ## Monte carlo integral of quantity
  for (ni in 1 : nSamp) {
    ## Which observation are we looking at first
    Xneg = X[samp[ni],groupNegJ]

    ## Sample nMC values from their conditional distribution
    XjSamp = t(matrix(rep(predJ[samp[ni],], nMC), ncol=nMC)) +
      t(t(cholVarJ) %*% matrix(rnorm(nMC*length(groupJ)), ncol=nMC))

    ## Now make predictions there
    newX = matrix(NA, nrow = nMC, ncol = p)
    newX[,groupJ] = XjSamp
    newX[,groupNegJ] = t(matrix(rep(Xneg, nMC), nrow=length(groupNegJ)))
    pred <- predict_MVmix(obj, newX = newX, include_intercept=TRUE, allx=TRUE)

    ## Take average across MC samples
    for (kk in 1 : K) {
      predStore[ni,kk] = mean(pred$summary[[kk]]$mean)
    }
  }

  ## Final variance estimates
  integratedVariances = apply(predStore, 2, var)

  ## Variance estimates without integrating out X_j
  pred <- predict_MVmix(obj, newX = X, include_intercept=TRUE, allx=TRUE)
  origVariances = rep(NA, K)

  for (kk in 1 : K) {
    origVariances[kk] = var(pred$summary[[kk]]$mean)
  }

  ## Final value of VIM
  VIMs = 1 - integratedVariances/origVariances

  return(VIMs)
}
