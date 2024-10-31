require(randomForest)
require(tidyverse)

compute_integratedVariances <- function(obj,exps,X,n,p,K,nSamp,nMC){

  ## Create indices for group of interest and for remaining exposures
  groupJ = exps
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

  ## Monte carlo integral of quantity
  ## first define X mat
  newX = matrix(NA, nrow = nMC*nSamp, ncol = p)
  ## add in Xneg which changes across nSamp, but not across nMC
  for(jj in groupNegJ){
    newX[,jj] = rep(X[samp,jj],each=nMC)
  }
  ## Sample nMC values from their conditional distribution (Normal approx)
  for (ni in 1 : nSamp) {
    newX[(ni-1)*nMC+1:nMC,groupJ] = t(matrix(rep(predJ[samp[ni],], nMC), ncol=nMC)) +
      t(t(cholVarJ) %*% matrix(rnorm(nMC*length(groupJ)), ncol=nMC))
  }

  ## Now make predictions on stacked newX
  pred <- predict_MVmix(obj, newX = newX, include_intercept=TRUE, allx=TRUE)

  ## Take average across MC samples
  integratedVariances <- rep(NA,K)
  for (kk in 1 : K) {
    ## take mean across nMC
    dfpred <- data.frame(nsamp=rep(1:nSamp,each=nMC),
                         pred=pred$summary[[kk]]$mean)
    predmean <- dfpred %>% group_by(nsamp) %>% summarize_all(mean)

    ## Final variance estimates across nSamp
    integratedVariances[kk] <- var(predmean$pred)

  }

  return(integratedVariances)
}

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
  if(!is.list(exposures)){
    exposures <- list(exposures)
  }

  ## numerators
  integratedVariances <- lapply(exposures,function(exp_set){compute_integratedVariances(obj,exp_set,X,n,p,K,nSamp,nMC)})

  ## Variance estimates without integrating out X_j
  pred <- predict_MVmix(obj, newX = X, include_intercept=TRUE, allx=TRUE)
  origVariances = sapply(1:K,function(kk){var(pred$summary[[kk]]$mean)})

  ## Final value of VIM
  VIMs = sapply(integratedVariances,function(x) 1 - x/origVariances)

  return(list(VIM=VIMs,
              varF=origVariances))
}
