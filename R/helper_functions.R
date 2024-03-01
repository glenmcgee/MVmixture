## Notes:
### When we add more exposures, we'll need absorb.cons=TRUE and then manually add an intercept. for now intercept is included in B implicitly
### For betastar: currently using N(0,1) prior on unpenalized components
##### Fix this

## function to compute pi
compute_pi <- function(params,const){

  ## compute marginals under independence
  pi_beta <- params$Vbeta*c(1,cumprod(1-params$Vbeta)[-const$C])
  pi_theta <- params$Vtheta*c(1,cumprod(1-params$Vtheta)[-const$C])

  ## joint probs
  pimat <- (pi_beta)%*%t(pi_theta) ## pi^I (under independence)
  diag(pimat) <- (1+exp(params$logrho))*diag(pimat) ## multiply diagonal by 1+rho

  ## standardize
  pistarmat <- pimat/sum(pimat)

  params$pimat <- pimat
  params$pistarmat <- pistarmat
  return(params)
}

compute_logpi <- function(params,const){

  ## compute marginals under independence
  ## making it add up to 1... I guess we need not sample this parameter
  logpi_beta <- log(params$Vbeta)+c(0,cumsum(log(1-params$Vbeta))[-const$C])
  logpi_theta <- log(params$Vtheta)+c(0,cumsum(log(1-params$Vtheta))[-const$C])

  ## joint probs
  logpimat <- matrix(logpi_beta,nrow=const$C,ncol=const$C,byrow=FALSE)+matrix(logpi_theta,nrow=const$C,ncol=const$C,byrow=TRUE) ## pi^I (under independence)
  diag(logpimat) <- log(1+exp(params$logrho))+diag(logpimat) ## multiply diagonal by 1+rho

  ## standardize
  logpistarmat <- logpimat-log(sum(exp(logpimat)))

  params$logpimat <- logpimat
  params$logpistarmat <- logpistarmat
  return(params)
}

## get basis functions
get_Btheta <- function(Xtheta,const){
  return(mgcv::PredictMat(const$SS,data=data.frame(Xtheta)))
}

## get derivatives of Basis functions
get_DerivBtheta <- function(Xtheta,const){
  return(mgcv::PredictMat(const$SSderiv,data=data.frame(Xtheta)))
}

## compute B times the corresponding beta for all obs
get_B_beta <- function(params,const){

  B_beta <- rep(NA,const$n*const$K)
  for(kk in 1:const$K){
    B_beta[const$k_index==kk] <- params$Btheta[const$k_index==kk,]%*%params$betastar[(params$Zbeta[kk]-1)*const$d+(1:const$d)]
  }
  return(B_beta)
}


## computes density of V_c ## up to proportionality constant
## called iteratively by update_V and update_V_MH
get_Vlogdensity <- function(vv,cc,params,
                         Vbeta=TRUE, ## for Vbeta, vs Vtheta
                         const){

  params$Vbeta[cc] <- vv
  params <- compute_pi(params,const)

  if(Vbeta==TRUE){
    alpha=params$alpha[1]
  }else{
    alpha=params$alpha[2]
  }

  ## log density
  return((alpha-1)*log(1-vv) + sum(sapply(1:const$K,function(kk){log(params$pistarmat[params$Zbeta[kk],params$Ztheta[kk]])})))
}

assign_betas <- function(params,const){
  for(kk in 1:const$K){
    params$beta[(kk-1)*const$d+(1:const$d)] <- params$betastar[(params$Zbeta[kk]-1)*const$d+(1:const$d)]
  }
  return(params)
}

assign_thetas <- function(params,const){
  for(kk in 1:const$K){
    params$theta[(kk-1)*const$L+(1:const$L)] <- params$thetastar[(params$Ztheta[kk]-1)*const$L+(1:const$L)]
  }
  return(params)
}


## initialize constants
initialize_const <- function(Y, ## response
                             X, ## list of exposures/mixture components
                             Z, ## confounders to be adjusted
                             ## MCMC specs
                             niter, ## number of iterations
                             nburn, ## burn-in fraction
                             nthin, ## thinning number
                             nchains, ## number of chains ## not yet implemented
                             ncores, ## number of cores for mclapply (set to 1 for non-parallel) ## only used if nchains>1
                             ## prior hyperparameters
                             maxClusters,
                             prior_alpha_beta,
                             prior_alpha_theta,
                             prior_rho,
                             prior_tau_theta,
                             prior_lambda_beta,
                             prior_lambda_theta,
                             prior_sigma2_u,
                             prior_sigma2,
                             sharedlambda,
                             DLM,
                             lagOrder,
                             ## MH tuning
                             stepsize_logrho,
                             stepsize_loglambda_theta,
                             Vgridsearch,
                             gridsize,
                             rfbtries,
                             wls_steps,
                             MHwls,
                             stepsize_theta){

  const <- list(y=c(Y), ## convert nxK matrix to a single vector of outcomes,
                X=X,
                Z=Z,
                ## MCMC specs
                niter=niter,
                nburn=nburn,
                nthin=nthin,
                nchains=nchains,
                ncores=ncores,
                ## prior hyperparameters
                maxClusters=maxClusters,
                prior_alpha_beta=prior_alpha_beta,
                prior_alpha_theta=prior_alpha_theta,
                prior_rho=prior_rho,
                prior_tau_theta=prior_tau_theta,
                prior_lambda_beta=prior_lambda_beta,
                prior_lambda_theta=prior_lambda_theta,
                prior_sigma2_u=prior_sigma2_u,
                prior_sigma2=prior_sigma2,
                sharedlambda=sharedlambda,
                DLM=DLM,
                lagOrder=lagOrder,
                ## MH tuning
                stepsize_logrho=stepsize_logrho,
                stepsize_loglambda_theta=stepsize_loglambda_theta,
                Vgridsearch=Vgridsearch,
                gridsize=gridsize,
                rfbtries=rfbtries,
                wls_steps=wls_steps,
                MHwls=MHwls,
                stepsize_theta=stepsize_theta)

  ## indices
  const$n <- nrow(X)
  const$K <- ncol(Y)
  if(is.null(maxClusters)){
    const$C <- const$K ## default number of clusters
  }else{
    const$C <- maxClusters
  }
  const$L <- ncol(X)
  const$i_index <- rep(1:const$n,const$K) ## 1....1,2...2,....
  const$k_index <- rep(1:const$K,each=const$n)      ## 1:K,1:K,...

  ##
  if(DLM==TRUE){ ## create penalty matrix
    D <- diag(const$L)
    for(jj in 1:(lagOrder)){
      D <- diff(D)
    }
    const$PEN <- t(D)%*%D
  }else{ ## no penalties
    const$PEN <- matrix(0,ncol=const$L,nrow=const$L)
  }


  ## basis functions
  temptheta <- rep(1,const$L)/sqrt(const$L) ## just a standard theta value in order to compute basis functions
  Xtheta <- matrix(sapply(1:const$K,function(kk){X%*%temptheta}))
  const$SS <- mgcv::smoothCon(s(Xtheta,bs="bs"),data=data.frame(Xtheta),absorb.cons = FALSE)[[1]]
  const$SSderiv <- const$SS ## make a smooth object for computing first order derivatives
  const$SSderiv$deriv <- 1 ## first order derivatives
  Btheta <- const$SS$X
  const$d <- ncol(Btheta)
  const$mu0 <- rep(0,const$d) ## probably no need to change this from 0
  const$invSig0 <- const$SS$S[[1]]
 ## accounting for the unpenalized columns (not invertible)
  unpenalized_params <- which(rowSums(const$invSig0)==0)
  diag(const$invSig0)[unpenalized_params] <- 1 ## currently using N(0,1) prior on unpenalized components

  ## grid for
  const$grid <- (1:(const$gridsize-1))/const$gridsize ## exclude 0s and 1s

  return(const)
}


## get starting parameter values
get_starting_vals <- function(const){
  params <- list()

  ## cluster memberships (completely random)
  # params$Zbeta <-  sample(const$C,const$K,replace=TRUE)
  # params$Ztheta <- sample(const$C,const$K,replace=TRUE)

  ## other params from prior
  params$alpha <- c(rgamma(1,shape=const$prior_alpha_beta[1],rate=const$prior_alpha_beta[2]),
                    rgamma(1,shape=const$prior_alpha_theta[1],rate=const$prior_alpha_theta[2]) )

  params$Vbeta <- c(rbeta(const$C-1,1,params$alpha[1]),1) ## final value is 1 (truncated)
  params$Vtheta <- c(rbeta(const$C-1,1,params$alpha[2]),1) ## final value is 1 (truncated)

  params$logrho <- log(rgamma(1,shape=const$prior_rho[1],rate=const$prior_rho[2]))

  params <- compute_pi(params,const) ## update weights pimat and pimatstar

  ## draw Zbeta and Ztheta according to pistarmat
  for(kk in 1:const$K){

        ## sample 1 of C^2 with correct probabilities
        ab <- sample(1:(const$C^2),1,prob=c(params$pistarmat))
        ## check which a and b this corresponds to
        Zk <- which(matrix(1:(const$C^2),ncol=const$C)==ab,arr.ind = TRUE) ## which element of CxC matrix
        params$Zbeta[kk] <- Zk[1]  ## a=corresponding row
        params$Ztheta[kk] <- Zk[2] ## b=corresponding column

  }


  if(const$sharedlambda==TRUE){
    params$lambda_beta <- rgamma(1,shape=1,rate=1)
  }else{
    params$lambda_beta <- rgamma(const$C,shape=1,rate=1)
  }


  if(const$DLM==TRUE){
    params$loglambda_theta <- log(rgamma(1,shape=1,rate=1))
  }else{
    params$loglambda_theta <- -Inf
  }

  # params$sigma2 <- 1/rgamma(1,shape=const$prior_sigma2[1],rate=const$prior_sigma2[2])
  params$sigma2 <- 1/rgamma(1,shape=10,rate=10)

  # params$sigma2_u <- 1/rgamma(1,shape=const$prior_sigma2_u[1],rate=const$prior_sigma2_u[2])
  params$sigma2_u <- 1/rgamma(1,shape=10,rate=10)

  params$u <- rnorm(n,0,sqrt(params$sigma2_u))

  ## update thetastar
  params$thetastar <- c(t(rFisherBingham(const$C,mu = const$prior_tau_theta*rep(1,const$L), Aplus = 0)))##rep(rep(1,const$L)/sqrt(const$L),const$C)# ## 0 for vMF
  params$gammastar <-  params$thetastar
  Xtheta <- matrix(sapply(1:const$K,function(kk){const$X%*%params$thetastar[(params$Zbeta[kk]-1)*const$L+(1:const$L)]}))
  params$Btheta <- get_Btheta(Xtheta,const)
  params$DerivBtheta <- get_DerivBtheta(Xtheta,const)
  params$theta <- rep(NA,const$L*const$K)
  params <- assign_thetas(params, const)

  # testing
  # params$Zbeta <- params$Ztheta <- rep(1,length(params$Zbeta))
  # params$thetastar <- c(t(rFisherBingham(const$C,mu = const$prior_tau_theta*runif(const$L,-1,1), Aplus = 0)))##rep(rep(1,const$L)/sqrt(const$L),const$C)# ## 0 for vMF
  # params$gammastar <-  params$thetastar
  # Xtheta <- matrix(sapply(1:const$K,function(kk){const$X%*%params$thetastar[(params$Zbeta[kk]-1)*const$L+(1:const$L)]}))
  # params$Btheta <- get_Btheta(Xtheta,const)
  # params$DerivBtheta <- get_DerivBtheta(Xtheta,const)
  # params$theta <- rep(NA,const$L*const$K)
  # params <- assign_thetas(params, const)
  # params$thetastar <- rep(1/sqrt(const$L),length(params$thetastar))
  # params$theta <- rep(1/sqrt(const$L),length(params$theta))
  # Xtheta <- matrix(sapply(1:const$K,function(kk){const$X%*%params$thetastar[(params$Zbeta[kk]-1)*const$L+(1:const$L)]}))
  # params$Btheta <- get_Btheta(Xtheta,const)
  # params$DerivBtheta <- get_DerivBtheta(Xtheta,const)
  # params$theta <- rep(NA,const$L*const$K)
  # params <- assign_thetas(params, const)
  params$lambda_beta <- 1#0.0001
  params$sigma2 <- 1#sqrt(0.01)
  params$sigma2_u <- 0.000000000001
  params$u <- rnorm(n,0,sqrt(params$sigma2_u))
  ## end testing

  # print(params$theta[1:4])

  ## betastar
  # params$betastar <- c(t(rmvnorm(const$C,const$mu0,solve(const$invSig0))))
  # params$betastar <- c(t(rmvnorm(const$C,const$mu0,MASS::ginv(const$invSig0))))
  params$betastar <- c(t(rmvnorm(const$C,const$mu0,diag(const$d))))
  params$beta <- rep(NA,const$d*const$K)
  params <- assign_betas(params, const)
  # params$betastar <- truebetastar
  # params$beta <- truebeta

  # Sig0 <- const$invSig0 ## accounting for the unpenalized columns (not invertible)
  # unpenalized_params <- which(rowSums(Sig0)==0)
  # diag(Sig0)[unpenalized_params] <- 1
  # Sig0 <-  solve(Sig0)
  # params$betastar <- c(t(rmvnorm(const$C,const$mu0,Sig0)))


  return(params)
}
