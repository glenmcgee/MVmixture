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

## get W and Wytilde
get_W <- function(cc,whichk,params,const){
  Wytilde <- c(Reduce('+',## summing over the k vectors
                      lapply(whichk, ## only k in cluster cc
                             function(kk){ ## Wk*ytildek
                               ((params$DerivBtheta[const$k_index==kk,]%*%params$betastar[(cc-1)*const$d+(1:const$d)])^2)* ## Wk
                                 (const$X%*%params$thetastar[(cc-1)*const$L+(1:const$L)]+((const$y[const$k_index==kk]-params$u)-params$Btheta[const$k_index==kk,]%*%params$betastar[(cc-1)*const$d+(1:const$d)])/(params$DerivBtheta[const$k_index==kk,]%*%params$betastar[(cc-1)*const$d+(1:const$d)]))
                             })))
  W <- diag(c(## vectorizing then turning into single diagonal matrix at the end
    Reduce('+',## summing over the k
           lapply(whichk, ## only k in cluster cc
                  function(kk){ ## Wk
                    ((params$DerivBtheta[const$k_index==kk,]%*%params$betastar[(cc-1)*const$d+(1:const$d)])^2)
                  }))))
  return(list(Wytilde=Wytilde,
              W=W))
}

## obtain (standardized) WLS estimate
get_WLS <- function(cc,whichk,params,const){
  newparams <- params

  for(ss in 1:const$wls_steps){
    oldparams <- newparams

    ## get W and Wytilde
    Wcomps <- get_W(cc,whichk,oldparams,const)
    Wytilde <- Wcomps$Wytilde
    W <- Wcomps$W

    # wls <- solve(t(const$X)%*%W%*%const$X,t(const$X)%*%c(Wytilde))
    wls <- params$sigma2*solve(t(const$X)%*%W%*%const$X+exp(params$loglambda_theta)*const$PEN)%*%(const$prior_tau_theta*as.matrix(rep(1,const$L)) + (1/params$sigma2)*t(const$X)%*%c(Wytilde))

    newparams$thetastar[(cc-1)*const$L+(1:const$L)] <- c(wls/sqrt(sum(wls^2))) ## standardize
    for(kk in whichk){
      newparams$Btheta[const$k_index==kk,] <- get_Btheta(const$X%*%newparams$thetastar[(cc-1)*const$L+(1:const$L)],const)
      newparams$DerivBtheta[const$k_index==kk,] <- get_DerivBtheta(const$X%*%newparams$thetastar[(cc-1)*const$L+(1:const$L)],const)
    }

  }

  return(newparams)

}




assign_betas <- function(params,const){
  for(kk in 1:const$K){
    params$beta[(kk-1)*const$d+(1:const$d)] <- params$betastar[(params$Zbeta[kk]-1)*const$d+(1:const$d)]
  }
  return(params)
}


## compute thetastar from omegastar
get_theta <- function(omegastar){
  sin_w <- c(sin(omegastar),1)
  cos_w <- c(1,cos(omegastar))
  return(sin_w*cumprod(cos_w))
}

## compute omegastar from thetastar
get_omega <- function(thetastar){
  omegastar <- c()
  omegastar <- c(omegastar,asin(thetastar[1]))
  for(jj in 2:(length(thetastar)-1)){
    omegastar <- c(omegastar,asin(thetastar[jj]/(prod(cos(omegastar)))))
  }
  return(omegastar)
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
                             prior_omega_a,
                             sharedlambda,
                             DLM,
                             lagOrder,
                             ## MH tuning
                             stepsize_logrho,
                             stepsize_loglambda_theta,
                             Vgridsearch,
                             gridsize,
                             rfbtries,
                             thetaMethod,
                             wls_steps,
                             MHwls,
                             stepsize_theta,
                             thetagridsize){

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
                prior_omega_a=prior_omega_a,
                sharedlambda=sharedlambda,
                DLM=DLM,
                lagOrder=lagOrder,
                ## MH tuning
                stepsize_logrho=stepsize_logrho,
                stepsize_loglambda_theta=stepsize_loglambda_theta,
                Vgridsearch=Vgridsearch,
                gridsize=gridsize,
                rfbtries=rfbtries,
                thetaMethod=thetaMethod,
                wls_steps=wls_steps,
                MHwls=MHwls,
                stepsize_theta=stepsize_theta,
                thetagridsize=thetagridsize)

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
  ## grid for theta
  const$thetagrid <- (1:(const$thetagridsize-1))/const$thetagridsize ## exclude 0s and 1s


  return(const)
}


## get starting parameter values
get_starting_vals <- function(const){
  params <- list()

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

  ## initialize thetastar
  if(const$thetaMethod=="MH_beta"){ ## use omegastar parameterization
    params$omegastar <- rbeta((const$L-1)*const$C,1,1)
    params$thetastar <- sapply(1:const$C,function(cc){
        return(get_theta(params$omegastar[(cc-1)*(const$L-1)+1:(const$L-1)]))
      })
  }else{ ## otherwise just draw from fb
    params$omegastar <- NULL
    params$thetastar <- c(t(rFisherBingham(const$C,mu = const$prior_tau_theta*rep(1,const$L), Aplus = 0)))##rep(rep(1,const$L)/sqrt(const$L),const$C)# ## 0 for vMF
  }

  params$gammastar <-  params$thetastar ## only used for MH_mvn
  Xtheta <- matrix(sapply(1:const$K,function(kk){const$X%*%params$thetastar[(params$Zbeta[kk]-1)*const$L+(1:const$L)]}))
  params$Btheta <- get_Btheta(Xtheta,const)
  params$DerivBtheta <- get_DerivBtheta(Xtheta,const)
  params$theta <- rep(NA,const$L*const$K)
  params <- assign_thetas(params, const)


  ## betastar
  # params$betastar <- c(t(rmvnorm(const$C,const$mu0,solve(const$invSig0))))
  # params$betastar <- c(t(rmvnorm(const$C,const$mu0,MASS::ginv(const$invSig0))))
  params$betastar <- c(t(rmvnorm(const$C,const$mu0,diag(const$d))))
  params$beta <- rep(NA,const$d*const$K)
  params <- assign_betas(params, const)

  # Sig0 <- const$invSig0 ## accounting for the unpenalized columns (not invertible)
  # unpenalized_params <- which(rowSums(Sig0)==0)
  # diag(Sig0)[unpenalized_params] <- 1
  # Sig0 <-  solve(Sig0)
  # params$betastar <- c(t(rmvnorm(const$C,const$mu0,Sig0)))


  return(params)
}


