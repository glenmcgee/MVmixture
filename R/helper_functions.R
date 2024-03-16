## Notes:
### For betastar: currently using N(0,1) prior on unpenalized components
##### Fix this

## function to compute pi
compute_pi <- function(params,const){

  ## compute "marginals" under independence
  pi_beta <- params$Vbeta*c(1,cumprod(1-params$Vbeta)[-const$C])
  pi_theta <- params$Vtheta*c(1,cumprod(1-params$Vtheta)[-const$C])

  ## joint probs
  pimat <- (pi_beta)%*%t(pi_theta) + diag(exp(params$logrho)*(pi_beta*pi_theta))
  # pimat <- (pi_beta)%*%t(pi_theta) ## pi^I (under independence)
  # diag(pimat) <- (1+exp(params$logrho))*diag(pimat) ## multiply diagonal by 1+rho

  ## standardize
  pistarmat <- pimat/sum(pimat)

  params$pimat <- pimat
  params$pistarmat <- pistarmat
  return(params)
}


## get basis functions
get_Btheta <- function(Xtheta,const){
  if(ncol(as.matrix(Xtheta))>1){
    return(sapply(1:ncol(Xtheta),function(jj){mgcv::PredictMat(const$SS,data=data.frame(Xtheta[,jj]))}))
  }else{
    return(mgcv::PredictMat(const$SS,data=data.frame(Xtheta)))
  }
}

## get derivatives of Basis functions
get_DerivBtheta <- function(Xtheta,const){
  if(ncol(as.matrix(Xtheta))>1){
    return(sapply(1:ncol(Xtheta),function(jj){mgcv::PredictMat(const$SSderiv,data=data.frame(Xtheta[,jj]))}))
  }else{
    return(mgcv::PredictMat(const$SSderiv,data=data.frame(Xtheta)))
  }
}


## compute B times the corresponding beta for single cluster
get_B_beta_k <- function(params,const,kk){

  return(sapply(1:const$p,function(jj){
    get_Btheta(const$X[[jj]]%*%params$thetastar[(params$Ztheta[kk,jj]-1)*const$L+(1:const$L)],const)%*%params$betastar[(params$Zbeta[kk,jj]-1)*const$d+(1:const$d)]
  }))

}

## compute B times the corresponding beta for all obs (d)
get_B_beta <- function(params,const){

  return(Reduce("rbind",lapply(1:const$K,function(kk){
    get_B_beta_k(params,const,kk)
  })))

}


## compute DerivB times the corresponding beta for given k,j
get_DerivB_beta <- function(params,const,kk,jj){

  return( get_DerivBtheta(const$X[[jj]]%*%params$thetastar[(params$Ztheta[kk,jj]-1)*const$L+(1:const$L)],const)%*%params$betastar[(params$Zbeta[kk,jj]-1)*const$d+(1:const$d)] )

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
  return((alpha-1)*log(1-vv) + sum(sapply(1:const$p,function(jj){sapply(1:const$K,function(kk){log(params$pistarmat[params$Zbeta[kk,jj],params$Ztheta[kk,jj]])})})))
}


## get components for thetastar sampling
get_XTytilde <- function(cc,whichk,whichkj,params,const){

  Xtilde_k <- lapply(whichk,function(kk){
    Reduce("+",
      lapply(whichkj[[kk]],function(jj){
        c(get_DerivB_beta(params,const,kk,jj))*const$X[[jj]]
        })
      )
    })

  ytilde_k <- lapply(whichk,function(kk){
    const$y[const$k_index==kk]-params$b0[kk]-apply(get_B_beta_k(params,const,kk),1,sum)-params$u+
      Reduce("+",
             lapply(whichkj[[kk]],function(jj){
               (c(get_DerivB_beta(params,const,kk,jj))*const$X[[jj]])%*%params$thetastar[(cc-1)*const$L+(1:const$L)]
             })
       )
  })

  return(list(XTX=Reduce("+",lapply(Xtilde_k,function(XX){t(XX)%*%XX})),
              XTy=Reduce("+",lapply(1:length(Xtilde_k),function(kk){t(Xtilde_k[[kk]])%*%ytilde_k[[kk]]}))))
}

# ## obtain (standardized) WLS estimate
# get_WLS <- function(cc,whichk,params,const){
#   newparams <- params
#
#   for(ss in 1:const$wls_steps){
#     oldparams <- newparams
#
#     ## get W and Wytilde
#     Wcomps <- get_W(cc,whichk,oldparams,const)
#     Wytilde <- Wcomps$Wytilde
#     W <- Wcomps$W
#
#     # wls <- solve(t(const$X)%*%W%*%const$X,t(const$X)%*%c(Wytilde))
#     wls <- params$sigma2*solve(t(const$X)%*%W%*%const$X+exp(params$loglambda_theta)*const$PEN)%*%(const$prior_tau_theta*as.matrix(rep(1,const$L)) + (1/params$sigma2)*t(const$X)%*%c(Wytilde))
#
#     newparams$thetastar[(cc-1)*const$L+(1:const$L)] <- c(wls/sqrt(sum(wls^2))) ## standardize
#     for(kk in whichk){
#       newparams$Btheta[const$k_index==kk,] <- get_Btheta(const$X%*%newparams$thetastar[(cc-1)*const$L+(1:const$L)],const)
#       newparams$DerivBtheta[const$k_index==kk,] <- get_DerivBtheta(const$X%*%newparams$thetastar[(cc-1)*const$L+(1:const$L)],const)
#     }
#
#   }
#
#   return(newparams)
#
# }



#
# assign_betas <- function(params,const){
#   for(kk in 1:const$K){
#     params$beta[(kk-1)*const$d+(1:const$d)] <- params$betastar[(params$Zbeta[kk]-1)*const$d+(1:const$d)]
#   }
#   return(params)
# }
#
# assign_thetas <- function(params,const){
#   for(kk in 1:const$K){
#     params$theta[(kk-1)*const$L+(1:const$L)] <- params$thetastar[(params$Ztheta[kk]-1)*const$L+(1:const$L)]
#   }
#   return(params)
# }


# assign_betas <- function(params,const){
#   for(kk in 1:const$K){
#     params$beta[(kk-1)*const$d+(1:const$d)] <- params$betastar[(params$Zbeta[kk]-1)*const$d+(1:const$d)]
#   }
#   return(params)
# }

## compute thetastar from omegastar
get_theta <- function(omegastar){
  return(c(sin(omegastar),1) * cumprod(c(1,cos(omegastar))))
}

## compute omegastar from thetastar
get_omega <- function(thetastar){

  omegastar <- c(asin(thetastar[1]))
  for(jj in 2:(length(thetastar)-1)){
    omegastar <- c(omegastar,asin(thetastar[jj]/(prod(cos(omegastar)))))
  }
  return(omegastar)
}



assign_betas <- function(obj){

  betas <- lapply(1:obj$const$p,function(jj){
    Reduce("cbind",lapply(1:obj$const$K,function(kk){
      Reduce("rbind",lapply(1:nrow(obj$Zbeta),function(rr){
        return(obj$betastar[rr, (obj$Zbeta[rr,obj$const$K*(jj-1)+kk]-1)*obj$const$d+(1:obj$const$d)])
      }))
    }))
  })

  return(betas)
}

assign_thetas <- function(obj){

  thetas <- lapply(1:obj$const$p,function(jj){
    Reduce("cbind",lapply(1:obj$const$K,function(kk){
      Reduce("rbind",lapply(1:nrow(obj$Ztheta),function(rr){
        return(obj$thetastar[rr, (obj$Ztheta[rr,obj$const$K*(jj-1)+kk]-1)*obj$const$L+(1:obj$const$L)])
      }))
    }))
  })

  return(thetas)
}



summarize_clusters <- function(obj){

  Zbeta <- round(100*t(apply(obj$Zbeta,2,function(x) table(factor(x,levels=1:obj$const$C))))/nrow(obj$Zbeta))
  Ztheta <- round(100*t(apply(obj$Ztheta,2,function(x) table(factor(x,levels=1:obj$const$C))))/nrow(obj$Ztheta))
  outcome <- rep(1:obj$const$K,obj$const$p)
  exposure <- rep(1:obj$const$p,each=obj$const$K)
  summ <- data.frame(outcome,exposure,Zbeta=Zbeta,Ztheta=Ztheta)
  return(summ[order(summ$outcome),])
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
                thetagridsize=thetagridsize)

  ## indices
  if(is.list(X)==FALSE){
    const$X <- list(X)
  }
  const$p <- length(const$X)
  const$n <- nrow(const$X[[1]])
  const$K <- ncol(Y)
  const$L <- ncol(const$X[[1]])
  if(is.null(maxClusters)){
    const$C <- const$K ## default number of clusters
  }else{
    const$C <- maxClusters
  }
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
  ### Uses same basis functions for all exposures.
  temptheta <- rep(1,const$L)/sqrt(const$L) ## just a standard theta value in order to compute basis functions
  ## combining all exposures to get basis functions
  Xtheta <- c(sapply(1:const$p,function(jj){const$X[[jj]]%*%temptheta}))
  const$SS <- mgcv::smoothCon(s(Xtheta,bs="bs"),data=data.frame(Xtheta),absorb.cons = TRUE)[[1]] ## should be true for multiple exposures
  const$SSderiv <- const$SS ## make a smooth object for computing first order derivatives
  const$SSderiv$deriv <- 1 ## first order derivatives
  const$d <- ncol(const$SS$X)
  ## now a list of n x d matrices B
  Btheta <- lapply(1:const$p,function(jj){mgcv::PredictMat(const$SS,data=data.frame(Xtheta=const$X[[jj]]%*%temptheta))})
  const$mu0 <- rep(0,const$d) ## probably no need to change this from 0
  const$invSig0 <- const$SS$S[[1]]
  # Xtheta <- matrix(sapply(1:const$K,function(kk){X%*%temptheta}))
  # const$SS <- mgcv::smoothCon(s(Xtheta,bs="bs"),data=data.frame(Xtheta),absorb.cons = FALSE)[[1]]
  # const$SSderiv <- const$SS ## make a smooth object for computing first order derivatives
  # const$SSderiv$deriv <- 1 ## first order derivatives
  # Btheta <- const$SS$X
  # const$d <- ncol(Btheta)
  # const$mu0 <- rep(0,const$d) ## probably no need to change this from 0
  # const$invSig0 <- const$SS$S[[1]]

 ## accounting for the unpenalized columns (not invertible)
  # unpenalized_params <- which(rowSums(const$invSig0)==0)
  # diag(const$invSig0)[unpenalized_params] <- 1 ## currently using N(0,1) prior on unpenalized components

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
  params$Zbeta <- params$Ztheta <- matrix(NA,nrow=const$K,ncol=const$p)
  for(kk in 1:const$K){
    for(jj in 1:const$p){

        ## sample 1 of C^2 with correct probabilities
        ab <- sample(1:(const$C^2),1,prob=c(params$pistarmat))
        ## check which a and b this corresponds to
        Zk <- which(matrix(1:(const$C^2),ncol=const$C)==ab,arr.ind = TRUE) ## which element of CxC matrix
        params$Zbeta[kk,jj] <- Zk[1]  ## a=corresponding row
        params$Ztheta[kk,jj] <- Zk[2] ## b=corresponding column
    }
  }

  # # ## TESTING
  # params$Zbeta <- matrix(c(rep(c(1,1),each=2),rep(c(3,3),each=2)),ncol=2)
  # params$Ztheta <- matrix(c(rep(c(1,2),each=2),rep(c(3,4),each=2)),ncol=2)


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


  ## Get Xtheta and basis functions
  Xtheta <- sapply(1:const$p,function(jj){matrix(sapply(1:const$K,function(kk){const$X[[jj]]%*%params$thetastar[(params$Zbeta[kk,jj]-1)*const$L+(1:const$L)]}))})
  params$Btheta <- lapply(1:const$p,function(jj){get_Btheta(Xtheta[,jj],const)})
  params$DerivBtheta <- lapply(1:const$p,function(jj){get_DerivBtheta(Xtheta[,jj],const)})
  params$theta <- matrix(NA,ncol=const$L*const$p,nrow=const$K)
  # params <- assign_thetas(params, const)

  ## betastar
  params$betastar <- c(t(rmvnorm(const$C,const$mu0,diag(const$d))))
  params$beta <- matrix(NA,ncol=const$d*const$p,nrow=const$K)
  # params <- assign_betas(params, const)

  # intercept
  params$b0 <- rnorm(const$K)


  return(params)
}

