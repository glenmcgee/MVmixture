### MCMC for MVmixture


library(parallel)
library(mvtnorm)
library(simdd)
library(mgcv)
library(MASS)
setwd("~/GitHub/MVmixture/R")
source("helper_functions.R")
source("fb.saddle.R")
source("updates.R")
source("build_sampler.R")
source("predict_MVmix.R")
##
MVmix <- function(Y, ## n x K matrix of responses
                  X, ## n x L matrix of exposures
                  Z, ## confounders to be adjusted
                  ## MCMC specs
                  niter=10000, ## number of iterations
                  nburn=0.5*niter, ## burn-in fraction
                  nthin=10, ## thinning number
                  nchains=1, ## number of chains ## not yet implemented
                  ncores=1, ## number of cores for mclapply (set to 1 for non-parallel) ## only used if nchains>1
                  ## prior hyperparameters
                  maxClusters=NULL,
                  prior_alpha_beta=c(1,1),
                  prior_alpha_theta=c(1,1),
                  prior_rho=c(1,1),
                  prior_tau_theta=1,
                  prior_lambda_beta=c(1,1),
                  prior_lambda_theta=c(1,1),
                  prior_sigma2_u=c(0.01,0.01),
                  prior_sigma2=c(0.01,0.01),
                  sharedlambda=TRUE,
                  DLM=FALSE,
                  lagOrder=2, ## min 1. order for lag penalty (only if DLM=TRUE)
                  ## MH tuning
                  stepsize_logrho=1, ## sd for random walk
                  stepsize_loglambda_theta=1, ## sd for random walk
                  Vgridsearch=TRUE, ## use grid search for approximate sampling of V_c
                  gridsize=100, ## size of grid. Not used it Vgridsearch==FALSE
                  rfbtries=1000, ## mtop for rFisherBingham (default 1000)
                  wls_steps=4, ## number of steps in WLS approximation
                  MHwls=FALSE, ## center theta proposals around WLS
                  stepsize_theta=1){ ## k parameter for theta proposals

  ## set up constants
  const <- initialize_const(Y, ## response
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
                            stepsize_theta)


  ## set up MCMC sampler with options
  update_params <- build_sampler(const)

  ## mcmc function (called by mclapply for multiple chains)
  run_sampler <- function(ind){

    ## get starting values
    params_ss <- get_starting_vals(const)


    ## initialize results storage
    nkeep <- floor((niter-nburn)/nthin)# round

    keep_Zbeta <- matrix(0,ncol=const$K,nrow=nkeep)
    keep_Ztheta <- matrix(0,ncol=const$K,nrow=nkeep)
    keep_Vbeta <- matrix(0,ncol=const$C,nrow=nkeep)
    keep_Vtheta <- matrix(0,ncol=const$C,nrow=nkeep)
    keep_betastar <- matrix(0,ncol=const$d*const$C,nrow=nkeep)
    keep_beta <- matrix(0,ncol=const$d*const$K,nrow=nkeep)
    keep_gammastar <- matrix(0,ncol=const$L*const$C,nrow=nkeep)
    keep_thetastar <- matrix(0,ncol=const$L*const$C,nrow=nkeep)
    keep_theta <- matrix(0,ncol=const$L*const$K,nrow=nkeep)
    keep_alpha <- matrix(0,ncol=2,nrow=nkeep)
    keep_logrho <- matrix(0,ncol=1,nrow=nkeep)
    keep_lambda_beta <- matrix(0,ncol=const$C,nrow=nkeep)
    keep_loglambda_theta <- matrix(0,ncol=1,nrow=nkeep)
    keep_u <- matrix(0,ncol=const$n,nrow=nkeep)
    keep_sigma2_u <- matrix(0,ncol=1,nrow=nkeep)
    keep_sigma2 <- matrix(0,ncol=1,nrow=nkeep)


    ## MCMC
    for(ss in 1:niter){

      ## update parameters
      params_ss <- update_params(params_ss)

      ## retain samples after burn-in and thinning
      if(ss>nburn & ss%%nthin==0){
        skeep <- (ss-nburn)/nthin

        keep_Zbeta[skeep,] <- params_ss$Zbeta
        keep_Ztheta[skeep,] <- params_ss$Ztheta
        keep_Vbeta[skeep,] <- params_ss$Vbeta
        keep_Vtheta[skeep,] <- params_ss$Vtheta
        keep_betastar[skeep,] <- params_ss$betastar
        keep_beta[skeep,] <- params_ss$beta
        keep_gammastar[skeep,] <- params_ss$gammastar
        keep_thetastar[skeep,] <- params_ss$thetastar
        keep_theta[skeep,] <- params_ss$theta
        keep_alpha[skeep,] <- params_ss$alpha
        keep_logrho[skeep,] <- params_ss$logrho
        keep_lambda_beta[skeep,] <- params_ss$lambda_beta
        keep_loglambda_theta[skeep,] <- params_ss$loglambda_theta
        keep_u[skeep,] <- params_ss$u
        keep_sigma2_u[skeep,] <- params_ss$sigma2_u
        keep_sigma2[skeep,] <- params_ss$sigma2

      }
    }

    return(list(Zbeta=keep_Zbeta,
                Ztheta=keep_Ztheta,
                Vbeta=keep_Vbeta,
                Vtheta=keep_Vtheta,
                betastar=keep_betastar,
                beta=keep_beta,
                gammastar=keep_gammastar,
                thetastar=keep_thetastar,
                theta=keep_theta,
                alpha=keep_alpha,
                logrho=keep_logrho,
                lambda_beta=keep_lambda_beta,
                loglambda_theta=keep_loglambda_theta,
                u=keep_u,
                sigma2_u=keep_sigma2_u,
                sigma2=keep_sigma2,
                const=const) )
  }


  ## run chains
  if(nchains==1){ ## run single chain
    samples <- run_sampler()
    # samples$PSR <- NULL
  }else{ ## run multiple chains
    if(ncores==1){ ## run multiple chains (not parallel)
      chains <- vector(mode = "list", length = nchains)
      for(cc in 1:nchains){
        chains[[cc]] <- run_sampler()
      }
    }else{ ## run in parallel with nchains>1 and ncores>1
      chains <- mclapply(1:nchains,run_sampler,mc.cosamples=ncores)
    }

    ## append chains
    samples <- chains[[1]] ## use first for names/formatting
    samples <- lapply(names(samples), ## looping over names of variables to append
                  function(v){ Reduce('rbind',lapply(chains,function(chn) chn[[v]]))} )## for each variable name, extract element from each chain, then collapse them via rbind


    # ## compute PSR
    # samples$PSR <- vector(mode = "list", length = length(names(samples)))
    # samples$PSR <- lapply(names(samples),function(v){
    #   compute_PSR(lapply(chains,function(chn){
    #     chn[[v]]
    #   }))
    # })
    # names(samples$PSR) <- names(samples)
  }

  return(samples)
}
