


# require(parallel)
# require(mvtnorm)
# require(simdd)
# require(mgcv)
# require(MASS)
# require(reshape)
# require(randomForest)
# require(tidyverse)
# setwd("~/GitHub/MVmixture/R")
# source("helper_functions.R")
# source("fb.saddle.R")
# source("updates.R")
# source("build_sampler.R")
# source("predict_MVmix.R")
# source("pairwise_clusters.R")
# source("ExposureImportance.R")
##

#' Fit multivariate mixtures model with adaptive clustering
#'
#' This function will take in the observed data and fit a possibly multivariate
#' additive mixtures model with or without adaptive clustering over outcomes/exposures.
#'
#' @param Y   n x K matrix of responses
#' @param X  p-list of n x L matrix of exposures
#' @param Z confounders to be adjusted
#' @param niter number of iterations
#' @param nburn  burn-in fraction
#' @param nthin  thinning number
#' @param nchains  number of chains ## not yet implemented
#' @param ncores  number of cores for mclapply (set to 1 for non-parallel) ## only used if nchains>1
#' @param clustering  one of "both","theta","beta" or ,"neither";
#' @param fixedZbeta  optional KxP matrix of fixed cluster IDs for beta; only used if clustering="theta" or "neither"
#' @param fixedZtheta  optional KxP matrix of fixed cluster IDs for theta; only used if clustering="theta" or "neither"
#' @param maxClusters max number of distinct clusters
#' @param prior_alpha_beta  gamma prior hyperparameters for alpha_beta
#' @param prior_alpha_theta  gamma prior hyperparameters for alpha_theta
#' @param prior_rho  gamma prior hyperparameters for rho
#' @param prior_tau_theta Do not change. prior hyperparameter for thetastar direction
#' @param prior_lambda_beta gamma prior hyperparameters for lambda_beta controlling smoothness of exposure response functions
#' @param prior_lambda_theta gamma prior hyperparameters for lambda_theta controlling smoothness of weight functions (if DLMpenalty=TRUE)
#' @param prior_xi inverse gamma prior hyperparameters for xi
#' @param prior_sigma2 inverse gamma prior hyperparameters for sigma^2
#' @param prop_phi_a  hyperparameter a for the beta(a,b) proposal on phistar. Higher value means smaller steps.
#' @param sharedlambda Should a single lambda value controlling smoothness be shared across all exposure-response functions? T/F
#' @param DLM Use B-spline approximation to impose smoothness in weights over time
#' @param DLMpenalty  Include smoothness penalty in weights over time; only if DLM=TRUE
#' @param lagOrder Number of basis functions for weights (if DLM=TRUE); NULL indicates no dimension reduction
#' @param diff  Degree of difference penalty matrix (if DLMpenalty=TRUE)
#' @param MIM Fit MIM version with unknown index structure (T/F)
#' @param MIMorder Maximum order of MIM (i.e. number of indices); ignored if MIM=FALSE
#' @param LM Force linear exposure response relationships (T/F)
#' @param betaOrder B-spline basis dimension for exposure response functions. default -1 allows mgcv to choose automatically
#' @param stepsize_logrho SD for random walk updates on log(rho)
#' @param stepsize_loglambda_theta SD for random walk updates on log(lambda_theta)
#' @param stepsize_logxi SD for random walk updates on log(xi)
#' @param stepsize_logsigma2 SD for random walk updates on log(sigma2)
#' @param Vgridsearch Use grid search for approximate sampling of V_c
#' @param gridsize Size of grid. Not used if Vgridsearch=FALSE
#' @param rfbtries Max number of tries for drawing from Fisher-Bingham; equiv.to "mtop" for rFisherBingham (default 1000)
#' @param approx Should rFB approximation be used? Otherwise use polar transformation and MH updates
#' @param appendchain Option to use last values of a previous fit as starting values for a new chain.
#'
#' @return Posterior samples for all parameters.
#'
#' @export
MVmix <- function(Y, ## n x K matrix of responses
                  X, ## p-list of n x L matrix of exposures
                  Z, ## confounders to be adjusted
                  ## MCMC specs
                  niter=10000, ## number of iterations
                  nburn=0.5*niter, ## burn-in fraction
                  nthin=10, ## thinning number
                  nchains=1, ## number of chains ## not yet implemented
                  ncores=1, ## number of cores for mclapply (set to 1 for non-parallel) ## only used if nchains>1
                  ## prior hyperparameters
                  clustering="both", ## one of "both","theta","beta" or ,"neither";
                  fixedZbeta=NULL,  ## optional KxP matrix of fixed cluster IDs for beta; only used if clustering="theta" or "neither"
                  fixedZtheta=NULL, ## optional KxP matrix of fixed cluster IDs for theta; only used if clustering="theta" or "neither"
                  maxClusters=NULL,
                  prior_alpha_beta=c(1,1),
                  prior_alpha_theta=c(1,1),
                  prior_rho=c(1,1),
                  prior_tau_theta=0,# do not change#1,
                  prior_lambda_beta=c(1,1),
                  prior_lambda_theta=c(1,0.001),
                  prior_xi=c(0.01,0.01),
                  prior_sigma2=c(0.01,0.01),
                  prop_phi_a=200, ## hyperparameter a for the beta(a,b) proposal on phistar. higher value means smaller steps
                  sharedlambda=TRUE,
                  DLM=FALSE, ## use b-spline approximation to impose smoothness over time
                  DLMpenalty=FALSE, ## include smoothness penalty over time, only if DLM=TRUE
                  lagOrder=6, ## no. of bases for omega weight function (if DLM=TRUE); NULL indicates no dimension reduction
                  diff=2, ## degree of difference penalty matrix (if DLM=TRUE)
                  MIM=FALSE, ## fit MIM version
                  MIMorder=4, ## maximum order of MIM (ignored if MIM=FALSE)
                  LM=FALSE, ## force linear effects
                  betaOrder=-1, ## b spline basis dimension for exposure response functions. default -1 allows mgcv to choose automatically
                  ## MH tuning
                  stepsize_logrho=1, ## sd for random walk
                  stepsize_loglambda_theta=1, ## sd for random walk
                  stepsize_logxi=1,
                  stepsize_logsigma2=1,
                  Vgridsearch=TRUE, ## use grid search for approximate sampling of V_c
                  gridsize=10, ## size of grid. Not used it Vgridsearch==FALSE
                  rfbtries=1000, ## mtop for rFisherBingham (default 1000)
                  approx=TRUE,  ## TRUE=MVN/rFB sampling. FALSE=MH_Beta sampling
                  appendchain=NULL){ ## set to previous_fit to keep running chain

  ## set up constants
  if(!is.null(appendchain)){
    if(!is.null(appendchain$const)){
      const <- appendchain$const
    }else{
      stop("Pass previous fitted object to appendchain argument")
    }
  }else{
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
                              clustering,
                              fixedZbeta,
                              fixedZtheta,
                              maxClusters,
                              prior_alpha_beta,
                              prior_alpha_theta,
                              prior_rho,
                              prior_tau_theta,
                              prior_lambda_beta,
                              prior_lambda_theta,
                              prior_xi,
                              prior_sigma2,
                              prop_phi_a,
                              sharedlambda,
                              DLM,
                              DLMpenalty,
                              lagOrder,
                              diff,
                              MIM,
                              MIMorder,
                              LM,
                              betaOrder,
                              ## MH tuning
                              stepsize_logrho,
                              stepsize_loglambda_theta,
                              stepsize_logxi,
                              stepsize_logsigma2,
                              Vgridsearch,
                              gridsize,
                              rfbtries,
                              approx)
  }


  ## set up MCMC sampler with options
  update_params <- build_sampler(const)

  ## mcmc function (called by mclapply for multiple chains)
  run_sampler <- function(ind){

    ## get starting values
    if(!is.null(appendchain)){
      params_ss <- appendchain$lastvals ## start where last chain ended
    }else{
      params_ss <- get_starting_vals(const)
    }

    ## initialize results storage
    nkeep <- floor((niter-nburn)/nthin)# round

    keep_Zbeta <- matrix(0,ncol=length(params_ss$Zbeta),nrow=nkeep)
    keep_Ztheta <- matrix(0,ncol=length(params_ss$Ztheta),nrow=nkeep)
    keep_Vbeta <- matrix(0,ncol=const$Cbeta,nrow=nkeep)
    keep_Vtheta <- matrix(0,ncol=const$Ctheta,nrow=nkeep)
    keep_betastar <- matrix(0,ncol=const$d*const$Cbeta,nrow=nkeep)
    keep_thetastar <- matrix(0,ncol=const$Lq*const$Ctheta,nrow=nkeep)
    keep_omegastar <- matrix(0,ncol=const$L*const$Ctheta,nrow=nkeep)
    keep_alpha <- matrix(0,ncol=2,nrow=nkeep)
    keep_logrho <- matrix(0,ncol=1,nrow=nkeep)
    keep_lambda_beta <- matrix(0,ncol=const$Cbeta,nrow=nkeep)
    keep_loglambda_theta <- matrix(0,ncol=1,nrow=nkeep)
    keep_u <- matrix(0,ncol=const$n,nrow=nkeep)
    keep_xi <- matrix(0,ncol=1,nrow=nkeep)
    keep_sigma2 <- matrix(0,ncol=const$K,nrow=nkeep)
    keep_b0 <- matrix(0,ncol=const$K,nrow=nkeep)
    keep_betaZk <- matrix(0,ncol=const$K*const$pz,nrow=nkeep)

    ## MCMC
    maxerr <- 0 ## tracking errors
    for(ss in 1:niter){

      ## update parameters
      params_ss <- update_params(params_ss)
      if(params_ss$err==1){
        maxerr <- ss ## max iteration with clustering error
      }

      ## retain samples after burn-in and thinning
      if(ss>nburn & (ss-nburn)%%nthin==0){
        skeep <- (ss-nburn)/nthin

        keep_Zbeta[skeep,] <- c(params_ss$Zbeta)
        keep_Ztheta[skeep,] <- c(params_ss$Ztheta)
        keep_Vbeta[skeep,] <- params_ss$Vbeta
        keep_Vtheta[skeep,] <- params_ss$Vtheta
        keep_betastar[skeep,] <- params_ss$betastar
        keep_thetastar[skeep,] <- params_ss$thetastar
        keep_omegastar[skeep,] <- params_ss$omegastar
        keep_alpha[skeep,] <- params_ss$alpha
        keep_logrho[skeep,] <- params_ss$logrho
        keep_lambda_beta[skeep,] <- params_ss$lambda_beta
        keep_loglambda_theta[skeep,] <- params_ss$loglambda_theta
        keep_u[skeep,] <- params_ss$u
        keep_xi[skeep,] <- params_ss$xi
        keep_sigma2[skeep,] <- params_ss$sigma2
        keep_b0[skeep,] <- params_ss$b0
        keep_betaZk[skeep,] <- c(t(params_ss$betaZk)) #(each row for each k)

      }

    }

    res <- list(Zbeta=keep_Zbeta,
                Ztheta=keep_Ztheta,
                Vbeta=keep_Vbeta,
                Vtheta=keep_Vtheta,
                betastar=keep_betastar,
                thetastar=keep_thetastar,
                omegastar=keep_omegastar,
                alpha=keep_alpha,
                logrho=keep_logrho,
                lambda_beta=keep_lambda_beta,
                loglambda_theta=keep_loglambda_theta,
                u=keep_u,
                xi=keep_xi,
                sigma2=keep_sigma2,
                b0=keep_b0,
                betaZk=keep_betaZk,
                const=const,
                lastvals=params_ss)
    attr(res,"err") <- maxerr
    return(res)
  }


  ## run chains
  if(nchains==1){ ## run single chain
    samples <- run_sampler()
    maxerr <- attr(samples,"err")
    # samples$PSR <- NULL
  }else{ ## run multiple chains
    if(ncores==1){ ## run multiple chains (not parallel)
      chains <- vector(mode = "list", length = nchains)
      for(cc in 1:nchains){
        chains[[cc]] <- run_sampler()
      }
    }else{ ## run in parallel with nchains>1 and ncores>1
      chains <- parallel::mclapply(1:nchains,run_sampler,mc.cosamples=ncores)
    }


    ## append chains
    samples <- chains[[1]] ## use first for names/formatting
    samples <- lapply(names(samples)[!(names(samples)%in%c("const","lastvals"))], ## looping over names of variables to append
                  function(v){ Reduce('rbind',lapply(chains,function(chn) chn[[v]]))} )## for each variable name, extract element from each chain, then collapse them via rbind

    samples$const <- chains[[1]]$const
    samples$lastvals <- chains[[1]]$lastvals
    names(samples) <- names(chains[[1]])
    ## errs
    maxerr <- sapply(chains,function(x)attr(x,"err"))

  }

  ## summarize clusters
  Zbeta <- round(100*t(apply(samples$Zbeta,2,function(x) table(factor(x,levels=1:samples$const$Cbeta))))/nrow(samples$Zbeta))
  Ztheta <- round(100*t(apply(samples$Ztheta,2,function(x) table(factor(x,levels=1:samples$const$Ctheta))))/nrow(samples$Ztheta))
  outcome <- rep(1:samples$const$K,samples$const$p)
  exposure <- rep(1:samples$const$p,each=samples$const$K)
  summ <- data.frame(outcome,exposure,Zbeta=Zbeta,Ztheta=Ztheta)
  samples$cluster_summary <- summ[order(summ$outcome),]

  ## convert big matrix of Z to array (KxPxR)
  samples$Zbeta <- array(t(samples$Zbeta),dim=c(samples$const$K,samples$const$p,nrow(samples$Zbeta)))
  samples$Ztheta <- array(t(samples$Ztheta),dim=c(samples$const$K,samples$const$p,nrow(samples$Ztheta)))

  ## convert big coefficient matrices to lists of length C
  samples$betastar <- lapply(1:samples$const$Cbeta,function(cc)
    samples$betastar[,(cc-1)*samples$const$d+(1:samples$const$d)])
  samples$thetastar <- lapply(1:samples$const$Ctheta,function(cc)
    samples$thetastar[,(cc-1)*samples$const$Lq+(1:samples$const$Lq)])
  samples$omegastar <- lapply(1:samples$const$Ctheta,function(cc)
    samples$omegastar[,(cc-1)*samples$const$L+(1:samples$const$L)])

  ## assign beta and theta from betastar and thetastar
  samples$beta <- lapply(1:samples$const$p,function(jj){
                    lapply(1:samples$const$K,function(kk){
                      Reduce("rbind",lapply(1:nrow(as.matrix(samples$betastar[[1]])),function(rr){
                        return(as.matrix(samples$betastar[[samples$Zbeta[kk,jj,rr]]])[rr,,drop=F])
                      }))
                    })
                   })
  samples$theta <- lapply(1:samples$const$p,function(jj){
                      lapply(1:samples$const$K,function(kk){
                        Reduce("rbind",lapply(1:nrow(as.matrix(samples$thetastar[[1]])),function(rr){
                          return(as.matrix(samples$thetastar[[samples$Ztheta[kk,jj,rr]]])[rr,,drop=F])
                        }))
                      })
                    })
  samples$omega <- lapply(1:samples$const$p,function(jj){
    lapply(1:samples$const$K,function(kk){
      Reduce("rbind",lapply(1:nrow(as.matrix(samples$omegastar[[1]])),function(rr){
        return(as.matrix(samples$omegastar[[samples$Ztheta[kk,jj,rr]]])[rr,,drop=F])
      }))
    })
  })

  ## track errors
  samples$maxerr <- paste0(round(100*maxerr/niter),"%")

  return(samples)
}
