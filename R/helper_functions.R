

#' Compute pi matrix
#' @keywords internal
compute_pi <- function(params,const){

  if(const$clustering=="both"){
    ## compute "marginals" under independence
    pi_beta <- params$Vbeta*c(1,cumprod(1-params$Vbeta)[-const$C])
    pi_theta <- params$Vtheta*c(1,cumprod(1-params$Vtheta)[-const$C])

    ## joint probs
    pimat <- (pi_beta)%*%t(pi_theta) + diag(exp(params$logrho)*(pi_beta*pi_theta))

  }else if(const$clustering=="beta"){ ## just repeat for the theta columns
    pi_beta <- params$Vbeta*c(1,cumprod(1-params$Vbeta)[-const$Cbeta])
    pimat <- matrix(pi_beta,nrow=length(pi_beta),ncol=max(const$fixedZtheta))
  }else if(const$clustering=="theta"){ ## just repeat for the beta rows
    pi_theta <- params$Vtheta*c(1,cumprod(1-params$Vtheta)[-const$Ctheta])
    pimat <- matrix(pi_theta,byrow=TRUE,nrow=max(const$fixedZbeta),ncol=length(pi_theta))
  }

  ## standardize
  pistarmat <- pimat/sum(pimat)

  params$pimat <- pimat
  params$pistarmat <- pistarmat
  return(params)
}


#' Get basis functions
#' @keywords internal
get_Btheta <- function(Xomega,const,params=NULL,k,j){

    if(const$LM==TRUE | const$NONSEP==TRUE){ ## if forcing linearity
      return(Xomega)
    }else{
      return(mgcv::PredictMat(const$SS,data=data.frame(Xomega)))
    }

}


#' Get derivatives of Basis functions
#' @keywords internal
get_DerivBtheta <- function(Xomega,const,params,k,j){

  if(const$LM==TRUE){ ## deriv is 1 if forcing linear effects
    return(((Xomega)^0))
  }else{
    return(mgcv::PredictMat(const$SSderiv,data=data.frame(Xomega)))
  }

}


#' Compute B times the corresponding beta for single cluster
#' @keywords internal
get_B_beta_k <- function(params,const,kk){

  if(const$NONSEP==FALSE){
    return(sapply(1:const$p,function(jj){
      get_Btheta(const$X[[jj]]%*%params$omegastar[(params$Ztheta[kk,jj]-1)*const$L+(1:const$L)],const,params,kk,jj)%*%params$betastar[(params$Zbeta[kk,jj]-1)*const$d+(1:const$d)]
    }))
  }else{
    return(sapply(1:const$p,function(jj){
      get_Btheta(const$X[[jj]],const,params,kk,jj)%*%params$betastar[(params$Zbeta[kk,jj]-1)*const$d+(1:const$d)]
    }))
  }


}

#' Compute B times the corresponding beta for all obs (d)
#' @keywords internal
get_B_beta <- function(params,const){

  return(Reduce("rbind",lapply(1:const$K,function(kk){
    get_B_beta_k(params,const,kk)
  })))

}


#' Compute DerivB times the corresponding beta for given k,j
#' @keywords internal
get_DerivB_beta <- function(params,const,kk,jj){

  return( get_DerivBtheta(const$X[[jj]]%*%params$omegastar[(params$Ztheta[kk,jj]-1)*const$L+(1:const$L)],const,params,kk,jj)%*%params$betastar[(params$Zbeta[kk,jj]-1)*const$d+(1:const$d)] )

}


#' Get Btheta (alternate version for updating Ztheta=b in update_clustMemb)
#' @keywords internal
get_Btheta_b <- function(Xomega,const,Ztheta=NULL,k,j){

  if(const$LM==TRUE | const$NONSEP==TRUE){ ## if forcing linear effects
    return(Xomega)
  }else{
    return(mgcv::PredictMat(const$SS,data=data.frame(Xomega)))
  }

}

#' Get Btheta times beta (alternate version for updating Ztheta=b in update_clustMemb)
#' @keywords internal
get_B_beta_k_b <- function(params,const,kk,Ztheta){

  if(const$NONSEP==FALSE){
    return(sapply(1:const$p,function(jj){
      get_Btheta_b(const$X[[jj]]%*%params$omegastar[(params$Ztheta[kk,jj]-1)*const$L+(1:const$L)],const,Ztheta,kk,jj)%*%params$betastar[(params$Zbeta[kk,jj]-1)*const$d+(1:const$d)]
    }))
  }else{
    return(sapply(1:const$p,function(jj){
      get_Btheta_b(const$X[[jj]],const,Ztheta,kk,jj)%*%params$betastar[(params$Zbeta[kk,jj]-1)*const$d+(1:const$d)]
    }))
  }


}



#' Computes density of V_c ## up to proportionality constant
#' @keywords internal
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


#' Get components for omegastar sampling
#' @keywords internal
get_XTyhat <- function(cc,whichk,whichkj,params,const){

  Xhat_k <- lapply(whichk,function(kk){
    Reduce("+",
      lapply(whichkj[[kk]],function(jj){
        c(get_DerivB_beta(params,const,kk,jj))*const$XPsi[[jj]]
        })
      )/sqrt(params$sigma2[kk]) ## sqrt because it gets multiplied below
    })

  yhat_k <- lapply(whichk,function(kk){
    (const$y[const$k_index==kk]-params$b0[kk]-apply(get_B_beta_k(params,const,kk),1,sum)-params$xi*sqrt(params$sigma2[kk])*params$u-const$Zcovariates%*%params$betaZk[kk,]+
      Reduce("+",
             lapply(whichkj[[kk]],function(jj){
               (c(get_DerivB_beta(params,const,kk,jj))*const$XPsi[[jj]])%*%params$thetastar[(cc-1)*const$Lq+(1:const$Lq)]
             })
       ))/sqrt(params$sigma2[kk]) ## sqrt because it gets multiplied below
  })

  return(list(XTX=Reduce("+",lapply(Xhat_k,function(XX){t(XX)%*%XX})),
              XTy=Reduce("+",lapply(1:length(Xhat_k),function(kk){t(Xhat_k[[kk]])%*%yhat_k[[kk]]}))))
}



#' Compute omegastar from phistar
#' @keywords internal
get_theta <- function(phistar){
  return(c(sin(phistar),1) * cumprod(c(1,cos(phistar))))
}

#' Compute phistar from omegastar
#' @keywords internal
get_phi <- function(omegastar){

  phistar <- c(asin(omegastar[1]))
  if(length(omegastar)>2){
    for(jj in 2:(length(omegastar)-1)){
      phistar <- c(phistar,asin(omegastar[jj]/(prod(cos(phistar)))))
    }
  }
  return(phistar)
}




#' Initialize constants
#' @keywords internal
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
                             # prior_sigma2_u,
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
                             NONSEP,
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
                             approx){


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
                clustering=clustering,
                fixedZbeta=fixedZbeta,
                fixedZtheta=fixedZtheta,
                maxClusters=maxClusters,
                prior_alpha_beta=prior_alpha_beta,
                prior_alpha_theta=prior_alpha_theta,
                prior_rho=prior_rho,
                prior_tau_theta=prior_tau_theta,
                prior_lambda_beta=prior_lambda_beta,
                prior_lambda_theta=prior_lambda_theta,
                # prior_sigma2_u=prior_sigma2_u,
                prior_xi=prior_xi,
                prior_sigma2=prior_sigma2,
                prop_phi_a=prop_phi_a,
                sharedlambda=sharedlambda,
                DLM=DLM,
                DLMpenalty=DLMpenalty,
                lagOrder=lagOrder,
                diff=diff,
                MIM=MIM,
                MIMorder=MIMorder,
                NONSEP=NONSEP,
                LM=LM,
                betaOrder=betaOrder,
                ## MH tuning
                stepsize_logrho=stepsize_logrho,
                stepsize_loglambda_theta=stepsize_loglambda_theta,
                stepsize_logxi=stepsize_logxi,
                stepsize_logsigma2=stepsize_logsigma2,
                Vgridsearch=Vgridsearch,
                gridsize=gridsize,
                rfbtries=rfbtries,
                approx=approx)

  if(is.list(X)==FALSE){
    const$X <- list(X)
  }
  if(const$MIM==TRUE){ ## no longer using the MIM IDprod paramaterization
    const$DLM <- FALSE ## mutually exclusive options
    const$DLMpenalty <- FALSE ## mutually exclusive options
    const$X <- rep(const$X,const$MIMorder)
    const$NONSEP <- FALSE ## NONSEP only used for DLM
  }


  const$p <- length(const$X)
  const$n <- nrow(const$X[[1]])
  const$K <- ncol(Y)
  const$L <- ncol(const$X[[1]])
  if(is.null(const$Z)){
    const$Zcovariates <- matrix(0,ncol=1,nrow=const$n)
    const$pz <- 1
  }else{
    const$Zcovariates <- as.matrix(const$Z)
    const$pz <- ncol(const$Z)
  }
  if(is.null(maxClusters)){
    const$C <- const$K ## default number of clusters
  }else{
    const$C <- maxClusters
  }
  const$i_index <- rep(1:const$n,const$K) ## 1....1,2...2,....
  const$k_index <- rep(1:const$K,each=const$n)      ## 1:K,1:K,...

  ## max number of clusters
  ### different ONLY when clustering!="both"
  const$Cbeta <- const$Ctheta <- const$C

  ## no theta updates if single exposure
  if(const$L==1){
    if(const$clustering=="both"){
      const$clustering <- "beta"
    }else if(const$clustering=="theta"){
      const$clustering <- "neither"
    }
  }

  ## NONSEP version of DLM only has beta
  if(const$NONSEP==TRUE){
    ## only cluster beta
    if(const$clustering=="both"){
      const$clustering <- "beta"
    }else if(const$clustering=="theta"){
      const$clustering <- "neither"
    }

    ## force DLM=TRUE for nonseparable version
    const$DLM=TRUE
    const$DLMpenalty=TRUE
    const$LM=FALSE
    # const$sharedlambda=TRUE

    # force dim<=20
    if(const$betaOrder*const$lagOrder>20){
      print("Setting dim 5x4 for NONSEP version")
      const$betaOrder <- 5
      const$lagOrder <- 4
    }


  }

  ## cannot provide fixed IDs if clustering
  if(const$clustering=="both" | const$clustering=="beta"){
    const$fixedZbeta <- NULL
  }
  if(const$clustering=="both" | const$clustering=="theta"){
    const$fixedZtheta <- NULL
  }

  ## if not clustering and no fixed IDs, set to 1:KP
  if(const$clustering=="neither" | const$clustering=="theta"){
    if(is.null(const$fixedZbeta)){
      const$fixedZbeta <- matrix(1:(const$K*const$p),nrow=const$K,ncol=const$p,byrow=TRUE)
    }else if(nrow(as.matrix(const$fixedZbeta))!=const$K | ncol(as.matrix(const$fixedZbeta))!=const$p){
      print("Fixed Z of wrong dimension.")
      const$fixedZbeta <- matrix(1:(const$K*const$p),nrow=const$K,ncol=const$p,byrow=TRUE)
    }
    const$Cbeta <- max(const$fixedZbeta)
  }
  if(const$clustering=="neither" | const$clustering=="beta"){
    if(is.null(const$fixedZtheta)){
      const$fixedZtheta <- matrix(1:(const$K*const$p),nrow=const$K,ncol=const$p,byrow=TRUE)
    }else if(nrow(as.matrix(const$fixedZtheta))!=const$K | ncol(as.matrix(const$fixedZtheta))!=const$p){
      print("Fixed Z of wrong dimension.")
      const$fixedZtheta <- matrix(1:(const$K*const$p),nrow=const$K,ncol=const$p,byrow=TRUE)
    }
    const$Ctheta <- max(const$fixedZtheta)
  }


  ##
  if(NONSEP==FALSE){ ## standard parameterization
    if(DLM==TRUE){ ## create penalty matrix

      if(!is.null(const$lagOrder)){
        #########################################
        ## B splines (with difference penalty) ##
        QQ <- dlnm::ps(1:const$L,
                       diff=const$diff, ##
                       df=const$lagOrder, ##
                       intercept=TRUE)
        B1 <- as.matrix(data.frame(QQ)) # basis matrix for smoothed l
        qrB <- qr(B1)
        Q <- qr.Q(qrB)
        const$Psi <- Q
        const$Lq <- ncol(const$Psi) ## =lagOrder
        const$XPsi <- lapply(1:length(const$X),function(jj){ return(const$X[[jj]]%*%const$Psi) })
        # DDtemp = getDtf(L, ord = 12) ## same as below without package dependency
        DDtemp <- diag(const$L)
        for(jj in 1:const$diff){DDtemp <- diff(DDtemp)}
        PEN <- MASS::ginv(t(Q) %*% MASS::ginv(t(DDtemp) %*% DDtemp) %*% Q)
        const$PEN = (PEN + t(PEN))/2 ## ensuring symmetry

      }else{## if is.null(lagOrder) --> dont do basis expansion (no dimension reduction), but still do penalty

        const$Psi <- diag(const$L)
        const$Lq <- const$L
        const$XPsi <- const$X #lapply(1:length(const$X),function(jj){ return(const$X[[jj]]%*%const$Psi) })
        D <- diag(const$L)
        for(jj in 1:const$diff){D <- diff(D)}
        diffpen <- t(D)%*%D ## gives same as dlnm::ps(1:L,diff=diff,df=L,intercept=TRUE)
        const$PEN <- diffpen

      }



    }else{ ## no penalties
      const$PEN <- matrix(0,ncol=const$L,nrow=const$L)
      const$Lq <- const$L ##
      const$Psi <- diag( const$L) ## basis matrix
      const$XPsi <- const$X #
    }


    ## basis functions
    ### Uses same basis functions for all exposures.
    if(const$LM==TRUE){
      const$sharedlambda <- TRUE ## single ridge penalty
      const$SS <- NULL
      const$SSderiv <- NULL ## make a smooth object for computing first order derivatives
      const$d <- 1
      const$mu0 <- rep(0,const$d) ## probably no need to change this from 0
      const$invSig0 <- as.matrix(1)
    }else{
      tempomega <- rep(1,const$L)/sqrt(const$L) ## just a standard theta value in order to compute basis functions
      ## combining all exposures to get basis functions
      Xomega <- c(sapply(1:const$p,function(jj){const$X[[jj]]%*%tempomega}))
      const$SS <- mgcv::smoothCon(mgcv::s(Xomega,bs="bs",k=betaOrder),data=data.frame(Xomega),
                                  absorb.cons = TRUE)[[1]] ## should be true for multiple exposures
      const$SSderiv <- const$SS ## make a smooth object for computing first order derivatives
      const$SSderiv$deriv <- 1 ## first order derivatives
      const$d <- ncol(const$SS$X)
      const$mu0 <- rep(0,const$d) ## probably no need to change this from 0
      const$invSig0 <- const$SS$S[[1]]
    }

  }else if (const$NONSEP==TRUE & const$DLM==TRUE){
    ## non-separable parametrization requires tensor product basis

    # ## DLNM version --didn't apply same constraints. new version below
    # Xlong <- Reduce("rbind",const$X) ## stack all exposure matrices to get single basis representation.
    # cb.obj <- dlnm::crossbasis(x=Xlong, lag=c(1,const$L),
    #                      argvar=list(fun="ps",df=const$betaOrder,degree=2,intercept=FALSE),
    #                      arglag=list(fun="ps",df=const$lagOrder,degree=2,intercept=TRUE))
    #
    # ## penalty matrices
    # penMatlist <- dlnm::cbPen(cb.obj)
    # const$invSig0 <- penMatlist$Svar ## penalty matrix for X
    # const$PEN <- penMatlist$Slag ## penalty matrix for L
    #
    # ## define B matrices a priori since they wont change (not a fuction of theta)
    # const$X <- lapply(1:const$p,function(pp){cb.obj[(pp-1)*const$n+(1:const$n),]})
    #
    # ## dimension
    # const$d <- ncol(cb.obj) ## both now combined into single nonseparable beta
    # const$Lq <- 1
    #
    # ## unused components carried over from X basis
    # const$SS <- cb.obj ## now a cross-basis object
    # const$SSderiv <- NULL
    # const$mu0 <- rep(0,const$d)

    # Xlong <- Reduce("rbind",const$X) ## stack all exposure matrices to get single basis representation.
    # lag <- matrix(1:const$L,byrow=TRUE,ncol=ncol(Xlong),nrow=nrow(Xlong)) ## matrix of lag times to construct tensor
    #
    # const$SS <- mgcv::smoothCon(mgcv::te(x,lag,bs="ps",k=c(const$betaOrder,const$lagOrder)),data=list(x=Xlong,lag=lag),n=nrow(Xlong),absorb.cons=TRUE)[[1]]
    # const$invSig0 <- const$SS$S[[1]]
    # const$PEN <- const$SS$S[[2]]
    #
    # ## define B matrices a priori since they wont change (not a fuction of theta)
    # const$X <- lapply(1:const$p,function(pp){const$SS$X[(pp-1)*const$n+(1:const$n),]})
    # const$lagmat <- lag ## for PredictMat ## EDIT. might be long?

    Xlong <- Reduce("rbind",const$X) ## stack all exposure matrices to get single basis representation.
    laglong <- matrix(1:const$L,byrow=TRUE,ncol=ncol(Xlong),nrow=nrow(Xlong)) ## matrix of lag times to construct tensor
    lag <- matrix(1:const$L,byrow=TRUE,ncol=ncol(Xlong),nrow=nrow(const$X[[1]]))

    const$SS <- mgcv::gam(rep(Y[,1],length(const$X))~s(x,bs="ps",k=const$betaOrder)+
                            ti(x,lag,bs="ps",k=c(const$betaOrder,const$lagOrder)),
                                data=list(x=Xlong,lag=laglong),
                         sp=rep(1,3*length(const$X)))# just using this for basis #,fit=TRUE

    # predict(obj$const$SS,data=list(x=Xlong,lag=laglong),type="lpmatrix"))[,-1]

    const$invSig0 <- as.matrix(Matrix::bdiag(const$SS$smooth[[1]]$S[[1]],
                                             const$SS$smooth[[2]]$S[[1]]))
    const$PEN <- as.matrix(Matrix::bdiag(matrix(0,nrow=nrow(const$SS$smooth[[1]]$S[[1]]),ncol=ncol(const$SS$smooth[[1]]$S[[1]])),
                               const$SS$smooth[[2]]$S[[2]]))

    ## define B matrices a priori since they wont change (not a fuction of theta)
    SS_X <- model.matrix(const$SS)[,-1]
    const$X <- lapply(1:const$p,function(pp){SS_X[(pp-1)*const$n+(1:const$n),]})
    const$lagmat <- lag ## for PredictMat ## EDIT.

    ## dimension
    const$d <- ncol(SS_X) ## both now combined into single nonseparable beta
    const$Lq <- 1

    ## unused components carried over from X basis
    const$SSderiv <- NULL
    const$mu0 <- rep(0,const$d)

    ## unused components carried over from L basis
    const$Psi <- NULL
    const$XPsi <- NULL

  }



  ## grid for
  const$grid <- (1:(const$gridsize-1))/const$gridsize ## exclude 0s and 1s

  if(const$approx==FALSE & const$prior_tau_theta!=0){
    print("Forcing prior_tau_theta=0 for polar method.")
    const$prior_tau_theta <- 0 ## force to be 0 for polar sampling, otherwise normalizing constant wont work
  }

  return(const)
}


#' Get starting parameter values
#' @keywords internal
get_starting_vals <- function(const){
  params <- list()

  ## other params from prior
  params$alpha <- c(stats::rgamma(1,shape=const$prior_alpha_beta[1],rate=const$prior_alpha_beta[2]),
                    stats::rgamma(1,shape=const$prior_alpha_theta[1],rate=const$prior_alpha_theta[2]) )

  params$Vbeta <- c(stats::rbeta(const$Cbeta-1,1,params$alpha[1]),1) ## final value is 1 (truncated)
  params$Vtheta <- c(stats::rbeta(const$Ctheta-1,1,params$alpha[2]),1) ## final value is 1 (truncated)

  if(const$clustering=="both"){ ## only defined if clustering jointly
    params$logrho <- log(stats::rgamma(1,shape=const$prior_rho[1],rate=const$prior_rho[2]))
  }else{
    params$logrho <- -999
  }

  if(const$clustering!="neither"){
    params <- compute_pi(params,const) ## update weights pimat and pimatstar

    ## draw Zbeta and Ztheta according to pistarmat
    params$Zbeta <- params$Ztheta <- matrix(NA,nrow=const$K,ncol=const$p)
    for(kk in 1:const$K){
      for(jj in 1:const$p){
        ## sample 1 of C^2 with correct probabilities
        ab <- sample(1:(const$Cbeta*const$Ctheta),1,prob=c(params$pistarmat))
        ## check which a and b this corresponds to
        Zk <- which(matrix(1:(const$Cbeta*const$Ctheta),ncol=const$Ctheta)==ab,arr.ind = TRUE) ## which element of CxC matrix
        params$Zbeta[kk,jj] <- Zk[1]  ## a=corresponding row
        params$Ztheta[kk,jj] <- Zk[2] ## b=corresponding column
      }
    }
  }
  ## use fixed IDs if not clustering
  if(!is.null(const$fixedZbeta)){
    params$Zbeta <- const$fixedZbeta
  }
  if(!is.null(const$fixedZtheta)){
    params$Ztheta <- const$fixedZtheta
  }


  if(const$sharedlambda==TRUE){
    params$lambda_beta <- stats::rgamma(1,shape=1,rate=1)
  }else{
    params$lambda_beta <- stats::rgamma(const$Cbeta,shape=1,rate=1)
  }
  if(const$LM==TRUE){
    params$lambda_beta <- 5
  }

  if(const$DLM==TRUE & const$DLMpenalty==TRUE){
    params$loglambda_theta <- log(stats::rgamma(1,shape=const$prior_lambda_theta[1],rate=const$prior_lambda_theta[2]))
  }else{
    params$loglambda_theta <- -Inf
  }

  params$sigma2 <- 1/stats::rgamma(const$K,shape=5,rate=1)

  #
  if(const$K>1){
    params$xi <- 1/stats::rgamma(1,shape=10,rate=1)
    params$u <- stats::rnorm(const$n,0,1)

  }else{
    params$xi <- 0
    params$u <- rep(0,const$n)
  }


  if(const$L==1 | const$Lq==1){ ## dont need theta if L==1
    params$phistar <- NULL
    params$omegastar <- params$thetastar <- rep(1,const$Ctheta)

  }else{
    ## initialize omegastar
    if(const$approx==FALSE){ ## use polar coordinate parameterization
      ## starting with thetastar from fisherbingham, flipping sign if necessary, then converting to phistar
      thetadraw <- t(simdd::rFisherBingham(const$Ctheta,mu = const$prior_tau_theta*rep(1,const$Lq), Aplus = 0))
      signid <- thetadraw[const$Lq,]<0
      thetadraw[,signid] <- -thetadraw[,signid]
      params$thetastar <- c(thetadraw)
      params$phistar <- c(sapply(1:const$Ctheta,function(cc){
        return(get_phi(params$thetastar[(cc-1)*(const$Lq)+1:(const$Lq)]))
      }))

    }else{ ## otherwise just draw from fb
      params$phistar <- NULL
      params$thetastar <- c(t(simdd::rFisherBingham(const$Ctheta,mu = const$prior_tau_theta*rep(1,const$Lq), Aplus = 0)))
    }
    ## convert thetastar to omegastar if necessary (Psi is just Identity if not a DLM, then omegastar=thetastar)
    params$omegastar <- sapply(1:const$Ctheta,function(cc){
      return(c(const$Psi%*%params$thetastar[(cc-1)*(const$Lq)+1:(const$Lq)]))
    })
  }

  ## Get Xomega and basis functions
  # Xomega <- sapply(1:const$p,function(jj){matrix(sapply(1:const$K,function(kk){const$X[[jj]]%*%params$omegastar[(params$Zbeta[kk,jj]-1)*const$L+(1:const$L)]}))})

  ## betastar
  params$betastar <- c(t(mvtnorm::rmvnorm(const$Cbeta,const$mu0,diag(const$d))))

  ## intercepts
  params$b0 <- stats::rnorm(const$K)

  # linear coefficients
  if(is.null(const$Z)){
    params$betaZk <- matrix(0,ncol=1,nrow=const$K)
  }else{
    # matrix(stats::rnorm(const$pz*const$K),ncol=const$pz,nrow=const$K)
    # alternatively could just start with LS:
    ZZ <- cbind(1,const$Zcovariates)
    LS <- t(solve(t(ZZ)%*%ZZ)%*%t(ZZ)%*%matrix(const$y,ncol=const$K))
    params$b0 <- c(LS[,1])
    params$betaZk <- LS[,-1,drop=F]
  }

  params$err <- 0
  return(params)
}

