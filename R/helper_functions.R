
## function to compute pi
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


## get basis functions
get_Btheta <- function(Xomega,const,params=NULL,k,j){

    if(const$MIM==FALSE | j==1 ){
      IDprod <- 1
    }else{  ## identifiability product for MIM
      IDprod <- prod(params$Ztheta[k,1:(j-1)]!=params$Ztheta[k,j])
    }

    if(const$LM==TRUE){ ## if forcing linearity
      return(IDprod*Xomega)
    }else{
      return(IDprod*mgcv::PredictMat(const$SS,data=data.frame(Xomega)))
    }



}


## get derivatives of Basis functions
get_DerivBtheta <- function(Xomega,const,params,k,j){
  # if(const$MIM==TRUE){
  if(const$MIM==FALSE | j==1 ){
    IDprod <- 1
  }else{  ## identifiability product for MIM
    IDprod <- prod(params$Ztheta[k,1:(j-1)]!=params$Ztheta[k,j])
  }
  if(const$LM==TRUE){ ## deriv is 1 if forcing linear effects
    return(IDprod*((Xomega)^0))
  }else{
    return(IDprod*mgcv::PredictMat(const$SSderiv,data=data.frame(Xomega)))
  }


}


## compute B times the corresponding beta for single cluster
get_B_beta_k <- function(params,const,kk){

  return(sapply(1:const$p,function(jj){
    get_Btheta(const$X[[jj]]%*%params$omegastar[(params$Ztheta[kk,jj]-1)*const$L+(1:const$L)],const,params,kk,jj)%*%params$betastar[(params$Zbeta[kk,jj]-1)*const$d+(1:const$d)]
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

  return( get_DerivBtheta(const$X[[jj]]%*%params$omegastar[(params$Ztheta[kk,jj]-1)*const$L+(1:const$L)],const,params,kk,jj)%*%params$betastar[(params$Zbeta[kk,jj]-1)*const$d+(1:const$d)] )

}


## alternate version for updating Ztheta=b in update_clustMemb
get_Btheta_b <- function(Xomega,const,Ztheta=NULL,k,j){
  # if(const$MIM==TRUE){
  if(const$MIM==FALSE | j==1 ){
    IDprod <- 1
  }else{  ## identifiability product for MIM
    IDprod <- prod(Ztheta[k,1:(j-1)]!=Ztheta[k,j])
  }
  if(const$LM==TRUE){ ## if forcing linear effects
    return(IDprod*Xomega)
  }else{
    return(IDprod*mgcv::PredictMat(const$SS,data=data.frame(Xomega)))
  }

}
## alternate version for updating Ztheta=b in update_clustMemb
get_B_beta_k_b <- function(params,const,kk,Ztheta){

  return(sapply(1:const$p,function(jj){
    get_Btheta_b(const$X[[jj]]%*%params$omegastar[(params$Ztheta[kk,jj]-1)*const$L+(1:const$L)],const,Ztheta,kk,jj)%*%params$betastar[(params$Zbeta[kk,jj]-1)*const$d+(1:const$d)]
  }))

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


## get components for omegastar sampling
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




## compute omegastar from phistar
get_theta <- function(phistar){
  return(c(sin(phistar),1) * cumprod(c(1,cos(phistar))))
}

## compute phistar from omegastar
get_phi <- function(omegastar){

  phistar <- c(asin(omegastar[1]))
  if(length(omegastar)>2){
    for(jj in 2:(length(omegastar)-1)){
      phistar <- c(phistar,asin(omegastar[jj]/(prod(cos(phistar)))))
    }
  }
  return(phistar)
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
    const$DLM <- TRUE ##
    const$DLMpenalty <- FALSE ## mutually exclusive options
    const$X <- rep(const$X,const$MIMorder)
    const$MIM <- FALSE ## shut off the MIM IDprod parameterization . can delete the IDprod components later
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
      PEN <- ginv(t(Q) %*% ginv(t(DDtemp) %*% DDtemp) %*% Q)
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
    const$SS <- mgcv::smoothCon(s(Xomega,bs="bs",k=betaOrder),data=data.frame(Xomega),
                                absorb.cons = TRUE)[[1]] ## should be true for multiple exposures
    const$SSderiv <- const$SS ## make a smooth object for computing first order derivatives
    const$SSderiv$deriv <- 1 ## first order derivatives
    const$d <- ncol(const$SS$X)
    const$mu0 <- rep(0,const$d) ## probably no need to change this from 0
    const$invSig0 <- const$SS$S[[1]]
  }

  ## grid for
  const$grid <- (1:(const$gridsize-1))/const$gridsize ## exclude 0s and 1s

  if(const$approx==FALSE & const$prior_tau_theta!=0){
    print("Forcing prior_tau_theta=0 for polar method.")
    const$prior_tau_theta <- 0 ## force to be 0 for polar sampling, otherwise normalizing constant wont work
  }

  return(const)
}


## get starting parameter values
get_starting_vals <- function(const){
  params <- list()

  ## other params from prior
  params$alpha <- c(rgamma(1,shape=const$prior_alpha_beta[1],rate=const$prior_alpha_beta[2]),
                    rgamma(1,shape=const$prior_alpha_theta[1],rate=const$prior_alpha_theta[2]) )

  params$Vbeta <- c(rbeta(const$Cbeta-1,1,params$alpha[1]),1) ## final value is 1 (truncated)
  params$Vtheta <- c(rbeta(const$Ctheta-1,1,params$alpha[2]),1) ## final value is 1 (truncated)

  if(const$clustering=="both"){ ## only defined if clustering jointly
    params$logrho <- log(rgamma(1,shape=const$prior_rho[1],rate=const$prior_rho[2]))
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
    params$lambda_beta <- rgamma(1,shape=1,rate=1)
  }else{
    params$lambda_beta <- rgamma(const$Cbeta,shape=1,rate=1)
  }
  if(const$LM==TRUE){
    params$lambda_beta <- 5
  }


  if(const$DLM==TRUE & const$DLMpenalty==TRUE){
    params$loglambda_theta <- log(rgamma(1,shape=const$prior_lambda_theta[1],rate=const$prior_lambda_theta[2]))
  }else{
    params$loglambda_theta <- -Inf
  }

  params$sigma2 <- 1/rgamma(const$K,shape=5,rate=1)

  #
  if(const$K>1){
    params$xi <- 1/rgamma(1,shape=10,rate=1)
    params$u <- rnorm(const$n,0,1)

  }else{
    params$xi <- 0
    params$u <- rep(0,const$n)
  }

  if(const$L==1){ ## dont need theta if L==1
    params$phistar <- NULL
    params$omegastar <- params$thetastar <- rep(1,const$Ctheta)

  }else{
    ## initialize omegastar
    if(const$approx==FALSE){ ## use polar coordinate parameterization
      ## starting with thetastar from fisherbingham, flipping sign if necessary, then converting to phistar
      thetadraw <- t(rFisherBingham(const$Ctheta,mu = const$prior_tau_theta*rep(1,const$Lq), Aplus = 0))
      signid <- thetadraw[const$Lq,]<0
      thetadraw[,signid] <- -thetadraw[,signid]
      params$thetastar <- c(thetadraw)
      params$phistar <- sapply(1:const$Ctheta,function(cc){
        return(get_phi(params$thetastar[(cc-1)*(const$Lq)+1:(const$Lq)]))
      })

    }else{ ## otherwise just draw from fb
      params$phistar <- NULL
      params$thetastar <- c(t(rFisherBingham(const$Ctheta,mu = const$prior_tau_theta*rep(1,const$Lq), Aplus = 0)))
    }
    ## convert thetastar to omegastar if necessary (Psi is just Identity if not a DLM, then omegastar=thetastar)
    params$omegastar <- sapply(1:const$Ctheta,function(cc){
      return(c(const$Psi%*%params$thetastar[(cc-1)*(const$Lq)+1:(const$Lq)]))
    })
  }

  ## Get Xomega and basis functions
  Xomega <- sapply(1:const$p,function(jj){matrix(sapply(1:const$K,function(kk){const$X[[jj]]%*%params$omegastar[(params$Zbeta[kk,jj]-1)*const$L+(1:const$L)]}))})

  ## betastar
  params$betastar <- c(t(rmvnorm(const$Cbeta,const$mu0,diag(const$d))))

  ## intercepts
  params$b0 <- rnorm(const$K)

  # linear coefficients
  if(is.null(const$Z)){
    params$betaZk <- matrix(0,ncol=1,nrow=const$K)
  }else{
    # matrix(rnorm(const$pz*const$K),ncol=const$pz,nrow=const$K)
    # alternatively could just start with LS:
    ZZ <- cbind(1,const$Zcovariates)
    LS <- t(solve(t(ZZ)%*%ZZ)%*%t(ZZ)%*%matrix(const$y,ncol=const$K))
    params$b0 <- c(LS[,1])
    params$betaZk <- LS[,-1,drop=F]
  }



  params$err <- 0
  return(params)
}

