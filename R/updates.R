## update functions



#' update cluster Membership
#' @keywords internal
update_clustMemb <- function(params,const){

  params$err <- 0 ## for error tracking

  for(kk in 1:const$K){## loop over outcomes
    for(jj in 1:const$p){## loop over exposures

      ## Zbetakj ##
      if(const$clustering=="both" | const$clustering=="beta"){

        ## necessary components
        B_beta <- get_B_beta_k(params,const,kk)
        y_B_u <- const$y[const$k_index==kk]-params$b0[kk]-(apply(B_beta[,-jj,drop=F],1,sum) +params$xi*sqrt(params$sigma2[kk])*params$u+const$Zcovariates%*%params$betaZk[kk,])
        if(const$NONSEP==FALSE){ ## allow clustering beta in non-separable parameterization
          Bth_kj <- get_Btheta(const$X[[jj]]%*%params$omegastar[(params$Ztheta[kk,jj]-1)*const$L+(1:const$L)],const,params,kk,jj)
        }else{
          Bth_kj <- get_Btheta(const$X[[jj]],const,params,kk,jj)
        }

        ## compute probabilities for all possible a
        logprobs <- c(sapply(1:const$Cbeta, function(a){  ## loop over a (rows; beta clusters)
          (log(params$pimat[a,params$Ztheta[kk,jj]]) -(0.5/params$sigma2[kk])*sum((y_B_u-Bth_kj%*%params$betastar[(a-1)*const$d+(1:const$d)])^2))
        })) ## to be standardized below

        ## error handling for very small values
        if(!is.finite(sum(exp(logprobs)/sum(exp(logprobs))))){
          stbfctr <- -500-max(logprobs) ## max largest value -500
          logprobs <- logprobs+stbfctr ## multiply all probs by common factor
        }

        ## sample 1 of C with correct probabilities
        newZbeta <- tryCatch(sample(1:const$Cbeta,1,prob=exp(logprobs)/sum(exp(logprobs))), ## standardized probs
                 error=function(err){NULL})

        if(!is.null(newZbeta)){
          params$Zbeta[kk,jj] <- newZbeta
        }else{## track errors
          params$err <- 1
        }

      }


      ## Zthetakj ##
      if(const$clustering=="both" | const$clustering=="theta"){

        ## compute probabilities for all possible b

        ## now have changed
        B_beta <- get_B_beta_k(params,const,kk)
        y_B_u <- const$y[const$k_index==kk]-params$b0[kk]-(apply(B_beta[,-jj,drop=F],1,sum) +params$xi*sqrt(params$sigma2[kk])*params$u+const$Zcovariates%*%params$betaZk[kk,])

        logprobs <- c(sapply(1:const$Ctheta, function(b){  ## loop over b (columns; theta clusters)
          (log(params$pimat[params$Zbeta[kk,jj],b]) -(0.5/params$sigma2[kk])*sum((y_B_u- get_Btheta(const$X[[jj]]%*%params$omegastar[(b-1)*const$L+(1:const$L)],const,params,kk,jj)%*%params$betastar[(params$Zbeta[kk,jj]-1)*const$d+(1:const$d)])^2))
        }))## to be standardized below

        ## error handling for very small values
        if(!is.finite(sum(exp(logprobs)/sum(exp(logprobs))))){
          stbfctr <- -500-max(logprobs) ## max largest value -500
          logprobs <- logprobs+stbfctr ## multiply all probs by common factor
        }

        ## sample 1 of C with correct probabilities
        newZtheta <- tryCatch(sample(1:const$Ctheta,1,prob=exp(logprobs)/sum(exp(logprobs))),
                 error=function(err){NULL})

        if(!is.null(newZtheta)){
          params$Ztheta[kk,jj] <- newZtheta
        }else{## track errors
          params$err <- 1
        }

      }

    }

  }
  return(params)
}




#' grid search for V
#' @keywords internal
update_V <- function(params,const){

  for(cc in 1:(const$C-1)){

    ## sample V_c^beta
    if(const$clustering=="both" | const$clustering=="beta"){
      logprobs <- sapply(const$grid,function(vv){get_Vlogdensity(vv,cc,params,Vbeta=TRUE,const) }) ## looping over grid, compute density
      ## error handling for very small values
      if(sum(exp(logprobs))==0){
        stbfctr <- -500-max(logprobs) ## max largest value -500
        logprobs <- logprobs+stbfctr ## multiply all probs by common factor
      }
      params$Vbeta[cc] <- sample(const$grid,1,prob=exp(logprobs)) ## sample from the grid
    }

    ## sample V_c^theta
    if(const$clustering=="both" | const$clustering=="theta"){
      logprobs <- sapply(const$grid,function(vv){get_Vlogdensity(vv,cc,params,Vbeta=FALSE,const) }) ## looping over grid, compute density
      ## error handling for very small values
      if(sum(exp(logprobs))==0){
        stbfctr <- -500-max(logprobs) ## max largest value -500
        logprobs <- logprobs+stbfctr ## multiply all probs by common factor
      }
      params$Vtheta[cc] <- sample(const$grid,1,prob=exp(logprobs)) ## sample from the grid
    }

  }
  ## update pimat and pistarmat with new V (not needed since get_Vlogdensity computes it anyway)
  params <- compute_pi(params,const)
  return(params)
}


#' MH step for V
#' @keywords internal
update_V_MH <- function(params,const){

  nbeta <- sapply(1:const$C, function(cc){sum(params$Zbeta==cc)})   ## n_c^beta
  ntheta <- sapply(1:const$C, function(cc){sum(params$Ztheta==cc)}) ## n_c^theta

  for(cc in 1:(const$C-1)){

    ## update Vbeta ##
    if(const$clustering=="both" | const$clustering=="beta"){
      prop_params <- params

      ## setting minimum of shape2 param to be 1 for proposals
      s1 <- 1+nbeta[cc]
      s2 <- max(1,prop_params$alpha[1]+sum(nbeta[(cc+1):const$C]))
      prop_params$Vbeta[cc] <- stats::rbeta(1,shape1=s1,shape2=s2 )

      ## compute log-acceptance ratio
      logPostRatio <- get_Vlogdensity(prop_params$Vbeta[cc],cc,prop_params,Vbeta=TRUE,const)-
        get_Vlogdensity(params$Vbeta[cc],cc,params,Vbeta=TRUE,const)

      logPropRatio <- stats::dbeta(prop_params$Vbeta[cc],shape1=s1,shape2=s2 ,log=TRUE)-
        stats::dbeta(params$Vbeta[cc],shape1=s1,shape2=s2 ,log=TRUE)

      logRatio <- logPostRatio-logPropRatio
      if(log(stats::runif(1,0,1)) < logRatio){ ## accept
        ## update pimat and pistarmat with new V (not needed since get_Vlogdensity computes it anyway)
        prop_params <- compute_pi(prop_params,const)
        params <- prop_params
      }
    }



    ## update Vtheta ##
    if(const$clustering=="both" | const$clustering=="theta"){
      prop_params <- params
      ## setting minimum of shape2 param to be 1 for proposals
      s1 <- 1+ntheta[cc]
      s2 <- max(1,prop_params$alpha[2]+sum(ntheta[(cc+1):const$C]))
      prop_params$Vtheta[cc] <- stats::rbeta(1,shape1=s1,shape2=s2 )

      ## compute log-acceptance ratio
      logPostRatio <- get_Vlogdensity(prop_params$Vtheta[cc],cc,prop_params,Vbeta=FALSE,const)-
        get_Vlogdensity(params$Vtheta[cc],cc,params,Vbeta=FALSE,const)

      logPropRatio <- stats::dbeta(prop_params$Vtheta[cc],shape1=s1,shape2=s2 ,log=TRUE)-
        stats::dbeta(params$Vtheta[cc],shape1=s1,shape2=s2 ,log=TRUE)

      logRatio <- logPostRatio-logPropRatio
      if(log(stats::runif(1,0,1)) < logRatio){ ## accept
        ## update pimat and pistarmat with new V (not needed since get_Vlogdensity computes it anyway)
        prop_params <- compute_pi(prop_params,const)
        params <- prop_params
      }
    }


  }

  return(params)
}

#' Update intercepts
#' @keywords internal
update_intercept <- function(params,const){

  for(kk in 1:const$K){
    params$b0[kk] <- stats::rnorm(1,
                       mean=sum((const$y[const$k_index==kk]-apply(get_B_beta_k(params,const,kk),1,sum)-params$xi*sqrt(params$sigma2[kk])*params$u-const$Zcovariates%*%params$betaZk[kk,]))/(const$n),
                       sd=sqrt(params$sigma2[kk]/(const$n)))
  }

  return(params)
}


#' Update betaZk
#' @keywords internal
update_betaZk <- function(params,const){

  ZTZinv <- solve(t(const$Zcovariates)%*%const$Zcovariates)

  for(kk in 1:const$K){
    params$betaZk[kk,] <- mvtnorm::rmvnorm(1,
                           mean=ZTZinv%*%t(const$Zcovariates)%*%(const$y[const$k_index==kk]-params$b0[kk]-apply(get_B_beta_k(params,const,kk),1,sum)-params$xi*sqrt(params$sigma2[kk])*params$u),
                           sigma=params$sigma2[kk]*ZTZinv)
  }

  return(params)
}


#' Update betastar
#' @keywords internal
update_betastar <- function(params,const){

  ## computing only once
  Btheta <- lapply(1:const$p,function(jj){
    Reduce("rbind",lapply(1:const$K,function(kk){
      get_Btheta(const$X[[jj]]%*%params$omegastar[(params$Ztheta[kk,jj]-1)*const$L+(1:const$L)],const,params,kk,jj)
    }))
  })


  for(cc in 1:const$Cbeta){
    n_c <- sum(params$Zbeta==cc)


    if(n_c>0){
      if(const$sharedlambda==TRUE){
        lambda_beta <- params$lambda_beta
      }else{
        lambda_beta <- params$lambda_beta[cc]
      }

      whichZ <- which(params$Zbeta==cc,arr.ind=TRUE)
      whichk <- sort(unique(whichZ[,1]))
      whichkj <- lapply(1:const$K,function(kk){sort(whichZ[whichZ[,1]==kk,2])})
      whichkNotj <- lapply(1:const$K,function(kk){which(!(1:const$p)%in%whichkj[[kk]])  })

      ## Btheta for relevant k,j pairs
      B_kc <- lapply(whichk,function(kk){
        Reduce("+",lapply(whichkj[[kk]],function(jj){
          Btheta[[jj]][const$k_index==kk,]/sqrt(params$sigma2[kk]) ## sqrt since it gets squared in BTB
        }))
      })

      ## sum of B^TB across relevant k
      BTB <- Reduce("+",lapply(B_kc,function(BB){t(BB)%*%BB}))

      ##
      y_u_B_k <- lapply(whichk,function(kk){
          y_u <- const$y[const$k_index==kk]-params$b0[kk]-params$xi*sqrt(params$sigma2[kk])*params$u-const$Zcovariates%*%params$betaZk[kk,]
            if(length(whichkNotj[[kk]])>0){
              y_u <- y_u - Reduce("+",lapply(whichkNotj[[kk]],function(jj){
                Btheta[[jj]][const$k_index==kk,,drop=F]%*%params$betastar[(params$Zbeta[kk,jj]-1)*const$d+(1:const$d),drop=F]
              }))
            }
          return(y_u/sqrt(params$sigma2[kk])) ## sqrt because it multiples B
      })

      ##
      yTB <- Reduce("+",lapply(1:length(whichk),function(kk){
        t(y_u_B_k[[kk]])%*%B_kc[[kk]]
      }))

      ## compute Vmat only once ## summing over all k in cluster cc
      Vmat <- solve(lambda_beta*const$invSig0+BTB)
      Vmat <- (Vmat+t(Vmat))/2

      params$betastar[(cc-1)*const$d+(1:const$d)] <- mvtnorm::rmvnorm(n=1,
                                                             mean=Vmat%*%t(lambda_beta*t(const$mu0)%*%const$invSig0+yTB  ),
                                                             sigma=Vmat)


    }else{ ## if n_c=0, draw from the prior
      if(const$sharedlambda==TRUE){
        lambda_beta <- params$lambda_beta
      }else{
        lambda_beta <- params$lambda_beta[cc]
      }
      params$betastar[(cc-1)*const$d+(1:const$d)] <- mvtnorm::rmvnorm(n=1,mean=const$mu0, sigma=(1/lambda_beta)*MASS::ginv(const$invSig0) )

    }

  }

  return(params)
}


#' Update betastar for NONSEP parameterization
#' @keywords internal
update_betastar_NONSEP <- function(params,const){

  ## computing only once
  Btheta <- lapply(1:const$p,function(jj){
    Reduce("rbind",lapply(1:const$K,function(kk){
      get_Btheta(const$X[[jj]],const,params,kk,jj)
    }))
  })


  for(cc in 1:const$Cbeta){
    n_c <- sum(params$Zbeta==cc)


    if(n_c>0){
      if(const$sharedlambda==TRUE){
        lambda_beta <- params$lambda_beta
      }else{
        lambda_beta <- params$lambda_beta[cc]
      }

      whichZ <- which(params$Zbeta==cc,arr.ind=TRUE)
      whichk <- sort(unique(whichZ[,1]))
      whichkj <- lapply(1:const$K,function(kk){sort(whichZ[whichZ[,1]==kk,2])})
      whichkNotj <- lapply(1:const$K,function(kk){which(!(1:const$p)%in%whichkj[[kk]])  })

      ## Btheta for relevant k,j pairs
      B_kc <- lapply(whichk,function(kk){
        Reduce("+",lapply(whichkj[[kk]],function(jj){
          Btheta[[jj]][const$k_index==kk,]/sqrt(params$sigma2[kk]) ## sqrt since it gets squared in BTB
        }))
      })

      ## sum of B^TB across relevant k
      BTB <- Reduce("+",lapply(B_kc,function(BB){t(BB)%*%BB}))

      ##
      y_u_B_k <- lapply(whichk,function(kk){
        y_u <- const$y[const$k_index==kk]-params$b0[kk]-params$xi*sqrt(params$sigma2[kk])*params$u-const$Zcovariates%*%params$betaZk[kk,]
        if(length(whichkNotj[[kk]])>0){
          y_u <- y_u - Reduce("+",lapply(whichkNotj[[kk]],function(jj){
            Btheta[[jj]][const$k_index==kk,,drop=F]%*%params$betastar[(params$Zbeta[kk,jj]-1)*const$d+(1:const$d),drop=F]
          }))
        }
        return(y_u/sqrt(params$sigma2[kk])) ## sqrt because it multiples B
      })

      ##
      yTB <- Reduce("+",lapply(1:length(whichk),function(kk){
        t(y_u_B_k[[kk]])%*%B_kc[[kk]]
      }))

      ## compute Vmat only once ## summing over all k in cluster cc
      Pmat <- lambda_beta*(const$invSig0+const$PEN)
      Vmat <- MASS::ginv(Pmat+BTB)
      Vmat <- (Vmat+t(Vmat))/2

      params$betastar[(cc-1)*const$d+(1:const$d)] <- mvtnorm::rmvnorm(n=1,
                                                                      mean=Vmat%*%t(t(const$mu0)%*%Pmat+yTB  ),
                                                                      sigma=Vmat)


    }else{ ## if n_c=0, draw from the prior
      if(const$sharedlambda==TRUE){
        lambda_beta <- params$lambda_beta
      }else{
        lambda_beta <- params$lambda_beta[cc]
      }
      Pmat <- lambda_beta*(const$invSig0+const$PEN)
      params$betastar[(cc-1)*const$d+(1:const$d)] <- mvtnorm::rmvnorm(n=1,
                                                                      mean=const$mu0,
                                                                      sigma=MASS::ginv(Pmat) )

    }

  }

  return(params)
}



#' Update thetastar via approx method
#' @keywords internal
update_thetastar <- function(params,const){

  for(cc in 1:const$Ctheta){
    if(sum(params$Ztheta==cc)>0){ ## n_c>0 (otherwise draw from prior)
      whichZ <- which(params$Ztheta==cc,arr.ind=TRUE)
      whichk <- sort(unique(whichZ[,1]))
      whichkj <- lapply(1:const$K,function(kk){sort(whichZ[whichZ[,1]==kk,2])})
      whichkNotj <- lapply(1:const$K,function(kk){which(!(1:const$p)%in%whichkj[[kk]])  })

      ## get components of posterior
      comps <- get_XTyhat(cc,whichk,whichkj,params,const)
      XTX <- comps$XTX
      XTy <- comps$XTy

      thetadraw <-tryCatch({
        c(simdd::rFisherBingham(1,
                       mu = const$prior_tau_theta*as.matrix(rep(1,const$Lq)) + XTy,#XtWkyhat,
                       Aplus = -0.5*XTX#XtWX
                       -0.5*exp(params$loglambda_theta)*const$PEN,
                       mtop=const$rfbtries))},
        error=function(err){ ## if rFB fails
          # print(paste0("Scaled normal approximation for omegastar in cluster ",cc))
          PRECmat <- XTX+exp(params$loglambda_theta)*const$PEN
          PRECmat <- (PRECmat+t(PRECmat))/2
          Vmat <- solve(PRECmat)
          Vmat <- (Vmat+t(Vmat))/2
          approxthetadraw <-c(mvtnorm::rmvnorm(1,
                                               mean=Vmat%*%(const$prior_tau_theta*as.matrix(rep(1,const$Lq)) + XTy),#mean=solve((1/params$sigma2)*XtWX+exp(params$loglambda_theta)*const$PEN)%*%(const$prior_tau_theta*as.matrix(rep(1,const$L)) + (1/params$sigma2)*XtWkyhat),
                                               sigma=Vmat))
          return(approxthetadraw/sqrt(sum(approxthetadraw^2)))
        }
      )
      params$thetastar[(cc-1)*const$Lq+(1:const$Lq)] <- thetadraw
      ## basis transformation to get back to omega scale
      params$omegastar[(cc-1)*const$L+(1:const$L)] <- c(const$Psi%*%thetadraw)

    }else{## if n_c=0, draw from the prior
      thetadraw <- c(simdd::rFisherBingham(1,
                                  mu = const$prior_tau_theta*rep(1,const$Lq),
                                  Aplus = -0.5*exp(params$loglambda_theta)*const$PEN))
      params$thetastar[(cc-1)*const$Lq+(1:const$Lq)] <- thetadraw
      ## basis transformation to get back to omega scale
      params$omegastar[(cc-1)*const$L+(1:const$L)] <- c(const$Psi%*%thetadraw)
    }

  }

  return(params)
}

#' Update thetastar via MH step on polar coordinates
#' @keywords internal
update_thetastar_MH_beta <- function(params, const){

  for(cc in 1:const$Cbeta){
    if(sum(params$Ztheta==cc)>0){ ## n_c>0 (otherwise draw from prior)

      whichZ <- which(params$Ztheta==cc,arr.ind=TRUE)
      whichk <- sort(unique(whichZ[,1]))
      whichkj <- lapply(1:const$K,function(kk){sort(whichZ[whichZ[,1]==kk,2])})
      whichkNotj <- lapply(1:const$K,function(kk){which(!(1:const$p)%in%whichkj[[kk]])  })

      for(jj in sample(1:(const$Lq-1))){
        prop_params <- params

        phibeta <- (prop_params$phistar[(cc-1)*(const$Lq-1)+jj]+pi/2)/pi

        ## propose new phi_j from pi*beta distribution (using previous value as mode)
        phibetamode <- phibeta
        b_mode <- ((1-(phibetamode))*const$prop_phi_a+2*(phibetamode)-1)/(phibetamode) ## using previous value as mode
        prop_phibeta <- stats::rbeta(1,const$prop_phi_a,b_mode)
        b_mode_reverse <- ((1-(prop_phibeta))*const$prop_phi_a+2*(prop_phibeta)-1)/(prop_phibeta) ## using previous value as mode

        prop_params$phistar[(cc-1)*(const$Lq-1)+jj] <- pi*prop_phibeta-(pi/2)

        ## compute corresponding omegastar proposal and Btheta
        ## basis transformation to get back to omega scale
        thetadraw <- get_theta(prop_params$phistar[(cc-1)*(const$Lq-1)+(1:(const$Lq-1))])
        prop_params$thetastar[(cc-1)*const$Lq+(1:const$Lq)] <- thetadraw
        prop_params$omegastar[(cc-1)*const$L+(1:const$L)] <- c(const$Psi%*%thetadraw)

        ## components of likelihood
        prop_y_B_u2 <- sapply(whichk,function(kk){sum((const$y[const$k_index==kk]-prop_params$b0[kk]-apply(get_B_beta_k(prop_params,const,kk),1,sum)-prop_params$xi*sqrt(prop_params$sigma2[kk])*prop_params$u-const$Zcovariates%*%prop_params$betaZk[kk,])^2)/prop_params$sigma2[kk]})
        y_B_u2 <- sapply(whichk,function(kk){sum((const$y[const$k_index==kk]-params$b0[kk]-apply(get_B_beta_k(params,const,kk),1,sum)-params$xi*sqrt(params$sigma2[kk])*params$u-const$Zcovariates%*%params$betaZk[kk,])^2)/params$sigma2[kk]})

        ## calculate acceptance ratio
        logLikRatio <- (-0.5*sum(prop_y_B_u2)) -
          (-0.5*sum(y_B_u2))

        logPriorRatio <- (const$prior_tau_theta*rep(1,const$Lq)%*%prop_params$thetastar[(cc-1)*const$Lq+(1:const$Lq)]-0.5*exp(prop_params$loglambda_theta)*c(t(prop_params$thetastar[(cc-1)*const$Lq+(1:const$Lq)])%*%const$PEN%*%prop_params$thetastar[(cc-1)*const$Lq+(1:const$Lq)])) -
          (const$prior_tau_theta*rep(1,const$Lq)%*%params$thetastar[(cc-1)*const$Lq+(1:const$Lq)]-0.5*exp(params$loglambda_theta)*c(t(params$thetastar[(cc-1)*const$Lq+(1:const$Lq)])%*%const$PEN%*%params$thetastar[(cc-1)*const$Lq+(1:const$Lq)]))+
          ## from change of variable
          (const$Lq-jj)*log(abs(cos(prop_params$phistar[(cc-1)*(const$Lq-1)+jj]))) -  ## log(abs(cos(prop_phibeta)^(const$Lq-jj))) -
          (const$Lq-jj)*log(abs(cos(params$phistar[(cc-1)*(const$Lq-1)+jj])))  ## log(abs(cos(phibeta)^(const$Lq-jj)))

        ## compare on beta scales
        logPropRatio <- stats::dbeta(prop_phibeta,const$prop_phi_a,b_mode,log=TRUE) -
          stats::dbeta(phibeta,const$prop_phi_a,b_mode_reverse,log=TRUE)

        logRatio <- logLikRatio+logPriorRatio-logPropRatio

        if(log(stats::runif(1,0,1)) < logRatio){ ## accept
          params <- prop_params
        }

      }

    }else{## if n_c=0, draw from the prior
      thetadraw <- c(simdd::rFisherBingham(1,
                                  mu = const$prior_tau_theta*rep(1,const$Lq),
                                  Aplus = -0.5*exp(params$loglambda_theta)*const$PEN))
      ## fix later: change of variable on vMF (or set tau=0)
      if(thetadraw[length(thetadraw)]<0){
        thetadraw <- -thetadraw
      }
      params$thetastar[(cc-1)*const$Lq+(1:const$Lq)] <- thetadraw
      params$omegastar[(cc-1)*const$L+(1:const$L)] <- c(const$Psi%*%thetadraw)
      params$phistar[(cc-1)*(const$Lq-1)+(1:(const$Lq-1))] <- get_phi(params$thetastar[(cc-1)*const$Lq+(1:const$Lq)])
    }

  }

  return(params)
}

#' Update alpha
#' @keywords internal
update_alpha <- function(params,const){


  if(const$clustering=="both" | const$clustering=="beta"){
    params$alpha[1] <- stats::rgamma(1,shape=const$prior_alpha_beta[1]+const$Cbeta-1,rate=const$prior_alpha_beta[2]-sum(log(1-params$Vbeta[-const$Cbeta])))
  }
  if(const$clustering=="both" | const$clustering=="theta"){
    params$alpha[2] <- stats::rgamma(1,shape=const$prior_alpha_theta[1]+const$Ctheta-1,rate=const$prior_alpha_theta[2]-sum(log(1-params$Vtheta[-const$Ctheta])))
  }


  return(params) ## alpha_beta then alpha_theta
}

#' Update logrho
#' @keywords internal
update_logrho <- function(params,const){

  ## proposal
  prop_params <- params
  prop_params$logrho <- params$logrho+stats::rnorm(1,0,const$stepsize_logrho)
  prop_params <- compute_pi(prop_params,const) ## update proposed pistar and pi accordingly

  ## compute log-acceptance ratio
  logLikRatio <- sum(sapply(1:const$p,function(jj){sapply(1:const$K,function(kk){log(prop_params$pistarmat[params$Zbeta[kk,jj],params$Ztheta[kk,jj]])})}))-
    sum(sapply(1:const$p,function(jj){sapply(1:const$K,function(kk){log(params$pistarmat[params$Zbeta[kk],params$Ztheta[kk]])})}))

  logPriorRatio <- stats::dgamma(exp(prop_params$logrho),shape=const$prior_rho[1],rate=const$prior_rho[2],log=TRUE)-
    stats::dgamma(exp(params$logrho),shape=const$prior_rho[1],rate=const$prior_rho[2],log=TRUE)+
    prop_params$logrho-params$logrho ## from change of variable

  logRatio <- logLikRatio+logPriorRatio-0
  if(log(stats::runif(1,0,1)) < logRatio){ ## accept
    params <- prop_params
  }

  return(params)
}


#' Update smoothness penalty on beta
#' @keywords internal
update_lambda_beta <- function(params,const){

  if(const$sharedlambda==TRUE){
    gamrate <- const$prior_lambda_beta[2]+0.5*sum(sapply(1:const$Cbeta, ## summing over all clusters cc
                                                         function(cc){c(t(params$betastar[(cc-1)*const$d+(1:const$d)])%*%const$invSig0%*%params$betastar[(cc-1)*const$d+(1:const$d)])}))
    if(is.finite(gamrate)){
      params$lambda_beta <- stats::rgamma(1,
                                   shape=const$prior_lambda_beta[1]+0.5*const$Cbeta*const$d,
                                   rate=gamrate)
    }
  }else{
    for(cc in 1:const$Cbeta){

      params$lambda_beta[cc] <- stats::rgamma(1,
                                       shape=const$prior_lambda_beta[1]+0.5*const$d,
                                       rate=const$prior_lambda_beta[2]+0.5*c(t(params$betastar[(cc-1)*const$d+(1:const$d)])%*%const$invSig0%*%params$betastar[(cc-1)*const$d+(1:const$d)]))
    }

  }


  return(params)
}


#' Update smoothness penalty on theta
#' @keywords internal
update_loglambda_theta <- function(params,const){

  ## proposal
  prop_params <- params
  prop_params$loglambda_theta <- params$loglambda_theta+stats::rnorm(1,0,const$stepsize_loglambda_theta)

  ## compute log-acceptance ratio
  logLikRatio <- (-0.5*exp(prop_params$loglambda_theta)*sum(sapply(1:const$Ctheta, ## summing over all clusters cc
                                                                  function(cc){c(t(prop_params$thetastar[(cc-1)*const$Lq+(1:const$Lq)])%*%const$PEN%*%prop_params$thetastar[(cc-1)*const$Lq+(1:const$Lq)])}))
                        - const$Ctheta*(fb.saddle(gam=const$prior_tau_theta*rep(1,const$Lq),lam=0.5*exp(prop_params$loglambda_theta)*eigen(const$PEN)$values)[3]) ) -  ## third element is third order approximation
    (-0.5*exp(params$loglambda_theta)*sum(sapply(1:const$Ctheta, ## summing over all clusters cc
                                                function(cc){c(t(params$thetastar[(cc-1)*const$Lq+(1:const$Lq)])%*%const$PEN%*%params$thetastar[(cc-1)*const$Lq+(1:const$Lq)])}))
            - const$Ctheta*(fb.saddle(gam=const$prior_tau_theta*rep(1,const$Lq),lam=0.5*exp(params$loglambda_theta)*eigen(const$PEN)$values)[3]) ) ## third element is third order approximation
  ## NOTE: fb.saddle uses the parameterization form Kume and wood 2005; rFisherBingham above uses a different one. Hence we have no negative in the lam term here; but we had negative in the Aplus term above


  logPriorRatio <- stats::dgamma(exp(prop_params$loglambda_theta),shape=const$prior_lambda_theta[1],rate=const$prior_lambda_theta[2],log=TRUE)-
    stats::dgamma(exp(params$loglambda_theta),shape=const$prior_lambda_theta[1],rate=const$prior_lambda_theta[2],log=TRUE)+
    prop_params$loglambda_theta-params$loglambda_theta ## from change of variable

  logRatio <- logLikRatio+logPriorRatio-0
  tryCatch({if(log(stats::runif(1,0,1)) < logRatio){ ## accept
                params <- prop_params}
            },
           error=function(err){print("Skipping loglambda_theta")})


  return(params)
}


#' Update both smoothness penalties for NONSEP parameterization
#' @keywords internal
update_lambda_NONSEP <- function(params,const){

  Pmat <- const$invSig0+const$PEN ## elsewhere this includes the penalty params (for later use possibly); ## here contains only constant penalty matrix



  if(const$sharedlambda==TRUE){
    gamrate <- const$prior_lambda_beta[2]+0.5*sum(sapply(1:const$Cbeta, ## summing over all clusters cc
                                                         function(cc){c(t(params$betastar[(cc-1)*const$d+(1:const$d)])%*%Pmat%*%params$betastar[(cc-1)*const$d+(1:const$d)])}))
    if(is.finite(gamrate)){
      params$lambda_beta <- stats::rgamma(1,
                                          shape=const$prior_lambda_beta[1]+0.5*const$Cbeta*Matrix::rankMatrix(Pmat),
                                          rate=gamrate)

    }
  }else{
    for(cc in 1:const$Cbeta){

      params$lambda_beta[cc] <- stats::rgamma(1,
                                              shape=const$prior_lambda_beta[1]+0.5*Matrix::rankMatrix(Pmat),
                                              rate=const$prior_lambda_beta[2]+0.5*c(t(params$betastar[(cc-1)*const$d+(1:const$d)])%*%Pmat%*%params$betastar[(cc-1)*const$d+(1:const$d)]))
    }

  }



  return(params)
}



#' Update u (ranef)
#' @keywords internal
update_u <- function(params,const){

  ## first compute B times the corresponding beta for all obs
  sumB_beta <- apply(get_B_beta(params,const),1,sum)

  for(ii in 1:const$n){## loop over all subjects
    params$u[ii] <- stats::rnorm(1,
                          mean=(1/(1+(const$K*(params$xi^2))))*sum((params$xi/sqrt(params$sigma2))*(const$y[const$i_index==ii]-params$b0-sumB_beta[const$i_index==ii]-t(const$Zcovariates[ii,]%*%t(params$betaZk)))), ## DEBUG: previously had xi*sigma, not xi/sigma
                          sd=sqrt(1/(1+(const$K*(params$xi^2))))  ) ## DEBUG:  previously had sum(params$xi^2) instead of K*params$xi^2
  }

  ## if forcing orthogonality
  if(const$condkrig==TRUE){

    ## first construct Btheta1 (design mat from first outcome)
    if(const$NONSEP==FALSE & const$LM==FALSE){
      if(const$fixorthog==TRUE){ ## constant theta for orthogonalization
        fixtheta <- c(rep(1,const$L)/sqrt(const$L))
        B1 <- Reduce("cbind",lapply(1:const$p,function(jj){
          get_Btheta(const$X[[jj]]%*%fixtheta,const,params,1,jj)
        }))
      }else{
        B1 <- Reduce("cbind",lapply(1:const$p,function(jj){
          get_Btheta(const$X[[jj]]%*%params$omegastar[(params$Ztheta[1,jj]-1)*const$L+(1:const$L)],const,params,1,jj)
        }))
      }
    }else{ ## nonseparable version is exact (no theta)
      B1 <- Reduce("cbind",lapply(1:const$p,function(jj){
        get_Btheta(const$X[[jj]],const,params,1,jj)
      }))
    }

    params$u <- (diag(const$n)-B1%*%ginv(t(B1)%*%B1)%*%t(B1))%*%params$u ### CONDITIONING BY KRIGING
  }

  return(params)
}

#' #' Update u (ranef) when fixorthog=TRUE (joint update)
#' #' @keywords internal
#' update_u_fixorthog <- function(params,const){
#'
#'   ## mean vector (except for the SigmaU which needs to be premultiplied)
#'   mu <- Reduce("+",lapply(1:const$K,function(kk){
#'     return((params$xi/sqrt(params$sigma2[kk]))*(const$y[const$k_index==kk]-params$b0[kk]-apply(get_B_beta_k(params,const,kk),1,sum)-const$Zcovariates%*%params$betaZk[kk,]))
#'   }))
#'
#'   ## variance of U from woodbury formula, with orthogonalized/RSR parameterization
#'   SigmaU <- (1/(const$K*params$xi^2))*diag(const$n)+(1/(const$K*params$xi^2))*(1/(const$K*params$xi^2))*const$fixPmat
#'
#'   ## draw u
#'   params$u <- c(mvtnorm::rmvnorm(n=1,
#'                                  mean=SigmaU%*%mu,
#'                                  sigma=SigmaU))
#'
#'
#'   return(params)
#' }



#' Update xi
#' @keywords internal
update_xi <- function(params,const){

  ## first compute B times the corresponding beta for all obs
  sumB_beta <- apply(get_B_beta(params,const),1,sum)

  ## random walk proposal (on log scale)
  prop_params <- params
  prop_params$xi <- exp(log(params$xi)+stats::rnorm(1,0,const$stepsize_logxi))

  ## compute log-acceptance ratio
  logLikRatio <-  (-0.5*sum(sapply(1:const$K,function(kk) {(1/prop_params$sigma2[kk])*sum((const$y[const$k_index==kk]-prop_params$b0[kk]-sumB_beta[const$k_index==kk]-prop_params$xi*sqrt(prop_params$sigma2[kk])*prop_params$u-const$Zcovariates%*%prop_params$betaZk[kk,])^2)})))   - # prop
    (-0.5*sum(sapply(1:const$K,function(kk) {(1/params$sigma2[kk])*sum((const$y[const$k_index==kk]-params$b0[kk]-sumB_beta[const$k_index==kk]-params$xi*sqrt(params$sigma2[kk])*params$u-const$Zcovariates%*%params$betaZk[kk,])^2)})))    ## current

  logPriorRatio <- ((-const$prior_xi[1]-1)*log(prop_params$xi) -(const$prior_xi[2]/prop_params$xi))- # log of inverse gamma density for proposal
    ((-const$prior_xi[1]-1)*log(params$xi) -(const$prior_xi[2]/params$xi))+  # log of inverse gamma density for current
    log(prop_params$xi)-log(params$xi) ## from change of variable

  logRatio <- logLikRatio+logPriorRatio-0
  if(log(stats::runif(1,0,1)) < logRatio){ ## accept
    params <- prop_params
  }

  return(params)
}


## old version (where var(u)=sigma2_u)
#
# update_u <- function(params,const){
#
#   ## first compute B times the corresponding beta for all obs
#   sumB_beta <- apply(get_B_beta(params,const),1,sum)
#
#   for(ii in 1:const$n){## loop over all subjects
#      params$u[ii] <- stats::rnorm(1,
#                           mean=(params$sigma2_u/(params$sigma2+const$K*params$sigma2_u))*sum((const$y[const$i_index==ii]-params$b0-sumB_beta[const$i_index==ii])),
#                           sd=sqrt(params$sigma2_u*params$sigma2/(params$sigma2+const$K*params$sigma2_u))  )
#   }
#   return(params)
# }
#
# update_sigma2_u <- function(params,const){
#
#   params$sigma2_u <- 1/stats::rgamma(1,
#                             shape=const$prior_sigma2_u[1]+0.5*const$n,
#                             rate=const$prior_sigma2_u[2]+0.5*sum(params$u^2) )
#   return(params)
# }




#' Update sigma2 (from conjugate when K=1)
#' @keywords internal
update_sigma2 <- function(params,const){

  ## first compute B times the corresponding beta for all obs
  sumB_beta <- apply(get_B_beta(params,const),1,sum)

  params$sigma2 <- 1/stats::rgamma(1,
                            shape=const$prior_sigma2[1]+0.5*const$n*const$K,
                            rate=const$prior_sigma2[2]+0.5*(sum((const$y-(rep(params$b0,each=const$n)+sumB_beta+rep(params$u,const$K)+const$Zcovariates%*%params$betaZk[1,]))^2)) )

  return(params)
}


#' Update variances (via MH when K>1)
#' @keywords internal
update_sigma2_k <- function(params,const){

  ## first compute B times the corresponding beta for all obs
  sumB_beta <- apply(get_B_beta(params,const),1,sum)

  for(kk in 1:const$K){## loop over all outcomes

    ## random walk proposal (on log scale)
    prop_params <- params
    prop_params$sigma2[kk] <- exp(log(params$sigma2[kk])+stats::rnorm(1,0,const$stepsize_logsigma2))

    ## compute log-acceptance ratio
    logLikRatio <-  (-0.5*const$n*log(prop_params$sigma2[kk])-0.5*(1/prop_params$sigma2[kk])*sum((const$y[const$k_index==kk]-prop_params$b0[kk]-sumB_beta[const$k_index==kk]-prop_params$xi*sqrt(prop_params$sigma2[kk])*prop_params$u-const$Zcovariates%*%prop_params$betaZk[kk,])^2))   - # prop
      (-0.5*const$n*log(params$sigma2[kk])-0.5*(1/params$sigma2[kk])*sum((const$y[const$k_index==kk]-params$b0[kk]-sumB_beta[const$k_index==kk]-params$xi*sqrt(params$sigma2[kk])*params$u-const$Zcovariates%*%params$betaZk[kk,])^2))    ## current

    logPriorRatio <- ((-const$prior_sigma2[1]-1)*log(prop_params$sigma2[kk]) -(const$prior_sigma2[2]/prop_params$sigma2[kk]))- # log of inverse gamma density for proposal
      ((-const$prior_sigma2[1]-1)*log(params$sigma2[kk]) -(const$prior_sigma2[2]/params$sigma2[kk]))+  # log of inverse gamma density for current
       log(prop_params$sigma2[kk])-log(params$sigma2[kk]) ## from change of variable

    logRatio <- logLikRatio+logPriorRatio-0
    if(log(stats::runif(1,0,1)) < logRatio){ ## accept
      params <- prop_params
    }

  }

  return(params)
}






