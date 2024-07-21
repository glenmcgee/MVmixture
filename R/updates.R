## update functions


##
update_clustMemb <- function(params,const){

  ## used for the components that are not being summed over
  B_beta <- get_B_beta(params,const)

  for(kk in 1:const$K){## loop over outcomes
    for(jj in 1:const$p){## loop over exposures

      ## components that get re-used repeatedly
      y_B_u <- const$y[const$k_index==kk]-params$b0[kk]-(apply(B_beta[const$k_index==kk,-jj,drop=F],1,sum) +params$u)
      Bth_kj <- get_Btheta(const$X[[jj]]%*%params$omegastar[(params$Ztheta[kk,jj]-1)*const$L+(1:const$L)],const,params,kk,jj)

      ## Zbetakj ##
      if(const$clustering=="both" | const$clustering=="beta"){

        ## compute probabilities for all possible a
        probs <- c(sapply(1:const$Cbeta, function(a){  ## loop over a (rows; beta clusters)
          exp(log(params$pimat[a,params$Ztheta[kk,jj]]) -(0.5/params$sigma2)*sum((y_B_u-Bth_kj%*%params$betastar[(a-1)*const$d+(1:const$d)])^2))
        })) ## to be standardized below

        ## sample 1 of C with correct probabilities
        tryCatch({params$Zbeta[kk,jj] <- sample(1:const$Cbeta,1,prob=probs/sum(probs))}, ## standardized probs
                 error=function(err){print(paste0("Skipping cluster member update for (k,j)=",kk,",",jj))})
      }


      ## Zthetakj ##
      if(const$clustering=="both" | const$clustering=="theta"){
        ## compute probabilities for all possible b
        probs <- c(sapply(1:const$Ctheta, function(b){  ## loop over a (rows; beta clusters)
          exp(log(params$pimat[params$Zbeta[kk,jj],b]) -(0.5/params$sigma2)*sum((y_B_u- get_Btheta(const$X[[jj]]%*%params$omegastar[(b-1)*const$L+(1:const$L)],const,params,kk,jj)%*%params$betastar[(params$Zbeta[kk,jj]-1)*const$d+(1:const$d)])^2))
        }))## to be standardized below

        ## sample 1 of C with correct probabilities
        tryCatch({params$Ztheta[kk,jj] <- sample(1:const$Ctheta,1,prob=probs/sum(probs))},
                 error=function(err){print(paste0("Skipping cluster member update for (k,j)=",kk,",",jj))})
      }

    }

  }
  return(params)
}




## grid search for V
update_V <- function(params,const){

  for(cc in 1:(const$C-1)){

    ## sample V_c^beta
    if(const$clustering=="both" | const$clustering=="beta"){
      probs <- exp(sapply(const$grid,function(vv){get_Vlogdensity(vv,cc,params,Vbeta=TRUE,const) })) ## looping over grid, compute density
      params$Vbeta[cc] <- sample(const$grid,1,prob=probs) ## sample from the grid
    }

    ## sample V_c^theta
    if(const$clustering=="both" | const$clustering=="theta"){
      probs <- exp(sapply(const$grid,function(vv){get_Vlogdensity(vv,cc,params,Vbeta=FALSE,const) })) ## looping over grid, compute density
      params$Vtheta[cc] <- sample(const$grid,1,prob=probs) ## sample from the grid
    }

  }
  ## update pimat and pistarmat with new V (not needed since get_Vlogdensity computes it anyway)
  params <- compute_pi(params,const)
  return(params)
}

## MH step for V
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
      prop_params$Vbeta[cc] <- rbeta(1,shape1=s1,shape2=s2 )

      ## compute log-acceptance ratio
      logPostRatio <- get_Vlogdensity(prop_params$Vbeta[cc],cc,prop_params,Vbeta=TRUE,const)-
        get_Vlogdensity(params$Vbeta[cc],cc,params,Vbeta=TRUE,const)

      logPropRatio <- dbeta(prop_params$Vbeta[cc],shape1=s1,shape2=s2 ,log=TRUE)-
        dbeta(params$Vbeta[cc],shape1=s1,shape2=s2 ,log=TRUE)

      logRatio <- logPostRatio-logPropRatio
      if(log(runif(1,0,1)) < logRatio){ ## accept
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
      prop_params$Vtheta[cc] <- rbeta(1,shape1=s1,shape2=s2 )

      ## compute log-acceptance ratio
      logPostRatio <- get_Vlogdensity(prop_params$Vtheta[cc],cc,prop_params,Vbeta=FALSE,const)-
        get_Vlogdensity(params$Vtheta[cc],cc,params,Vbeta=FALSE,const)

      logPropRatio <- dbeta(prop_params$Vtheta[cc],shape1=s1,shape2=s2 ,log=TRUE)-
        dbeta(params$Vtheta[cc],shape1=s1,shape2=s2 ,log=TRUE)

      logRatio <- logPostRatio-logPropRatio
      if(log(runif(1,0,1)) < logRatio){ ## accept
        ## update pimat and pistarmat with new V (not needed since get_Vlogdensity computes it anyway)
        prop_params <- compute_pi(prop_params,const)
        params <- prop_params
      }
    }


  }

  return(params)
}

update_intercept <- function(params,const){

  for(kk in 1:const$K){
    params$b0[kk] <- rnorm(1,
                       mean=sum((const$y[const$k_index==kk]-apply(get_B_beta_k(params,const,kk),1,sum)-params$u))/(const$n),
                       sd=sqrt(params$sigma2/(const$n)))
  }

  return(params)
}


update_betastar <- function(params,const){
  ## computing only once
  Btheta <- lapply(1:const$p,function(jj){
    Reduce("rbind",lapply(1:const$K,function(kk){
      get_Btheta(const$X[[jj]]%*%params$omegastar[(params$Ztheta[kk,jj]-1)*const$L+(1:const$L)],const,params,kk,jj)
    }))
  })

  if(const$MIM==TRUE){ ## for identifiability product in MIM
    IDprodmat <- matrix(1,nrow=const$K,ncol=const$p)
    for(kk in 1:const$K){
      for(jj in 2:const$p){
        IDprodmat[kk,jj] <- prod(params$Ztheta[kk,1:(jj-1)]!=params$Ztheta[kk,jj])
      }
    }
  }

  for(cc in 1:const$Cbeta){
    if(const$MIM==TRUE){
      n_c <- sum((params$Zbeta==cc)*IDprodmat)
    }else{
      n_c <- sum(params$Zbeta==cc)
    }



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
          Btheta[[jj]][const$k_index==kk,]
        }))
      })

      ## sum of B^TB across relevant k
      BTB <- Reduce("+",lapply(B_kc,function(BB){t(BB)%*%BB}))

      ##
      y_u_B_k <- lapply(whichk,function(kk){
          y_u <- const$y[const$k_index==kk]-params$b0[kk]-params$u
            if(length(whichkNotj[[kk]])>0){
              y_u <- y_u - Reduce("+",lapply(whichkNotj[[kk]],function(jj){
                Btheta[[jj]][const$k_index==kk,]%*%params$betastar[(params$Zbeta[kk,jj]-1)*const$d+(1:const$d)]
              }))
            }
          return(y_u)
      })

      ##
      yTB <- Reduce("+",lapply(1:length(whichk),function(kk){
        t(y_u_B_k[[kk]])%*%B_kc[[kk]]
      }))

      ## compute Vmat only once ## summing over all k in cluster cc
      Vmat <- solve(lambda_beta*const$invSig0+(1/params$sigma2)*BTB)

      params$betastar[(cc-1)*const$d+(1:const$d)] <- mvtnorm::rmvnorm(n=1,
                                                             mean=Vmat%*%t(lambda_beta*t(const$mu0)%*%const$invSig0+(1/params$sigma2)*yTB  ),
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
        c(rFisherBingham(1,
                       mu = const$prior_tau_theta*as.matrix(rep(1,const$Lq)) + (1/params$sigma2)*XTy,#XtWkyhat,
                       Aplus = -(0.5/params$sigma2)*XTX#XtWX
                       -0.5*exp(params$loglambda_theta)*const$PEN,
                       mtop=const$rfbtries))},
        error=function(err){ ## if rFB fails
          # print(paste0("Scaled normal approximation for omegastar in cluster ",cc))
          approxthetadraw <-c(mvtnorm::rmvnorm(1,
                                               mean=solve((1/params$sigma2)*XTX+exp(params$loglambda_theta)*const$PEN)%*%(const$prior_tau_theta*as.matrix(rep(1,const$Lq)) + (1/params$sigma2)*XTy),#mean=solve((1/params$sigma2)*XtWX+exp(params$loglambda_theta)*const$PEN)%*%(const$prior_tau_theta*as.matrix(rep(1,const$L)) + (1/params$sigma2)*XtWkyhat),
                                               sigma=solve((1/params$sigma2)*XTX+exp(params$loglambda_theta)*const$PEN)))
          return(approxthetadraw/sqrt(sum(approxthetadraw^2)))
        }
      )
      params$thetastar[(cc-1)*const$Lq+(1:const$Lq)] <- thetadraw
      ## basis transformation to get back to omega scale
      params$omegastar[(cc-1)*const$L+(1:const$L)] <- c(const$Psi%*%thetadraw)
      # for(kk in whichk){
      #   for(jj in whichkj[[kk]]){
      #     params$Btheta[[jj]][const$k_index==kk,] <- get_Btheta(const$X[[jj]]%*%params$omegastar[(cc-1)*const$L+(1:const$L)],const,params,kk,jj)
      #     params$DerivBtheta[[jj]][const$k_index==kk,] <- get_DerivBtheta(const$X[[jj]]%*%params$omegastar[(cc-1)*const$L+(1:const$L)],const,params,kk,jj)
      #   }
      # }

    }else{## if n_c=0, draw from the prior
      thetadraw <- c(rFisherBingham(1,
                                  mu = const$prior_tau_theta*rep(1,const$Lq),
                                  Aplus = -0.5*exp(params$loglambda_theta)*const$PEN))
      params$thetastar[(cc-1)*const$Lq+(1:const$Lq)] <- thetadraw
      ## basis transformation to get back to omega scale
      params$omegastar[(cc-1)*const$L+(1:const$L)] <- c(const$Psi%*%thetadraw)
    }

  }

  return(params)
}


## do MH on reparameterized index with beta proposals
update_thetastar_MH_beta <- function(params, const){

  for(cc in 1:const$Cbeta){
    if(sum(params$Ztheta==cc)>0){ ## n_c>0 (otherwise draw from prior)

      whichZ <- which(params$Ztheta==cc,arr.ind=TRUE)
      whichk <- sort(unique(whichZ[,1]))
      whichkj <- lapply(1:const$K,function(kk){sort(whichZ[whichZ[,1]==kk,2])})
      whichkNotj <- lapply(1:const$K,function(kk){which(!(1:const$p)%in%whichkj[[kk]])  })

      for(jj in 1:(const$Lq-1)){
        prop_params <- params

        ## phi on beta scale
        if(jj==1){ ## phi_1 defined on [0,pi/2]
          phibeta <- prop_params$phistar[(cc-1)*(const$Lq-1)+jj]/(pi/2)
        }else{ ## phi_j defined on [-pi/2,pi/2]
          phibeta <- (prop_params$phistar[(cc-1)*(const$Lq-1)+jj]+(pi/2))/pi
        }

        ## propose new phi_j from pi*beta distribution (using previous value as mode)
        b_mode <- ((1-(phibeta))*const$prior_phi_a+2*(phibeta)-1)/(phibeta) ## using previous value as mode
        prop_phibeta <- rbeta(1,const$prior_phi_a,b_mode)
        if(jj==1){ ## phi_1 defined on [0,pi/2]
          prop_params$phistar[(cc-1)*(const$Lq-1)+jj] <- (pi/2)*prop_phibeta
        }else{ ## phi_j defined on [-pi/2,pi/2]
          prop_params$phistar[(cc-1)*(const$Lq-1)+jj] <- pi*prop_phibeta-(pi/2)
        }

        ## compute corresponding omegastar proposal and Btheta
        # prop_params$omegastar[(cc-1)*const$L+(1:const$L)] <- get_theta(prop_params$phistar[(cc-1)*(const$L-1)+(1:(const$L-1))])
        ## basis transformation to get back to omega scale
        thetadraw <- get_theta(prop_params$phistar[(cc-1)*(const$Lq-1)+(1:(const$Lq-1))])
        prop_params$thetastar[(cc-1)*const$Lq+(1:const$Lq)] <- thetadraw
        prop_params$omegastar[(cc-1)*const$L+(1:const$L)] <- c(const$Psi%*%thetadraw)
        # for(kk in whichk){
        #   for(jj in whichkj[[kk]]){
        #     prop_params$Btheta[[jj]][const$k_index==kk,] <- get_Btheta(const$X[[jj]]%*%prop_params$omegastar[(cc-1)*const$L+(1:const$L)],const,params,kk,jj)
        #     prop_params$DerivBtheta[[jj]][const$k_index==kk,] <- get_DerivBtheta(const$X[[jj]]%*%prop_params$omegastar[(cc-1)*const$L+(1:const$L)],const,params,kk,jj)
        #   }
        # }

        ## components of likelihood
        prop_y_B_u2 <- sapply(whichk,function(kk){sum((const$y[const$k_index==kk]-prop_params$b0[kk]-apply(get_B_beta_k(prop_params,const,kk),1,sum)-prop_params$u)^2)})
        y_B_u2 <- sapply(whichk,function(kk){sum((const$y[const$k_index==kk]-params$b0[kk]-apply(get_B_beta_k(params,const,kk),1,sum)-params$u)^2)})

        ## calculate acceptance ratio
        logLikRatio <- -0.5*(1/prop_params$sigma2)*sum(prop_y_B_u2) -
          -0.5*(1/params$sigma2)*sum(y_B_u2)

        logPriorRatio <- const$prior_tau_theta*rep(1,const$Lq)%*%prop_params$thetastar[(cc-1)*const$Lq+(1:const$Lq)]-0.5*exp(prop_params$loglambda_theta)*c(t(prop_params$thetastar[(cc-1)*const$Lq+(1:const$Lq)])%*%const$PEN%*%prop_params$thetastar[(cc-1)*const$Lq+(1:const$Lq)]) -
          const$prior_tau_theta*rep(1,const$Lq)%*%params$thetastar[(cc-1)*const$Lq+(1:const$Lq)]-0.5*exp(params$loglambda_theta)*c(t(params$thetastar[(cc-1)*const$Lq+(1:const$Lq)])%*%const$PEN%*%params$thetastar[(cc-1)*const$Lq+(1:const$Lq)])

        ## compare on beta scales
        logPropRatio <- dbeta(prop_phibeta,const$prior_phi_a,b_mode,log=TRUE) -
          dbeta(phibeta,const$prior_phi_a,b_mode,log=TRUE) +
          ## from change of variable
          log(abs(cos(prop_params$phistar[(cc-1)*(const$Lq-1)+jj])^(const$Lq-jj))) -  ## log(abs(cos(prop_phibeta)^(const$Lq-jj))) -
          log(abs(cos(params$phistar[(cc-1)*(const$Lq-1)+jj])^(const$Lq-jj)))  ## log(abs(cos(phibeta)^(const$Lq-jj)))

        logRatio <- logLikRatio+logPriorRatio-logPropRatio
        if(log(runif(1,0,1)) < logRatio){ ## accept
          params <- prop_params
        }

      }

    }else{## if n_c=0, draw from the prior
      thetadraw <- c(rFisherBingham(1,
                                  mu = const$prior_tau_theta*rep(1,const$Lq),
                                  Aplus = -0.5*exp(params$loglambda_theta)*const$PEN))
      ## fix later: change of variable on vMF (or set tau=0)
      if(thetadraw[1]<0){
        thetadraw <- -thetadraw
      }
      params$thetastar[(cc-1)*const$Lq+(1:const$Lq)] <- thetadraw
      params$omegastar[(cc-1)*const$L+(1:const$L)] <- c(const$Psi%*%thetadraw)
      params$phistar[(cc-1)*(const$Lq-1)+(1:(const$Lq-1))] <- get_phi(params$thetastar[(cc-1)*const$Lq+(1:const$Lq)])
    }

  }

  return(params)
}


update_alpha <- function(params,const){


  if(const$clustering=="both" | const$clustering=="beta"){
    params$alpha[1] <- rgamma(1,shape=const$prior_alpha_beta[1]+const$Cbeta-1,rate=const$prior_alpha_beta[2]-sum(log(1-params$Vbeta[-const$Cbeta])))
  }
  if(const$clustering=="both" | const$clustering=="theta"){
    params$alpha[2] <- rgamma(1,shape=const$prior_alpha_theta[1]+const$Ctheta-1,rate=const$prior_alpha_theta[2]-sum(log(1-params$Vtheta[-const$Ctheta])))
  }

  # params$alpha <- c(
  #   rgamma(1,shape=const$prior_alpha_beta[1]+const$C-1,rate=const$prior_alpha_beta[2]-sum(log(1-params$Vbeta[-const$C]))),
  #   rgamma(1,shape=const$prior_alpha_theta[1]+const$C-1,rate=const$prior_alpha_theta[2]-sum(log(1-params$Vtheta[-const$C])))
  # )

  return(params) ## alpha_beta then alpha_theta
}


update_logrho <- function(params,const){

  ## proposal
  prop_params <- params
  prop_params$logrho <- params$logrho+rnorm(1,0,const$stepsize_logrho)
  prop_params <- compute_pi(prop_params,const) ## update proposed pistar and pi accordingly

  ## compute log-acceptance ratio
  logLikRatio <- sum(sapply(1:const$p,function(jj){sapply(1:const$K,function(kk){log(prop_params$pistarmat[params$Zbeta[kk,jj],params$Ztheta[kk,jj]])})}))-
    sum(sapply(1:const$p,function(jj){sapply(1:const$K,function(kk){log(params$pistarmat[params$Zbeta[kk],params$Ztheta[kk]])})}))

  logPriorRatio <- dgamma(exp(prop_params$logrho),shape=const$prior_rho[1],rate=const$prior_rho[2],log=TRUE)-
    dgamma(exp(params$logrho),shape=const$prior_rho[1],rate=const$prior_rho[2],log=TRUE)+
    prop_params$logrho-params$logrho ## from change of variable

  logRatio <- logLikRatio+logPriorRatio-0
  if(log(runif(1,0,1)) < logRatio){ ## accept
    params <- prop_params
  }

  return(params)
}

## smoothness penalty on beta
update_lambda_beta <- function(params,const){

  if(const$sharedlambda==TRUE){
    gamrate <- const$prior_lambda_beta[2]+0.5*sum(sapply(1:const$Cbeta, ## summing over all clusters cc
                                                         function(cc){c(t(params$betastar[(cc-1)*const$d+(1:const$d)])%*%const$invSig0%*%params$betastar[(cc-1)*const$d+(1:const$d)])}))
    if(is.finite(gamrate)){
      params$lambda_beta <- rgamma(1,
                                   shape=const$prior_lambda_beta[1]+0.5*const$Cbeta*const$d,
                                   rate=gamrate)
    }
  }else{
    for(cc in 1:const$Cbeta){

      params$lambda_beta[cc] <- rgamma(1,
                                       shape=const$prior_lambda_beta[1]+0.5*const$d,
                                       rate=const$prior_lambda_beta[2]+0.5*c(t(params$betastar[(cc-1)*const$d+(1:const$d)])%*%const$invSig0%*%params$betastar[(cc-1)*const$d+(1:const$d)]))
    }

  }


  return(params)
}

update_loglambda_theta <- function(params,const){

  ## proposal
  prop_params <- params
  prop_params$loglambda_theta <- params$loglambda_theta+rnorm(1,0,const$stepsize_loglambda_theta)

  ## compute log-acceptance ratio
  logLikRatio <- (-0.5*exp(prop_params$loglambda_theta)*sum(sapply(1:const$Ctheta, ## summing over all clusters cc
                                                                  function(cc){c(t(prop_params$thetastar[(cc-1)*const$Lq+(1:const$Lq)])%*%const$PEN%*%prop_params$thetastar[(cc-1)*const$Lq+(1:const$Lq)])}))
                        - const$Ctheta*(fb.saddle(gam=const$prior_tau_theta*rep(1,const$Lq),lam=-0.5*exp(prop_params$loglambda_theta)*eigen(const$PEN)$values)[3]) ) -  ## third element is third order approximation
    (-0.5*exp(params$loglambda_theta)*sum(sapply(1:const$Ctheta, ## summing over all clusters cc
                                                      function(cc){c(t(params$thetastar[(cc-1)*const$Lq+(1:const$Lq)])%*%const$PEN%*%params$thetastar[(cc-1)*const$Lq+(1:const$Lq)])}))
            - const$Ctheta*(fb.saddle(gam=const$prior_tau_theta*rep(1,const$Lq),lam=-0.5*exp(params$loglambda_theta)*eigen(const$PEN)$values)[3]) ) ## third element is third order approximation

  logPriorRatio <- dgamma(exp(prop_params$loglambda_theta),shape=const$prior_lambda_theta[1],rate=const$prior_lambda_theta[2],log=TRUE)-
    dgamma(exp(params$loglambda_theta),shape=const$prior_lambda_theta[1],rate=const$prior_lambda_theta[2],log=TRUE)+
    prop_params$loglambda_theta-params$loglambda_theta ## from change of variable

  logRatio <- logLikRatio+logPriorRatio-0
  tryCatch({if(log(runif(1,0,1)) < logRatio){ ## accept
                params <- prop_params}
            },
           error=function(err){print("Skipping loglambda_theta")})


  return(params)
}


update_u <- function(params,const){

  ## first compute B times the corresponding beta for all obs
  sumB_beta <- apply(get_B_beta(params,const),1,sum)

  for(ii in 1:const$n){## loop over all subjects
     params$u[ii] <- rnorm(1,
                          mean=(params$sigma2_u/(params$sigma2+const$K*params$sigma2_u))*sum((const$y[const$i_index==ii]-params$b0-sumB_beta[const$i_index==ii])),
                          sd=sqrt(params$sigma2_u*params$sigma2/(params$sigma2+const$K*params$sigma2_u))  )
  }
  return(params)
}


update_sigma2_u <- function(params,const){

  params$sigma2_u <- 1/rgamma(1,
                            shape=const$prior_sigma2_u[1]+0.5*const$n,
                            rate=const$prior_sigma2_u[2]+0.5*sum(params$u^2) )
  return(params)
}


update_sigma2 <- function(params,const){

  ## first compute B times the corresponding beta for all obs
  sumB_beta <- apply(get_B_beta(params,const),1,sum)

  params$sigma2 <- 1/rgamma(1,
                            shape=const$prior_sigma2[1]+0.5*const$n*const$K,
                            rate=const$prior_sigma2[2]+0.5*(sum((const$y-(rep(params$b0,each=const$n)+sumB_beta+rep(params$u,const$K)))^2)) )

  return(params)
}







