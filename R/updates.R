## update functions
## Debug rFisherBingham
##### it sometimes tries 1000s of times, at which point the function stops and errors out
##### maybe 1-2% of the time
##### currently skipping when it errors
## Debug loglambda_theta
##### errors out

#### EDIT: using separate lambda_beta for each cluster now


# ## not currently being used. numerical issues
# update_JOINTclustMemb <- function(params,const){
#
#   for(kk in 1:const$K){
#
#     ## compute probabilities for all possible a,b
#     probs <-sapply(1:const$C,function(b){     ## loop over b (columns; theta clusters)
#               sapply(1:const$C, function(a){  ## loop over a (rows; beta clusters)
#                 # exp(log(params$pimat[a,b]) -(0.5/params$sigma2)*sum((const$y[const$k_index==kk]-(get_Btheta(const$X%*%params$theta[(b-1)*const$L+(1:const$L)],const)%*%params$betastar[(a-1)*const$d+(1:const$d)]+params$u))^2))
#                 exp(log(params$pimat[a,b]) -(0.5/params$sigma2)*sum((const$y[const$k_index==kk]-(get_Btheta(const$X%*%params$theta[(b-1)*const$L+(1:const$L)],const)%*%params$betastar[(a-1)*const$d+(1:const$d)]+params$u))^2))
#               })
#             })
#     probs <- probs/sum(probs)  ## standardize them
#
#     ## sample 1 of C^2 with correct probabilities
#     ab <- sample(1:(const$C^2),1,prob=c(probs))
#     ## check which a and b this corresponds to
#     Zk <- which(matrix(1:(const$C^2),ncol=const$C)==ab,arr.ind = TRUE) ## which element of CxC matrix
#     params$Zbeta[kk] <- Zk[1]  ## a=corresponding row
#     params$Ztheta[kk] <- Zk[2] ## b=corresponding column
#   }
#   return(params)
# }

## not currently being used. numerical issues
update_clustMemb <- function(params,const){

  for(kk in 1:const$K){

    ## first Zbetak ##

    ## compute probabilities for all possible a
    probs <- sapply(1:const$C, function(a){  ## loop over a (rows; beta clusters)
                exp(log(params$pimat[a,params$Ztheta[kk]]) -(0.5/params$sigma2)*sum((const$y[const$k_index==kk]-(get_Btheta(const$X%*%params$thetastar[(params$Ztheta[kk]-1)*const$L+(1:const$L)],const)%*%params$betastar[(a-1)*const$d+(1:const$d)]+params$u))^2))
             })
    probs <- probs/sum(probs)  ## standardize them

    ## sample 1 of C with correct probabilities
    tryCatch({params$Zbeta[kk] <- sample(1:const$C,1,prob=c(probs))},
             error=function(err){print(paste0("Skipping cluster member update for k=",kk))})



    ## now Zthetak ##

    ## compute probabilities for all possible b
    probs <- sapply(1:const$C, function(b){  ## loop over a (rows; beta clusters)
      exp(log(params$pimat[params$Zbeta[kk],b]) -(0.5/params$sigma2)*sum((const$y[const$k_index==kk]-(get_Btheta(const$X%*%params$thetastar[(b-1)*const$L+(1:const$L)],const)%*%params$betastar[(params$Zbeta[kk]-1)*const$d+(1:const$d)]+params$u))^2))
    })
    probs <- probs/sum(probs)  ## standardize them

    ## sample 1 of C with correct probabilities
    tryCatch({params$Ztheta[kk] <- sample(1:const$C,1,prob=c(probs))},
             error=function(err){print(paste0("Skipping cluster member update for k=",kk))})


  }
  return(params)
}




## grid search for V
update_V <- function(params,const){

  for(cc in 1:(const$C-1)){

    ## sample V_c^beta
    probs <- exp(sapply(const$grid,function(vv){get_Vlogdensity(vv,cc,params,Vbeta=TRUE,const) })) ## looping over grid, compute density
    params$Vbeta[cc] <- sample(const$grid,1,prob=probs) ## sample from the grid

    ## sample V_c^theta
    probs <- exp(sapply(const$grid,function(vv){get_Vlogdensity(vv,cc,params,Vbeta=FALSE,const) })) ## looping over grid, compute density
    params$Vtheta[cc] <- sample(const$grid,1,prob=probs) ## sample from the grid

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
    prop_params <- params
    ## old version often returned 1 when shape2 is too small.
    ## setting minimum of shape2 param to be 1 for proposals
    s1 <- 1+nbeta[cc]
    s2 <- max(1,prop_params$alpha[1]+sum(nbeta[(cc+1):const$C]))
    # prop_params$Vbeta[cc] <- rbeta(1,shape1=1+nbeta[cc],shape2=prop_params$alpha[1]+sum(nbeta[(cc+1):const$C]) )
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


    ## update Vtheta ##

    prop_params <- params
    ## old version often returned 1 when shape2 is too small.
    ## setting minimum of shape2 param to be 1 for proposals
    s1 <- 1+ntheta[cc]
    s2 <- max(1,prop_params$alpha[2]+sum(ntheta[(cc+1):const$C]))
    # prop_params$Vtheta[cc] <- rbeta(1,shape1=1+ntheta[cc],shape2=prop_params$alpha[2]+sum(ntheta[(cc+1):const$C]) )
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

  return(params)
}


update_betastar <- function(params,const){

  for(cc in 1:const$C){

    if(sum(params$Zbeta==cc)>0){
      whichk <- which(params$Zbeta==cc) ## which k are in cluster cc
      if(const$sharedlambda==TRUE){
        lambda_beta <- params$lambda_beta
      }else{
        lambda_beta <- params$lambda_beta[cc]
      }


      ## compute Vmat only once ## summing over all k in cluster cc
      Vmat <- solve(lambda_beta*const$invSig0+(1/params$sigma2)*t(params$Btheta[const$k_index%in%whichk,])%*%params$Btheta[const$k_index%in%whichk,])

      params$betastar[(cc-1)*const$d+(1:const$d)] <- mvtnorm::rmvnorm(n=1,
                                                             mean=Vmat%*%t(lambda_beta*t(const$mu0)%*%const$invSig0+(1/params$sigma2)*(t((const$y-rep(params$u,each=const$K))[const$k_index%in%whichk])%*%params$Btheta[const$k_index%in%whichk,])  ),
                                                             sigma=Vmat)
    }else{ ## if n_c=0, draw from the prior
      if(const$sharedlambda==TRUE){
        lambda_beta <- params$lambda_beta
      }else{
        lambda_beta <- params$lambda_beta[cc]
      }
      params$betastar[(cc-1)*const$d+(1:const$d)] <- mvtnorm::rmvnorm(n=1,
                                                             mean=const$mu0 ,
                                                             sigma=(1/lambda_beta)*solve(const$invSig0) )

    }

  }
  params <- assign_betas(params,const)
  return(params)
}


update_thetastar <- function(params,const){

  for(cc in 1:const$C){
    if(sum(params$Ztheta==cc)>0){ ## n_c>0 (otherwise draw from prior)

      tryCatch({params$thetastar[(cc-1)*const$L+(1:const$L)] <-
        rFisherBingham(1,
                       mu = const$prior_tau_theta*as.matrix(rep(1,const$L)) + (1/params$sigma2)*t(const$X)%*%c(Reduce('+',## summing over the k vectors
                                                                                                                      lapply(which(params$Ztheta==cc), ## only k in cluster cc
                                                                                                                             function(kk){ ## Wk*ytildek
                                                                                                                               ((params$DerivBtheta[const$k_index==kk,]%*%params$betastar[(cc-1)*const$d+(1:const$d)])^2)* ## Wk
                                                                                                                                 (const$X%*%params$thetastar[(cc-1)*const$L+(1:const$L)]+((const$y[const$k_index==kk]-params$u)-params$Btheta[const$k_index==kk,]%*%params$betastar[(cc-1)*const$d+(1:const$d)])/(params$DerivBtheta[const$k_index==kk,]%*%params$betastar[(cc-1)*const$d+(1:const$d)]))
                                                                                                                             }))),
                       Aplus = -(0.5/params$sigma2)*t(const$X)%*%diag(c(## vectorizing then turning into single diagonal matrix at the end
                         Reduce('+',## summing over the k
                                lapply(which(params$Ztheta==cc), ## only k in cluster cc
                                       function(kk){ ## Wk
                                         ((params$DerivBtheta[const$k_index==kk,]%*%params$betastar[(cc-1)*const$d+(1:const$d)])^2)
                                       }))))%*%const$X
                       -0.5*exp(params$loglambda_theta)*const$PEN,
                       mtop=const$rfbtries)},
               error=function(err){print(paste0("Skipping thetastar for cluster c=",cc))} )


      params$Btheta[const$k_index%in%const$Ztheta[which(const$Ztheta==cc)],] <- get_Btheta(const$X%*%params$thetastar[(cc-1)*const$L+(1:const$L)],const)
      params$DerivBtheta[const$k_index%in%const$Ztheta[which(const$Ztheta==cc)],] <- get_DerivBtheta(const$X%*%params$thetastar[(cc-1)*const$L+(1:const$L)],const)

    }else{## if n_c=0, draw from the prior
      params$thetastar[(cc-1)*const$L+(1:const$L)] <- rFisherBingham(1,mu = const$prior_tau_theta*rep(1,const$L), Aplus = 0)
    }

  }
  params <- assign_thetas(params,const)
  return(params)
}


update_alpha <- function(params,const){

  params$alpha <- c(
    rgamma(1,shape=const$prior_alpha_beta[1]+const$C-1,rate=const$prior_alpha_beta[2]-sum(log(1-params$Vbeta[-const$C]))),
    rgamma(1,shape=const$prior_alpha_theta[1]+const$C-1,rate=const$prior_alpha_theta[2]-sum(log(1-params$Vtheta[-const$C])))
  )
  return(params) ## alpha_beta then alpha_theta
}


update_logrho <- function(params,const){

  ## proposal
  prop_params <- params
  prop_params$logrho <- params$logrho+rnorm(1,0,const$stepsize_logrho)
  prop_params <- compute_pi(prop_params,const) ## update proposed pistar and pi accordingly

  ## compute log-acceptance ratio
  logLikRatio <- sum(sapply(1:const$K,function(kk){log(prop_params$pistarmat[params$Zbeta[kk],params$Ztheta[kk]])}))-
    sum(sapply(1:const$K,function(kk){log(params$pistarmat[params$Zbeta[kk],params$Ztheta[kk]])}))

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
    gamrate <- const$prior_lambda_beta[2]+0.5*sum(sapply(1:const$C, ## summing over all clusters cc
                                                         function(cc){c(t(params$betastar[(cc-1)*const$d+(1:const$d)])%*%const$invSig0%*%params$betastar[(cc-1)*const$d+(1:const$d)])}))
    if(is.finite(gamrate)){
      params$lambda_beta <- rgamma(1,
                                   shape=const$prior_lambda_beta[1]+0.5*const$C*const$d,
                                   rate=gamrate)
    }
  }else{
    for(cc in 1:const$C){

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
  logLikRatio <- (-0.5*exp(prop_params$loglambda_theta)*sum(sapply(1:const$C, ## summing over all clusters cc
                                                                  function(cc){c(t(prop_params$thetastar[(cc-1)*const$L+(1:const$L)])%*%const$PEN%*%prop_params$thetastar[(cc-1)*const$L+(1:const$L)])}))
                        - const$C*(fb.saddle(gam=const$prior_tau_theta*rep(1,const$L),lam=-0.5*exp(prop_params$loglambda_theta)*eigen(const$PEN)$values)[3]) ) -  ## third element is third order approximation
    (-0.5*exp(params$loglambda_theta)*sum(sapply(1:const$C, ## summing over all clusters cc
                                                      function(cc){c(t(params$thetastar[(cc-1)*const$L+(1:const$L)])%*%const$PEN%*%params$thetastar[(cc-1)*const$L+(1:const$L)])}))
            - const$C*(fb.saddle(gam=const$prior_tau_theta*rep(1,const$L),lam=-0.5*exp(params$loglambda_theta)*eigen(const$PEN)$values)[3]) ) ## third element is third order approximation

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
  B_beta <- get_B_beta(params,const)

  for(ii in 1:const$n){## loop over all subjects
     params$u[ii] <- rnorm(1,
                          mean=(params$sigma2_u/(params$sigma2+const$K*params$sigma2_u))*sum((const$y[const$i_index==ii]-B_beta[const$i_index==ii])),
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
  B_beta <- get_B_beta(params,const)

  params$sigma2 <- 1/rgamma(1,
                            shape=const$prior_sigma2[1]+0.5*const$n*const$K,
                            rate=const$prior_sigma2[2]+0.5*(sum((const$y-(B_beta+rep(params$u,each=const$K)))^2)) )

  return(params)
}




