### Old versions of helper/update functions for MIM==TRUE
##### Not used in current version
### All used a previous parameterization, defining an identifiability product term to avoid double counting
### Current version does not use this identifiability product


# ### helper_functions
#
# ## get basis functions
# get_Btheta <- function(Xomega,const,params=NULL,k,j){
#
#   if(const$MIM==FALSE | j==1 ){
#     IDprod <- 1
#   }else{  ## identifiability product for MIM
#     IDprod <- prod(params$Ztheta[k,1:(j-1)]!=params$Ztheta[k,j])
#   }
#
#   if(const$LM==TRUE){ ## if forcing linearity
#     return(IDprod*Xomega)
#   }else{
#     return(IDprod*mgcv::PredictMat(const$SS,data=data.frame(Xomega)))
#   }
#
# }
#
#
# ## get derivatives of Basis functions
# get_DerivBtheta <- function(Xomega,const,params,k,j){
#   # if(const$MIM==TRUE){
#   if(const$MIM==FALSE | j==1 ){
#     IDprod <- 1
#   }else{  ## identifiability product for MIM
#     IDprod <- prod(params$Ztheta[k,1:(j-1)]!=params$Ztheta[k,j])
#   }
#   if(const$LM==TRUE){ ## deriv is 1 if forcing linear effects
#     return(IDprod*((Xomega)^0))
#   }else{
#     return(IDprod*mgcv::PredictMat(const$SSderiv,data=data.frame(Xomega)))
#   }
#
# }
#
#
#
#
# ## alternate version for updating Ztheta=b in update_clustMemb
# get_Btheta_b <- function(Xomega,const,Ztheta=NULL,k,j){
#   if(const$MIM==FALSE | j==1 ){
#     IDprod <- 1
#   }else{  ## identifiability product for MIM
#     IDprod <- prod(Ztheta[k,1:(j-1)]!=Ztheta[k,j])
#   }
#   if(const$LM==TRUE){ ## if forcing linear effects
#     return(IDprod*Xomega)
#   }else{
#     return(IDprod*mgcv::PredictMat(const$SS,data=data.frame(Xomega)))
#   }
#
# }
#
#
# ### updates
#
# ##
# update_clustMemb <- function(params,const){
#
#   params$err <- 0 ## for error tracking
#
#   for(kk in 1:const$K){## loop over outcomes
#     for(jj in 1:const$p){## loop over exposures
#
#       ## Zbetakj ##
#       if(const$clustering=="both" | const$clustering=="beta"){
#
#         ## necessary components
#         B_beta <- get_B_beta_k(params,const,kk)
#         y_B_u <- const$y[const$k_index==kk]-params$b0[kk]-(apply(B_beta[,-jj,drop=F],1,sum) +params$xi*sqrt(params$sigma2[kk])*params$u+const$Zcovariates%*%params$betaZk[kk,])
#         Bth_kj <- get_Btheta(const$X[[jj]]%*%params$omegastar[(params$Ztheta[kk,jj]-1)*const$L+(1:const$L)],const,params,kk,jj)
#
#         ## compute probabilities for all possible a
#         logprobs <- c(sapply(1:const$Cbeta, function(a){  ## loop over a (rows; beta clusters)
#           (log(params$pimat[a,params$Ztheta[kk,jj]]) -(0.5/params$sigma2[kk])*sum((y_B_u-Bth_kj%*%params$betastar[(a-1)*const$d+(1:const$d)])^2))
#         })) ## to be standardized below
#
#         ## error handling for very small values
#         if(!is.finite(sum(exp(logprobs)/sum(exp(logprobs))))){
#           stbfctr <- -500-max(logprobs) ## max largest value -500
#           logprobs <- logprobs+stbfctr ## multiply all probs by common factor
#         }
#
#         ## sample 1 of C with correct probabilities
#         newZbeta <- tryCatch(sample(1:const$Cbeta,1,prob=exp(logprobs)/sum(exp(logprobs))), ## standardized probs
#                              error=function(err){NULL})
#
#         if(!is.null(newZbeta)){
#           params$Zbeta[kk,jj] <- newZbeta
#         }else{## track errors
#           params$err <- 1
#         }
#
#       }
#
#
#       ## Zthetakj ##
#       if(const$clustering=="both" | const$clustering=="theta"){
#
#         ## compute probabilities for all possible b
#         if(const$MIM==TRUE){ ## changing Ztheta may change IDmat
#           y_u <- const$y[const$k_index==kk]-params$b0[kk]-(params$xi*sqrt(params$sigma2[kk])*params$u+const$Zcovariates%*%params$betaZk[kk,])
#           tempZtheta <- params$Ztheta
#           logprobs <- c(sapply(1:const$Ctheta, function(b){  ## loop over b (columns; theta clusters)
#             tempZtheta[kk,jj] <- b
#             testconst <- const; testconst$MIM=FALSE
#             (log(params$pimat[params$Zbeta[kk,jj],b]) -(0.5/params$sigma2[kk])*sum((y_u-apply(get_B_beta_k_b(params,const,kk,tempZtheta),1,sum))^2))
#
#           }))## to be standardized below
#
#         }else{ ## dont need to worry about ID mat
#           ## now have changed
#           B_beta <- get_B_beta_k(params,const,kk)
#           y_B_u <- const$y[const$k_index==kk]-params$b0[kk]-(apply(B_beta[,-jj,drop=F],1,sum) +params$xi*sqrt(params$sigma2[kk])*params$u+const$Zcovariates%*%params$betaZk[kk,])
#
#           logprobs <- c(sapply(1:const$Ctheta, function(b){  ## loop over b (columns; theta clusters)
#             (log(params$pimat[params$Zbeta[kk,jj],b]) -(0.5/params$sigma2[kk])*sum((y_B_u- get_Btheta(const$X[[jj]]%*%params$omegastar[(b-1)*const$L+(1:const$L)],const,params,kk,jj)%*%params$betastar[(params$Zbeta[kk,jj]-1)*const$d+(1:const$d)])^2))
#           }))## to be standardized below
#         }
#
#         ## error handling for very small values
#         if(!is.finite(sum(exp(logprobs)/sum(exp(logprobs))))){
#           stbfctr <- -500-max(logprobs) ## max largest value -500
#           logprobs <- logprobs+stbfctr ## multiply all probs by common factor
#         }
#
#         ## sample 1 of C with correct probabilities
#         newZtheta <- tryCatch(sample(1:const$Ctheta,1,prob=exp(logprobs)/sum(exp(logprobs))),
#                               error=function(err){NULL})
#
#         if(!is.null(newZtheta)){
#           params$Ztheta[kk,jj] <- newZtheta
#         }else{## track errors
#           params$err <- 1
#         }
#
#       }
#
#     }
#
#   }
#   return(params)
# }
#
# ##
# update_betastar <- function(params,const){
#
#   ## computing only once
#   Btheta <- lapply(1:const$p,function(jj){
#     Reduce("rbind",lapply(1:const$K,function(kk){
#       get_Btheta(const$X[[jj]]%*%params$omegastar[(params$Ztheta[kk,jj]-1)*const$L+(1:const$L)],const,params,kk,jj)
#     }))
#   })
#
#   if(const$MIM==TRUE){ ## for identifiability product in MIM
#     IDprodmat <- matrix(1,nrow=const$K,ncol=const$p)
#     for(kk in 1:const$K){
#       if(const$p>1){
#         for(jj in 2:const$p){
#           IDprodmat[kk,jj] <- prod(params$Ztheta[kk,1:(jj-1)]!=params$Ztheta[kk,jj])
#         }
#       }
#     }
#   }
#
#   for(cc in 1:const$Cbeta){
#     if(const$MIM==TRUE){
#       n_c <- sum((params$Zbeta==cc)*IDprodmat)
#     }else{
#       n_c <- sum(params$Zbeta==cc)
#     }
#
#
#     if(n_c>0){
#       if(const$sharedlambda==TRUE){
#         lambda_beta <- params$lambda_beta
#       }else{
#         lambda_beta <- params$lambda_beta[cc]
#       }
#
#       whichZ <- which(params$Zbeta==cc,arr.ind=TRUE)
#       whichk <- sort(unique(whichZ[,1]))
#       whichkj <- lapply(1:const$K,function(kk){sort(whichZ[whichZ[,1]==kk,2])})
#       whichkNotj <- lapply(1:const$K,function(kk){which(!(1:const$p)%in%whichkj[[kk]])  })
#
#       ## Btheta for relevant k,j pairs
#       B_kc <- lapply(whichk,function(kk){
#         Reduce("+",lapply(whichkj[[kk]],function(jj){
#           Btheta[[jj]][const$k_index==kk,]/sqrt(params$sigma2[kk]) ## sqrt since it gets squared in BTB
#         }))
#       })
#
#       ## sum of B^TB across relevant k
#       BTB <- Reduce("+",lapply(B_kc,function(BB){t(BB)%*%BB}))
#
#       ##
#       y_u_B_k <- lapply(whichk,function(kk){
#         y_u <- const$y[const$k_index==kk]-params$b0[kk]-params$xi*sqrt(params$sigma2[kk])*params$u-const$Zcovariates%*%params$betaZk[kk,]
#         if(length(whichkNotj[[kk]])>0){
#           y_u <- y_u - Reduce("+",lapply(whichkNotj[[kk]],function(jj){
#             Btheta[[jj]][const$k_index==kk,,drop=F]%*%params$betastar[(params$Zbeta[kk,jj]-1)*const$d+(1:const$d),drop=F]
#           }))
#         }
#         return(y_u/sqrt(params$sigma2[kk])) ## sqrt because it multiples B
#       })
#
#       ##
#       yTB <- Reduce("+",lapply(1:length(whichk),function(kk){
#         t(y_u_B_k[[kk]])%*%B_kc[[kk]]
#       }))
#
#       ## compute Vmat only once ## summing over all k in cluster cc
#       Vmat <- solve(lambda_beta*const$invSig0+BTB)
#       Vmat <- (Vmat+t(Vmat))/2
#
#       params$betastar[(cc-1)*const$d+(1:const$d)] <- mvtnorm::rmvnorm(n=1,
#                                                                       mean=Vmat%*%t(lambda_beta*t(const$mu0)%*%const$invSig0+yTB  ),
#                                                                       sigma=Vmat)
#
#
#     }else{ ## if n_c=0, draw from the prior
#       if(const$sharedlambda==TRUE){
#         lambda_beta <- params$lambda_beta
#       }else{
#         lambda_beta <- params$lambda_beta[cc]
#       }
#       params$betastar[(cc-1)*const$d+(1:const$d)] <- mvtnorm::rmvnorm(n=1,mean=const$mu0, sigma=(1/lambda_beta)*MASS::ginv(const$invSig0) )
#
#     }
#
#   }
#
#   return(params)
# }





