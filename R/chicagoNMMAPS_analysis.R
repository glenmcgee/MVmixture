
library(dlnm)
library(ggplot2)
library(patchwork)
library(reshape)
library(mgcv)
data(chicagoNMMAPS)
maxlag <- 14
covlag <- 3 ## for temperature
nit <- 30000
nburn <- 0.5*nit
nthin = 5


## carry forward exposure values if missing
for(ii in 2:nrow(chicagoNMMAPS)){
  if(is.na(chicagoNMMAPS$pm10[ii])){chicagoNMMAPS$pm10[ii] <- chicagoNMMAPS$pm10[ii-1]}
  if(is.na(chicagoNMMAPS$rhum[ii])){chicagoNMMAPS$rhum[ii] <- chicagoNMMAPS$rhum[ii-1]}
}

## standardize
mn_pm10 <- mean(chicagoNMMAPS$pm10)
sd_pm10 <- sd(chicagoNMMAPS$pm10)
mn_o3 <- mean(chicagoNMMAPS$o3)
sd_o3 <- sd(chicagoNMMAPS$o3)
chicagoNMMAPS$pm10 <- (chicagoNMMAPS$pm10-mn_pm10)/sd_pm10
chicagoNMMAPS$o3 <- (chicagoNMMAPS$o3-mn_o3)/sd_o3

## make lagged exposures
pm10 <- o3 <-  matrix(NA,ncol=maxlag,nrow=nrow(chicagoNMMAPS))
for(ii in (maxlag):nrow(chicagoNMMAPS)){
  pm10[ii,] <- c(chicagoNMMAPS$pm10[(ii-maxlag+1):ii])
  o3[ii,] <- c(chicagoNMMAPS$o3[(ii-maxlag+1):ii])
}

## remove a pm10 outlier 15 SDs away from the mean
outlier_IDs <- which(apply(pm10,1,max)>15)
pm10 <- pm10[-outlier_IDs,]
o3 <- o3[-outlier_IDs,]
chicagoNMMAPS <- chicagoNMMAPS[-outlier_IDs,]

Xlist <- list(pm10[maxlag:nrow(chicagoNMMAPS),,drop=F],
              o3[maxlag:nrow(chicagoNMMAPS),,drop=F])

chicagoNMMAPS$avgtemp <- chicagoNMMAPS$avgdptp <- rep(0,nrow(chicagoNMMAPS))
for(ii in (covlag):nrow(chicagoNMMAPS)){
  chicagoNMMAPS$avgtemp[ii] <- mean(chicagoNMMAPS$temp[(ii-covlag+1):ii])
  chicagoNMMAPS$avgdptp[ii] <- mean(chicagoNMMAPS$dptp[(ii-covlag+1):ii])
}

## get spline terms for time trends (and covariates)
prelim_gam <- gam(log(death) ~ s(year,k=4)+
                    s(month, bs = "cc",k=4)+
                    # s(doy,k=5)+
                    as.factor(dow)+
                    s(temp,k=6)+
                    s(avgtemp,k=3)+
                    s(dptp,k=6)+
                    s(avgdptp,k=3),
                    method = "REML",data=chicagoNMMAPS)

## prep covariates
Ztime <- model.matrix(prelim_gam)[,-1]
Ztime <- Ztime[maxlag:nrow(chicagoNMMAPS),]

## prep outcomes
Ycvd <- log(chicagoNMMAPS$cvd[maxlag:nrow(chicagoNMMAPS)]+0.5)
Yresp <- log(chicagoNMMAPS$resp[maxlag:nrow(chicagoNMMAPS)]+0.5)
Yother <- log((chicagoNMMAPS$death-chicagoNMMAPS$cvd-chicagoNMMAPS$resp)[maxlag:nrow(chicagoNMMAPS)]+0.5)
Y <- as.matrix(cbind(Ycvd,Yresp,Yother))

## prep new X for predictions
npred <- 50
seq1 <- seq(quantile(Xlist[[1]],0.1),quantile(Xlist[[1]],0.9),length=npred)
seq2 <- seq(quantile(Xlist[[2]],0.1),quantile(Xlist[[2]],0.9),length=npred)
Xnew <- list(matrix(seq1,ncol=maxlag,nrow=npred,byrow=F),
             matrix(seq2,ncol=maxlag,nrow=npred,byrow=F))


#################################
## Plotting functions
#################################
spaghettiplot <- function(x,constraint=TRUE){
  if(constraint==TRUE){ ## post hoc identifiability constraint for signflipping
    # id <- apply(x,1,sum)<0
    id <- x[,ncol(x)]<0
    x <- x*((-1)^id)
  }


  df <- data.frame(Weight=c(t(x)),
                   Lag=rep((maxlag-1):0,nrow(x)),
                   ID=rep(1:nrow(x),each=maxlag))

  pp <- ggplot(df,aes(y=Weight, x=Lag, group=ID))+
    geom_line(alpha=0.02,col="#4E84C4")+
    ylim(-1,1)+
    geom_hline(yintercept=0,linetype="dashed")+
    # ylab(bquote(theta))+
    scale_x_reverse()+
    theme_bw() +
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

  return(pp)
}

erfplot <- function(pred,
                    jj=1, ## which exposure
                    kk=1, ## which outcome
                    ylims=NULL){

    grid <- list(seq1,seq2)[[jj]]


  df <- data.frame(
    gridx=grid,
    fitted=pred$summary[[kk]][[jj]]$mean,#predict(loess(cbmim_pred_ind[[jj]]$mean~seq(-2,2,length=21))),
    uci=pred$summary[[kk]][[jj]]$upper,#predict(loess(cbmim_pred_ind[[jj]]$uci~seq(-2,2,length=21))),
    lci=pred$summary[[kk]][[jj]]$lower)#predict(loess(cbmim_pred_ind[[jj]]$lci~seq(-2,2,length=21))) )

  if(is.null(ylims)){
    ylims <-c( min(df$lci)-0.002,max(df$uci)+0.002)
  }
  pp <- ggplot(df,aes(x=gridx,y=fitted))+
    geom_hline(yintercept=0,colour="red",linetype=2)+
    geom_ribbon(aes_string(ymin="lci",ymax="uci"),fill="lightgray")+
    geom_line(linetype=3)+ ## 3 is dotted
    ylim(ylims[1],ylims[2])+
    ylab(bquote(f[.(jj)](X[.(jj)])))+
    xlab(bquote(X[.(jj)]))+
    scale_x_continuous(expand = c(0, 0)) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  return(pp)

}



lagplot <- function(pred,
                    ylims=NULL){

  df <- data.frame(
    Lag=c((maxlag-1):0),
    fitted=pred$mean,#predict(loess(cbmim_pred_ind[[jj]]$mean~seq(-2,2,length=21))),
    uci=pred$upper,#predict(loess(cbmim_pred_ind[[jj]]$uci~seq(-2,2,length=21))),
    lci=pred$lower)#predict(loess(cbmim_pred_ind[[jj]]$lci~seq(-2,2,length=21))) )

  if(is.null(ylims)){
    ylims <-c( min(df$lci)-0.002,max(df$uci)+0.002)
  }

  pp <- ggplot(df,aes(y=fitted, x=Lag))+
    geom_point()+
    geom_errorbar(aes(ymin=lci, ymax=uci), width=.2)+
    ylim(ylims[1],ylims[2])+
    geom_hline(yintercept=0,linetype="dashed")+
    scale_x_reverse()+
    ylab("Contrast")+
    theme_bw() +
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())

  return(pp)

}


# #######################################
# ## Fit single XY models
# #######################################
# set.seed(1)
# nit <- 5000
# nburn <- 0.5*nit
# nthin = 5
#
# chicagoNMMAPS_DLAG11 <- MVmix(Y[,1,drop=F],list(Xlist[[1]]),Z=Ztime,
#                             niter=nit,nburn=nburn,nthin=nthin,
#                             Vgridsearch = TRUE,gridsize=10,
#                             DLM=TRUE,DLMpenalty=TRUE,lagOrder=4,diff=2,
#                             cluster="neither",maxClusters=1,approx=TRUE)
# save(chicagoNMMAPS_DLAG11, file = paste0("Results/chicagoNMMAPS_DLAG11",maxlag,".RData"))
#
# pred_DLAG11 <- predict_MVmix(chicagoNMMAPS_DLAG11,
#                            newX = Xnew[[1]],# newZ = Ztime,
#                            include_intercept=FALSE,fixomega = FALSE,
#                            allx=FALSE,contrast=TRUE)
#
# erfplot(pred_DLAG11,1)
#
# spaghettiplot(chicagoNMMAPS_DLAG11$omega[[1]][[1]])
#
#
# chicagoNMMAPS_DLAG12 <- MVmix(Y[,1,drop=F],list(Xlist[[2]]),Z=Ztime,
#                               niter=nit,nburn=nburn,nthin=nthin,
#                               Vgridsearch = TRUE,gridsize=10,
#                               DLM=TRUE,DLMpenalty=TRUE,lagOrder=4,diff=2,
#                               cluster="neither",maxClusters=1,approx=TRUE)
# save(chicagoNMMAPS_DLAG12, file = paste0("Results/chicagoNMMAPS_DLAG12",maxlag,".RData"))
#
# pred_DLAG12 <- predict_MVmix(chicagoNMMAPS_DLAG12,
#                              newX = Xnew[[2]],# newZ = Ztime,
#                              include_intercept=FALSE,fixomega = FALSE,
#                              allx=FALSE,contrast=TRUE)
#
# erfplot(pred_DLAG12,1)
#
# spaghettiplot(chicagoNMMAPS_DLAG12$omega[[1]][[1]])
#
#
#
#
# chicagoNMMAPS_DLAG21 <- MVmix(Y[,2,drop=F],list(Xlist[[1]]),Z=Ztime,
#                               niter=nit,nburn=nburn,nthin=nthin,
#                               Vgridsearch = TRUE,gridsize=10,
#                               DLM=TRUE,DLMpenalty=TRUE,lagOrder=4,diff=2,
#                               cluster="neither",maxClusters=1,approx=TRUE)
# save(chicagoNMMAPS_DLAG21, file = paste0("Results/chicagoNMMAPS_DLAG21",maxlag,".RData"))
#
# pred_DLAG21 <- predict_MVmix(chicagoNMMAPS_DLAG21,
#                              newX = Xnew[[1]],# newZ = Ztime,
#                              include_intercept=FALSE,fixomega = FALSE,
#                              allx=FALSE,contrast=TRUE)
#
# erfplot(pred_DLAG21,1)
#
# spaghettiplot(chicagoNMMAPS_DLAG21$omega[[1]][[1]])
#
#
# chicagoNMMAPS_DLAG22 <- MVmix(Y[,2,drop=F],list(Xlist[[2]]),Z=Ztime,
#                               niter=nit,nburn=nburn,nthin=nthin,
#                               Vgridsearch = TRUE,gridsize=10,
#                               DLM=TRUE,DLMpenalty=TRUE,lagOrder=4,diff=2,
#                               cluster="neither",maxClusters=1,approx=TRUE)
# save(chicagoNMMAPS_DLAG22, file = paste0("Results/chicagoNMMAPS_DLAG22",maxlag,".RData"))
#
# pred_DLAG22 <- predict_MVmix(chicagoNMMAPS_DLAG22,
#                              newX = Xnew[[2]],# newZ = Ztime,
#                              include_intercept=FALSE,fixomega = FALSE,
#                              allx=FALSE,contrast=TRUE)
#
# erfplot(pred_DLAG22,1)
#
# spaghettiplot(chicagoNMMAPS_DLAG22$omega[[1]][[1]])
#
#
#
#
# chicagoNMMAPS_DLAG31 <- MVmix(Y[,3,drop=F],list(Xlist[[1]]),Z=Ztime,
#                               niter=nit,nburn=nburn,nthin=nthin,
#                               Vgridsearch = TRUE,gridsize=10,
#                               DLM=TRUE,DLMpenalty=TRUE,lagOrder=4,diff=2,
#                               cluster="neither",maxClusters=1,approx=TRUE)
# save(chicagoNMMAPS_DLAG31, file = paste0("Results/chicagoNMMAPS_DLAG31",maxlag,".RData"))
#
# pred_DLAG31 <- predict_MVmix(chicagoNMMAPS_DLAG31,
#                              newX = Xnew[[1]],# newZ = Ztime,
#                              include_intercept=FALSE,fixomega = FALSE,
#                              allx=FALSE,contrast=TRUE)
#
# erfplot(pred_DLAG31,1)
#
# spaghettiplot(chicagoNMMAPS_DLAG31$omega[[1]][[1]])
#
#
# chicagoNMMAPS_DLAG32 <- MVmix(Y[,3,drop=F],list(Xlist[[2]]),Z=Ztime,
#                               niter=nit,nburn=nburn,nthin=nthin,
#                               Vgridsearch = TRUE,gridsize=10,
#                               DLM=TRUE,DLMpenalty=TRUE,lagOrder=4,diff=2,
#                               cluster="neither",maxClusters=1,approx=TRUE)
# save(chicagoNMMAPS_DLAG32, file = paste0("Results/chicagoNMMAPS_DLAG32",maxlag,".RData"))
#
# pred_DLAG32 <- predict_MVmix(chicagoNMMAPS_DLAG32,
#                              newX = Xnew[[2]],# newZ = Ztime,
#                              include_intercept=FALSE,fixomega = FALSE,
#                              allx=FALSE,contrast=TRUE)
#
# erfplot(pred_DLAG32,1)
#
# spaghettiplot(chicagoNMMAPS_DLAG32$omega[[1]][[1]])
#
#
#
# ## combine results of separate fits
#
# erfplot(pred_DLAG11,1)/
# erfplot(pred_DLAG21,1)/
# erfplot(pred_DLAG31,1)
#
# erfplot(pred_DLAG12,1)/
# erfplot(pred_DLAG22,1)/
# erfplot(pred_DLAG32,1)
#
#
# spaghettiplot(chicagoNMMAPS_DLAG11$omega[[1]][[1]])/
# spaghettiplot(chicagoNMMAPS_DLAG21$omega[[1]][[1]])/
# spaghettiplot(chicagoNMMAPS_DLAG31$omega[[1]][[1]])
#
# spaghettiplot(chicagoNMMAPS_DLAG12$omega[[1]][[1]])/
# spaghettiplot(chicagoNMMAPS_DLAG22$omega[[1]][[1]])/
# spaghettiplot(chicagoNMMAPS_DLAG32$omega[[1]][[1]])
#
# lagplot(est_lag(chicagoNMMAPS_DLAG11,Xhold=-1)$summary[[1]][[1]])/
# lagplot(est_lag(chicagoNMMAPS_DLAG21,Xhold=-1)$summary[[1]][[1]])/
# lagplot(est_lag(chicagoNMMAPS_DLAG31,Xhold=-1)$summary[[1]][[1]])
#
# lagplot(est_lag(chicagoNMMAPS_DLAG12,Xhold=-1)$summary[[1]][[1]])/
# lagplot(est_lag(chicagoNMMAPS_DLAG22,Xhold=-1)$summary[[1]][[1]])/
# lagplot(est_lag(chicagoNMMAPS_DLAG32,Xhold=-1)$summary[[1]][[1]])
#
#
#
#
#
#
# #######################################
# ## Fit LINEAR single XY models
# #######################################
# set.seed(1)
# nit <- 5000
# nburn <- 0.5*nit
# nthin = 5
#
# chicagoNMMAPS_DLAGLIN11 <- MVmix(Y[,1,drop=F],list(Xlist[[1]]),Z=Ztime,
#                               niter=nit,nburn=nburn,nthin=nthin,
#                               Vgridsearch = TRUE,gridsize=10,LM=TRUE,
#                               DLM=TRUE,DLMpenalty=TRUE,lagOrder=4,diff=2,
#                               cluster="neither",maxClusters=1,approx=TRUE)
# save(chicagoNMMAPS_DLAGLIN11, file = paste0("Results/chicagoNMMAPS_DLAGLIN11",maxlag,".RData"))
#
# pred_DLAGLIN11 <- predict_MVmix(chicagoNMMAPS_DLAGLIN11,
#                              newX = Xnew[[1]],# newZ = Ztime,
#                              include_intercept=FALSE,fixomega = FALSE,
#                              allx=FALSE,contrast=TRUE)
#
# erfplot(pred_DLAGLIN11,1)
#
# spaghettiplot(chicagoNMMAPS_DLAGLIN11$omega[[1]][[1]])
#
#
# chicagoNMMAPS_DLAGLIN12 <- MVmix(Y[,1,drop=F],list(Xlist[[2]]),Z=Ztime,
#                               niter=nit,nburn=nburn,nthin=nthin,
#                               Vgridsearch = TRUE,gridsize=10,LM=TRUE,
#                               DLM=TRUE,DLMpenalty=TRUE,lagOrder=4,diff=2,
#                               cluster="neither",maxClusters=1,approx=TRUE)
# save(chicagoNMMAPS_DLAGLIN12, file = paste0("Results/chicagoNMMAPS_DLAGLIN12",maxlag,".RData"))
#
# pred_DLAGLIN12 <- predict_MVmix(chicagoNMMAPS_DLAGLIN12,
#                              newX = Xnew[[2]],# newZ = Ztime,
#                              include_intercept=FALSE,fixomega = FALSE,
#                              allx=FALSE,contrast=TRUE)
#
# erfplot(pred_DLAGLIN12,1)
#
# spaghettiplot(chicagoNMMAPS_DLAGLIN12$omega[[1]][[1]])
#
#
#
#
# chicagoNMMAPS_DLAGLIN21 <- MVmix(Y[,2,drop=F],list(Xlist[[1]]),Z=Ztime,
#                               niter=nit,nburn=nburn,nthin=nthin,
#                               Vgridsearch = TRUE,gridsize=10,LM=TRUE,
#                               DLM=TRUE,DLMpenalty=TRUE,lagOrder=4,diff=2,
#                               cluster="neither",maxClusters=1,approx=TRUE)
# save(chicagoNMMAPS_DLAGLIN21, file = paste0("Results/chicagoNMMAPS_DLAGLIN21",maxlag,".RData"))
#
# pred_DLAGLIN21 <- predict_MVmix(chicagoNMMAPS_DLAGLIN21,
#                              newX = Xnew[[1]],# newZ = Ztime,
#                              include_intercept=FALSE,fixomega = FALSE,
#                              allx=FALSE,contrast=TRUE)
#
# erfplot(pred_DLAGLIN21,1)
#
# spaghettiplot(chicagoNMMAPS_DLAGLIN21$omega[[1]][[1]])
#
#
# chicagoNMMAPS_DLAGLIN22 <- MVmix(Y[,2,drop=F],list(Xlist[[2]]),Z=Ztime,
#                               niter=nit,nburn=nburn,nthin=nthin,
#                               Vgridsearch = TRUE,gridsize=10,LM=TRUE,
#                               DLM=TRUE,DLMpenalty=TRUE,lagOrder=4,diff=2,
#                               cluster="neither",maxClusters=1,approx=TRUE)
# save(chicagoNMMAPS_DLAGLIN22, file = paste0("Results/chicagoNMMAPS_DLAGLIN22",maxlag,".RData"))
#
# pred_DLAGLIN22 <- predict_MVmix(chicagoNMMAPS_DLAGLIN22,
#                              newX = Xnew[[2]],# newZ = Ztime,
#                              include_intercept=FALSE,fixomega = FALSE,
#                              allx=FALSE,contrast=TRUE)
#
# erfplot(pred_DLAGLIN22,1)
#
# spaghettiplot(chicagoNMMAPS_DLAGLIN22$omega[[1]][[1]])
#
#
#
#
# chicagoNMMAPS_DLAGLIN31 <- MVmix(Y[,3,drop=F],list(Xlist[[1]]),Z=Ztime,
#                               niter=nit,nburn=nburn,nthin=nthin,
#                               Vgridsearch = TRUE,gridsize=10,LM=TRUE,
#                               DLM=TRUE,DLMpenalty=TRUE,lagOrder=4,diff=2,
#                               cluster="neither",maxClusters=1,approx=TRUE)
# save(chicagoNMMAPS_DLAGLIN31, file = paste0("Results/chicagoNMMAPS_DLAGLIN31",maxlag,".RData"))
#
# pred_DLAGLIN31 <- predict_MVmix(chicagoNMMAPS_DLAGLIN31,
#                              newX = Xnew[[1]],# newZ = Ztime,
#                              include_intercept=FALSE,fixomega = FALSE,
#                              allx=FALSE,contrast=TRUE)
#
# erfplot(pred_DLAGLIN31,1)
#
# spaghettiplot(chicagoNMMAPS_DLAGLIN31$omega[[1]][[1]])
#
#
# chicagoNMMAPS_DLAGLIN32 <- MVmix(Y[,3,drop=F],list(Xlist[[2]]),Z=Ztime,
#                               niter=nit,nburn=nburn,nthin=nthin,
#                               Vgridsearch = TRUE,gridsize=10,LM=TRUE,
#                               DLM=TRUE,DLMpenalty=TRUE,lagOrder=4,diff=2,
#                               cluster="neither",maxClusters=1,approx=TRUE)
# save(chicagoNMMAPS_DLAGLIN32, file = paste0("Results/chicagoNMMAPS_DLAGLIN32",maxlag,".RData"))
#
# pred_DLAGLIN32 <- predict_MVmix(chicagoNMMAPS_DLAGLIN32,
#                              newX = Xnew[[2]],# newZ = Ztime,
#                              include_intercept=FALSE,fixomega = FALSE,
#                              allx=FALSE,contrast=TRUE)
#
# erfplot(pred_DLAGLIN32,1)
#
# spaghettiplot(chicagoNMMAPS_DLAGLIN32$omega[[1]][[1]])
#
#
#
# ## combine results of separate fits
#
# erfplot(pred_DLAGLIN11,1)/
#   erfplot(pred_DLAGLIN21,1)/
#   erfplot(pred_DLAGLIN31,1)
#
# erfplot(pred_DLAGLIN12,1)/
#   erfplot(pred_DLAGLIN22,1)/
#   erfplot(pred_DLAGLIN32,1)
#
#
# spaghettiplot(chicagoNMMAPS_DLAGLIN11$omega[[1]][[1]])/
#   spaghettiplot(chicagoNMMAPS_DLAGLIN21$omega[[1]][[1]])/
#   spaghettiplot(chicagoNMMAPS_DLAGLIN31$omega[[1]][[1]])
#
# spaghettiplot(chicagoNMMAPS_DLAGLIN12$omega[[1]][[1]])/
#   spaghettiplot(chicagoNMMAPS_DLAGLIN22$omega[[1]][[1]])/
#   spaghettiplot(chicagoNMMAPS_DLAGLIN32$omega[[1]][[1]])
#
# lagplot(est_lag(chicagoNMMAPS_DLAGLIN11)$summary[[1]][[1]])/
#   lagplot(est_lag(chicagoNMMAPS_DLAGLIN21)$summary[[1]][[1]])/
#   lagplot(est_lag(chicagoNMMAPS_DLAGLIN31)$summary[[1]][[1]])
#
# lagplot(est_lag(chicagoNMMAPS_DLAGLIN12)$summary[[1]][[1]])/
#   lagplot(est_lag(chicagoNMMAPS_DLAGLIN22)$summary[[1]][[1]])/
#   lagplot(est_lag(chicagoNMMAPS_DLAGLIN32)$summary[[1]][[1]])
#
#
# #######################################
# ## Fit separate outcome models
# #######################################
# set.seed(1)
# nit <- 10000
# nburn <- 0.5*nit
# nthin = 5
#
# chicagoNMMAPS_DLAG1 <- MVmix(Y[,1,drop=F],Xlist,Z=Ztime,
#                               niter=nit,nburn=nburn,nthin=nthin,
#                               Vgridsearch = TRUE,gridsize=10,
#                               DLM=TRUE,DLMpenalty=TRUE,lagOrder=4,diff=2,
#                               cluster="neither",maxClusters=2,approx=TRUE)
# save(chicagoNMMAPS_DLAG1, file = paste0("Results/chicagoNMMAPS_DLAG1",maxlag,".RData"))
#
# pred_DLAG1 <- predict_MVmix(chicagoNMMAPS_DLAG1,
#                              newX = Xnew,# newZ = Ztime,
#                              include_intercept=FALSE,fixomega = FALSE,
#                              allx=FALSE,contrast=TRUE)
#
#
#
# chicagoNMMAPS_DLAG2 <- MVmix(Y[,2,drop=F],Xlist,Z=Ztime,
#                               niter=nit,nburn=nburn,nthin=nthin,
#                               Vgridsearch = TRUE,gridsize=10,
#                               DLM=TRUE,DLMpenalty=TRUE,lagOrder=4,diff=2,
#                               cluster="neither",maxClusters=2,approx=TRUE)
# save(chicagoNMMAPS_DLAG2, file = paste0("Results/chicagoNMMAPS_DLAG2",maxlag,".RData"))
#
# pred_DLAG2 <- predict_MVmix(chicagoNMMAPS_DLAG2,
#                              newX = Xnew,# newZ = Ztime,
#                              include_intercept=FALSE,fixomega = FALSE,
#                              allx=FALSE,contrast=TRUE)
#
#
# chicagoNMMAPS_DLAG3 <- MVmix(Y[,3,drop=F],Xlist,Z=Ztime,
#                               niter=nit,nburn=nburn,nthin=nthin,
#                               Vgridsearch = TRUE,gridsize=10,
#                               DLM=TRUE,DLMpenalty=TRUE,lagOrder=4,diff=2,
#                               cluster="neither",maxClusters=2,approx=TRUE)
# save(chicagoNMMAPS_DLAG3, file = paste0("Results/chicagoNMMAPS_DLAG3",maxlag,".RData"))
#
# pred_DLAG3 <- predict_MVmix(chicagoNMMAPS_DLAG3,
#                              newX = Xnew,# newZ = Ztime,
#                              include_intercept=FALSE,fixomega = FALSE,
#                              allx=FALSE,contrast=TRUE)
#
#
# ## combine results of separate fits
#
# erfplot(pred_DLAG1,1)/
#   erfplot(pred_DLAG2,1)/
#   erfplot(pred_DLAG3,1)
#
# erfplot(pred_DLAG1,2)/
#   erfplot(pred_DLAG2,2)/
#   erfplot(pred_DLAG3,2)
#
#
# spaghettiplot(chicagoNMMAPS_DLAG1$omega[[1]][[1]])/
#   spaghettiplot(chicagoNMMAPS_DLAG2$omega[[2]][[1]])/
#   spaghettiplot(chicagoNMMAPS_DLAG3$omega[[3]][[1]])
#
# spaghettiplot(chicagoNMMAPS_DLAG1$omega[[1]][[2]])/
#   spaghettiplot(chicagoNMMAPS_DLAG2$omega[[2]][[2]])/
#   spaghettiplot(chicagoNMMAPS_DLAG3$omega[[3]][[2]])
#
# lagplot(est_lag(chicagoNMMAPS_DLAG1,Xhold=-1)$summary[[1]][[1]])/
#   lagplot(est_lag(chicagoNMMAPS_DLAG2,Xhold=-1)$summary[[1]][[1]])/
#   lagplot(est_lag(chicagoNMMAPS_DLAG3,Xhold=-1)$summary[[1]][[1]])
#
# lagplot(est_lag(chicagoNMMAPS_DLAG1,Xhold=-1)$summary[[1]][[1]])/
#   lagplot(est_lag(chicagoNMMAPS_DLAG2,Xhold=-1)$summary[[1]][[1]])/
#   lagplot(est_lag(chicagoNMMAPS_DLAG3,Xhold=-1)$summary[[1]][[1]])
#



#######################################
## Fit DLNM
#######################################
set.seed(1)
nit <- 30000
nburn <- 0.5*nit
nthin = 5

chicagoNMMAPS_DLAG <- MVmix(scale(Y),Xlist,Z=Ztime,
                            niter=nit,nburn=nburn,nthin=nthin,
                            Vgridsearch = TRUE,gridsize=10,
                            DLM=TRUE,DLMpenalty=TRUE,lagOrder=NULL,diff=2,
                            cluster="both",maxClusters=6,approx=FALSE)
save(chicagoNMMAPS_DLAG, file = paste0("Results/chicagoNMMAPSpolar_DLAG",maxlag,".RData"))

pred_DLAG <- predict_MVmix(chicagoNMMAPS_DLAG,
                           newX = Xnew,
                           # newZ = Ztime,
                           include_intercept=FALSE,
                           fixomega = TRUE,
                           allx=FALSE,contrast=TRUE)

erfplot(pred_DLAG,jj=1,kk=1)
erfplot(pred_DLAG,jj=1,kk=2)
erfplot(pred_DLAG,jj=1,kk=3)

erfplot(pred_DLAG,jj=2,kk=1)
erfplot(pred_DLAG,jj=2,kk=2)
erfplot(pred_DLAG,jj=2,kk=3)


spaghettiplot(chicagoNMMAPS_DLAG$omega[[1]][[1]])
spaghettiplot(chicagoNMMAPS_DLAG$omega[[2]][[1]])

spaghettiplot(chicagoNMMAPS_DLAG$omega[[1]][[2]])
spaghettiplot(chicagoNMMAPS_DLAG$omega[[2]][[2]])

spaghettiplot(chicagoNMMAPS_DLAG$omega[[1]][[3]])
spaghettiplot(chicagoNMMAPS_DLAG$omega[[2]][[3]])


pc_DLAG <- pairwise_clusters(chicagoNMMAPS_DLAG)
make_heatplot(pc_DLAG$beta_y)+
  make_heatplot(pc_DLAG$theta_y)
ggsave(paste0("Results/Plots/chicagoNMMAPSpolarHeat_lag",maxlag,".pdf"),width=8,height=5)



lagplot(est_lag(chicagoNMMAPS_DLAG,Xhold=-1)$summary[[1]][[1]])/
  lagplot(est_lag(chicagoNMMAPS_DLAG,Xhold=-1)$summary[[2]][[1]])/
  lagplot(est_lag(chicagoNMMAPS_DLAG,Xhold=-1)$summary[[3]][[1]])

lagplot(est_lag(chicagoNMMAPS_DLAG,Xhold=-1)$summary[[1]][[2]])/
  lagplot(est_lag(chicagoNMMAPS_DLAG,Xhold=-1)$summary[[2]][[2]])/
  lagplot(est_lag(chicagoNMMAPS_DLAG,Xhold=-1)$summary[[3]][[2]])



#######################################
## Fit DLNM -- No clustering
#######################################
set.seed(1234)
nit <- 100000
nburn <- 0.5*nit
nthin = 5

chicagoNMMAPS_DLAGnoclust <- MVmix(scale(Y),Xlist,Z=Ztime,
                            niter=nit,nburn=nburn,nthin=nthin,
                            Vgridsearch = TRUE,gridsize=10,
                            DLM=TRUE,DLMpenalty=TRUE,lagOrder=6,diff=2,
                            cluster="neither",maxClusters=6,approx=FALSE) ## DOING POLAR
save(chicagoNMMAPS_DLAGnoclust, file = paste0("Results/chicagoNMMAPSpolar_DLAGnoclust",maxlag,".RData"))

pred_DLAGnoclust <- predict_MVmix(chicagoNMMAPS_DLAGnoclust,
                           newX = Xnew,
                           # newZ = Ztime,
                           include_intercept=FALSE,
                           fixomega = TRUE,
                           allx=FALSE,contrast=TRUE)


erfplot(pred_DLAGnoclust,jj=1,kk=1)
erfplot(pred_DLAGnoclust,jj=1,kk=2)
erfplot(pred_DLAGnoclust,jj=1,kk=3)

erfplot(pred_DLAGnoclust,jj=2,kk=1)
erfplot(pred_DLAGnoclust,jj=2,kk=2)
erfplot(pred_DLAGnoclust,jj=2,kk=3)


spaghettiplot(chicagoNMMAPS_DLAGnoclust$omega[[1]][[1]])
spaghettiplot(chicagoNMMAPS_DLAGnoclust$omega[[1]][[2]])
spaghettiplot(chicagoNMMAPS_DLAGnoclust$omega[[2]][[1]])
spaghettiplot(chicagoNMMAPS_DLAGnoclust$omega[[2]][[2]])


pc_DLAGnoclust <- pairwise_clusters(chicagoNMMAPS_DLAGnoclust)
make_heatplot(pc_DLAGnoclust$beta_y)+
  make_heatplot(pc_DLAGnoclust$theta_y)
# ggsave(paste0("Results/Plots/chicagoNMMAPSpolarHeat_lagnoclust",maxlag,".pdf"),width=8,height=5)


lagplot(est_lag(chicagoNMMAPS_DLAGnoclust,Xhold=-1)$summary[[1]][[1]])/
  lagplot(est_lag(chicagoNMMAPS_DLAGnoclust,Xhold=-1)$summary[[2]][[1]])/
  lagplot(est_lag(chicagoNMMAPS_DLAGnoclust,Xhold=-1)$summary[[3]][[1]])

lagplot(est_lag(chicagoNMMAPS_DLAGnoclust,Xhold=-1)$summary[[1]][[2]])/
  lagplot(est_lag(chicagoNMMAPS_DLAGnoclust,Xhold=-1)$summary[[2]][[2]])/
  lagplot(est_lag(chicagoNMMAPS_DLAGnoclust,Xhold=-1)$summary[[3]][[2]])

# #######################################
# ## Fit DLNM -- only cluster theta
# #######################################
# set.seed(1)
# nit <- 20000
# nburn <- 0.5*nit
# nthin = 5
#
# chicagoNMMAPS_DLAGthetaonly <- MVmix(scale(Y),Xlist,Z=Ztime,
#                             niter=nit,nburn=nburn,nthin=nthin,
#                             Vgridsearch = TRUE,gridsize=10,
#                             DLM=TRUE,DLMpenalty=TRUE,lagOrder=4,diff=2,
#                             cluster="theta",maxClusters=6,approx=TRUE)
# save(chicagoNMMAPS_DLAGthetaonly, file = paste0("Results/chicagoNMMAPS_DLAGthetaonly",maxlag,".RData"))
#
# pred_DLAGthetaonly <- predict_MVmix(chicagoNMMAPS_DLAGthetaonly,
#                            newX = Xnew,
#                            # newZ = Ztime,
#                            include_intercept=FALSE,
#                            fixomega = FALSE,
#                            allx=FALSE,contrast=TRUE)
#
# erfplot(pred_DLAGthetaonly,1)
# erfplot(pred_DLAGthetaonly,2)
#
#
# spaghettiplot(chicagoNMMAPS_DLAGthetaonly$omega[[1]][[1]])
# spaghettiplot(chicagoNMMAPS_DLAGthetaonly$omega[[2]][[1]])
# spaghettiplot(chicagoNMMAPS_DLAGthetaonly$omega[[3]][[1]])
#
# pc_DLAGthetaonly <- pairwise_clusters(chicagoNMMAPS_DLAGthetaonly)
# make_heatplot(pc_DLAGthetaonly$beta_y)+
#   make_heatplot(pc_DLAGthetaonly$theta_y)
# ggsave(paste0("Results/Plots/chicagoNMMAPSHeat_lag",maxlag,".pdf"),width=8,height=5)
#
#
#
# lagplot(est_lag(chicagoNMMAPS_DLAGthetaonly,Xhold=-1)$summary[[1]][[1]])/
#   lagplot(est_lag(chicagoNMMAPS_DLAGthetaonly,Xhold=-1)$summary[[2]][[1]])/
#   lagplot(est_lag(chicagoNMMAPS_DLAGthetaonly,Xhold=-1)$summary[[3]][[1]])
#
# lagplot(est_lag(chicagoNMMAPS_DLAGthetaonly,Xhold=-1)$summary[[1]][[2]])/
#   lagplot(est_lag(chicagoNMMAPS_DLAGthetaonly,Xhold=-1)$summary[[2]][[2]])/
#   lagplot(est_lag(chicagoNMMAPS_DLAGthetaonly,Xhold=-1)$summary[[3]][[2]])
#
#
#
# #######################################
# ## Fit DLNM Linear
# #######################################
set.seed(0)
nit <- 10000
nburn <- 0.5*nit
nthin = 5

chicagoNMMAPS_DLAGLIN <- MVmix(Y,Xlist,Z=Ztime,
                            niter=nit,nburn=nburn,nthin=nthin,
                            Vgridsearch = TRUE,gridsize=10,LM=TRUE,
                            DLM=TRUE,DLMpenalty=TRUE,lagOrder=4,diff=2,
                            cluster="both",maxClusters=6,approx=TRUE)
save(chicagoNMMAPS_DLAGLIN, file = paste0("Results/chicagoNMMAPS_DLAGLIN",maxlag,".RData"))
#
# pred_DLAGLIN <- predict_MVmix(chicagoNMMAPS_DLAGLIN,
#                            newX = Xnew,
#                            # newZ = Ztime,
#                            include_intercept=FALSE,
#                            fixomega = FALSE,
#                            allx=FALSE,contrast=TRUE)
#
# erfplot(pred_DLAGLIN,1)
# erfplot(pred_DLAGLIN,2)
#
#
# pc_DLAGLIN <- pairwise_clusters(chicagoNMMAPS_DLAGLIN)
# make_heatplot(pc_DLAGLIN$beta_y)+
#   make_heatplot(pc_DLAGLIN$theta_y)
# # ggsave(paste0("Results/Plots/chicagoNMMAPSHeat_lagLIN",maxlag,".pdf"),width=8,height=5)
#
#
# lagplot(est_lag(chicagoNMMAPS_DLAGLIN)$summary[[1]][[1]])/
#   lagplot(est_lag(chicagoNMMAPS_DLAGLIN)$summary[[2]][[1]])/
#   lagplot(est_lag(chicagoNMMAPS_DLAGLIN)$summary[[3]][[1]])
#
# lagplot(est_lag(chicagoNMMAPS_DLAGLIN)$summary[[1]][[2]])/
#   lagplot(est_lag(chicagoNMMAPS_DLAGLIN)$summary[[2]][[2]])/
#   lagplot(est_lag(chicagoNMMAPS_DLAGLIN)$summary[[3]][[2]])


#######################################
## Final Plots
#######################################

(lagplot(est_lag(chicagoNMMAPS_DLAGnoclust,Xhold=-1,Xshift=2)$summary[[1]][[1]],ylims=c(-0.15,0.3))+lagplot(est_lag(chicagoNMMAPS_DLAG,Xhold=-1,Xshift=2)$summary[[1]][[1]],ylims=c(-0.15,0.3)))/
  (lagplot(est_lag(chicagoNMMAPS_DLAGnoclust,Xhold=-1,Xshift=2)$summary[[2]][[1]],ylims=c(-0.15,0.3))+lagplot(est_lag(chicagoNMMAPS_DLAG,Xhold=-1,Xshift=2)$summary[[2]][[1]],ylims=c(-0.15,0.3)))/
  (lagplot(est_lag(chicagoNMMAPS_DLAGnoclust,Xhold=-1,Xshift=2)$summary[[3]][[1]],ylims=c(-0.15,0.3))+lagplot(est_lag(chicagoNMMAPS_DLAG,Xhold=-1,Xshift=2)$summary[[3]][[1]],ylims=c(-0.15,0.3)))
ggsave(paste0("Results/Plots/chicagoNMMAPSpolar_contrast1_lag",maxlag,".pdf"),width=6,height=8)

(lagplot(est_lag(chicagoNMMAPS_DLAGnoclust,Xhold=-1,Xshift=2)$summary[[1]][[2]],ylims=c(-0.1,0.15))+lagplot(est_lag(chicagoNMMAPS_DLAG,Xhold=-1,Xshift=2)$summary[[1]][[2]],ylims=c(-0.1,0.15)))/
  (lagplot(est_lag(chicagoNMMAPS_DLAGnoclust,Xhold=-1,Xshift=2)$summary[[2]][[2]],ylims=c(-0.1,0.15))+lagplot(est_lag(chicagoNMMAPS_DLAG,Xhold=-1,Xshift=2)$summary[[2]][[2]],ylims=c(-0.1,0.15)))/
  (lagplot(est_lag(chicagoNMMAPS_DLAGnoclust,Xhold=-1,Xshift=2)$summary[[3]][[2]],ylims=c(-0.1,0.15))+lagplot(est_lag(chicagoNMMAPS_DLAG,Xhold=-1,Xshift=2)$summary[[3]][[2]],ylims=c(-0.1,0.15)))
ggsave(paste0("Results/Plots/chicagoNMMAPSpolar_contrast2_lag",maxlag,".pdf"),width=6,height=8)



(spaghettiplot(chicagoNMMAPS_DLAGnoclust$omega[[1]][[1]])+spaghettiplot(chicagoNMMAPS_DLAG$omega[[1]][[1]]))/
(spaghettiplot(chicagoNMMAPS_DLAGnoclust$omega[[1]][[2]])+spaghettiplot(chicagoNMMAPS_DLAG$omega[[1]][[2]]))/
(spaghettiplot(chicagoNMMAPS_DLAGnoclust$omega[[1]][[3]])+spaghettiplot(chicagoNMMAPS_DLAG$omega[[1]][[3]]))
ggsave(paste0("Results/Plots/chicagoNMMAPSconstrained_weight1_lag",maxlag,".pdf"),width=6,height=8)

(spaghettiplot(chicagoNMMAPS_DLAGnoclust$omega[[2]][[1]])+spaghettiplot(chicagoNMMAPS_DLAG$omega[[2]][[1]]))/
(spaghettiplot(chicagoNMMAPS_DLAGnoclust$omega[[2]][[2]])+spaghettiplot(chicagoNMMAPS_DLAG$omega[[2]][[2]]))/
(spaghettiplot(chicagoNMMAPS_DLAGnoclust$omega[[2]][[3]])+spaghettiplot(chicagoNMMAPS_DLAG$omega[[2]][[3]]))
ggsave(paste0("Results/Plots/chicagoNMMAPSconstrained_weight2_lag",maxlag,".pdf"),width=6,height=8)


(erfplot(pred_DLAGnoclust,jj=1,kk=1)+erfplot(pred_DLAG,jj=1,kk=1))/
(erfplot(pred_DLAGnoclust,jj=1,kk=2)+erfplot(pred_DLAG,jj=1,kk=2))/
(erfplot(pred_DLAGnoclust,jj=1,kk=3)+erfplot(pred_DLAG,jj=1,kk=3))
ggsave(paste0("Results/Plots/chicagoNMMAPSpolar_erf1_lag",maxlag,".pdf"),width=6,height=8)

(erfplot(pred_DLAGnoclust,jj=2,kk=1)+erfplot(pred_DLAG,jj=2,kk=1))/
(erfplot(pred_DLAGnoclust,jj=2,kk=2)+erfplot(pred_DLAG,jj=2,kk=2))/
(erfplot(pred_DLAGnoclust,jj=2,kk=3)+erfplot(pred_DLAG,jj=2,kk=3))
ggsave(paste0("Results/Plots/chicagoNMMAPSpolar_erf2_lag",maxlag,".pdf"),width=6,height=8)

