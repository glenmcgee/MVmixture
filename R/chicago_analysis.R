library(gamair)
library(ggplot2)
library(patchwork)
library(reshape)
library(mgcv)
data(chicago)
maxlag <- 7
nit <- 30000
nburn <- 0.5*nit
nthin = 5


## month and year
n_years <- ceiling(nrow(chicago)/365)
month <- rep(rep(1:12,each=30), times = n_years)
chicago$month <- month[1:nrow(chicago)]
year <- rep(1:(n_years), each = 12*30)
chicago$year <- year[1:nrow(chicago)]

# ## remove pm25median (mostly missing)
# chicago <- chicago[,-3]

## carry forward exposure values if missing
chicago$pm25median[1] <- 0 ## mostly missing until 2000
for(ii in 2:nrow(chicago)){
  if(is.na(chicago$pm25median[ii])){chicago$pm25median[ii] <- chicago$pm25median[ii-1]}
  if(is.na(chicago$pm10median[ii])){chicago$pm10median[ii] <- chicago$pm10median[ii-1]}
  if(is.na(chicago$o3median[ii])){chicago$o3median[ii] <- chicago$o3median[ii-1]}
  if(is.na(chicago$so2median[ii])){chicago$so2median[ii] <- chicago$so2median[ii-1]}
}

## make lagged exposures
pm25 <- pm10 <- o3 <- so2 <- tmpd <- matrix(NA,ncol=maxlag,nrow=nrow(chicago))
for(ii in (maxlag):nrow(chicago)){
  pm25[ii,] <- c(chicago$pm25median[(ii-maxlag+1):ii])
  pm10[ii,] <- c(chicago$pm10median[(ii-maxlag+1):ii])
  o3[ii,] <- c(chicago$o3median[(ii-maxlag+1):ii])
  so2[ii,] <- c(chicago$so2median[(ii-maxlag+1):ii])
  tmpd[ii,] <- c(chicago$tmpd[(ii-maxlag+1):ii])
}
## exclude outliers!
outliers <- which(chicago$death>200)
chicago <- chicago[-outliers,,drop=F]
pm25 <- scale(pm25[-outliers,,drop=F])
pm10 <- scale(pm10[-outliers,,drop=F])#+40)
o3 <- log(o3[-outliers,,drop=F]+30)
so2 <- scale(so2[-outliers,,drop=F])#+10)
tmpd <- tmpd[-outliers,,drop=F]

Xlist <- list(pm10[maxlag:nrow(chicago),,drop=F],
              o3[maxlag:nrow(chicago),,drop=F],
              so2[maxlag:nrow(chicago),,drop=F])
XTemplist <- list(pm10[maxlag:nrow(chicago),,drop=F],
                  o3[maxlag:nrow(chicago),,drop=F],
                  so2[maxlag:nrow(chicago),,drop=F],
                  tmpd[maxlag:nrow(chicago),,drop=F])
X2000list <- list(tail(pm25,366-maxlag),
                  tail(pm10,366-maxlag),
                  tail(o3,366-maxlag),
                  tail(so2,366-maxlag))


## get spline terms for time trends (and tmpd)
prelim_gam <- gam(log(death) ~ s(year,k=4)+
                                s(month, bs = "cc",k=4)+
                                s(time,k=5)+
                                s(tmpd,k=5),
                  method = "REML",data=chicago)

## prep covariates
Ztime <- model.matrix(prelim_gam)[,-1]
Ztime <- Ztime[maxlag:nrow(chicago),]
ZTemptime <- Ztime[,1:9] ## excluding temp that is already in XTemp

prelim2000_gam <- gam(log(death) ~s(month, bs = "cc",k=4)+
                    s(time,k=5)+
                    s(tmpd,k=5),
                  method = "REML",data=tail(chicago,366-maxlag))
Z2000time <- model.matrix(prelim2000_gam)[,-1]

## prep outcome
Y <- log(chicago$death[maxlag:nrow(chicago)])
Y2000 <- tail(Y,366-maxlag)

## prep new X for predictions
npred <- 50
seq0 <- seq(quantile(X2000list[[1]],0.1),quantile(X2000list[[1]],0.9),length=npred)
seq1 <- seq(quantile(Xlist[[1]],0.1),quantile(Xlist[[1]],0.9),length=npred)
seq2 <- seq(quantile(Xlist[[2]],0.1),quantile(Xlist[[2]],0.9),length=npred)
seq3 <- seq(quantile(Xlist[[3]],0.1),quantile(Xlist[[3]],0.9),length=npred)
seq4 <- seq(quantile(XTemplist[[4]],0.1),quantile(XTemplist[[4]],0.9),length=npred)
Xnew <- list(matrix(seq1,ncol=maxlag,nrow=npred,byrow=F),
             matrix(seq2,ncol=maxlag,nrow=npred,byrow=F),
             matrix(seq3,ncol=maxlag,nrow=npred,byrow=F))
XTempnew <- list(matrix(seq1,ncol=maxlag,nrow=npred,byrow=F),
                 matrix(seq2,ncol=maxlag,nrow=npred,byrow=F),
                 matrix(seq3,ncol=maxlag,nrow=npred,byrow=F),
                 matrix(seq4,ncol=maxlag,nrow=npred,byrow=F))
X2000new <- list(matrix(seq0,ncol=maxlag,nrow=npred,byrow=F),
                 matrix(seq1,ncol=maxlag,nrow=npred,byrow=F),
                 matrix(seq2,ncol=maxlag,nrow=npred,byrow=F),
                 matrix(seq3,ncol=maxlag,nrow=npred,byrow=F))

#################################
## Plotting functions
#################################
spaghettiplot <- function(x){

  # plot(x[1,],type="l",col=alpha("black", .1),ylim=c(-1,1),ylab="w")
  # for(ii in 2:nrow(x)){
  #   lines(x[ii,],col=alpha("black", .1))
  # }

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

erfplot <- function(pred,jj,
                    y2k=FALSE,## y2k analysis with pm2.5
                    ylims=NULL){
  if(y2k==TRUE){
    grid <- list(seq0,seq1,seq2,seq3)[[jj]]
  }else{
    grid <- list(seq1,seq2,seq3,seq4)[[jj]]
  }

  # plot(pred$summary[[1]][[jj]]$mean~grid,type='l',
  #      ylim=c(min(pred$summary[[1]][[jj]]),max(pred$summary[[1]][[jj]])))
  # lines(pred$summary[[1]][[jj]]$lower~grid,lty=2)
  # lines(pred$summary[[1]][[jj]]$upper~grid,lty=2)

  df <- data.frame(
    gridx=grid,
    fitted=pred$summary[[1]][[jj]]$mean,#predict(loess(cbmim_pred_ind[[jj]]$mean~seq(-2,2,length=21))),
    uci=pred$summary[[1]][[jj]]$upper,#predict(loess(cbmim_pred_ind[[jj]]$uci~seq(-2,2,length=21))),
    lci=pred$summary[[1]][[jj]]$lower)#predict(loess(cbmim_pred_ind[[jj]]$lci~seq(-2,2,length=21))) )

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







#######################################
## Fit DLNM
#######################################
set.seed(0)

chicago_DLAG <- MVmix(as.matrix(Y),Xlist,Z=Ztime,
              niter=nit,nburn=nburn,nthin=nthin,
              Vgridsearch = TRUE,gridsize=10,
              DLM=TRUE,DLMpenalty=TRUE,lagOrder=NULL,diff=2,
              cluster="both",maxClusters=3,approx=TRUE)
save(chicago_DLAG, file = paste0("Results/chicago_DLAG",maxlag,".RData"))

pred_DLAG <- predict_MVmix(chicago_DLAG,
                           newX = Xnew,
                           # newZ = Ztime,
                           include_intercept=FALSE,
                           fixomega = FALSE,
                           allx=FALSE,contrast=TRUE)

erfplot(pred_DLAG,1)
erfplot(pred_DLAG,2)
erfplot(pred_DLAG,3)

spaghettiplot(chicago_DLAG$omega[[1]][[1]])
spaghettiplot(chicago_DLAG$omega[[2]][[1]])
spaghettiplot(chicago_DLAG$omega[[3]][[1]])

pc_DLAG <- pairwise_clusters(chicago_DLAG)
make_heatplot(pc_DLAG$beta_y)+
make_heatplot(pc_DLAG$theta_y)
ggsave(paste0("Results/Plots/chicagoHeat_lag",maxlag,".pdf"),width=8,height=5)


#######################################
## Fit DLNM --- No clustering
#######################################
set.seed(0)

chicago_DLAG_noclust <- MVmix(as.matrix(Y),Xlist,Z=Ztime,
                      niter=nit,nburn=nburn,nthin=nthin,
                      Vgridsearch = TRUE,gridsize=10,
                      DLM=TRUE,DLMpenalty=TRUE,lagOrder=NULL,diff=2,
                      cluster="neither",maxClusters=3,approx=TRUE)
save(chicago_DLAG_noclust, file = paste0("Results/chicago_DLAG_noclust",maxlag,".RData"))

pred_DLAG_noclust <- predict_MVmix(chicago_DLAG_noclust,
                           newX = Xnew,
                           include_intercept=FALSE,
                           fixomega = FALSE,
                           allx=FALSE,contrast=TRUE)

erfplot(pred_DLAG_noclust,1)
erfplot(pred_DLAG_noclust,2)
erfplot(pred_DLAG_noclust,3)

spaghettiplot(chicago_DLAG_noclust$omega[[1]][[1]])
spaghettiplot(chicago_DLAG_noclust$omega[[2]][[1]])
spaghettiplot(chicago_DLAG_noclust$omega[[3]][[1]])


### PLOTS
## maxlag=14
if(maxlag==1){
  lims=list(c(-0.025,0.021),c(-0.01,0.026),c(-0.006,0.006))
}else if(maxlag==7){
  lims=list(c(-0.03,0.05),c(-0.04,0.07),c(-0.02,0.02))
}else if(maxlag==14){
  lims=list(c(-0.05,0.05),c(-0.04,0.07),c(-0.02,0.02))
}else if(maxlag==28){
  lims=list(c(-0.05,0.16),c(-0.07,0.07),c(-0.02,0.035))
}
(erfplot(pred_DLAG,1,ylims=lims[[1]])+
erfplot(pred_DLAG_noclust,1,ylims=lims[[1]]))/
(erfplot(pred_DLAG,2,ylims=lims[[2]])+
erfplot(pred_DLAG_noclust,2,ylims=lims[[2]]))/
(erfplot(pred_DLAG,3,ylims=lims[[3]])+
erfplot(pred_DLAG_noclust,3,ylims=lims[[3]]))
ggsave(paste0("Results/Plots/chicago_lag",maxlag,".pdf"),width=8,height=10)

(spaghettiplot(chicago_DLAG$omega[[1]][[1]])+
  spaghettiplot(chicago_DLAG_noclust$omega[[1]][[1]]))/
(spaghettiplot(chicago_DLAG$omega[[2]][[1]])+
  spaghettiplot(chicago_DLAG_noclust$omega[[2]][[1]]))/
(spaghettiplot(chicago_DLAG$omega[[3]][[1]])+
spaghettiplot(chicago_DLAG_noclust$omega[[3]][[1]]))
ggsave(paste0("Results/Plots/chicagoWeights_lag",maxlag,".pdf"),width=8,height=10)



#######################################
## Fit DLNM -- YEAR 2000 (with pm2.5)
#######################################
set.seed(0)

chicago2000_DLAG <- MVmix(as.matrix(Y2000),X2000list,Z=Z2000time,
                      niter=nit,nburn=nburn,nthin=nthin,
                      Vgridsearch = TRUE,gridsize=10,
                      DLM=TRUE,DLMpenalty=TRUE,lagOrder=NULL,diff=2,
                      cluster="both",maxClusters=4,approx=TRUE)
save(chicago2000_DLAG, file = paste0("Results/chicago2000_DLAG",maxlag,".RData"))

pred2000_DLAG <- predict_MVmix(chicago2000_DLAG,
                           newX = X2000new,
                           # newZ = Ztime,
                           include_intercept=FALSE,
                           fixomega = FALSE,
                           allx=FALSE,contrast=TRUE)

erfplot(pred2000_DLAG,1,y2k=TRUE)
erfplot(pred2000_DLAG,2,y2k=TRUE)
erfplot(pred2000_DLAG,3,y2k=TRUE)
erfplot(pred2000_DLAG,4,y2k=TRUE)

spaghettiplot(chicago2000_DLAG$omega[[1]][[1]])
spaghettiplot(chicago2000_DLAG$omega[[2]][[1]])
spaghettiplot(chicago2000_DLAG$omega[[3]][[1]])
spaghettiplot(chicago2000_DLAG$omega[[4]][[1]])

pc2000_DLAG <- pairwise_clusters(chicago2000_DLAG)
make_heatplot(pc2000_DLAG$beta_y)+
make_heatplot(pc2000_DLAG$theta_y)
ggsave(paste0("Results/Plots/chicago2000Heat_lag",maxlag,".pdf"),width=8,height=5)


#############################################################
## Fit DLNM --- No clustering -- YEAR 2000 (with pm2.5)
#############################################################
set.seed(0)

chicago2000_DLAG_noclust <- MVmix(as.matrix(Y2000),X2000list,Z=Z2000time,
                              niter=nit,nburn=nburn,nthin=nthin,
                              Vgridsearch = TRUE,gridsize=10,
                              DLM=TRUE,DLMpenalty=TRUE,lagOrder=NULL,diff=2,
                              cluster="neither",maxClusters=4,approx=TRUE)
save(chicago2000_DLAG_noclust, file = paste0("Results/chicago2000_DLAG_noclust",maxlag,".RData"))

pred2000_DLAG_noclust <- predict_MVmix(chicago2000_DLAG_noclust,
                                   newX = X2000new,
                                   include_intercept=FALSE,
                                   fixomega = FALSE,
                                   allx=FALSE,contrast=TRUE)

erfplot(pred2000_DLAG_noclust,1,y2k=TRUE)
erfplot(pred2000_DLAG_noclust,2,y2k=TRUE)
erfplot(pred2000_DLAG_noclust,3,y2k=TRUE)
erfplot(pred2000_DLAG_noclust,4,y2k=TRUE)

spaghettiplot(chicago2000_DLAG_noclust$omega[[1]][[1]])
spaghettiplot(chicago2000_DLAG_noclust$omega[[2]][[1]])
spaghettiplot(chicago2000_DLAG_noclust$omega[[3]][[1]])
spaghettiplot(chicago2000_DLAG_noclust$omega[[4]][[1]])


## make plots
##
if(maxlag==7){
  lims <- list(c(-4,4),c(-1,1),c(-1,1),c(-1,1))
}else if(maxlag==14){
  lims <- list(c(-0.25,0.4),c(-0.6,0.6),c(-0.2,0.2),c(-0.1,0.15))
}else if(maxlag==28){
  lims <- list(c(-1.3,0.9),c(-2.2,2.4),c(-1.4,1.2),c(-0.2,0.5))
}
(erfplot(pred2000_DLAG,1,ylims=lims[[1]])+
erfplot(pred2000_DLAG_noclust,1,ylims=lims[[1]]))/
(erfplot(pred2000_DLAG,2,ylims=lims[[2]])+
erfplot(pred2000_DLAG_noclust,2,ylims=lims[[2]]))/
(erfplot(pred2000_DLAG,3,ylims=lims[[3]])+
erfplot(pred2000_DLAG_noclust,3,ylims=lims[[3]]))/
(erfplot(pred2000_DLAG,4,ylims=lims[[4]])+
erfplot(pred2000_DLAG_noclust,4,ylims=lims[[4]]))
ggsave(paste0("Results/Plots/chicago2000_lag",maxlag,".pdf"),width=8,height=12)

(spaghettiplot(chicago2000_DLAG$omega[[1]][[1]])+
    spaghettiplot(chicago2000_DLAG_noclust$omega[[1]][[1]]))/
  (spaghettiplot(chicago2000_DLAG$omega[[2]][[1]])+
     spaghettiplot(chicago2000_DLAG_noclust$omega[[2]][[1]]))/
  (spaghettiplot(chicago2000_DLAG$omega[[3]][[1]])+
     spaghettiplot(chicago2000_DLAG_noclust$omega[[3]][[1]]))/
  (spaghettiplot(chicago2000_DLAG$omega[[4]][[1]])+
     spaghettiplot(chicago2000_DLAG_noclust$omega[[4]][[1]]))
ggsave(paste0("Results/Plots/chicago2000Weights_lag",maxlag,".pdf"),width=8,height=12)

#
# #######################################
# ## Fit DLNM + TEMP
# #######################################
# set.seed(0)
#
chicagoTemp_DLAG <- MVmix(as.matrix(Y),list(XTemplist[[1]],XTemplist[[3]],XTemplist[[4]]),,Z=ZTemptime,
                      niter=nit,nburn=nburn,nthin=nthin,
                      Vgridsearch = TRUE,gridsize=10,
                      DLM=TRUE,DLMpenalty=TRUE,lagOrder=NULL,diff=2,
                      cluster="both",maxClusters=3,approx=TRUE)
save(chicagoTemp_DLAG, file = paste0("Results/chicagoTemp_DLAG",maxlag,".RData"))
#
predTemp_DLAG <- predict_MVmix(chicagoTemp_DLAG,
                           newX = list(XTempnew[[1]],XTempnew[[3]],XTempnew[[4]]),
                           include_intercept=FALSE,
                           fixomega = FALSE,
                           allx=FALSE,contrast=TRUE)


erfplot(predTemp_DLAG,1)
erfplot(predTemp_DLAG,2)
erfplot(predTemp_DLAG,3)

spaghettiplot(chicagoTemp_DLAG$omega[[1]][[1]])
spaghettiplot(chicagoTemp_DLAG$omega[[2]][[1]])
spaghettiplot(chicagoTemp_DLAG$omega[[3]][[1]])

pcTemp_DLAG <- pairwise_clusters(chicagoTemp_DLAG)
make_heatplot(pcTemp_DLAG$beta_y)
make_heatplot(pcTemp_DLAG$theta_y)


# #######################################
# ## Fit DLNM + TEMP --- No clustering
# #######################################
# set.seed(0)
#
# chicagoTemp_DLAG_noclust <- MVmix(as.matrix(Y),XTemplist,Z=ZTemptime,
#                               niter=nit,nburn=nburn,nthin=nthin,
#                               Vgridsearch = TRUE,gridsize=10,
#                               DLM=TRUE,DLMpenalty=TRUE,lagOrder=8,diff=4,
#                               cluster="neither",maxClusters=3,approx=TRUE)
# save(chicagoTemp_DLAG_noclust, file = paste0("Results/chicagoTemp_DLAG_noclust",maxlag,".RData"))
#
# predTemp_DLAG_noclust <- predict_MVmix(chicagoTemp_DLAG_noclust,
#                                    newX = XTempnew,
#                                    include_intercept=FALSE,
#                                    fixomega = TRUE,
#                                    allx=FALSE)
#
# erfplot(predTemp_DLAG_noclust,1)
# erfplot(predTemp_DLAG_noclust,2)
# erfplot(predTemp_DLAG_noclust,3)
#
# spaghettiplot(chicagoTemp_DLAG_noclust$omega[[1]][[1]])
# spaghettiplot(chicagoTemp_DLAG_noclust$omega[[2]][[1]])
# spaghettiplot(chicagoTemp_DLAG_noclust$omega[[3]][[1]])
#
# pcTemp_DLAG_noclust <- pairwise_clusters(chicagoTemp_DLAG_noclust)
# make_heatplot(pcTemp_DLAG_noclust$beta_y)
# make_heatplot(pcTemp_DLAG_noclust$theta_y)



#######################################
## Fit GAMs as sanity check
#######################################


test_gam <- gam(log(death) ~
                    s(pm10median)+
                    s(o3median)+
                    s(so2median)+
                    s(year,k=4)+
                    s(month, bs = "cc",k=4)+
                    s(time,k=5)+
                    s(tmpd,k=5),
                  method = "REML",data=chicago)
# plot(test_gam,xlim=c(-20,25),ylim=c(-0.03,0.07)) ## 1
# plot(test_gam,xlim=c(-15,12),ylim=c(-0.02,0.04)) ## 2
# plot(test_gam,xlim=c(-4,3),ylim=c(-0.04,0.06)) ## 3

## unpenalized time and tmpd
tempdat <- data.frame(Y,
                           pm10median=Xlist[[1]][,maxlag],
                           o3median=Xlist[[2]][,maxlag],
                           so2median=Xlist[[3]][,maxlag],
                           Ztime)
newdat <- data.frame(pm10median=seq1,
           o3median=seq2,
           so2median=seq3,
           zz=matrix(0,nrow=length(seq1),ncol=13)
)
colnames(tempdat) <- c("Y","pm10median","o3median","so2median",
              paste0("Z",1:ncol(Ztime)))
colnames(newdat) <- c("pm10median","o3median","so2median",
                      paste0("Z",1:ncol(Ztime)))
simple_gam <- gam(Y ~ s(pm10median)+
                  s(o3median)+
                  s(so2median)+
                  Z1+Z2+Z3+Z4+Z5+Z6+Z7+Z8+Z9+Z10+Z11+Z12+Z13,
                method = "REML",data=tempdat)

pred_simple <- predict.gam(simple_gam,newdata=newdat,type="iterms",se=TRUE)
plot(pred_simple$fit[,"s(pm10median)"]~seq1,type="l")
lines(pred_simple$fit[,"s(pm10median)"]+1.96*pred_simple$se[,"s(pm10median)"]~seq1,lty=2)
lines(pred_simple$fit[,"s(pm10median)"]-1.96*pred_simple$se[,"s(pm10median)"]~seq1,lty=2)

plot(pred_simple$fit[,"s(o3median)"]~seq1,type="l")
lines(pred_simple$fit[,"s(o3median)"]+1.96*pred_simple$se[,"s(o3median)"]~seq1,lty=2)
lines(pred_simple$fit[,"s(o3median)"]-1.96*pred_simple$se[,"s(o3median)"]~seq1,lty=2)

plot(pred_simple$fit[,"s(so2median)"]~seq1,type="l",ylim=c(-0.014,0.006),xlim=c(-4,3))
lines(pred_simple$fit[,"s(so2median)"]+1.96*pred_simple$se[,"s(so2median)"]~seq1,lty=2)
lines(pred_simple$fit[,"s(so2median)"]-1.96*pred_simple$se[,"s(so2median)"]~seq1,lty=2)


test2000_gam <- gam(log(death) ~
                      s(pm25median)+
                      s(pm10median)+
                      s(o3median)+
                      s(so2median)+
                      s(month, bs = "cc",k=4)+
                      s(time,k=5)+
                      s(tmpd,k=5),
                method = "REML",data=tail(chicago,365))



# #############################################################
# ## Fit DLNM --- No clustering -- YEAR 2000 (with pm2.5)
# #############################################################
# set.seed(0)
# ## clearly it works
#
# chicago2000_1day_noclust <- MVmix(as.matrix(Y2000),lapply(X2000list,function(x)as.matrix(x[,maxlag])),Z=Z2000time,
#                                   niter=nit,nburn=nburn,nthin=nthin,
#                                   Vgridsearch = TRUE,gridsize=10,
#                                   DLM=TRUE,DLMpenalty=TRUE,lagOrder=8,diff=4,
#                                   cluster="neither",maxClusters=4,approx=TRUE)
# save(chicago2000_1day_noclust, file = paste0("Results/chicago2000_1day_noclust",maxlag,".RData"))
#
# pred2000_1day_noclust <- predict_MVmix(chicago2000_1day_noclust,
#                                        newX = lapply(X2000new,function(x)as.matrix(x[,maxlag])),
#                                        include_intercept=FALSE,
#                                        fixomega = TRUE,
#                                        allx=FALSE)
#
# erfplot(pred2000_1day_noclust,1,y2k=TRUE)
# erfplot(pred2000_1day_noclust,2,y2k=TRUE)
# erfplot(pred2000_1day_noclust,3,y2k=TRUE)
# erfplot(pred2000_1day_noclust,4,y2k=TRUE)
#
# spaghettiplot(chicago2000_1day_noclust$omega[[1]][[1]])
# spaghettiplot(chicago2000_1day_noclust$omega[[2]][[1]])
# spaghettiplot(chicago2000_1day_noclust$omega[[3]][[1]])
#
# pc2000_1day_noclust <- pairwise_clusters(chicago2000_1day_noclust)
# make_heatplot(pc2000_1day_noclust$beta_y)
# make_heatplot(pc2000_1day_noclust$theta_y)
#
#
