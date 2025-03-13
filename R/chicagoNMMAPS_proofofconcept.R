library(tsModel)
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

## make lagged exposures
pm10 <- o3 <-  matrix(NA,ncol=maxlag,nrow=nrow(chicagoNMMAPS))
for(ii in (maxlag):nrow(chicagoNMMAPS)){
  pm10[ii,] <- c(chicagoNMMAPS$pm10[(ii-maxlag+1):ii])
  o3[ii,] <- c(chicagoNMMAPS$o3[(ii-maxlag+1):ii])
}

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
Ztime <- model.matrix(prelim_gam)[,-1]

df1 <- df2 <- 4
L1 <- L2 <- 14
nlag1 <- L1
nlag2 <- L2
C1 <- onebasis(0:(L1-1), fun = "bs", degree = df1, intercept = FALSE)
C2 <- onebasis(0:(L2-1), fun = "bs", degree = df2, intercept = FALSE)
class(C1) <- "matrix"
class(C2) <- "matrix"
pm10_con <- pm10%*%C1
o3_con <- o3%*%C1
C <- as.matrix(bdiag(C1, C2))

dat <- cbind(chicagoNMMAPS,pm10_con,o3_con,Ztime)

## test CGLM
mod_CDLM <- glm(death ~  Ztime +pm10_con+o3_con, 
                family = poisson(link = "log"), 
                control = glm.control(epsilon = 1e-10, maxit = 1000),
                na.action = na.omit,data=dat)

est <- C%*%coef(mod_CDLM)[-(1:(ncol(Ztime)+1))]
var <- C%*%vcov(mod_CDLM)[-(1:(ncol(Ztime)+1)), -(1:(ncol(Ztime)+1))]%*%t(C)

## lags
plot(100*(exp(10*est[14:1])-1),type="l");abline(h=0,lty=2)
lines(100*(exp(10*(est[14:1]-1.96*sqrt(diag(var)[14:1])))-1),col="gray")
lines(100*(exp(10*(est[14:1]+1.96*sqrt(diag(var)[14:1])))-1),col="gray")
plot(100*(exp(10*est[28:15])-1),type="l");abline(h=0,lty=2)
lines(100*(exp(10*(est[28:15]-1.96*sqrt(diag(var)[28:15])))-1),col="gray")
lines(100*(exp(10*(est[28:15]+1.96*sqrt(diag(var)[28:15])))-1),col="gray")

#### linear version

lm_CDLM <- lm(log(death) ~ Ztime +pm10_con+o3_con, 
                na.action = na.omit,data=dat)

est <- C%*%coef(lm_CDLM)[-(1:(ncol(Ztime)+1))]
var <- C%*%vcov(lm_CDLM)[-(1:(ncol(Ztime)+1)), -(1:(ncol(Ztime)+1))]%*%t(C)

## lags
plot(100*(exp(10*est[14:1])-1),type="l");abline(h=0,lty=2)
lines(100*(exp(10*(est[14:1]-1.96*sqrt(diag(var)[14:1])))-1),col="gray")
lines(100*(exp(10*(est[14:1]+1.96*sqrt(diag(var)[14:1])))-1),col="gray")
plot(100*(exp(10*est[28:15])-1),type="l");abline(h=0,lty=2)
lines(100*(exp(10*(est[28:15]-1.96*sqrt(diag(var)[28:15])))-1),col="gray")
lines(100*(exp(10*(est[28:15]+1.96*sqrt(diag(var)[28:15])))-1),col="gray")


pm10_lag <- Lag(chicagoNMMAPS$pm10,0:(maxlag-1))
o3_lag <- Lag(chicagoNMMAPS$o3,0:(maxlag-1))
L <- matrix(0:(maxlag-1),nrow(pm10_lag),ncol(pm10_lag),byrow=TRUE)
gasp <- gam(death~s(pm10_lag,L,bs="cb",k=10)+
              s(o3_lag,L,bs="cb",k=10)+
              s(year,k=4)+
              s(month, bs = "cc",k=4)+
              as.factor(dow)+
              s(temp,k=6)+
              s(avgtemp,k=3)+
              s(dptp,k=6)+
              s(avgdptp,k=3),family=poisson(),
            data=chicagoNMMAPS,method='REML')

pm10gasp <- crosspred("pm10_lag",gasp,cen=30)
plot(pm10gasp,ptype="overall")
o3gasp <- crosspred("o3_lag",gasp,cen=20)
plot(o3gasp,ptype="overall")
# predslgasp <- crosspred("Q",gasp,by=0.2,bylag=0.2,cen=20)


## get spline terms for time trends (and covariates)
gam1_pm10 <- gam(log(death) ~ s(avgpm10)+s(year,k=4)+
                    s(month, bs = "cc",k=4)+
                    # s(doy,k=5)+
                    as.factor(dow)+
                    s(temp,k=6)+
                    s(avgtemp,k=3)+
                    s(dptp,k=6)+
                    s(avgdptp,k=3),
                  method = "REML",data=chicagoNMMAPS)

gam1_o3 <- gam(log(death) ~ s(avgo3)+s(year,k=4)+
                   s(month, bs = "cc",k=4)+
                   # s(doy,k=5)+
                   as.factor(dow)+
                   s(temp,k=6)+
                   s(avgtemp,k=3)+
                   s(dptp,k=6)+
                   s(avgdptp,k=3),
                 method = "REML",data=chicagoNMMAPS)

gam2_pm10 <- gam(log(cvd) ~ s(avgpm10)+s(year,k=4)+
                   s(month, bs = "cc",k=4)+
                   # s(doy,k=5)+
                   as.factor(dow)+
                   s(temp,k=6)+
                   s(avgtemp,k=3)+
                   s(dptp,k=6)+
                   s(avgdptp,k=3),
                 method = "REML",data=chicagoNMMAPS)

gam2_o3 <- gam(log(cvd) ~ s(avgo3)+s(year,k=4)+
                   s(month, bs = "cc",k=4)+
                   # s(doy,k=5)+
                   as.factor(dow)+
                   s(temp,k=6)+
                   s(avgtemp,k=3)+
                   s(dptp,k=6)+
                   s(avgdptp,k=3),
                 method = "REML",data=chicagoNMMAPS)


gam3_pm10 <- gam(log(resp+0.5) ~ s(avgpm10)+s(year,k=4)+
                   s(month, bs = "cc",k=4)+
                   # s(doy,k=5)+
                   as.factor(dow)+
                   s(temp,k=6)+
                   s(avgtemp,k=3)+
                   s(dptp,k=6)+
                   s(avgdptp,k=3),
                 method = "REML",data=chicagoNMMAPS)

gam3_o3 <- gam(log(resp+0.5) ~ s(avgo3)+s(year,k=4)+
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
