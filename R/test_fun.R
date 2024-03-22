###########################################
### Testing MVmix for multiple exposures
###########################################
set.seed(100)
source("MVmix.R")
n <- 100
K <- 4
L <- 4
X1 <- matrix(rnorm(n*L),ncol=L)
X2 <- matrix(rnorm(n*L),ncol=L)
w1 <- rep(1,L); w1 <- w1/sqrt(sum(w1^2))
w2 <- c((1:L)^1.5); w2 <- w2/sqrt(sum(w2^2))
w3 <- c((L:1)^1.5); w3 <- w3/sqrt(sum(w3^2))
w4 <- rep(c(0,1),L/2); w4 <- w4/sqrt(sum(w4^2))
x11th <- X1%*%w1 ## k=1:2 j=1
x12th <- X2%*%w2 ## k=1:2 j=2
x21th <- X1%*%w3 ## k=3:4 j=1
x22th <- X2%*%w4 ## k=3:4 j=2
f11 <- sin(x11th) ## k=1:2 j=1
f12 <- sin(x12th) ## k=1:2 j=2
f21 <- 2*cos(x21th) ## k=3:4 j=1
f22 <- 2*cos(x22th) ## k=3:4 j=2
fmat <- matrix(c(rep(f11,K/2),rep(f21,K/2))+c(rep(f12,K/2),rep(f22,K/2)),ncol=K,nrow=n,byrow=FALSE)
ranef <- matrix(rnorm(n,0,sqrt(0.00001)),ncol=K,nrow=n,byrow=FALSE)
eps <- matrix(rnorm(n*K,0,sqrt(0.5)),ncol=K)
Y <- fmat+ranef+eps


set.seed(1)
RR <- 4000
test1 <- MVmix(Y,X=list(X1,X2),Z=NULL,niter=RR,nburn=0.5*RR,nthin=2,
               Vgridsearch = TRUE,gridsize=10,
               maxClusters=4,DLM=FALSE,
               approx=TRUE) ## approx=TRUE means do the taylor approx version

## cluster membership
summarize_clusters(test1)

## plot weights
test1$beta <- assign_betas(test1)
test1$omega <- assign_omegas(test1)
for(kk in 1:test1$const$K){
  for(jj in 1:test1$const$p){
    boxplot(test1$omega[[jj]][,(kk-1)*L+1:L],
            ylim=c(-1,1),main=paste0("k=",kk,",  j=",jj))
  }
}


## fitted vals
pred <- predict_MVmix(test1,include_intercept=FALSE)
## intercept=FALSE--will add back in manually since we are summing two different fns

lm(c(mean(test1$b0[,1])+pred$summary[[1]][[1]]$mean+pred$summary[[1]][[2]]$mean)~c(f11+f12))
plot(c(mean(test1$b0[,1])+pred$summary[[1]][[1]]$mean+pred$summary[[1]][[2]]$mean)~c(f11+f12))
abline(0,1,col="red")

lm(c(mean(test1$b0[,(K/2+1)])+pred$summary[[(K/2+1)]][[1]]$mean+pred$summary[[(K/2+1)]][[2]]$mean)~c(f21+f22))
plot(c(mean(test1$b0[,(K/2+1)])+pred$summary[[(K/2+1)]][[1]]$mean+pred$summary[[(K/2+1)]][[2]]$mean)~c(f21+f22))
abline(0,1,col="red")



###########################################
### Testing MVmix for single exposure DLM
###########################################
set.seed(2)
source("MVmix.R")
n <- 100
K <- 8
L <- 52
# X1 <- matrix(rnorm(n*L),ncol=L)
X1 <- matrix(0,ncol=L,nrow=n)
X1[,1] <- rnorm(n)
for(ii in 2:L){X1[,ii] <- X1[,ii-1]+rnorm(n)}
w1 <- pnorm(seq(-6,3,length=L)); w1 <- w1/sqrt(sum(w1^2)) ##
w2 <- pnorm(-seq(-6,3,length=L)); w2 <- w2/sqrt(sum(w2^2)) ## DLM weights
x11th <- X1%*%w1 ## k=1:4 j=1
x21th <- X1%*%w2 ## k=5:8 j=1
f11 <- sin(x11th/(0.5*L)) ## k=1:4 j=1
f21 <- cos(x21th/(0.5*L)) ## k=5:8 j=1
fmat <- matrix(c(rep(f11,K/2),rep(f21,K/2)),ncol=K,nrow=n,byrow=FALSE)
ranef <- matrix(rnorm(n,0,sqrt(0.01)),ncol=K,nrow=n,byrow=FALSE)
eps <- matrix(rnorm(n*K,0,sqrt(0.02)),ncol=K)
Y <- fmat+ranef+eps


RR <- 4000
set.seed(100)
testDLM1 <- MVmix(Y,list(X1),Z=NULL,niter=RR,nburn=0.5*RR,nthin=5,
                  Vgridsearch = TRUE,gridsize=10,
                  maxClusters=5,DLM=TRUE,sharedlambda = FALSE,lagOrder=3,
                  approx=FALSE) ## doing MH

summarize_clusters(testDLM1)

testDLM1$beta <- assign_betas(testDLM1)
testDLM1$omega <- assign_omegas(testDLM1)

for(kk in 1:testDLM1$const$K){
  for(jj in 1:testDLM1$const$p){
    boxplot(testDLM1$omega[[jj]][,(kk-1)*L+1:L],
            ylim=c(-1,1),main=paste0("k=",kk,",  j=",jj))
  }
}


pred <- predict_MVmix(testDLM1,include_intercept =TRUE)

par(mfrow=c(1,2))
lm(pred$summary[[1]][[1]]$mean~f11)
plot(pred$summary[[1]][[1]]$mean~f11)
abline(0,1,col="red")
lm(pred$summary[[(K/2+1)]][[1]]$mean~f21)
plot(pred$summary[[(K/2+1)]][[1]]$mean~f21)
abline(0,1,col="red")
par(mfrow=c(1,1))





#################################################
### Testing MVmix for multiple DLM -- not working
#################################################
set.seed(2)
source("MVmix.R")
n <- 100
K <- 4
L <- 20
# X1 <- matrix(rnorm(n*L),ncol=L)
X1 <- matrix(0,ncol=L,nrow=n)
X1[,1] <- rnorm(n)
for(ii in 2:L){X1[,ii] <- X1[,ii-1]+rnorm(n)}
# X2 <- matrix(rnorm(n*L),ncol=L)
X2 <- matrix(0,ncol=L,nrow=n)
X2[,1] <- rnorm(n)
for(ii in 2:L){X2[,ii] <- X2[,ii-1]+rnorm(n)}
w1 <- rep(1,L); w1 <- w1/sqrt(sum(w1^2)) ##
w2 <- pnorm(-seq(-3,3,length=L)); w2 <- w2/sqrt(sum(w2^2)) ## DLM weights
x11th <- X1%*%w1 ## k=1 j=1
x12th <- X2%*%w1 ## k=1 j=2
x21th <- X1%*%w2 ## k=2 j=1
x22th <- X2%*%w2 ## k=2 j=2
f11 <- sin(x11th/(0.5*L)) ## k=1 j=1
f12 <- sin(x12th/(0.5*L)) ## k=1 j=2
f21 <- 2*cos(x21th/(0.5*L)) ## k=2 j=1
f22 <- 2*cos(x22th/(0.5*L)) ## k=2 j=2
fmat <- matrix(c(rep(f11,K/2),rep(f21,K/2),rep(f12,K/2),rep(f22,K/2)),ncol=K,nrow=n,byrow=FALSE)
ranef <- matrix(rnorm(n,0,sqrt(0.01)),ncol=K,nrow=n,byrow=FALSE)
eps <- matrix(rnorm(n*K,0,sqrt(0.02)),ncol=K)
Y <- fmat+ranef+eps

##
RR <- 4000
set.seed(100)
testDLM2 <- MVmix(Y,list(X1,X2),Z=NULL,niter=RR,nburn=0.5*RR,nthin=5,
              Vgridsearch = TRUE,gridsize=10,
              maxClusters=5,DLM=TRUE,sharedlambda = FALSE,lagOrder=3,
              approx=F)

summarize_clusters(testDLM2)

testDLM2$beta <- assign_betas(testDLM2)
testDLM2$omega <- assign_omegas(testDLM2)

for(kk in 1:testDLM2$const$K){
  for(jj in 1:testDLM2$const$p){
    boxplot(testDLM2$omega[[jj]][,(kk-1)*L+1:L],
            ylim=c(-1,1),main=paste0("k=",kk,",  j=",jj))
  }
}


pred <- predict_MVmix(testDLM2,include_intercept =FALSE)

par(mfrow=c(1,2))
lm(c(mean(testDLM2$b0[,1])+pred$summary[[1]][[1]]$mean+pred$summary[[1]][[2]]$mean)~c(f11+f12))
plot(c(mean(testDLM2$b0[,1])+pred$summary[[1]][[1]]$mean+pred$summary[[1]][[2]]$mean)~c(f11+f12))
abline(0,1,col="red")
lm(c(mean(testDLM2$b0[,(K/2+1)])+pred$summary[[(K/2+1)]][[1]]$mean+pred$summary[[(K/2+1)]][[2]]$mean)~c(f21+f22))
plot(c(mean(testDLM2$b0[,(K/2+1)])+pred$summary[[(K/2+1)]][[1]]$mean+pred$summary[[(K/2+1)]][[2]]$mean)~c(f21+f22))
abline(0,1,col="red")
par(mfrow=c(1,1))

