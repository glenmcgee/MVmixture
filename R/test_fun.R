### Testing MVmix for DLM
set.seed(2)
source("MVmix.R")
n <- 100
K <- 8
p <- 20
X <- matrix(rnorm(n*p),ncol=p)
theta1 <- rep(1,p); theta1 <- theta1/sqrt(sum(theta1^2)) ##
theta2 <- 0.2+dnorm(seq(-3,3,length=p)); theta2 <- theta2/sqrt(sum(theta2^2)) ## DLM weights
x1th <- X%*%theta1
x2th <- X%*%theta2
f1 <- 3*sin(x1th)
f2 <- 5*cos(x2th)
fmat <- matrix(c(rep(f1,K/2),rep(f2,K/2)),ncol=K,nrow=n,byrow=FALSE)
ranef <- matrix(rnorm(n,0,sqrt(0.01)),ncol=K,nrow=n,byrow=FALSE)
eps <- matrix(rnorm(n*K,0,sqrt(0.5)),ncol=K)
Y <- fmat+ranef+eps

## works well with normal approximation
RR <- 2000
set.seed(1)
test1 <- MVmix(Y,X,Z=NULL,niter=RR,nburn=0.5*RR,nthin=5,
              Vgridsearch = TRUE,gridsize=10,
              maxClusters=5,DLM=TRUE,lagOrder = 3)

## doing this with beta MH sampling is tough to converge
# RR <- 8000
# set.seed(1)
# test1 <- MVmix(Y,X,Z=NULL,niter=RR,nburn=0.5*RR,nthin=5,
#                Vgridsearch = TRUE,gridsize=10,
#                maxClusters=5,DLM=TRUE,lagOrder = 3,
#                thetaMethod="MH_beta",prior_omega_a=200)

pred <- predict_MVmix(test1)
# pred <- predict_MVmix(test1,fixtheta = TRUE)

par(mfrow=c(1,2))
par(mar = c(1,2,1,0.5))
plot(f1~x1th)
points(pred$summary[[1]]$mean~x1th,col="blue")
plot(f2~x2th)
points(pred$summary[[5]]$mean~x2th,col="blue")

boxplot(test1$Zbeta,ylim=c(1,5))
boxplot(test1$Ztheta,ylim=c(1,5))
par(mfrow=c(1,1))
# dev.off()

# for(jj in 1:K){
#   boxplot(test1$beta[,(jj-1)*10+1:10])
# }

for(jj in 1:K){
  boxplot(test1$theta[,(jj-1)*p+1:p],ylim=c(min(test1$theta),max(test1$theta)))
  if(jj<=(K/2)){
    points(theta1,col="red",pch="x")
  }else{
    points(theta2,col="red",pch="x")
  }
}


### Testing MVmix for non-DLM
set.seed(1)
source("MVmix.R")
n <- 100
K <- 8
p <- 4
X <- matrix(rnorm(n*p),ncol=p)
theta1 <- rep(1,p); theta1 <- theta1/sqrt(sum(theta1^2))
theta2 <- c(rep(c(1,0),each=p/2)); theta2 <- theta2/sqrt(sum(theta2^2)) #theta2 <- theta1 #
x1th <- X%*%theta1 ## TESTING
x2th <- X%*%theta2
f1 <- 3*sin(x1th)
f2 <- 5*cos(x2th)
fmat <- matrix(c(rep(f1,K/2),rep(f2,K/2)),ncol=K,nrow=n,byrow=FALSE)
ranef <- matrix(rnorm(n,0,sqrt(0.5)),ncol=K,nrow=n,byrow=FALSE)
eps <- matrix(rnorm(n*K,0,sqrt(0.5)),ncol=K)
Y <- fmat+ranef+eps


set.seed(1)
RR <- 4000
test2 <- MVmix(Y,X,Z=NULL,niter=RR,nburn=0.5*RR,nthin=2,
               Vgridsearch = TRUE,gridsize=10,
               maxClusters=5,DLM=FALSE,sharedlambda=TRUE,
               thetaMethod = "MH_beta",prior_omega_a = 1000)

pred <- predict_MVmix(test2)
# pred <- predict_MVmix(test2,fixtheta = TRUE)

par(mfrow=c(1,2))
par(mar = c(1,1,1,1))
plot(f1~x1th)
points(pred$summary[[1]]$mean~x1th,col="blue")
plot(f2~x2th)
points(pred$summary[[5]]$mean~x2th,col="blue")

boxplot(test2$Zbeta,ylim=c(1,5))
boxplot(test2$Ztheta,ylim=c(1,5))
par(mfrow=c(1,1))
# dev.off()

# for(jj in 1:K){
#   boxplot(test2$beta[,(jj-1)*10+1:10])
# }

for(jj in 1:K){
  boxplot(test2$theta[,(jj-1)*p+1:p],ylim=c(min(test2$theta),max(test2$theta)))
  if(jj<=(K/2)){
    points(theta1,col="red",pch="x")
  }else{
    points(theta2,col="red",pch="x")
  }
}



### Testing MVmix for multiple exposures
set.seed(100)
source("MVmix.R")
n <- 100
K <- 4
p <- 4
X1 <- matrix(rnorm(n*p),ncol=p)
X2 <- matrix(rnorm(n*p),ncol=p)
theta1 <- rep(1,p); theta1 <- theta1/sqrt(sum(theta1^2))
theta2 <- c((1:p)^1.5); theta2 <- theta2/sqrt(sum(theta2^2))
theta3 <- c((p:1)^1.5); theta3 <- theta3/sqrt(sum(theta3^2))
theta4 <- rep(c(0,1),p/2); theta4 <- theta4/sqrt(sum(theta4^2))
x11th <- X1%*%theta1 ## k=1 j=1
x12th <- X2%*%theta2 ## k=1 j=2
x21th <- X1%*%theta3 ## k=2 j=1
x22th <- X2%*%theta4 ## k=2 j=2
f11 <- sin(x11th) ## k=1 j=1
f12 <- 2*cos(x12th) ## k=1 j=2
f21 <- sin(x21th) ## k=2 j=1
f22 <- 2*cos(x22th) ## k=2 j=2
fmat <- matrix(c(rep(f11,K/2),rep(f21,K/2))+c(rep(f12,K/2),rep(f22,K/2)),ncol=K,nrow=n,byrow=FALSE)
ranef <- matrix(rnorm(n,0,sqrt(0.00001)),ncol=K,nrow=n,byrow=FALSE)
eps <- matrix(rnorm(n*K,0,sqrt(0.05)),ncol=K)
Y <- fmat+ranef+eps


set.seed(1)
RR <- 4000
test3 <- MVmix(Y,X=list(X1,X2),Z=NULL,niter=RR,nburn=0.5*RR,nthin=2,
               Vgridsearch = TRUE,gridsize=10,
               maxClusters=4,DLM=FALSE,sharedlambda=FALSE)

test3$beta <- assign_betas(test3)
test3$theta <- assign_thetas(test3)
summarize_clusters(test3)


pred <- predict_MVmix(test3)


lm(c(mean(-test3$b0[,1])+pred$summary[[1]][[1]]$mean+pred$summary[[1]][[2]]$mean)~c(f11+f12))
plot(c(mean(-test3$b0[,1])+pred$summary[[1]][[1]]$mean+pred$summary[[1]][[2]]$mean)~c(f11+f12))
abline(0,1,col="red")
lm(c(mean(-test3$b0[,3])+pred$summary[[3]][[1]]$mean+pred$summary[[3]][[2]]$mean)~c(f21+f22))
plot(c(mean(-test3$b0[,3])+pred$summary[[3]][[1]]$mean+pred$summary[[3]][[2]]$mean)~c(f21+f22))
abline(0,1,col="red")



for(kk in 1:test3$const$K){
  for(jj in 1:test3$const$p){
    boxplot(test3$theta[[jj]][,(kk-1)*p+1:p],
            ylim=c(-1,1),main=paste0("k=",kk,",  j=",jj))
  }
}
