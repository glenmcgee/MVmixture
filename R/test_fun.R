### Testing MVmix for DLM
set.seed(1)
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
ranef <- matrix(rnorm(n,0,sqrt(0.25)),ncol=K,nrow=n,byrow=FALSE)
eps <- matrix(rnorm(n*K,0,sqrt(0.5)),ncol=K)
Y <- fmat+ranef+eps


set.seed(1)
test2 <- MVmix(Y,X,Z=NULL,niter=2000,nburn=0.5,nthin=5,
              Vgridsearch = TRUE,gridsize=10,
              maxClusters=5,DLM=TRUE,lagOrder = 3)

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
ranef <- matrix(rnorm(n,0,sqrt(0.25)),ncol=K,nrow=n,byrow=FALSE)
eps <- matrix(rnorm(n*K,0,sqrt(0.5)),ncol=K)
Y <- fmat+ranef+eps


set.seed(1)
RR <- 5000
test2 <- MVmix(Y,X,Z=NULL,niter=RR,nburn=0.5*RR,nthin=2,
               Vgridsearch = TRUE,gridsize=10,
               maxClusters=5,DLM=FALSE,sharedlambda=TRUE,thetaMethod = "MH_vmf",stepsize_theta = 0.05)

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
