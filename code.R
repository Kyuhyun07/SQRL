#### Required library ####
library(devtools)
install_github("Kyuhyun07/qris")

library(qris)
library(quantreg)
library(survival)
library(nleqslv)
library(emplik)
library(MASS)

#### Computational code for table 1 ####
#### 1. Scenario and necessart assumption ####
# data size = 200
# Only beta0 effective
# Quantile 50%
# simulation = 2000
# eta = 200

# True Beta
#      beta_0    beta_1
# t_0=0 1.609438   0
# t_0=1 1.410748   0
# t_0=2 1.219403   0
# t_0=3 1.040613   0

# Assumptions
# exp(beta_0) = 5
exp.beta.initial.0 <- 5
# kappa = 2 (Shape parameter of Weibull distribution)
k <- 2
# calculated rho_0 given exp(beta_0) = 5
r.initial.0 <- (log(10/5))^(1/k)/exp.beta.initial.0

#### 2. Data Generation function ####
data.gen <- function(samplesize, censor){
  sim <- matrix(NA,samplesize,5)
  colnames(sim) <- c("T","C","Z","X","censored")
  # Generate C_i
  sim[,2] <- runif(samplesize,0,censor)
  # Covariates (Control=0, Treatment=1)
  sim[,4] <- rbinom(samplesize,size=1,p=0.5)
  # Generate T_i (Given Condition r=rho_0, k=2, exp(beta_0)=5, exp(beta_1)=2))
  unif <- runif(n=samplesize ,min = 0,max = 1)
  for (q in 1:samplesize){
    sim[q,1] <- {{-log(1-unif[q])}^(1/k)}/r.initial.0
  }
  
  # Generate Y_i (min(T,C))
  sim[,3] <- apply(sim[,1:2], 1, FUN=min)
  # Censoring indicator (Censored=0, Not censored=1)
  sim[,5]<-I(sim[,1]<sim[,2])
  # Ordering
  sim <- sim[order(sim[,3]),]
  n <- nrow(sim)
  sim <- as.data.frame(sim)
  return(sim)
}

#### 3. Indicator function for measuring coverage proportion of 95% CI ####
ind <- function(a,b,c){
  if (a>=b&a<=c) {
    result <- 1
  } else {
    result <- 0
  }
  print(result)
}

#### 4. Make table for Beta estimation and variance estimation and Coverage ####
table0.isw1 <- matrix(NA,5,8)
rownames(table0.isw1) <- c(0,10,30,50,70)
colnames(table0.isw1) <- c("beta_0","SE of beta_0","SD of beta_0","Coverage of beta_0","beta_1","SE of beta_1","SD of beta_1","Coverage of beta_1")
table1.isw1<-matrix(NA,5,8)
rownames(table1.isw1) <- c(0,10,30,50,70)
colnames(table1.isw1) <- c("beta_0","SE of beta_0","SD of beta_0","Coverage of beta_0","beta_1","SE of beta_1","SD of beta_1","Coverage of beta_1")
table2.isw1<-matrix(NA,5,8)
rownames(table2.isw1) <- c(0,10,30,50,70)
colnames(table2.isw1) <- c("beta_0","SE of beta_0","SD of beta_0","Coverage of beta_0","beta_1","SE of beta_1","SD of beta_1","Coverage of beta_1")
table3.isw1<-matrix(NA,5,8)
rownames(table3.isw1) <- c(0,10,30,50,70)
colnames(table3.isw1) <- c("beta_0","SE of beta_0","SD of beta_0","Coverage of beta_0","beta_1","SE of beta_1","SD of beta_1","Coverage of beta_1")

#### 5. Execute simulation ####
#### censoring point at t_0=0 ####
c.0 <- 5000000
c.1 <- 52.15
c.3 <- 17.76
c.5 <- 10.5
c.7 <- 6.8

#### t_0=0 & c=0% ####
b0.isw1.00 <- c()
b0.isw1.sd.00 <- c()
b1.isw1.00 <- c()
b1.isw1.sd.00 <- c()
cover.isw1.00 <- matrix(NA,2000,8)
colnames(cover.isw1.00) <- c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

set.seed(1)
for (i in 1:2000){
  a <- data.gen(200,c.0)
  fm <- Surv(Z, censored) ~ X
  tryCatch({
    ismb.fit <- qris(fm, data = a, t0 = 0, Q = 0.5, ne = 200, "smooth", "pmb", c(1.609438, 0))
    b0.isw1.00[i] <- ismb.fit$coefficient[1]
    b0.isw1.sd.00[i] <- ismb.fit$stderr[1]
    b1.isw1.00[i] <- ismb.fit$coefficient[2]
    b1.isw1.sd.00[i] <- ismb.fit$stderr[2]
    # Coverage
    cover.isw1.00[i,1] <- ismb.fit$coefficient[1]-1.96*ismb.fit$stderr[1]
    cover.isw1.00[i,2] <- ismb.fit$coefficient[1]
    cover.isw1.00[i,3] <- ismb.fit$coefficient[1]+1.96*ismb.fit$stderr[1]
    cover.isw1.00[i,4] <- ind(1.609438, cover.isw1.00[i,1], cover.isw1.00[i,3])
    cover.isw1.00[i,5] <- ismb.fit$coefficient[2]-1.96*ismb.fit$stderr[2]
    cover.isw1.00[i,6] <- ismb.fit$coefficient[2]
    cover.isw1.00[i,7] <- ismb.fit$coefficient[2]+1.96*ismb.fit$stderr[2]
    cover.isw1.00[i,8] <- ind(0, cover.isw1.00[i,5], cover.isw1.00[i,7])}
    , error=function(e){
      b0.isw1.00[i] <- NA
      b0.isw1.sd.00[i] <- NA
      b1.isw1.00[i] <- NA
      b1.isw1.sd.00[i] <- NA
      # Coverage
      cover.isw1.00[i,1] <- NA
      cover.isw1.00[i,2] <- NA
      cover.isw1.00[i,3] <- NA
      cover.isw1.00[i,4] <- NA
      cover.isw1.00[i,5] <- NA
      cover.isw1.00[i,6] <- NA
      cover.isw1.00[i,7] <- NA
      cover.isw1.00[i,8] <- NA
    })
}
# IS beta table
table0.isw1[1,1] <- mean(b0.isw1.00,na.rm=TRUE)
table0.isw1[1,2] <- mean(b0.isw1.sd.00,na.rm=TRUE)
table0.isw1[1,3] <- sd(b0.isw1.00,na.rm=TRUE)
table0.isw1[1,4] <- mean(cover.isw1.00[,4],na.rm=TRUE)
table0.isw1[1,5] <- mean(b1.isw1.00,na.rm=TRUE)
table0.isw1[1,6] <- mean(b1.isw1.sd.00,na.rm=TRUE)
table0.isw1[1,7] <- sd(b1.isw1.00,na.rm=TRUE)
table0.isw1[1,8] <- mean(cover.isw1.00[,8],na.rm=TRUE)

#### t_0=0 & c=10% ####
b0.isw1.01 <- c()
b0.isw1.sd.01 <- c()
b1.isw1.01 <- c()
b1.isw1.sd.01 <- c()
cover.isw1.01 <- matrix(NA,2000,8)
colnames(cover.isw1.01) <- c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

set.seed(1)
for (i in 1:2000){
  a<-data.gen(200,c.1)
  fm <- Surv(Z, censored) ~ X
  tryCatch({
    ismb.fit <- qris(fm, data = a, t0 = 0, Q = 0.5, ne = 200, "smooth", "pmb", c(1.609438, 0))
    b0.isw1.01[i] <- ismb.fit$coefficient[1]
    b0.isw1.sd.01[i] <- ismb.fit$stderr[1]
    b1.isw1.01[i] <- ismb.fit$coefficient[2]
    b1.isw1.sd.01[i] <- ismb.fit$stderr[2]
    # Coverage
    cover.isw1.01[i,1] <- ismb.fit$coefficient[1]-1.96*ismb.fit$stderr[1]
    cover.isw1.01[i,2] <- ismb.fit$coefficient[1]
    cover.isw1.01[i,3] <- ismb.fit$coefficient[1]+1.96*ismb.fit$stderr[1]
    cover.isw1.01[i,4] <- ind(1.609438, cover.isw1.01[i,1], cover.isw1.01[i,3])
    cover.isw1.01[i,5] <- ismb.fit$coefficient[2]-1.96*ismb.fit$stderr[2]
    cover.isw1.01[i,6] <- ismb.fit$coefficient[2]
    cover.isw1.01[i,7] <- ismb.fit$coefficient[2]+1.96*ismb.fit$stderr[2]
    cover.isw1.01[i,8] <- ind(0, cover.isw1.01[i,5], cover.isw1.01[i,7])}
    , error=function(e){
      b0.isw1.01[i] <- NA
      b0.isw1.sd.01[i] <- NA
      b1.isw1.01[i] <- NA
      b1.isw1.sd.01[i] <- NA
      # Coverage
      cover.isw1.01[i,1] <- NA
      cover.isw1.01[i,2] <- NA
      cover.isw1.01[i,3] <- NA
      cover.isw1.01[i,4] <- NA
      cover.isw1.01[i,5] <- NA
      cover.isw1.01[i,6] <- NA
      cover.isw1.01[i,7] <- NA
      cover.isw1.01[i,8] <- NA
    })
}

# IS beta table
table0.isw1[2,1] <- mean(b0.isw1.01,na.rm=TRUE)
table0.isw1[2,2] <- mean(b0.isw1.sd.01,na.rm=TRUE)
table0.isw1[2,3] <- sd(b0.isw1.01,na.rm=TRUE)
table0.isw1[2,4] <- mean(cover.isw1.01[,4],na.rm=TRUE)
table0.isw1[2,5] <- mean(b1.isw1.01,na.rm=TRUE)
table0.isw1[2,6] <- mean(b1.isw1.sd.01,na.rm=TRUE)
table0.isw1[2,7] <- sd(b1.isw1.01,na.rm=TRUE)
table0.isw1[2,8] <- mean(cover.isw1.01[,8],na.rm=TRUE)

#### t_0=0 & c=30% ####
b0.isw1.03 <- c()
b0.isw1.sd.03 <- c()
b1.isw1.03 <- c()
b1.isw1.sd.03 <- c()
cover.isw1.03 <- matrix(NA,2000,8)
colnames(cover.isw1.03) <- c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

set.seed(1)
for (i in 1:2000){
  a <- data.gen(200,c.3)
  fm <- Surv(Z, censored) ~ X
  tryCatch({
    ismb.fit <- qris(fm, data = a, t0 = 0, Q = 0.5, ne = 200, "smooth", "pmb", c(1.609438, 0))
    b0.isw1.03[i] <- ismb.fit$coefficient[1]
    b0.isw1.sd.03[i] <- ismb.fit$stderr[1]
    b1.isw1.03[i] <- ismb.fit$coefficient[2]
    b1.isw1.sd.03[i] <- ismb.fit$stderr[2]
    # Coverage
    cover.isw1.03[i,1] <- ismb.fit$coefficient[1]-1.96*ismb.fit$stderr[1]
    cover.isw1.03[i,2] <- ismb.fit$coefficient[1]
    cover.isw1.03[i,3] <- ismb.fit$coefficient[1]+1.96*ismb.fit$stderr[1]
    cover.isw1.03[i,4] <- ind(1.609438, cover.isw1.03[i,1], cover.isw1.03[i,3])
    cover.isw1.03[i,5] <- ismb.fit$coefficient[2]-1.96*ismb.fit$stderr[2]
    cover.isw1.03[i,6] <- ismb.fit$coefficient[2]
    cover.isw1.03[i,7] <- ismb.fit$coefficient[2]+1.96*ismb.fit$stderr[2]
    cover.isw1.03[i,8] <- ind(0, cover.isw1.03[i,5], cover.isw1.03[i,7])}
    , error=function(e){
      b0.isw1.03[i] <- NA
      b0.isw1.sd.03[i] <- NA
      b1.isw1.03[i] <- NA
      b1.isw1.sd.03[i] <- NA
      # Coverage
      cover.isw1.03[i,1] <- NA
      cover.isw1.03[i,2] <- NA
      cover.isw1.03[i,3] <- NA
      cover.isw1.03[i,4] <- NA
      cover.isw1.03[i,5] <- NA
      cover.isw1.03[i,6] <- NA
      cover.isw1.03[i,7] <- NA
      cover.isw1.03[i,8] <- NA
    })
}

# IS beta table
table0.isw1[3,1] <- mean(b0.isw1.03,na.rm=TRUE)
table0.isw1[3,2] <- mean(b0.isw1.sd.03,na.rm=TRUE)
table0.isw1[3,3] <- sd(b0.isw1.03,na.rm=TRUE)
table0.isw1[3,4] <- mean(cover.isw1.03[,4],na.rm=TRUE)
table0.isw1[3,5] <- mean(b1.isw1.03,na.rm=TRUE)
table0.isw1[3,6] <- mean(b1.isw1.sd.03,na.rm=TRUE)
table0.isw1[3,7] <- sd(b1.isw1.03,na.rm=TRUE)
table0.isw1[3,8] <- mean(cover.isw1.03[,8],na.rm=TRUE)

#### t_0=0 & c=50% ####
b0.isw1.05 <- c()
b0.isw1.sd.05 <- c()
b1.isw1.05 <- c()
b1.isw1.sd.05 <- c()
cover.isw1.05 <- matrix(NA,2000,8)
colnames(cover.isw1.05) <- c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

set.seed(1)
for (i in 1:2000){
  a<-data.gen(200,c.5)
  fm <- Surv(Z, censored) ~ X
  tryCatch({
    ismb.fit <- qris(fm, data = a, t0 = 0, Q = 0.5, ne = 200, "smooth", "pmb", c(1.609438, 0))
    b0.isw1.05[i] <- ismb.fit$coefficient[1]
    b0.isw1.sd.05[i] <- ismb.fit$stderr[1]
    b1.isw1.05[i] <- ismb.fit$coefficient[2]
    b1.isw1.sd.05[i] <- ismb.fit$stderr[2]
    # Coverage
    cover.isw1.05[i,1] <- ismb.fit$coefficient[1]-1.96*ismb.fit$stderr[1]
    cover.isw1.05[i,2] <- ismb.fit$coefficient[1]
    cover.isw1.05[i,3] <- ismb.fit$coefficient[1]+1.96*ismb.fit$stderr[1]
    cover.isw1.05[i,4] <- ind(1.609438, cover.isw1.05[i,1], cover.isw1.05[i,3])
    cover.isw1.05[i,5] <- ismb.fit$coefficient[2]-1.96*ismb.fit$stderr[2]
    cover.isw1.05[i,6] <- ismb.fit$coefficient[2]
    cover.isw1.05[i,7] <- ismb.fit$coefficient[2]+1.96*ismb.fit$stderr[2]
    cover.isw1.05[i,8] <- ind(0, cover.isw1.05[i,5], cover.isw1.05[i,7])}
    , error=function(e){
      b0.isw1.05[i] <- NA
      b0.isw1.sd.05[i] <- NA
      b1.isw1.05[i] <- NA
      b1.isw1.sd.05[i] <- NA
      # Coverage
      cover.isw1.05[i,1] <- NA
      cover.isw1.05[i,2] <- NA
      cover.isw1.05[i,3] <- NA
      cover.isw1.05[i,4] <- NA
      cover.isw1.05[i,5] <- NA
      cover.isw1.05[i,6] <- NA
      cover.isw1.05[i,7] <- NA
      cover.isw1.05[i,8] <- NA
    })
}

# IS beta table
table0.isw1[4,1]<-mean(b0.isw1.05,na.rm=TRUE)
table0.isw1[4,2]<-mean(b0.isw1.sd.05,na.rm=TRUE)
table0.isw1[4,3]<-sd(b0.isw1.05,na.rm=TRUE)
table0.isw1[4,4]<-mean(cover.isw1.05[,4],na.rm=TRUE)
table0.isw1[4,5]<-mean(b1.isw1.05,na.rm=TRUE)
table0.isw1[4,6]<-mean(b1.isw1.sd.05,na.rm=TRUE)
table0.isw1[4,7]<-sd(b1.isw1.05,na.rm=TRUE)
table0.isw1[4,8]<-mean(cover.isw1.05[,8],na.rm=TRUE)

#### t_0=0 & c=70% ####
b0.isw1.07 <- c()
b0.isw1.sd.07 <- c()
b1.isw1.07 <- c()
b1.isw1.sd.07 <- c()
cover.isw1.07<-matrix(NA,2000,8)
colnames(cover.isw1.07)<-c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

set.seed(1)
for (i in 1:2000){
  a<-data.gen(200,c.7)
  fm <- Surv(Z, censored) ~ X
  tryCatch({
    ismb.fit <- qris(fm, data = a, t0 = 0, Q = 0.5, ne = 200, "smooth", "pmb", c(1.609438, 0))
    b0.isw1.07[i] <- ismb.fit$coefficient[1]
    b0.isw1.sd.07[i] <- ismb.fit$stderr[1]
    b1.isw1.07[i] <- ismb.fit$coefficient[2]
    b1.isw1.sd.07[i] <- ismb.fit$stderr[2]
    # Coverage
    cover.isw1.07[i,1] <- ismb.fit$coefficient[1]-1.96*ismb.fit$stderr[1]
    cover.isw1.07[i,2] <- ismb.fit$coefficient[1]
    cover.isw1.07[i,3] <- ismb.fit$coefficient[1]+1.96*ismb.fit$stderr[1]
    cover.isw1.07[i,4] <- ind(1.609438, cover.isw1.07[i,1], cover.isw1.07[i,3])
    cover.isw1.07[i,5] <- ismb.fit$coefficient[2]-1.96*ismb.fit$stderr[2]
    cover.isw1.07[i,6] <- ismb.fit$coefficient[2]
    cover.isw1.07[i,7] <- ismb.fit$coefficient[2]+1.96*ismb.fit$stderr[2]
    cover.isw1.07[i,8] <- ind(0, cover.isw1.07[i,5], cover.isw1.07[i,7])}
    , error=function(e){
      b0.isw1.07[i] <- NA
      b0.isw1.sd.07[i] <- NA
      b1.isw1.07[i] <- NA
      b1.isw1.sd.07[i] <- NA
      # Coverage
      cover.isw1.07[i,1] <- NA
      cover.isw1.07[i,2] <- NA
      cover.isw1.07[i,3] <- NA
      cover.isw1.07[i,4] <- NA
      cover.isw1.07[i,5] <- NA
      cover.isw1.07[i,6] <- NA
      cover.isw1.07[i,7] <- NA
      cover.isw1.07[i,8] <- NA
    })
}

# IS beta table
table0.isw1[5,1]<-mean(b0.isw1.07,na.rm=TRUE)
table0.isw1[5,2]<-mean(b0.isw1.sd.07,na.rm=TRUE)
table0.isw1[5,3]<-sd(b0.isw1.07,na.rm=TRUE)
table0.isw1[5,4]<-mean(cover.isw1.07[,4],na.rm=TRUE)
table0.isw1[5,5]<-mean(b1.isw1.07,na.rm=TRUE)
table0.isw1[5,6]<-mean(b1.isw1.sd.07,na.rm=TRUE)
table0.isw1[5,7]<-sd(b1.isw1.07,na.rm=TRUE)
table0.isw1[5,8]<-mean(cover.isw1.07[,8],na.rm=TRUE)

#### censoring point at t_0=1 ####
c.0<-5000000
c.1<-44.78
c.3<-15.85
c.5<-9.63
c.7<-6.26
#### t_0=1 & c=0% ####
b0.isw1.10 <- c()
b0.isw1.sd.10 <- c()
b1.isw1.10 <- c()
b1.isw1.sd.10 <- c()
cover.isw1.10<-matrix(NA,2000,8)
colnames(cover.isw1.10)<-c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

set.seed(1)
for (i in 1:2000){
  a <- data.gen(200,c.0)
  fm <- Surv(Z, censored) ~ X
  tryCatch({
    ismb.fit <- qris(fm, data = a, t0 = 1, Q = 0.5, ne = 200, "smooth", "pmb", c(1.410748, 0))
    b0.isw1.10[i] <- ismb.fit$coefficient[1]
    b0.isw1.sd.10[i] <- ismb.fit$stderr[1]
    b1.isw1.10[i] <- ismb.fit$coefficient[2]
    b1.isw1.sd.10[i] <- ismb.fit$stderr[2]
    # Coverage
    cover.isw1.10[i,1] <- ismb.fit$coefficient[1]-1.96*ismb.fit$stderr[1]
    cover.isw1.10[i,2] <- ismb.fit$coefficient[1]
    cover.isw1.10[i,3] <- ismb.fit$coefficient[1]+1.96*ismb.fit$stderr[1]
    cover.isw1.10[i,4] <- ind(1.410748, cover.isw1.10[i,1], cover.isw1.10[i,3])
    cover.isw1.10[i,5] <- ismb.fit$coefficient[2]-1.96*ismb.fit$stderr[2]
    cover.isw1.10[i,6] <- ismb.fit$coefficient[2]
    cover.isw1.10[i,7] <- ismb.fit$coefficient[2]+1.96*ismb.fit$stderr[2]
    cover.isw1.10[i,8] <- ind(0, cover.isw1.10[i,5], cover.isw1.10[i,7])}
    , error=function(e){
      b0.isw1.10[i] <- NA
      b0.isw1.sd.10[i] <- NA
      b1.isw1.10[i] <- NA
      b1.isw1.sd.10[i] <- NA
      # Coverage
      cover.isw1.10[i,1] <- NA
      cover.isw1.10[i,2] <- NA
      cover.isw1.10[i,3] <- NA
      cover.isw1.10[i,4] <- NA
      cover.isw1.10[i,5] <- NA
      cover.isw1.10[i,6] <- NA
      cover.isw1.10[i,7] <- NA
      cover.isw1.10[i,8] <- NA
    })
}

# IS beta table
table1.isw1[1,1]<-mean(b0.isw1.10,na.rm=TRUE)
table1.isw1[1,2]<-mean(b0.isw1.sd.10,na.rm=TRUE)
table1.isw1[1,3]<-sd(b0.isw1.10,na.rm=TRUE)
table1.isw1[1,4]<-mean(cover.isw1.10[,4],na.rm=TRUE)
table1.isw1[1,5]<-mean(b1.isw1.10,na.rm=TRUE)
table1.isw1[1,6]<-mean(b1.isw1.sd.10,na.rm=TRUE)
table1.isw1[1,7]<-sd(b1.isw1.10,na.rm=TRUE)
table1.isw1[1,8]<-mean(cover.isw1.10[,8],na.rm=TRUE)

#### t_0=1 & c=10% ####
b0.isw1.11 <- c()
b0.isw1.sd.11 <- c()
b1.isw1.11 <- c()
b1.isw1.sd.11 <- c()
cover.isw1.11<-matrix(NA,2000,8)
colnames(cover.isw1.11)<-c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

set.seed(1)
for (i in 1:2000){
  a<-data.gen(200,c.1)
  fm <- Surv(Z, censored) ~ X
  tryCatch({
    ismb.fit <- qris(fm, data = a, t0 = 1, Q = 0.5, ne = 200, "smooth", "pmb", c(1.410748, 0))
    b0.isw1.11[i] <- ismb.fit$coefficient[1]
    b0.isw1.sd.11[i] <- ismb.fit$stderr[1]
    b1.isw1.11[i] <- ismb.fit$coefficient[2]
    b1.isw1.sd.11[i] <- ismb.fit$stderr[2]
    # Coverage
    cover.isw1.11[i,1] <- ismb.fit$coefficient[1]-1.96*ismb.fit$stderr[1]
    cover.isw1.11[i,2] <- ismb.fit$coefficient[1]
    cover.isw1.11[i,3] <- ismb.fit$coefficient[1]+1.96*ismb.fit$stderr[1]
    cover.isw1.11[i,4] <- ind(1.410748, cover.isw1.11[i,1], cover.isw1.11[i,3])
    cover.isw1.11[i,5] <- ismb.fit$coefficient[2]-1.96*ismb.fit$stderr[2]
    cover.isw1.11[i,6] <- ismb.fit$coefficient[2]
    cover.isw1.11[i,7] <- ismb.fit$coefficient[2]+1.96*ismb.fit$stderr[2]
    cover.isw1.11[i,8] <- ind(0, cover.isw1.11[i,5], cover.isw1.11[i,7])}
    , error=function(e){
      b0.isw1.11[i] <- NA
      b0.isw1.sd.11[i] <- NA
      b1.isw1.11[i] <- NA
      b1.isw1.sd.11[i] <- NA
      # Coverage
      cover.isw1.11[i,1] <- NA
      cover.isw1.11[i,2] <- NA
      cover.isw1.11[i,3] <- NA
      cover.isw1.11[i,4] <- NA
      cover.isw1.11[i,5] <- NA
      cover.isw1.11[i,6] <- NA
      cover.isw1.11[i,7] <- NA
      cover.isw1.11[i,8] <- NA
    })
}

# IS beta table
table1.isw1[2,1]<-mean(b0.isw1.11,na.rm=TRUE)
table1.isw1[2,2]<-mean(b0.isw1.sd.11,na.rm=TRUE)
table1.isw1[2,3]<-sd(b0.isw1.11,na.rm=TRUE)
table1.isw1[2,4]<-mean(cover.isw1.11[,4],na.rm=TRUE)
table1.isw1[2,5]<-mean(b1.isw1.11,na.rm=TRUE)
table1.isw1[2,6]<-mean(b1.isw1.sd.11,na.rm=TRUE)
table1.isw1[2,7]<-sd(b1.isw1.11,na.rm=TRUE)
table1.isw1[2,8]<-mean(cover.isw1.11[,8],na.rm=TRUE)

#### t_0=1 & c=30% ####
b0.isw1.13 <- c()
b0.isw1.sd.13 <- c()
b1.isw1.13 <- c()
b1.isw1.sd.13 <- c()
cover.isw1.13<-matrix(NA,2000,8)
colnames(cover.isw1.13)<-c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

set.seed(1)
for (i in 1:2000){
  a<-data.gen(200,c.3)
  fm <- Surv(Z, censored) ~ X
  tryCatch({
    ismb.fit <- qris(fm, data = a, t0 = 1, Q = 0.5, ne = 200, "smooth", "pmb", c(1.410748, 0))
    b0.isw1.13[i] <- ismb.fit$coefficient[1]
    b0.isw1.sd.13[i] <- ismb.fit$stderr[1]
    b1.isw1.13[i] <- ismb.fit$coefficient[2]
    b1.isw1.sd.13[i] <- ismb.fit$stderr[2]
    # Coverage
    cover.isw1.13[i,1] <- ismb.fit$coefficient[1]-1.96*ismb.fit$stderr[1]
    cover.isw1.13[i,2] <- ismb.fit$coefficient[1]
    cover.isw1.13[i,3] <- ismb.fit$coefficient[1]+1.96*ismb.fit$stderr[1]
    cover.isw1.13[i,4] <- ind(1.410748, cover.isw1.13[i,1], cover.isw1.13[i,3])
    cover.isw1.13[i,5] <- ismb.fit$coefficient[2]-1.96*ismb.fit$stderr[2]
    cover.isw1.13[i,6] <- ismb.fit$coefficient[2]
    cover.isw1.13[i,7] <- ismb.fit$coefficient[2]+1.96*ismb.fit$stderr[2]
    cover.isw1.13[i,8] <- ind(0, cover.isw1.13[i,5], cover.isw1.13[i,7])}
    , error=function(e){
      b0.isw1.13[i] <- NA
      b0.isw1.sd.13[i] <- NA
      b1.isw1.13[i] <- NA
      b1.isw1.sd.13[i] <- NA
      # Coverage
      cover.isw1.13[i,1] <- NA
      cover.isw1.13[i,2] <- NA
      cover.isw1.13[i,3] <- NA
      cover.isw1.13[i,4] <- NA
      cover.isw1.13[i,5] <- NA
      cover.isw1.13[i,6] <- NA
      cover.isw1.13[i,7] <- NA
      cover.isw1.13[i,8] <- NA
    })
}

# IS beta table
table1.isw1[3,1]<-mean(b0.isw1.13,na.rm=TRUE)
table1.isw1[3,2]<-mean(b0.isw1.sd.13,na.rm=TRUE)
table1.isw1[3,3]<-sd(b0.isw1.13,na.rm=TRUE)
table1.isw1[3,4]<-mean(cover.isw1.13[,4],na.rm=TRUE)
table1.isw1[3,5]<-mean(b1.isw1.13,na.rm=TRUE)
table1.isw1[3,6]<-mean(b1.isw1.sd.13,na.rm=TRUE)
table1.isw1[3,7]<-sd(b1.isw1.13,na.rm=TRUE)
table1.isw1[3,8]<-mean(cover.isw1.13[,8],na.rm=TRUE)

#### t_0=1 & c=50% ####
b0.isw1.15 <- c()
b0.isw1.sd.15 <- c()
b1.isw1.15 <- c()
b1.isw1.sd.15 <- c()
cover.isw1.15<-matrix(NA,2000,8)
colnames(cover.isw1.15)<-c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

set.seed(1)
for (i in 1:2000){
  a<-data.gen(200,c.5)
  fm <- Surv(Z, censored) ~ X
  tryCatch({
    ismb.fit <- qris(fm, data = a, t0 = 1, Q = 0.5, ne = 200, "smooth", "pmb", c(1.410748, 0))
    b0.isw1.15[i] <- ismb.fit$coefficient[1]
    b0.isw1.sd.15[i] <- ismb.fit$stderr[1]
    b1.isw1.15[i] <- ismb.fit$coefficient[2]
    b1.isw1.sd.15[i] <- ismb.fit$stderr[2]
    # Coverage
    cover.isw1.15[i,1] <- ismb.fit$coefficient[1]-1.96*ismb.fit$stderr[1]
    cover.isw1.15[i,2] <- ismb.fit$coefficient[1]
    cover.isw1.15[i,3] <- ismb.fit$coefficient[1]+1.96*ismb.fit$stderr[1]
    cover.isw1.15[i,4] <- ind(1.410748, cover.isw1.15[i,1], cover.isw1.15[i,3])
    cover.isw1.15[i,5] <- ismb.fit$coefficient[2]-1.96*ismb.fit$stderr[2]
    cover.isw1.15[i,6] <- ismb.fit$coefficient[2]
    cover.isw1.15[i,7] <- ismb.fit$coefficient[2]+1.96*ismb.fit$stderr[2]
    cover.isw1.15[i,8] <- ind(0, cover.isw1.15[i,5], cover.isw1.15[i,7])}
    , error=function(e){
      b0.isw1.15[i] <- NA
      b0.isw1.sd.15[i] <- NA
      b1.isw1.15[i] <- NA
      b1.isw1.sd.15[i] <- NA
      # Coverage
      cover.isw1.15[i,1] <- NA
      cover.isw1.15[i,2] <- NA
      cover.isw1.15[i,3] <- NA
      cover.isw1.15[i,4] <- NA
      cover.isw1.15[i,5] <- NA
      cover.isw1.15[i,6] <- NA
      cover.isw1.15[i,7] <- NA
      cover.isw1.15[i,8] <- NA
    })
}

# IS beta table
table1.isw1[4,1]<-mean(b0.isw1.15,na.rm=TRUE)
table1.isw1[4,2]<-mean(b0.isw1.sd.15,na.rm=TRUE)
table1.isw1[4,3]<-sd(b0.isw1.15,na.rm=TRUE)
table1.isw1[4,4]<-mean(cover.isw1.15[,4],na.rm=TRUE)
table1.isw1[4,5]<-mean(b1.isw1.15,na.rm=TRUE)
table1.isw1[4,6]<-mean(b1.isw1.sd.15,na.rm=TRUE)
table1.isw1[4,7]<-sd(b1.isw1.15,na.rm=TRUE)
table1.isw1[4,8]<-mean(cover.isw1.15[,8],na.rm=TRUE)

#### t_0=1 & c=70% ####
b0.isw1.17 <- c()
b0.isw1.sd.17 <- c()
b1.isw1.17 <- c()
b1.isw1.sd.17 <- c()
cover.isw1.17<-matrix(NA,2000,8)
colnames(cover.isw1.17)<-c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

set.seed(1)
for (i in 1:2000){
  a<-data.gen(200,c.7)
  fm <- Surv(Z, censored) ~ X
  tryCatch({
    ismb.fit <- qris(fm, data = a, t0 = 1, Q = 0.5, ne = 200, "smooth", "pmb", c(1.410748, 0))
    b0.isw1.17[i] <- ismb.fit$coefficient[1]
    b0.isw1.sd.17[i] <- ismb.fit$stderr[1]
    b1.isw1.17[i] <- ismb.fit$coefficient[2]
    b1.isw1.sd.17[i] <- ismb.fit$stderr[2]
    # Coverage
    cover.isw1.17[i,1] <- ismb.fit$coefficient[1]-1.96*ismb.fit$stderr[1]
    cover.isw1.17[i,2] <- ismb.fit$coefficient[1]
    cover.isw1.17[i,3] <- ismb.fit$coefficient[1]+1.96*ismb.fit$stderr[1]
    cover.isw1.17[i,4] <- ind(1.410748, cover.isw1.17[i,1], cover.isw1.17[i,3])
    cover.isw1.17[i,5] <- ismb.fit$coefficient[2]-1.96*ismb.fit$stderr[2]
    cover.isw1.17[i,6] <- ismb.fit$coefficient[2]
    cover.isw1.17[i,7] <- ismb.fit$coefficient[2]+1.96*ismb.fit$stderr[2]
    cover.isw1.17[i,8] <- ind(0, cover.isw1.17[i,5], cover.isw1.17[i,7])}
    , error=function(e){
      b0.isw1.17[i] <- NA
      b0.isw1.sd.17[i] <- NA
      b1.isw1.17[i] <- NA
      b1.isw1.sd.17[i] <- NA
      # Coverage
      cover.isw1.17[i,1] <- NA
      cover.isw1.17[i,2] <- NA
      cover.isw1.17[i,3] <- NA
      cover.isw1.17[i,4] <- NA
      cover.isw1.17[i,5] <- NA
      cover.isw1.17[i,6] <- NA
      cover.isw1.17[i,7] <- NA
      cover.isw1.17[i,8] <- NA
    })
  
}

# IS beta table
table1.isw1[5,1]<-mean(b0.isw1.17,na.rm=TRUE)
table1.isw1[5,2]<-mean(b0.isw1.sd.17,na.rm=TRUE)
table1.isw1[5,3]<-sd(b0.isw1.17,na.rm=TRUE)
table1.isw1[5,4]<-mean(cover.isw1.17[,4],na.rm=TRUE)
table1.isw1[5,5]<-mean(b1.isw1.17,na.rm=TRUE)
table1.isw1[5,6]<-mean(b1.isw1.sd.17,na.rm=TRUE)
table1.isw1[5,7]<-sd(b1.isw1.17,na.rm=TRUE)
table1.isw1[5,8]<-mean(cover.isw1.17[,8],na.rm=TRUE)

#### censoring point at t_0=2 ####
c.0<-5000000
c.1<-39.03
c.3<-14.55
c.5<-9.22
c.7<-6.2
#### t_0=2 & c=0% ####
b0.isw1.20 <- c()
b0.isw1.sd.20 <- c()
b1.isw1.20 <- c()
b1.isw1.sd.20 <- c()
cover.isw1.20<-matrix(NA,2000,8)
colnames(cover.isw1.20)<-c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

set.seed(1)
for (i in 1:2000){
  a <- data.gen(200,c.0)
  fm <- Surv(Z, censored) ~ X
  tryCatch({
    ismb.fit <- qris(fm, data = a, t0 = 2, Q = 0.5, ne = 200, "smooth", "pmb", c(1.219403, 0))
    b0.isw1.20[i] <- ismb.fit$coefficient[1]
    b0.isw1.sd.20[i] <- ismb.fit$stderr[1]
    b1.isw1.20[i] <- ismb.fit$coefficient[2]
    b1.isw1.sd.20[i] <- ismb.fit$stderr[2]
    # Coverage
    cover.isw1.20[i,1] <- ismb.fit$coefficient[1]-1.96*ismb.fit$stderr[1]
    cover.isw1.20[i,2] <- ismb.fit$coefficient[1]
    cover.isw1.20[i,3] <- ismb.fit$coefficient[1]+1.96*ismb.fit$stderr[1]
    cover.isw1.20[i,4] <- ind(1.219403, cover.isw1.20[i,1], cover.isw1.20[i,3])
    cover.isw1.20[i,5] <- ismb.fit$coefficient[2]-1.96*ismb.fit$stderr[2]
    cover.isw1.20[i,6] <- ismb.fit$coefficient[2]
    cover.isw1.20[i,7] <- ismb.fit$coefficient[2]+1.96*ismb.fit$stderr[2]
    cover.isw1.20[i,8] <- ind(0, cover.isw1.20[i,5], cover.isw1.20[i,7])}
    , error=function(e){
      b0.isw1.20[i] <- NA
      b0.isw1.sd.20[i] <- NA
      b1.isw1.20[i] <- NA
      b1.isw1.sd.20[i] <- NA
      # Coverage
      cover.isw1.20[i,1] <- NA
      cover.isw1.20[i,2] <- NA
      cover.isw1.20[i,3] <- NA
      cover.isw1.20[i,4] <- NA
      cover.isw1.20[i,5] <- NA
      cover.isw1.20[i,6] <- NA
      cover.isw1.20[i,7] <- NA
      cover.isw1.20[i,8] <- NA
    })
}

# IS beta table
table2.isw1[1,1]<-mean(b0.isw1.20,na.rm=TRUE)
table2.isw1[1,2]<-mean(b0.isw1.sd.20,na.rm=TRUE)
table2.isw1[1,3]<-sd(b0.isw1.20,na.rm=TRUE)
table2.isw1[1,4]<-mean(cover.isw1.20[,4],na.rm=TRUE)
table2.isw1[1,5]<-mean(b1.isw1.20,na.rm=TRUE)
table2.isw1[1,6]<-mean(b1.isw1.sd.20,na.rm=TRUE)
table2.isw1[1,7]<-sd(b1.isw1.20,na.rm=TRUE)
table2.isw1[1,8]<-mean(cover.isw1.20[,8],na.rm=TRUE)

#### t_0=2 & c=10% ####
b0.isw1.21 <- c()
b0.isw1.sd.21 <- c()
b1.isw1.21 <- c()
b1.isw1.sd.21 <- c()
cover.isw1.21<-matrix(NA,2000,8)
colnames(cover.isw1.21)<-c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

set.seed(1)
for (i in 1:2000){
  a<-data.gen(200,c.1)
  fm <- Surv(Z, censored) ~ X
  tryCatch({
    ismb.fit <- qris(fm, data = a, t0 = 2, Q = 0.5, ne = 200, "smooth", "pmb", c(1.219403, 0))
    b0.isw1.21[i] <- ismb.fit$coefficient[1]
    b0.isw1.sd.21[i] <- ismb.fit$stderr[1]
    b1.isw1.21[i] <- ismb.fit$coefficient[2]
    b1.isw1.sd.21[i] <- ismb.fit$stderr[2]
    # Coverage
    cover.isw1.21[i,1] <- ismb.fit$coefficient[1]-1.96*ismb.fit$stderr[1]
    cover.isw1.21[i,2] <- ismb.fit$coefficient[1]
    cover.isw1.21[i,3] <- ismb.fit$coefficient[1]+1.96*ismb.fit$stderr[1]
    cover.isw1.21[i,4] <- ind(1.219403, cover.isw1.21[i,1], cover.isw1.21[i,3])
    cover.isw1.21[i,5] <- ismb.fit$coefficient[2]-1.96*ismb.fit$stderr[2]
    cover.isw1.21[i,6] <- ismb.fit$coefficient[2]
    cover.isw1.21[i,7] <- ismb.fit$coefficient[2]+1.96*ismb.fit$stderr[2]
    cover.isw1.21[i,8] <- ind(0, cover.isw1.21[i,5], cover.isw1.21[i,7])}
    , error=function(e){
      b0.isw1.21[i] <- NA
      b0.isw1.sd.21[i] <- NA
      b1.isw1.21[i] <- NA
      b1.isw1.sd.21[i] <- NA
      # Coverage
      cover.isw1.21[i,1] <- NA
      cover.isw1.21[i,2] <- NA
      cover.isw1.21[i,3] <- NA
      cover.isw1.21[i,4] <- NA
      cover.isw1.21[i,5] <- NA
      cover.isw1.21[i,6] <- NA
      cover.isw1.21[i,7] <- NA
      cover.isw1.21[i,8] <- NA
    })
}

# IS beta table
table2.isw1[2,1]<-mean(b0.isw1.21,na.rm=TRUE)
table2.isw1[2,2]<-mean(b0.isw1.sd.21,na.rm=TRUE)
table2.isw1[2,3]<-sd(b0.isw1.21,na.rm=TRUE)
table2.isw1[2,4]<-mean(cover.isw1.21[,4],na.rm=TRUE)
table2.isw1[2,5]<-mean(b1.isw1.21,na.rm=TRUE)
table2.isw1[2,6]<-mean(b1.isw1.sd.21,na.rm=TRUE)
table2.isw1[2,7]<-sd(b1.isw1.21,na.rm=TRUE)
table2.isw1[2,8]<-mean(cover.isw1.21[,8],na.rm=TRUE)

#### t_0=2 & c=30% ####
b0.isw1.23 <- c()
b0.isw1.sd.23 <- c()
b1.isw1.23 <- c()
b1.isw1.sd.23 <- c()
cover.isw1.23<-matrix(NA,2000,8)
colnames(cover.isw1.23)<-c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

set.seed(1)
for (i in 1:2000){
  a<-data.gen(200,c.3)
  fm <- Surv(Z, censored) ~ X
  tryCatch({
    ismb.fit <- qris(fm, data = a, t0 = 2, Q = 0.5, ne = 200, "smooth", "pmb", c(1.219403, 0))
    b0.isw1.23[i] <- ismb.fit$coefficient[1]
    b0.isw1.sd.23[i] <- ismb.fit$stderr[1]
    b1.isw1.23[i] <- ismb.fit$coefficient[2]
    b1.isw1.sd.23[i] <- ismb.fit$stderr[2]
    # Coverage
    cover.isw1.23[i,1] <- ismb.fit$coefficient[1]-1.96*ismb.fit$stderr[1]
    cover.isw1.23[i,2] <- ismb.fit$coefficient[1]
    cover.isw1.23[i,3] <- ismb.fit$coefficient[1]+1.96*ismb.fit$stderr[1]
    cover.isw1.23[i,4] <- ind(1.219403, cover.isw1.23[i,1], cover.isw1.23[i,3])
    cover.isw1.23[i,5] <- ismb.fit$coefficient[2]-1.96*ismb.fit$stderr[2]
    cover.isw1.23[i,6] <- ismb.fit$coefficient[2]
    cover.isw1.23[i,7] <- ismb.fit$coefficient[2]+1.96*ismb.fit$stderr[2]
    cover.isw1.23[i,8] <- ind(0, cover.isw1.23[i,5], cover.isw1.23[i,7])}
    , error=function(e){
      b0.isw1.23[i] <- NA
      b0.isw1.sd.23[i] <- NA
      b1.isw1.23[i] <- NA
      b1.isw1.sd.23[i] <- NA
      # Coverage
      cover.isw1.23[i,1] <- NA
      cover.isw1.23[i,2] <- NA
      cover.isw1.23[i,3] <- NA
      cover.isw1.23[i,4] <- NA
      cover.isw1.23[i,5] <- NA
      cover.isw1.23[i,6] <- NA
      cover.isw1.23[i,7] <- NA
      cover.isw1.23[i,8] <- NA
    })
}

# IS beta table
table2.isw1[3,1]<-mean(b0.isw1.23,na.rm=TRUE)
table2.isw1[3,2]<-mean(b0.isw1.sd.23,na.rm=TRUE)
table2.isw1[3,3]<-sd(b0.isw1.23,na.rm=TRUE)
table2.isw1[3,4]<-mean(cover.isw1.23[,4],na.rm=TRUE)
table2.isw1[3,5]<-mean(b1.isw1.23,na.rm=TRUE)
table2.isw1[3,6]<-mean(b1.isw1.sd.23,na.rm=TRUE)
table2.isw1[3,7]<-sd(b1.isw1.23,na.rm=TRUE)
table2.isw1[3,8]<-mean(cover.isw1.23[,8],na.rm=TRUE)

#### t_0=2 & c=50% ####
b0.isw1.25 <- c()
b0.isw1.sd.25 <- c()
b1.isw1.25 <- c()
b1.isw1.sd.25 <- c()
cover.isw1.25<-matrix(NA,2000,8)
colnames(cover.isw1.25)<-c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

set.seed(1)
for (i in 1:2000){
  a<-data.gen(200,c.5)
  fm <- Surv(Z, censored) ~ X
  tryCatch({
    ismb.fit <- qris(fm, data = a, t0 = 2, Q = 0.5, ne = 200, "smooth", "pmb", c(1.219403, 0))
    b0.isw1.25[i] <- ismb.fit$coefficient[1]
    b0.isw1.sd.25[i] <- ismb.fit$stderr[1]
    b1.isw1.25[i] <- ismb.fit$coefficient[2]
    b1.isw1.sd.25[i] <- ismb.fit$stderr[2]
    # Coverage
    cover.isw1.25[i,1] <- ismb.fit$coefficient[1]-1.96*ismb.fit$stderr[1]
    cover.isw1.25[i,2] <- ismb.fit$coefficient[1]
    cover.isw1.25[i,3] <- ismb.fit$coefficient[1]+1.96*ismb.fit$stderr[1]
    cover.isw1.25[i,4] <- ind(1.219403, cover.isw1.25[i,1], cover.isw1.25[i,3])
    cover.isw1.25[i,5] <- ismb.fit$coefficient[2]-1.96*ismb.fit$stderr[2]
    cover.isw1.25[i,6] <- ismb.fit$coefficient[2]
    cover.isw1.25[i,7] <- ismb.fit$coefficient[2]+1.96*ismb.fit$stderr[2]
    cover.isw1.25[i,8] <- ind(0, cover.isw1.25[i,5], cover.isw1.25[i,7])}
    , error=function(e){
      b0.isw1.25[i] <- NA
      b0.isw1.sd.25[i] <- NA
      b1.isw1.25[i] <- NA
      b1.isw1.sd.25[i] <- NA
      # Coverage
      cover.isw1.25[i,1] <- NA
      cover.isw1.25[i,2] <- NA
      cover.isw1.25[i,3] <- NA
      cover.isw1.25[i,4] <- NA
      cover.isw1.25[i,5] <- NA
      cover.isw1.25[i,6] <- NA
      cover.isw1.25[i,7] <- NA
      cover.isw1.25[i,8] <- NA
    })
}

# IS beta table
table2.isw1[4,1]<-mean(b0.isw1.25,na.rm=TRUE)
table2.isw1[4,2]<-mean(b0.isw1.sd.25,na.rm=TRUE)
table2.isw1[4,3]<-sd(b0.isw1.25,na.rm=TRUE)
table2.isw1[4,4]<-mean(cover.isw1.25[,4],na.rm=TRUE)
table2.isw1[4,5]<-mean(b1.isw1.25,na.rm=TRUE)
table2.isw1[4,6]<-mean(b1.isw1.sd.25,na.rm=TRUE)
table2.isw1[4,7]<-sd(b1.isw1.25,na.rm=TRUE)
table2.isw1[4,8]<-mean(cover.isw1.25[,8],na.rm=TRUE)

#### t_0=2 & c=70% ####
b0.isw1.27 <- c()
b0.isw1.sd.27 <- c()
b1.isw1.27 <- c()
b1.isw1.sd.27 <- c()
cover.isw1.27<-matrix(NA,2000,8)
colnames(cover.isw1.27)<-c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

set.seed(1)
for (i in 1:2000){
  a<-data.gen(200,c.7)
  fm <- Surv(Z, censored) ~ X
  tryCatch({
    ismb.fit <- qris(fm, data = a, t0 = 2, Q = 0.5, ne = 200, "smooth", "pmb", c(1.219403, 0))
    b0.isw1.27[i] <- ismb.fit$coefficient[1]
    b0.isw1.sd.27[i] <- ismb.fit$stderr[1]
    b1.isw1.27[i] <- ismb.fit$coefficient[2]
    b1.isw1.sd.27[i] <- ismb.fit$stderr[2]
    # Coverage
    cover.isw1.27[i,1] <- ismb.fit$coefficient[1]-1.96*ismb.fit$stderr[1]
    cover.isw1.27[i,2] <- ismb.fit$coefficient[1]
    cover.isw1.27[i,3] <- ismb.fit$coefficient[1]+1.96*ismb.fit$stderr[1]
    cover.isw1.27[i,4] <- ind(1.219403, cover.isw1.27[i,1], cover.isw1.27[i,3])
    cover.isw1.27[i,5] <- ismb.fit$coefficient[2]-1.96*ismb.fit$stderr[2]
    cover.isw1.27[i,6] <- ismb.fit$coefficient[2]
    cover.isw1.27[i,7] <- ismb.fit$coefficient[2]+1.96*ismb.fit$stderr[2]
    cover.isw1.27[i,8] <- ind(0, cover.isw1.27[i,5], cover.isw1.27[i,7])}
    , error=function(e){
      b0.isw1.27[i] <- NA
      b0.isw1.sd.27[i] <- NA
      b1.isw1.27[i] <- NA
      b1.isw1.sd.27[i] <- NA
      # Coverage
      cover.isw1.27[i,1] <- NA
      cover.isw1.27[i,2] <- NA
      cover.isw1.27[i,3] <- NA
      cover.isw1.27[i,4] <- NA
      cover.isw1.27[i,5] <- NA
      cover.isw1.27[i,6] <- NA
      cover.isw1.27[i,7] <- NA
      cover.isw1.27[i,8] <- NA
    })
}

# IS beta table
table2.isw1[5,1]<-mean(b0.isw1.27,na.rm=TRUE)
table2.isw1[5,2]<-mean(b0.isw1.sd.27,na.rm=TRUE)
table2.isw1[5,3]<-sd(b0.isw1.27,na.rm=TRUE)
table2.isw1[5,4]<-mean(cover.isw1.27[,4],na.rm=TRUE)
table2.isw1[5,5]<-mean(b1.isw1.27,na.rm=TRUE)
table2.isw1[5,6]<-mean(b1.isw1.sd.27,na.rm=TRUE)
table2.isw1[5,7]<-sd(b1.isw1.27,na.rm=TRUE)
table2.isw1[5,8]<-mean(cover.isw1.27[,8],na.rm=TRUE)

#### censoring point at t_0=3 ####
c.0<-5000000
c.1<-35
c.3<-13.78
c.5<-9.06
c.7<-6.36
#### t_0=3 & c=0% ####
b0.isw1.30 <- c()
b0.isw1.sd.30 <- c()
b1.isw1.30 <- c()
b1.isw1.sd.30 <- c()
cover.isw1.30<-matrix(NA,2000,8)
colnames(cover.isw1.30)<-c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

set.seed(1)
for (i in 1:2000){
  a <- data.gen(200,c.0)
  fm <- Surv(Z, censored) ~ X
  tryCatch({
    ismb.fit <- qris(fm, data = a, t0 = 3, Q = 0.5, ne = 200, "smooth", "pmb", c(1.040613, 0))
    b0.isw1.30[i] <- ismb.fit$coefficient[1]
    b0.isw1.sd.30[i] <- ismb.fit$stderr[1]
    b1.isw1.30[i] <- ismb.fit$coefficient[2]
    b1.isw1.sd.30[i] <- ismb.fit$stderr[2]
    # Coverage
    cover.isw1.30[i,1] <- ismb.fit$coefficient[1]-1.96*ismb.fit$stderr[1]
    cover.isw1.30[i,2] <- ismb.fit$coefficient[1]
    cover.isw1.30[i,3] <- ismb.fit$coefficient[1]+1.96*ismb.fit$stderr[1]
    cover.isw1.30[i,4] <- ind(1.040613, cover.isw1.30[i,1], cover.isw1.30[i,3])
    cover.isw1.30[i,5] <- ismb.fit$coefficient[2]-1.96*ismb.fit$stderr[2]
    cover.isw1.30[i,6] <- ismb.fit$coefficient[2]
    cover.isw1.30[i,7] <- ismb.fit$coefficient[2]+1.96*ismb.fit$stderr[2]
    cover.isw1.30[i,8] <- ind(0, cover.isw1.30[i,5], cover.isw1.30[i,7])}
    , error=function(e){
      b0.isw1.30[i] <- NA
      b0.isw1.sd.30[i] <- NA
      b1.isw1.30[i] <- NA
      b1.isw1.sd.30[i] <- NA
      # Coverage
      cover.isw1.30[i,1] <- NA
      cover.isw1.30[i,2] <- NA
      cover.isw1.30[i,3] <- NA
      cover.isw1.30[i,4] <- NA
      cover.isw1.30[i,5] <- NA
      cover.isw1.30[i,6] <- NA
      cover.isw1.30[i,7] <- NA
      cover.isw1.30[i,8] <- NA
    })
}

# IS beta table
table3.isw1[1,1]<-mean(b0.isw1.30,na.rm=TRUE)
table3.isw1[1,2]<-mean(b0.isw1.sd.30,na.rm=TRUE)
table3.isw1[1,3]<-sd(b0.isw1.30,na.rm=TRUE)
table3.isw1[1,4]<-mean(cover.isw1.30[,4],na.rm=TRUE)
table3.isw1[1,5]<-mean(b1.isw1.30,na.rm=TRUE)
table3.isw1[1,6]<-mean(b1.isw1.sd.30,na.rm=TRUE)
table3.isw1[1,7]<-sd(b1.isw1.30,na.rm=TRUE)
table3.isw1[1,8]<-mean(cover.isw1.30[,8],na.rm=TRUE)

#### t_0=3 & c=10% ####
b0.isw1.31 <- c()
b0.isw1.sd.31 <- c()
b1.isw1.31 <- c()
b1.isw1.sd.31 <- c()
cover.isw1.31<-matrix(NA,2000,8)
colnames(cover.isw1.31)<-c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

set.seed(1)
for (i in 1:2000){
  a<-data.gen(200,c.1)
  fm <- Surv(Z, censored) ~ X
  tryCatch({
    ismb.fit <- qris(fm, data = a, t0 = 3, Q = 0.5, ne = 200, "smooth", "pmb", c(1.040613, 0))
    b0.isw1.31[i] <- ismb.fit$coefficient[1]
    b0.isw1.sd.31[i] <- ismb.fit$stderr[1]
    b1.isw1.31[i] <- ismb.fit$coefficient[2]
    b1.isw1.sd.31[i] <- ismb.fit$stderr[2]
    # Coverage
    cover.isw1.31[i,1] <- ismb.fit$coefficient[1]-1.96*ismb.fit$stderr[1]
    cover.isw1.31[i,2] <- ismb.fit$coefficient[1]
    cover.isw1.31[i,3] <- ismb.fit$coefficient[1]+1.96*ismb.fit$stderr[1]
    cover.isw1.31[i,4] <- ind(1.040613, cover.isw1.31[i,1], cover.isw1.31[i,3])
    cover.isw1.31[i,5] <- ismb.fit$coefficient[2]-1.96*ismb.fit$stderr[2]
    cover.isw1.31[i,6] <- ismb.fit$coefficient[2]
    cover.isw1.31[i,7] <- ismb.fit$coefficient[2]+1.96*ismb.fit$stderr[2]
    cover.isw1.31[i,8] <- ind(0, cover.isw1.31[i,5], cover.isw1.31[i,7])}
    , error=function(e){
      b0.isw1.31[i] <- NA
      b0.isw1.sd.31[i] <- NA
      b1.isw1.31[i] <- NA
      b1.isw1.sd.31[i] <- NA
      # Coverage
      cover.isw1.31[i,1] <- NA
      cover.isw1.31[i,2] <- NA
      cover.isw1.31[i,3] <- NA
      cover.isw1.31[i,4] <- NA
      cover.isw1.31[i,5] <- NA
      cover.isw1.31[i,6] <- NA
      cover.isw1.31[i,7] <- NA
      cover.isw1.31[i,8] <- NA
    })
}

# IS beta table
table3.isw1[2,1]<-mean(b0.isw1.31,na.rm=TRUE)
table3.isw1[2,2]<-mean(b0.isw1.sd.31,na.rm=TRUE)
table3.isw1[2,3]<-sd(b0.isw1.31,na.rm=TRUE)
table3.isw1[2,4]<-mean(cover.isw1.31[,4],na.rm=TRUE)
table3.isw1[2,5]<-mean(b1.isw1.31,na.rm=TRUE)
table3.isw1[2,6]<-mean(b1.isw1.sd.31,na.rm=TRUE)
table3.isw1[2,7]<-sd(b1.isw1.31,na.rm=TRUE)
table3.isw1[2,8]<-mean(cover.isw1.31[,8],na.rm=TRUE)

#### t_0=3 & c=30% ####
b0.isw1.33 <- c()
b0.isw1.sd.33 <- c()
b1.isw1.33 <- c()
b1.isw1.sd.33 <- c()
cover.isw1.33<-matrix(NA,2000,8)
colnames(cover.isw1.33)<-c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

set.seed(1)
for (i in 1:2000){
  a<-data.gen(200,c.3)
  fm <- Surv(Z, censored) ~ X
  tryCatch({
    ismb.fit <- qris(fm, data = a, t0 = 3, Q = 0.5, ne = 200, "smooth", "pmb", c(1.040613, 0))
    b0.isw1.33[i] <- ismb.fit$coefficient[1]
    b0.isw1.sd.33[i] <- ismb.fit$stderr[1]
    b1.isw1.33[i] <- ismb.fit$coefficient[2]
    b1.isw1.sd.33[i] <- ismb.fit$stderr[2]
    # Coverage
    cover.isw1.33[i,1] <- ismb.fit$coefficient[1]-1.96*ismb.fit$stderr[1]
    cover.isw1.33[i,2] <- ismb.fit$coefficient[1]
    cover.isw1.33[i,3] <- ismb.fit$coefficient[1]+1.96*ismb.fit$stderr[1]
    cover.isw1.33[i,4] <- ind(1.040613, cover.isw1.33[i,1], cover.isw1.33[i,3])
    cover.isw1.33[i,5] <- ismb.fit$coefficient[2]-1.96*ismb.fit$stderr[2]
    cover.isw1.33[i,6] <- ismb.fit$coefficient[2]
    cover.isw1.33[i,7] <- ismb.fit$coefficient[2]+1.96*ismb.fit$stderr[2]
    cover.isw1.33[i,8] <- ind(0, cover.isw1.33[i,5], cover.isw1.33[i,7])}
    , error=function(e){
      b0.isw1.33[i] <- NA
      b0.isw1.sd.33[i] <- NA
      b1.isw1.33[i] <- NA
      b1.isw1.sd.33[i] <- NA
      # Coverage
      cover.isw1.33[i,1] <- NA
      cover.isw1.33[i,2] <- NA
      cover.isw1.33[i,3] <- NA
      cover.isw1.33[i,4] <- NA
      cover.isw1.33[i,5] <- NA
      cover.isw1.33[i,6] <- NA
      cover.isw1.33[i,7] <- NA
      cover.isw1.33[i,8] <- NA
    })
}

# IS beta table
table3.isw1[3,1]<-mean(b0.isw1.33,na.rm=TRUE)
table3.isw1[3,2]<-mean(b0.isw1.sd.33,na.rm=TRUE)
table3.isw1[3,3]<-sd(b0.isw1.33,na.rm=TRUE)
table3.isw1[3,4]<-mean(cover.isw1.33[,4],na.rm=TRUE)
table3.isw1[3,5]<-mean(b1.isw1.33,na.rm=TRUE)
table3.isw1[3,6]<-mean(b1.isw1.sd.33,na.rm=TRUE)
table3.isw1[3,7]<-sd(b1.isw1.33,na.rm=TRUE)
table3.isw1[3,8]<-mean(cover.isw1.33[,8],na.rm=TRUE)

#### t_0=3 & c=50% ####
b0.isw1.35 <- c()
b0.isw1.sd.35 <- c()
b1.isw1.35 <- c()
b1.isw1.sd.35 <- c()
cover.isw1.35<-matrix(NA,2000,8)
colnames(cover.isw1.35)<-c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
set.seed(1)
for (i in 1:2000){
  a<-data.gen(200,c.5)
  fm <- Surv(Z, censored) ~ X
  tryCatch({
    ismb.fit <- qris(fm, data = a, t0 = 3, Q = 0.5, ne = 200, "smooth", "pmb", c(1.040613, 0))
    b0.isw1.35[i] <- ismb.fit$coefficient[1]
    b0.isw1.sd.35[i] <- ismb.fit$stderr[1]
    b1.isw1.35[i] <- ismb.fit$coefficient[2]
    b1.isw1.sd.35[i] <- ismb.fit$stderr[2]
    # Coverage
    cover.isw1.35[i,1] <- ismb.fit$coefficient[1]-1.96*ismb.fit$stderr[1]
    cover.isw1.35[i,2] <- ismb.fit$coefficient[1]
    cover.isw1.35[i,3] <- ismb.fit$coefficient[1]+1.96*ismb.fit$stderr[1]
    cover.isw1.35[i,4] <- ind(1.040613, cover.isw1.35[i,1], cover.isw1.35[i,3])
    cover.isw1.35[i,5] <- ismb.fit$coefficient[2]-1.96*ismb.fit$stderr[2]
    cover.isw1.35[i,6] <- ismb.fit$coefficient[2]
    cover.isw1.35[i,7] <- ismb.fit$coefficient[2]+1.96*ismb.fit$stderr[2]
    cover.isw1.35[i,8] <- ind(0, cover.isw1.35[i,5], cover.isw1.35[i,7])}
    , error=function(e){
      b0.isw1.35[i] <- NA
      b0.isw1.sd.35[i] <- NA
      b1.isw1.35[i] <- NA
      b1.isw1.sd.35[i] <- NA
      # Coverage
      cover.isw1.35[i,1] <- NA
      cover.isw1.35[i,2] <- NA
      cover.isw1.35[i,3] <- NA
      cover.isw1.35[i,4] <- NA
      cover.isw1.35[i,5] <- NA
      cover.isw1.35[i,6] <- NA
      cover.isw1.35[i,7] <- NA
      cover.isw1.35[i,8] <- NA
    })
}

# IS beta table
table3.isw1[4,1]<-mean(b0.isw1.35,na.rm=TRUE)
table3.isw1[4,2]<-mean(b0.isw1.sd.35,na.rm=TRUE)
table3.isw1[4,3]<-sd(b0.isw1.35,na.rm=TRUE)
table3.isw1[4,4]<-mean(cover.isw1.35[,4],na.rm=TRUE)
table3.isw1[4,5]<-mean(b1.isw1.35,na.rm=TRUE)
table3.isw1[4,6]<-mean(b1.isw1.sd.35,na.rm=TRUE)
table3.isw1[4,7]<-sd(b1.isw1.35,na.rm=TRUE)
table3.isw1[4,8]<-mean(cover.isw1.35[,8],na.rm=TRUE)

#### t_0=3 & c=70% ####
b0.isw1.37 <- c()
b0.isw1.sd.37 <- c()
b1.isw1.37 <- c()
b1.isw1.sd.37 <- c()
cover.isw1.37<-matrix(NA,2000,8)
colnames(cover.isw1.37)<-c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

set.seed(1)
for (i in 1:2000){
  a<-data.gen(200,c.7)
  fm <- Surv(Z, censored) ~ X
  tryCatch({
    ismb.fit <- qris(fm, data = a, t0 = 3, Q = 0.5, ne = 200, "smooth", "pmb", c(1.040613, 0))
    b0.isw1.37[i] <- ismb.fit$coefficient[1]
    b0.isw1.sd.37[i] <- ismb.fit$stderr[1]
    b1.isw1.37[i] <- ismb.fit$coefficient[2]
    b1.isw1.sd.37[i] <- ismb.fit$stderr[2]
    # Coverage
    cover.isw1.37[i,1] <- ismb.fit$coefficient[1]-1.96*ismb.fit$stderr[1]
    cover.isw1.37[i,2] <- ismb.fit$coefficient[1]
    cover.isw1.37[i,3] <- ismb.fit$coefficient[1]+1.96*ismb.fit$stderr[1]
    cover.isw1.37[i,4] <- ind(1.040613, cover.isw1.37[i,1], cover.isw1.37[i,3])
    cover.isw1.37[i,5] <- ismb.fit$coefficient[2]-1.96*ismb.fit$stderr[2]
    cover.isw1.37[i,6] <- ismb.fit$coefficient[2]
    cover.isw1.37[i,7] <- ismb.fit$coefficient[2]+1.96*ismb.fit$stderr[2]
    cover.isw1.37[i,8] <- ind(0, cover.isw1.37[i,5], cover.isw1.37[i,7])}
    , error=function(e){
      b0.isw1.37[i] <- NA
      b0.isw1.sd.37[i] <- NA
      b1.isw1.37[i] <- NA
      b1.isw1.sd.37[i] <- NA
      # Coverage
      cover.isw1.37[i,1] <- NA
      cover.isw1.37[i,2] <- NA
      cover.isw1.37[i,3] <- NA
      cover.isw1.37[i,4] <- NA
      cover.isw1.37[i,5] <- NA
      cover.isw1.37[i,6] <- NA
      cover.isw1.37[i,7] <- NA
      cover.isw1.37[i,8] <- NA
    })
}

# IS beta table
table3.isw1[5,1]<-mean(b0.isw1.37,na.rm=TRUE)
table3.isw1[5,2]<-mean(b0.isw1.sd.37,na.rm=TRUE)
table3.isw1[5,3]<-sd(b0.isw1.37,na.rm=TRUE)
table3.isw1[5,4]<-mean(cover.isw1.37[,4],na.rm=TRUE)
table3.isw1[5,5]<-mean(b1.isw1.37,na.rm=TRUE)
table3.isw1[5,6]<-mean(b1.isw1.sd.37,na.rm=TRUE)
table3.isw1[5,7]<-sd(b1.isw1.37,na.rm=TRUE)
table3.isw1[5,8]<-mean(cover.isw1.37[,8],na.rm=TRUE)


#### Computational code for table 2 ####
#### 1. Scenario and necessart assumption ####
# data size = 200
# beta0 & beta1 effective
# Quantile 50%
# simulation = 2000
# eta = 200

# True Beta
#beta_0    beta_1
#t_0=0 1.609438 0.6931472
#t_0=1 1.410748 0.7974189
#t_0=2 1.219403 0.9070615
#t_0=3 1.040613 1.0174711

# Assumptions
# exp(beta_0) = 5
exp.beta.initial.0<-5
# exp(beta_0 + beta_1) = 10
exp.beta.initial.1<-10
# kappa = 2 (Shape parameter of Weibull distribution)
k<-2
# calculated rho_0 given exp(beta_0) = 5
r.initial.0<-(log(10/5))^(1/k)/exp.beta.initial.0
# calculated rho_1 given exp(beta_0+beta_1) = 10
r.initial.1<-(log(10/5))^(1/k)/exp.beta.initial.1

#### 2. Data Generation function ####
data.gen<-function(samplesize, censor){
  sim<-matrix(NA,samplesize,5)
  colnames(sim) <- c("T","C","Z","X","censored")
  # Generate C_i
  sim[,2] <- runif(samplesize,0,censor)
  # Covariates (Control=0, Treatment=1)
  sim[,4] <- rbinom(samplesize,size=1,p=0.5)
  # Generate T_i (Given Condition r=rho_0, k=2, exp(beta_0)=5, exp(beta_1)=2))
  unif <- runif(n=samplesize ,min = 0,max = 1)
  for (q in 1:samplesize){
    if (sim[q,4]==0){
      sim[q,1]<-{{-log(1-unif[q])}^(1/k)}/r.initial.0
    } else {
      sim[q,1]<-{{-log(1-unif[q])}^(1/k)}/r.initial.1
    }
  }
  # Generate Y_i (min(T,C))
  sim[,3] <- apply(sim[,1:2], 1, FUN=min)
  # Censoring indicator (Censored=0, Not censored=1)
  sim[,5]<-I(sim[,1]<sim[,2])
  # Ordering
  sim <- sim[order(sim[,3]),]
  n <- nrow(sim)
  sim <- as.data.frame(sim)
  return(sim)
}

#### 3. Indicator function for measuring coverage proportion of 95% CI ####
ind <- function(a,b,c){
  if (a>=b&a<=c) {
    result <- 1
  } else {
    result <- 0
  }
  print(result)
}

#### 4. Make table for Beta estimation and variance estimation and Coverage ####
table0.isw1<-matrix(NA,5,8)
rownames(table0.isw1)<-c(0,10,30,50,70)
colnames(table0.isw1)<-c("beta_0","SE of beta_0","SD of beta_0","Coverage of beta_0","beta_1","SE of beta_1","SD of beta_1","Coverage of beta_1")
table1.isw1<-matrix(NA,5,8)
rownames(table1.isw1)<-c(0,10,30,50,70)
colnames(table1.isw1)<-c("beta_0","SE of beta_0","SD of beta_0","Coverage of beta_0","beta_1","SE of beta_1","SD of beta_1","Coverage of beta_1")
table2.isw1<-matrix(NA,5,8)
rownames(table2.isw1)<-c(0,10,30,50,70)
colnames(table2.isw1)<-c("beta_0","SE of beta_0","SD of beta_0","Coverage of beta_0","beta_1","SE of beta_1","SD of beta_1","Coverage of beta_1")
table3.isw1<-matrix(NA,5,8)
rownames(table3.isw1)<-c(0,10,30,50,70)
colnames(table3.isw1)<-c("beta_0","SE of beta_0","SD of beta_0","Coverage of beta_0","beta_1","SE of beta_1","SD of beta_1","Coverage of beta_1")

#### 5. Execute simulation ####
#### censoring point at t_0=0 ####
c.0<-5000000
c.1<-78.11
c.3<-26.36
c.5<-15.08
c.7<-9.09

#### t_0=0 & c=0% ####
b0.isw1.00 <- c()
b0.isw1.sd.00 <- c()
b1.isw1.00 <- c()
b1.isw1.sd.00 <- c()
cover.isw1.00<-matrix(NA,2000,8)
colnames(cover.isw1.00)<-c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

set.seed(1)
for (i in 1:2000){
  a <- data.gen(200,c.0)
  fm <- Surv(Z, censored) ~ X
  tryCatch({
    ismb.fit <- qris(fm, data = a, t0 = 0, Q = 0.5, ne = 200, "smooth", "pmb", c(1.609438,0.6931472))
    b0.isw1.00[i] <- ismb.fit$coefficient[1]
    b0.isw1.sd.00[i] <- ismb.fit$stderr[1]
    b1.isw1.00[i] <- ismb.fit$coefficient[2]
    b1.isw1.sd.00[i] <- ismb.fit$stderr[2]
    # Coverage
    cover.isw1.00[i,1] <- ismb.fit$coefficient[1]-1.96*ismb.fit$stderr[1]
    cover.isw1.00[i,2] <- ismb.fit$coefficient[1]
    cover.isw1.00[i,3] <- ismb.fit$coefficient[1]+1.96*ismb.fit$stderr[1]
    cover.isw1.00[i,4] <- ind(1.609438, cover.isw1.00[i,1], cover.isw1.00[i,3])
    cover.isw1.00[i,5] <- ismb.fit$coefficient[2]-1.96*ismb.fit$stderr[2]
    cover.isw1.00[i,6] <- ismb.fit$coefficient[2]
    cover.isw1.00[i,7] <- ismb.fit$coefficient[2]+1.96*ismb.fit$stderr[2]
    cover.isw1.00[i,8] <- ind(0.6931472, cover.isw1.00[i,5], cover.isw1.00[i,7])}
    , error=function(e){
      b0.isw1.00[i] <- NA
      b0.isw1.sd.00[i] <- NA
      b1.isw1.00[i] <- NA
      b1.isw1.sd.00[i] <- NA
      # Coverage
      cover.isw1.00[i,1] <- NA
      cover.isw1.00[i,2] <- NA
      cover.isw1.00[i,3] <- NA
      cover.isw1.00[i,4] <- NA
      cover.isw1.00[i,5] <- NA
      cover.isw1.00[i,6] <- NA
      cover.isw1.00[i,7] <- NA
      cover.isw1.00[i,8] <- NA
    })
}
# IS beta table
table0.isw1[1,1]<-mean(b0.isw1.00,na.rm=TRUE)
table0.isw1[1,2]<-mean(b0.isw1.sd.00,na.rm=TRUE)
table0.isw1[1,3]<-sd(b0.isw1.00,na.rm=TRUE)
table0.isw1[1,4]<-mean(cover.isw1.00[,4],na.rm=TRUE)
table0.isw1[1,5]<-mean(b1.isw1.00,na.rm=TRUE)
table0.isw1[1,6]<-mean(b1.isw1.sd.00,na.rm=TRUE)
table0.isw1[1,7]<-sd(b1.isw1.00,na.rm=TRUE)
table0.isw1[1,8]<-mean(cover.isw1.00[,8],na.rm=TRUE)

#### t_0=0 & c=10% ####
b0.isw1.01 <- c()
b0.isw1.sd.01 <- c()
b1.isw1.01 <- c()
b1.isw1.sd.01 <- c()
cover.isw1.01<-matrix(NA,2000,8)
colnames(cover.isw1.01)<-c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

set.seed(1)
for (i in 1:2000){
  a<-data.gen(200,c.1)
  fm <- Surv(Z, censored) ~ X
  tryCatch({
    ismb.fit <- qris(fm, data = a, t0 = 0, Q = 0.5, ne = 200, "smooth", "pmb", c(1.609438,0.6931472))
    b0.isw1.01[i] <- ismb.fit$coefficient[1]
    b0.isw1.sd.01[i] <- ismb.fit$stderr[1]
    b1.isw1.01[i] <- ismb.fit$coefficient[2]
    b1.isw1.sd.01[i] <- ismb.fit$stderr[2]
    # Coverage
    cover.isw1.01[i,1] <- ismb.fit$coefficient[1]-1.96*ismb.fit$stderr[1]
    cover.isw1.01[i,2] <- ismb.fit$coefficient[1]
    cover.isw1.01[i,3] <- ismb.fit$coefficient[1]+1.96*ismb.fit$stderr[1]
    cover.isw1.01[i,4] <- ind(1.609438, cover.isw1.01[i,1], cover.isw1.01[i,3])
    cover.isw1.01[i,5] <- ismb.fit$coefficient[2]-1.96*ismb.fit$stderr[2]
    cover.isw1.01[i,6] <- ismb.fit$coefficient[2]
    cover.isw1.01[i,7] <- ismb.fit$coefficient[2]+1.96*ismb.fit$stderr[2]
    cover.isw1.01[i,8] <- ind(0.6931472, cover.isw1.01[i,5], cover.isw1.01[i,7])}
    , error=function(e){
      b0.isw1.01[i] <- NA
      b0.isw1.sd.01[i] <- NA
      b1.isw1.01[i] <- NA
      b1.isw1.sd.01[i] <- NA
      # Coverage
      cover.isw1.01[i,1] <- NA
      cover.isw1.01[i,2] <- NA
      cover.isw1.01[i,3] <- NA
      cover.isw1.01[i,4] <- NA
      cover.isw1.01[i,5] <- NA
      cover.isw1.01[i,6] <- NA
      cover.isw1.01[i,7] <- NA
      cover.isw1.01[i,8] <- NA
    })
}

# IS beta table
table0.isw1[2,1]<-mean(b0.isw1.01,na.rm=TRUE)
table0.isw1[2,2]<-mean(b0.isw1.sd.01,na.rm=TRUE)
table0.isw1[2,3]<-sd(b0.isw1.01,na.rm=TRUE)
table0.isw1[2,4]<-mean(cover.isw1.01[,4],na.rm=TRUE)
table0.isw1[2,5]<-mean(b1.isw1.01,na.rm=TRUE)
table0.isw1[2,6]<-mean(b1.isw1.sd.01,na.rm=TRUE)
table0.isw1[2,7]<-sd(b1.isw1.01,na.rm=TRUE)
table0.isw1[2,8]<-mean(cover.isw1.01[,8],na.rm=TRUE)

#### t_0=0 & c=30% ####
b0.isw1.03 <- c()
b0.isw1.sd.03 <- c()
b1.isw1.03 <- c()
b1.isw1.sd.03 <- c()
cover.isw1.03<-matrix(NA,2000,8)
colnames(cover.isw1.03)<-c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

set.seed(1)
for (i in 1:2000){
  a<-data.gen(200,c.3)
  fm <- Surv(Z, censored) ~ X
  tryCatch({
    ismb.fit <- qris(fm, data = a, t0 = 0, Q = 0.5, ne = 200, "smooth", "pmb", c(1.609438,0.6931472))
    b0.isw1.03[i] <- ismb.fit$coefficient[1]
    b0.isw1.sd.03[i] <- ismb.fit$stderr[1]
    b1.isw1.03[i] <- ismb.fit$coefficient[2]
    b1.isw1.sd.03[i] <- ismb.fit$stderr[2]
    # Coverage
    cover.isw1.03[i,1] <- ismb.fit$coefficient[1]-1.96*ismb.fit$stderr[1]
    cover.isw1.03[i,2] <- ismb.fit$coefficient[1]
    cover.isw1.03[i,3] <- ismb.fit$coefficient[1]+1.96*ismb.fit$stderr[1]
    cover.isw1.03[i,4] <- ind(1.609438, cover.isw1.03[i,1], cover.isw1.03[i,3])
    cover.isw1.03[i,5] <- ismb.fit$coefficient[2]-1.96*ismb.fit$stderr[2]
    cover.isw1.03[i,6] <- ismb.fit$coefficient[2]
    cover.isw1.03[i,7] <- ismb.fit$coefficient[2]+1.96*ismb.fit$stderr[2]
    cover.isw1.03[i,8] <- ind(0.6931472, cover.isw1.03[i,5], cover.isw1.03[i,7])}
    , error=function(e){
      b0.isw1.03[i] <- NA
      b0.isw1.sd.03[i] <- NA
      b1.isw1.03[i] <- NA
      b1.isw1.sd.03[i] <- NA
      # Coverage
      cover.isw1.03[i,1] <- NA
      cover.isw1.03[i,2] <- NA
      cover.isw1.03[i,3] <- NA
      cover.isw1.03[i,4] <- NA
      cover.isw1.03[i,5] <- NA
      cover.isw1.03[i,6] <- NA
      cover.isw1.03[i,7] <- NA
      cover.isw1.03[i,8] <- NA
    })
}

# IS beta table
table0.isw1[3,1]<-mean(b0.isw1.03,na.rm=TRUE)
table0.isw1[3,2]<-mean(b0.isw1.sd.03,na.rm=TRUE)
table0.isw1[3,3]<-sd(b0.isw1.03,na.rm=TRUE)
table0.isw1[3,4]<-mean(cover.isw1.03[,4],na.rm=TRUE)
table0.isw1[3,5]<-mean(b1.isw1.03,na.rm=TRUE)
table0.isw1[3,6]<-mean(b1.isw1.sd.03,na.rm=TRUE)
table0.isw1[3,7]<-sd(b1.isw1.03,na.rm=TRUE)
table0.isw1[3,8]<-mean(cover.isw1.03[,8],na.rm=TRUE)

#### t_0=0 & c=50% ####
b0.isw1.05 <- c()
b0.isw1.sd.05 <- c()
b1.isw1.05 <- c()
b1.isw1.sd.05 <- c()
cover.isw1.05<-matrix(NA,2000,8)
colnames(cover.isw1.05)<-c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

set.seed(1)
for (i in 1:2000){
  a<-data.gen(200,c.5)
  fm <- Surv(Z, censored) ~ X
  tryCatch({
    ismb.fit <- qris(fm, data = a, t0 = 0, Q = 0.5, ne = 200, "smooth", "pmb", c(1.609438,0.6931472))
    b0.isw1.05[i] <- ismb.fit$coefficient[1]
    b0.isw1.sd.05[i] <- ismb.fit$stderr[1]
    b1.isw1.05[i] <- ismb.fit$coefficient[2]
    b1.isw1.sd.05[i] <- ismb.fit$stderr[2]
    # Coverage
    cover.isw1.05[i,1] <- ismb.fit$coefficient[1]-1.96*ismb.fit$stderr[1]
    cover.isw1.05[i,2] <- ismb.fit$coefficient[1]
    cover.isw1.05[i,3] <- ismb.fit$coefficient[1]+1.96*ismb.fit$stderr[1]
    cover.isw1.05[i,4] <- ind(1.609438, cover.isw1.05[i,1], cover.isw1.05[i,3])
    cover.isw1.05[i,5] <- ismb.fit$coefficient[2]-1.96*ismb.fit$stderr[2]
    cover.isw1.05[i,6] <- ismb.fit$coefficient[2]
    cover.isw1.05[i,7] <- ismb.fit$coefficient[2]+1.96*ismb.fit$stderr[2]
    cover.isw1.05[i,8] <- ind(0.6931472, cover.isw1.05[i,5], cover.isw1.05[i,7])}
    , error=function(e){
      b0.isw1.05[i] <- NA
      b0.isw1.sd.05[i] <- NA
      b1.isw1.05[i] <- NA
      b1.isw1.sd.05[i] <- NA
      # Coverage
      cover.isw1.05[i,1] <- NA
      cover.isw1.05[i,2] <- NA
      cover.isw1.05[i,3] <- NA
      cover.isw1.05[i,4] <- NA
      cover.isw1.05[i,5] <- NA
      cover.isw1.05[i,6] <- NA
      cover.isw1.05[i,7] <- NA
      cover.isw1.05[i,8] <- NA
    })
}

# IS beta table
table0.isw1[4,1]<-mean(b0.isw1.05,na.rm=TRUE)
table0.isw1[4,2]<-mean(b0.isw1.sd.05,na.rm=TRUE)
table0.isw1[4,3]<-sd(b0.isw1.05,na.rm=TRUE)
table0.isw1[4,4]<-mean(cover.isw1.05[,4],na.rm=TRUE)
table0.isw1[4,5]<-mean(b1.isw1.05,na.rm=TRUE)
table0.isw1[4,6]<-mean(b1.isw1.sd.05,na.rm=TRUE)
table0.isw1[4,7]<-sd(b1.isw1.05,na.rm=TRUE)
table0.isw1[4,8]<-mean(cover.isw1.05[,8],na.rm=TRUE)

#### t_0=0 & c=70% ####
b0.isw1.07 <- c()
b0.isw1.sd.07 <- c()
b1.isw1.07 <- c()
b1.isw1.sd.07 <- c()
cover.isw1.07<-matrix(NA,2000,8)
colnames(cover.isw1.07)<-c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

set.seed(1)
for (i in 1:2000){
  a<-data.gen(200,c.7)
  fm <- Surv(Z, censored) ~ X
  tryCatch({
    ismb.fit <- qris(fm, data = a, t0 = 0, Q = 0.5, ne = 200, "smooth", "pmb", c(1.609438,0.6931472))
    b0.isw1.07[i] <- ismb.fit$coefficient[1]
    b0.isw1.sd.07[i] <- ismb.fit$stderr[1]
    b1.isw1.07[i] <- ismb.fit$coefficient[2]
    b1.isw1.sd.07[i] <- ismb.fit$stderr[2]
    # Coverage
    cover.isw1.07[i,1] <- ismb.fit$coefficient[1]-1.96*ismb.fit$stderr[1]
    cover.isw1.07[i,2] <- ismb.fit$coefficient[1]
    cover.isw1.07[i,3] <- ismb.fit$coefficient[1]+1.96*ismb.fit$stderr[1]
    cover.isw1.07[i,4] <- ind(1.609438, cover.isw1.07[i,1], cover.isw1.07[i,3])
    cover.isw1.07[i,5] <- ismb.fit$coefficient[2]-1.96*ismb.fit$stderr[2]
    cover.isw1.07[i,6] <- ismb.fit$coefficient[2]
    cover.isw1.07[i,7] <- ismb.fit$coefficient[2]+1.96*ismb.fit$stderr[2]
    cover.isw1.07[i,8] <- ind(0.6931472, cover.isw1.07[i,5], cover.isw1.07[i,7])}
    , error=function(e){
      b0.isw1.07[i] <- NA
      b0.isw1.sd.07[i] <- NA
      b1.isw1.07[i] <- NA
      b1.isw1.sd.07[i] <- NA
      # Coverage
      cover.isw1.07[i,1] <- NA
      cover.isw1.07[i,2] <- NA
      cover.isw1.07[i,3] <- NA
      cover.isw1.07[i,4] <- NA
      cover.isw1.07[i,5] <- NA
      cover.isw1.07[i,6] <- NA
      cover.isw1.07[i,7] <- NA
      cover.isw1.07[i,8] <- NA
    })
}

# IS beta table
table0.isw1[5,1]<-mean(b0.isw1.07,na.rm=TRUE)
table0.isw1[5,2]<-mean(b0.isw1.sd.07,na.rm=TRUE)
table0.isw1[5,3]<-sd(b0.isw1.07,na.rm=TRUE)
table0.isw1[5,4]<-mean(cover.isw1.07[,4],na.rm=TRUE)
table0.isw1[5,5]<-mean(b1.isw1.07,na.rm=TRUE)
table0.isw1[5,6]<-mean(b1.isw1.sd.07,na.rm=TRUE)
table0.isw1[5,7]<-sd(b1.isw1.07,na.rm=TRUE)
table0.isw1[5,8]<-mean(cover.isw1.07[,8],na.rm=TRUE)

#### censoring point at t_0=1 ####
c.0<-5000000
c.1<-70.39
c.3<-24.35
c.5<-14.07
c.7<-8.49
#### t_0=1 & c=0% ####
b0.isw1.10 <- c()
b0.isw1.sd.10 <- c()
b1.isw1.10 <- c()
b1.isw1.sd.10 <- c()
cover.isw1.10<-matrix(NA,2000,8)
colnames(cover.isw1.10)<-c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

set.seed(1)
for (i in 1:2000){
  a <- data.gen(200,c.0)
  fm <- Surv(Z, censored) ~ X
  tryCatch({
    ismb.fit <- qris(fm, data = a, t0 = 1, Q = 0.5, ne = 200, "smooth", "pmb", c(1.410748, 0.7974189))
    b0.isw1.10[i] <- ismb.fit$coefficient[1]
    b0.isw1.sd.10[i] <- ismb.fit$stderr[1]
    b1.isw1.10[i] <- ismb.fit$coefficient[2]
    b1.isw1.sd.10[i] <- ismb.fit$stderr[2]
    # Coverage
    cover.isw1.10[i,1] <- ismb.fit$coefficient[1]-1.96*ismb.fit$stderr[1]
    cover.isw1.10[i,2] <- ismb.fit$coefficient[1]
    cover.isw1.10[i,3] <- ismb.fit$coefficient[1]+1.96*ismb.fit$stderr[1]
    cover.isw1.10[i,4] <- ind(1.410748, cover.isw1.10[i,1], cover.isw1.10[i,3])
    cover.isw1.10[i,5] <- ismb.fit$coefficient[2]-1.96*ismb.fit$stderr[2]
    cover.isw1.10[i,6] <- ismb.fit$coefficient[2]
    cover.isw1.10[i,7] <- ismb.fit$coefficient[2]+1.96*ismb.fit$stderr[2]
    cover.isw1.10[i,8] <- ind(0.7974189, cover.isw1.10[i,5], cover.isw1.10[i,7])}
    , error=function(e){
      b0.isw1.10[i] <- NA
      b0.isw1.sd.10[i] <- NA
      b1.isw1.10[i] <- NA
      b1.isw1.sd.10[i] <- NA
      # Coverage
      cover.isw1.10[i,1] <- NA
      cover.isw1.10[i,2] <- NA
      cover.isw1.10[i,3] <- NA
      cover.isw1.10[i,4] <- NA
      cover.isw1.10[i,5] <- NA
      cover.isw1.10[i,6] <- NA
      cover.isw1.10[i,7] <- NA
      cover.isw1.10[i,8] <- NA
    })
}

# IS beta table
table1.isw1[1,1]<-mean(b0.isw1.10,na.rm=TRUE)
table1.isw1[1,2]<-mean(b0.isw1.sd.10,na.rm=TRUE)
table1.isw1[1,3]<-sd(b0.isw1.10,na.rm=TRUE)
table1.isw1[1,4]<-mean(cover.isw1.10[,4],na.rm=TRUE)
table1.isw1[1,5]<-mean(b1.isw1.10,na.rm=TRUE)
table1.isw1[1,6]<-mean(b1.isw1.sd.10,na.rm=TRUE)
table1.isw1[1,7]<-sd(b1.isw1.10,na.rm=TRUE)
table1.isw1[1,8]<-mean(cover.isw1.10[,8],na.rm=TRUE)

#### t_0=1 & c=10% ####
b0.isw1.11 <- c()
b0.isw1.sd.11 <- c()
b1.isw1.11 <- c()
b1.isw1.sd.11 <- c()
cover.isw1.11<-matrix(NA,2000,8)
colnames(cover.isw1.11)<-c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

set.seed(1)
for (i in 1:2000){
  a<-data.gen(200,c.1)
  fm <- Surv(Z, censored) ~ X
  tryCatch({
    ismb.fit <- qris(fm, data = a, t0 = 1, Q = 0.5, ne = 200, "smooth", "pmb", c(1.410748, 0.7974189))
    b0.isw1.11[i] <- ismb.fit$coefficient[1]
    b0.isw1.sd.11[i] <- ismb.fit$stderr[1]
    b1.isw1.11[i] <- ismb.fit$coefficient[2]
    b1.isw1.sd.11[i] <- ismb.fit$stderr[2]
    # Coverage
    cover.isw1.11[i,1] <- ismb.fit$coefficient[1]-1.96*ismb.fit$stderr[1]
    cover.isw1.11[i,2] <- ismb.fit$coefficient[1]
    cover.isw1.11[i,3] <- ismb.fit$coefficient[1]+1.96*ismb.fit$stderr[1]
    cover.isw1.11[i,4] <- ind(1.410748, cover.isw1.11[i,1], cover.isw1.11[i,3])
    cover.isw1.11[i,5] <- ismb.fit$coefficient[2]-1.96*ismb.fit$stderr[2]
    cover.isw1.11[i,6] <- ismb.fit$coefficient[2]
    cover.isw1.11[i,7] <- ismb.fit$coefficient[2]+1.96*ismb.fit$stderr[2]
    cover.isw1.11[i,8] <- ind(0.7974189, cover.isw1.11[i,5], cover.isw1.11[i,7])}
    , error=function(e){
      b0.isw1.11[i] <- NA
      b0.isw1.sd.11[i] <- NA
      b1.isw1.11[i] <- NA
      b1.isw1.sd.11[i] <- NA
      # Coverage
      cover.isw1.11[i,1] <- NA
      cover.isw1.11[i,2] <- NA
      cover.isw1.11[i,3] <- NA
      cover.isw1.11[i,4] <- NA
      cover.isw1.11[i,5] <- NA
      cover.isw1.11[i,6] <- NA
      cover.isw1.11[i,7] <- NA
      cover.isw1.11[i,8] <- NA
    })
}

# IS beta table
table1.isw1[2,1]<-mean(b0.isw1.11,na.rm=TRUE)
table1.isw1[2,2]<-mean(b0.isw1.sd.11,na.rm=TRUE)
table1.isw1[2,3]<-sd(b0.isw1.11,na.rm=TRUE)
table1.isw1[2,4]<-mean(cover.isw1.11[,4],na.rm=TRUE)
table1.isw1[2,5]<-mean(b1.isw1.11,na.rm=TRUE)
table1.isw1[2,6]<-mean(b1.isw1.sd.11,na.rm=TRUE)
table1.isw1[2,7]<-sd(b1.isw1.11,na.rm=TRUE)
table1.isw1[2,8]<-mean(cover.isw1.11[,8],na.rm=TRUE)

#### t_0=1 & c=30% ####
b0.isw1.13 <- c()
b0.isw1.sd.13 <- c()
b1.isw1.13 <- c()
b1.isw1.sd.13 <- c()
cover.isw1.13<-matrix(NA,2000,8)
colnames(cover.isw1.13)<-c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

set.seed(1)
for (i in 1:2000){
  a<-data.gen(200,c.3)
  fm <- Surv(Z, censored) ~ X
  tryCatch({
    ismb.fit <- qris(fm, data = a, t0 = 1, Q = 0.5, ne = 200, "smooth", "pmb", c(1.410748, 0.7974189))
    b0.isw1.13[i] <- ismb.fit$coefficient[1]
    b0.isw1.sd.13[i] <- ismb.fit$stderr[1]
    b1.isw1.13[i] <- ismb.fit$coefficient[2]
    b1.isw1.sd.13[i] <- ismb.fit$stderr[2]
    # Coverage
    cover.isw1.13[i,1] <- ismb.fit$coefficient[1]-1.96*ismb.fit$stderr[1]
    cover.isw1.13[i,2] <- ismb.fit$coefficient[1]
    cover.isw1.13[i,3] <- ismb.fit$coefficient[1]+1.96*ismb.fit$stderr[1]
    cover.isw1.13[i,4] <- ind(1.410748, cover.isw1.13[i,1], cover.isw1.13[i,3])
    cover.isw1.13[i,5] <- ismb.fit$coefficient[2]-1.96*ismb.fit$stderr[2]
    cover.isw1.13[i,6] <- ismb.fit$coefficient[2]
    cover.isw1.13[i,7] <- ismb.fit$coefficient[2]+1.96*ismb.fit$stderr[2]
    cover.isw1.13[i,8] <- ind(0.7974189, cover.isw1.13[i,5], cover.isw1.13[i,7])}
    , error=function(e){
      b0.isw1.13[i] <- NA
      b0.isw1.sd.13[i] <- NA
      b1.isw1.13[i] <- NA
      b1.isw1.sd.13[i] <- NA
      # Coverage
      cover.isw1.13[i,1] <- NA
      cover.isw1.13[i,2] <- NA
      cover.isw1.13[i,3] <- NA
      cover.isw1.13[i,4] <- NA
      cover.isw1.13[i,5] <- NA
      cover.isw1.13[i,6] <- NA
      cover.isw1.13[i,7] <- NA
      cover.isw1.13[i,8] <- NA
    })
}

# IS beta table
table1.isw1[3,1]<-mean(b0.isw1.13,na.rm=TRUE)
table1.isw1[3,2]<-mean(b0.isw1.sd.13,na.rm=TRUE)
table1.isw1[3,3]<-sd(b0.isw1.13,na.rm=TRUE)
table1.isw1[3,4]<-mean(cover.isw1.13[,4],na.rm=TRUE)
table1.isw1[3,5]<-mean(b1.isw1.13,na.rm=TRUE)
table1.isw1[3,6]<-mean(b1.isw1.sd.13,na.rm=TRUE)
table1.isw1[3,7]<-sd(b1.isw1.13,na.rm=TRUE)
table1.isw1[3,8]<-mean(cover.isw1.13[,8],na.rm=TRUE)

#### t_0=1 & c=50% ####
b0.isw1.15 <- c()
b0.isw1.sd.15 <- c()
b1.isw1.15 <- c()
b1.isw1.sd.15 <- c()
cover.isw1.15<-matrix(NA,2000,8)
colnames(cover.isw1.15)<-c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

set.seed(1)
for (i in 1:2000){
  a<-data.gen(200,c.5)
  fm <- Surv(Z, censored) ~ X
  tryCatch({
    ismb.fit <- qris(fm, data = a, t0 = 1, Q = 0.5, ne = 200, "smooth", "pmb", c(1.410748, 0.7974189))
    b0.isw1.15[i] <- ismb.fit$coefficient[1]
    b0.isw1.sd.15[i] <- ismb.fit$stderr[1]
    b1.isw1.15[i] <- ismb.fit$coefficient[2]
    b1.isw1.sd.15[i] <- ismb.fit$stderr[2]
    # Coverage
    cover.isw1.15[i,1] <- ismb.fit$coefficient[1]-1.96*ismb.fit$stderr[1]
    cover.isw1.15[i,2] <- ismb.fit$coefficient[1]
    cover.isw1.15[i,3] <- ismb.fit$coefficient[1]+1.96*ismb.fit$stderr[1]
    cover.isw1.15[i,4] <- ind(1.410748, cover.isw1.15[i,1], cover.isw1.15[i,3])
    cover.isw1.15[i,5] <- ismb.fit$coefficient[2]-1.96*ismb.fit$stderr[2]
    cover.isw1.15[i,6] <- ismb.fit$coefficient[2]
    cover.isw1.15[i,7] <- ismb.fit$coefficient[2]+1.96*ismb.fit$stderr[2]
    cover.isw1.15[i,8] <- ind(0.7974189, cover.isw1.15[i,5], cover.isw1.15[i,7])}
    , error=function(e){
      b0.isw1.15[i] <- NA
      b0.isw1.sd.15[i] <- NA
      b1.isw1.15[i] <- NA
      b1.isw1.sd.15[i] <- NA
      # Coverage
      cover.isw1.15[i,1] <- NA
      cover.isw1.15[i,2] <- NA
      cover.isw1.15[i,3] <- NA
      cover.isw1.15[i,4] <- NA
      cover.isw1.15[i,5] <- NA
      cover.isw1.15[i,6] <- NA
      cover.isw1.15[i,7] <- NA
      cover.isw1.15[i,8] <- NA
    })
}

# IS beta table
table1.isw1[4,1]<-mean(b0.isw1.15,na.rm=TRUE)
table1.isw1[4,2]<-mean(b0.isw1.sd.15,na.rm=TRUE)
table1.isw1[4,3]<-sd(b0.isw1.15,na.rm=TRUE)
table1.isw1[4,4]<-mean(cover.isw1.15[,4],na.rm=TRUE)
table1.isw1[4,5]<-mean(b1.isw1.15,na.rm=TRUE)
table1.isw1[4,6]<-mean(b1.isw1.sd.15,na.rm=TRUE)
table1.isw1[4,7]<-sd(b1.isw1.15,na.rm=TRUE)
table1.isw1[4,8]<-mean(cover.isw1.15[,8],na.rm=TRUE)

#### t_0=1 & c=70% ####
b0.isw1.17 <- c()
b0.isw1.sd.17 <- c()
b1.isw1.17 <- c()
b1.isw1.sd.17 <- c()
cover.isw1.17<-matrix(NA,2000,8)
colnames(cover.isw1.17)<-c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

set.seed(1)
for (i in 1:2000){
  a<-data.gen(200,c.7)
  fm <- Surv(Z, censored) ~ X
  tryCatch({
    ismb.fit <- qris(fm, data = a, t0 = 1, Q = 0.5, ne = 200, "smooth", "pmb", c(1.410748, 0.7974189))
    b0.isw1.17[i] <- ismb.fit$coefficient[1]
    b0.isw1.sd.17[i] <- ismb.fit$stderr[1]
    b1.isw1.17[i] <- ismb.fit$coefficient[2]
    b1.isw1.sd.17[i] <- ismb.fit$stderr[2]
    # Coverage
    cover.isw1.17[i,1] <- ismb.fit$coefficient[1]-1.96*ismb.fit$stderr[1]
    cover.isw1.17[i,2] <- ismb.fit$coefficient[1]
    cover.isw1.17[i,3] <- ismb.fit$coefficient[1]+1.96*ismb.fit$stderr[1]
    cover.isw1.17[i,4] <- ind(1.410748, cover.isw1.17[i,1], cover.isw1.17[i,3])
    cover.isw1.17[i,5] <- ismb.fit$coefficient[2]-1.96*ismb.fit$stderr[2]
    cover.isw1.17[i,6] <- ismb.fit$coefficient[2]
    cover.isw1.17[i,7] <- ismb.fit$coefficient[2]+1.96*ismb.fit$stderr[2]
    cover.isw1.17[i,8] <- ind(0.7974189, cover.isw1.17[i,5], cover.isw1.17[i,7])}
    , error=function(e){
      b0.isw1.17[i] <- NA
      b0.isw1.sd.17[i] <- NA
      b1.isw1.17[i] <- NA
      b1.isw1.sd.17[i] <- NA
      # Coverage
      cover.isw1.17[i,1] <- NA
      cover.isw1.17[i,2] <- NA
      cover.isw1.17[i,3] <- NA
      cover.isw1.17[i,4] <- NA
      cover.isw1.17[i,5] <- NA
      cover.isw1.17[i,6] <- NA
      cover.isw1.17[i,7] <- NA
      cover.isw1.17[i,8] <- NA
    })
  
}

# IS beta table
table1.isw1[5,1]<-mean(b0.isw1.17,na.rm=TRUE)
table1.isw1[5,2]<-mean(b0.isw1.sd.17,na.rm=TRUE)
table1.isw1[5,3]<-sd(b0.isw1.17,na.rm=TRUE)
table1.isw1[5,4]<-mean(cover.isw1.17[,4],na.rm=TRUE)
table1.isw1[5,5]<-mean(b1.isw1.17,na.rm=TRUE)
table1.isw1[5,6]<-mean(b1.isw1.sd.17,na.rm=TRUE)
table1.isw1[5,7]<-sd(b1.isw1.17,na.rm=TRUE)
table1.isw1[5,8]<-mean(cover.isw1.17[,8],na.rm=TRUE)

#### censoring point at t_0=2 ####
c.0<-5000000
c.1<-64.86
c.3<-23.34
c.5<-13.62
c.7<-8.36
#### t_0=2 & c=0% ####
b0.isw1.20 <- c()
b0.isw1.sd.20 <- c()
b1.isw1.20 <- c()
b1.isw1.sd.20 <- c()
cover.isw1.20<-matrix(NA,2000,8)
colnames(cover.isw1.20)<-c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

set.seed(1)
for (i in 1:2000){
  a <- data.gen(200,c.0)
  fm <- Surv(Z, censored) ~ X
  tryCatch({
    ismb.fit <- qris(fm, data = a, t0 = 2, Q = 0.5, ne = 200, "smooth", "pmb", c(1.219403, 0.9070615))
    b0.isw1.20[i] <- ismb.fit$coefficient[1]
    b0.isw1.sd.20[i] <- ismb.fit$stderr[1]
    b1.isw1.20[i] <- ismb.fit$coefficient[2]
    b1.isw1.sd.20[i] <- ismb.fit$stderr[2]
    # Coverage
    cover.isw1.20[i,1] <- ismb.fit$coefficient[1]-1.96*ismb.fit$stderr[1]
    cover.isw1.20[i,2] <- ismb.fit$coefficient[1]
    cover.isw1.20[i,3] <- ismb.fit$coefficient[1]+1.96*ismb.fit$stderr[1]
    cover.isw1.20[i,4] <- ind(1.219403, cover.isw1.20[i,1], cover.isw1.20[i,3])
    cover.isw1.20[i,5] <- ismb.fit$coefficient[2]-1.96*ismb.fit$stderr[2]
    cover.isw1.20[i,6] <- ismb.fit$coefficient[2]
    cover.isw1.20[i,7] <- ismb.fit$coefficient[2]+1.96*ismb.fit$stderr[2]
    cover.isw1.20[i,8] <- ind(0.9070615, cover.isw1.20[i,5], cover.isw1.20[i,7])}
    , error=function(e){
      b0.isw1.20[i] <- NA
      b0.isw1.sd.20[i] <- NA
      b1.isw1.20[i] <- NA
      b1.isw1.sd.20[i] <- NA
      # Coverage
      cover.isw1.20[i,1] <- NA
      cover.isw1.20[i,2] <- NA
      cover.isw1.20[i,3] <- NA
      cover.isw1.20[i,4] <- NA
      cover.isw1.20[i,5] <- NA
      cover.isw1.20[i,6] <- NA
      cover.isw1.20[i,7] <- NA
      cover.isw1.20[i,8] <- NA
    })
}

# IS beta table
table2.isw1[1,1]<-mean(b0.isw1.20,na.rm=TRUE)
table2.isw1[1,2]<-mean(b0.isw1.sd.20,na.rm=TRUE)
table2.isw1[1,3]<-sd(b0.isw1.20,na.rm=TRUE)
table2.isw1[1,4]<-mean(cover.isw1.20[,4],na.rm=TRUE)
table2.isw1[1,5]<-mean(b1.isw1.20,na.rm=TRUE)
table2.isw1[1,6]<-mean(b1.isw1.sd.20,na.rm=TRUE)
table2.isw1[1,7]<-sd(b1.isw1.20,na.rm=TRUE)
table2.isw1[1,8]<-mean(cover.isw1.20[,8],na.rm=TRUE)

#### t_0=2 & c=10% ####
b0.isw1.21 <- c()
b0.isw1.sd.21 <- c()
b1.isw1.21 <- c()
b1.isw1.sd.21 <- c()
cover.isw1.21<-matrix(NA,2000,8)
colnames(cover.isw1.21)<-c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

set.seed(1)
for (i in 1:2000){
  a<-data.gen(200,c.1)
  fm <- Surv(Z, censored) ~ X
  tryCatch({
    ismb.fit <- qris(fm, data = a, t0 = 2, Q = 0.5, ne = 200, "smooth", "pmb", c(1.219403, 0.9070615))
    b0.isw1.21[i] <- ismb.fit$coefficient[1]
    b0.isw1.sd.21[i] <- ismb.fit$stderr[1]
    b1.isw1.21[i] <- ismb.fit$coefficient[2]
    b1.isw1.sd.21[i] <- ismb.fit$stderr[2]
    # Coverage
    cover.isw1.21[i,1] <- ismb.fit$coefficient[1]-1.96*ismb.fit$stderr[1]
    cover.isw1.21[i,2] <- ismb.fit$coefficient[1]
    cover.isw1.21[i,3] <- ismb.fit$coefficient[1]+1.96*ismb.fit$stderr[1]
    cover.isw1.21[i,4] <- ind(1.219403, cover.isw1.21[i,1], cover.isw1.21[i,3])
    cover.isw1.21[i,5] <- ismb.fit$coefficient[2]-1.96*ismb.fit$stderr[2]
    cover.isw1.21[i,6] <- ismb.fit$coefficient[2]
    cover.isw1.21[i,7] <- ismb.fit$coefficient[2]+1.96*ismb.fit$stderr[2]
    cover.isw1.21[i,8] <- ind(0.9070615, cover.isw1.21[i,5], cover.isw1.21[i,7])}
    , error=function(e){
      b0.isw1.21[i] <- NA
      b0.isw1.sd.21[i] <- NA
      b1.isw1.21[i] <- NA
      b1.isw1.sd.21[i] <- NA
      # Coverage
      cover.isw1.21[i,1] <- NA
      cover.isw1.21[i,2] <- NA
      cover.isw1.21[i,3] <- NA
      cover.isw1.21[i,4] <- NA
      cover.isw1.21[i,5] <- NA
      cover.isw1.21[i,6] <- NA
      cover.isw1.21[i,7] <- NA
      cover.isw1.21[i,8] <- NA
    })
}

# IS beta table
table2.isw1[2,1]<-mean(b0.isw1.21,na.rm=TRUE)
table2.isw1[2,2]<-mean(b0.isw1.sd.21,na.rm=TRUE)
table2.isw1[2,3]<-sd(b0.isw1.21,na.rm=TRUE)
table2.isw1[2,4]<-mean(cover.isw1.21[,4],na.rm=TRUE)
table2.isw1[2,5]<-mean(b1.isw1.21,na.rm=TRUE)
table2.isw1[2,6]<-mean(b1.isw1.sd.21,na.rm=TRUE)
table2.isw1[2,7]<-sd(b1.isw1.21,na.rm=TRUE)
table2.isw1[2,8]<-mean(cover.isw1.21[,8],na.rm=TRUE)

#### t_0=2 & c=30% ####
b0.isw1.23 <- c()
b0.isw1.sd.23 <- c()
b1.isw1.23 <- c()
b1.isw1.sd.23 <- c()
cover.isw1.23<-matrix(NA,2000,8)
colnames(cover.isw1.23)<-c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

set.seed(1)
for (i in 1:2000){
  a<-data.gen(200,c.3)
  fm <- Surv(Z, censored) ~ X
  tryCatch({
    ismb.fit <- qris(fm, data = a, t0 = 2, Q = 0.5, ne = 200, "smooth", "pmb", c(1.219403, 0.9070615))
    b0.isw1.23[i] <- ismb.fit$coefficient[1]
    b0.isw1.sd.23[i] <- ismb.fit$stderr[1]
    b1.isw1.23[i] <- ismb.fit$coefficient[2]
    b1.isw1.sd.23[i] <- ismb.fit$stderr[2]
    # Coverage
    cover.isw1.23[i,1] <- ismb.fit$coefficient[1]-1.96*ismb.fit$stderr[1]
    cover.isw1.23[i,2] <- ismb.fit$coefficient[1]
    cover.isw1.23[i,3] <- ismb.fit$coefficient[1]+1.96*ismb.fit$stderr[1]
    cover.isw1.23[i,4] <- ind(1.219403, cover.isw1.23[i,1], cover.isw1.23[i,3])
    cover.isw1.23[i,5] <- ismb.fit$coefficient[2]-1.96*ismb.fit$stderr[2]
    cover.isw1.23[i,6] <- ismb.fit$coefficient[2]
    cover.isw1.23[i,7] <- ismb.fit$coefficient[2]+1.96*ismb.fit$stderr[2]
    cover.isw1.23[i,8] <- ind(0.9070615, cover.isw1.23[i,5], cover.isw1.23[i,7])}
    , error=function(e){
      b0.isw1.23[i] <- NA
      b0.isw1.sd.23[i] <- NA
      b1.isw1.23[i] <- NA
      b1.isw1.sd.23[i] <- NA
      # Coverage
      cover.isw1.23[i,1] <- NA
      cover.isw1.23[i,2] <- NA
      cover.isw1.23[i,3] <- NA
      cover.isw1.23[i,4] <- NA
      cover.isw1.23[i,5] <- NA
      cover.isw1.23[i,6] <- NA
      cover.isw1.23[i,7] <- NA
      cover.isw1.23[i,8] <- NA
    })
}

# IS beta table
table2.isw1[3,1]<-mean(b0.isw1.23,na.rm=TRUE)
table2.isw1[3,2]<-mean(b0.isw1.sd.23,na.rm=TRUE)
table2.isw1[3,3]<-sd(b0.isw1.23,na.rm=TRUE)
table2.isw1[3,4]<-mean(cover.isw1.23[,4],na.rm=TRUE)
table2.isw1[3,5]<-mean(b1.isw1.23,na.rm=TRUE)
table2.isw1[3,6]<-mean(b1.isw1.sd.23,na.rm=TRUE)
table2.isw1[3,7]<-sd(b1.isw1.23,na.rm=TRUE)
table2.isw1[3,8]<-mean(cover.isw1.23[,8],na.rm=TRUE)

#### t_0=2 & c=50% ####
b0.isw1.25 <- c()
b0.isw1.sd.25 <- c()
b1.isw1.25 <- c()
b1.isw1.sd.25 <- c()
cover.isw1.25<-matrix(NA,2000,8)
colnames(cover.isw1.25)<-c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

set.seed(1)
for (i in 1:2000){
  a<-data.gen(200,c.5)
  fm <- Surv(Z, censored) ~ X
  tryCatch({
    ismb.fit <- qris(fm, data = a, t0 = 2, Q = 0.5, ne = 200, "smooth", "pmb", c(1.219403, 0.9070615))
    b0.isw1.25[i] <- ismb.fit$coefficient[1]
    b0.isw1.sd.25[i] <- ismb.fit$stderr[1]
    b1.isw1.25[i] <- ismb.fit$coefficient[2]
    b1.isw1.sd.25[i] <- ismb.fit$stderr[2]
    # Coverage
    cover.isw1.25[i,1] <- ismb.fit$coefficient[1]-1.96*ismb.fit$stderr[1]
    cover.isw1.25[i,2] <- ismb.fit$coefficient[1]
    cover.isw1.25[i,3] <- ismb.fit$coefficient[1]+1.96*ismb.fit$stderr[1]
    cover.isw1.25[i,4] <- ind(1.219403, cover.isw1.25[i,1], cover.isw1.25[i,3])
    cover.isw1.25[i,5] <- ismb.fit$coefficient[2]-1.96*ismb.fit$stderr[2]
    cover.isw1.25[i,6] <- ismb.fit$coefficient[2]
    cover.isw1.25[i,7] <- ismb.fit$coefficient[2]+1.96*ismb.fit$stderr[2]
    cover.isw1.25[i,8] <- ind(0.9070615, cover.isw1.25[i,5], cover.isw1.25[i,7])}
    , error=function(e){
      b0.isw1.25[i] <- NA
      b0.isw1.sd.25[i] <- NA
      b1.isw1.25[i] <- NA
      b1.isw1.sd.25[i] <- NA
      # Coverage
      cover.isw1.25[i,1] <- NA
      cover.isw1.25[i,2] <- NA
      cover.isw1.25[i,3] <- NA
      cover.isw1.25[i,4] <- NA
      cover.isw1.25[i,5] <- NA
      cover.isw1.25[i,6] <- NA
      cover.isw1.25[i,7] <- NA
      cover.isw1.25[i,8] <- NA
    })
}

# IS beta table
table2.isw1[4,1]<-mean(b0.isw1.25,na.rm=TRUE)
table2.isw1[4,2]<-mean(b0.isw1.sd.25,na.rm=TRUE)
table2.isw1[4,3]<-sd(b0.isw1.25,na.rm=TRUE)
table2.isw1[4,4]<-mean(cover.isw1.25[,4],na.rm=TRUE)
table2.isw1[4,5]<-mean(b1.isw1.25,na.rm=TRUE)
table2.isw1[4,6]<-mean(b1.isw1.sd.25,na.rm=TRUE)
table2.isw1[4,7]<-sd(b1.isw1.25,na.rm=TRUE)
table2.isw1[4,8]<-mean(cover.isw1.25[,8],na.rm=TRUE)

#### t_0=2 & c=70% ####
b0.isw1.27 <- c()
b0.isw1.sd.27 <- c()
b1.isw1.27 <- c()
b1.isw1.sd.27 <- c()
cover.isw1.27<-matrix(NA,2000,8)
colnames(cover.isw1.27)<-c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

set.seed(1)
for (i in 1:2000){
  a<-data.gen(200,c.7)
  fm <- Surv(Z, censored) ~ X
  tryCatch({
    ismb.fit <- qris(fm, data = a, t0 = 2, Q = 0.5, ne = 200, "smooth", "pmb", c(1.219403, 0.9070615))
    b0.isw1.27[i] <- ismb.fit$coefficient[1]
    b0.isw1.sd.27[i] <- ismb.fit$stderr[1]
    b1.isw1.27[i] <- ismb.fit$coefficient[2]
    b1.isw1.sd.27[i] <- ismb.fit$stderr[2]
    # Coverage
    cover.isw1.27[i,1] <- ismb.fit$coefficient[1]-1.96*ismb.fit$stderr[1]
    cover.isw1.27[i,2] <- ismb.fit$coefficient[1]
    cover.isw1.27[i,3] <- ismb.fit$coefficient[1]+1.96*ismb.fit$stderr[1]
    cover.isw1.27[i,4] <- ind(1.219403, cover.isw1.27[i,1], cover.isw1.27[i,3])
    cover.isw1.27[i,5] <- ismb.fit$coefficient[2]-1.96*ismb.fit$stderr[2]
    cover.isw1.27[i,6] <- ismb.fit$coefficient[2]
    cover.isw1.27[i,7] <- ismb.fit$coefficient[2]+1.96*ismb.fit$stderr[2]
    cover.isw1.27[i,8] <- ind(0.9070615, cover.isw1.27[i,5], cover.isw1.27[i,7])}
    , error=function(e){
      b0.isw1.27[i] <- NA
      b0.isw1.sd.27[i] <- NA
      b1.isw1.27[i] <- NA
      b1.isw1.sd.27[i] <- NA
      # Coverage
      cover.isw1.27[i,1] <- NA
      cover.isw1.27[i,2] <- NA
      cover.isw1.27[i,3] <- NA
      cover.isw1.27[i,4] <- NA
      cover.isw1.27[i,5] <- NA
      cover.isw1.27[i,6] <- NA
      cover.isw1.27[i,7] <- NA
      cover.isw1.27[i,8] <- NA
    })
}

# IS beta table
table2.isw1[5,1]<-mean(b0.isw1.27,na.rm=TRUE)
table2.isw1[5,2]<-mean(b0.isw1.sd.27,na.rm=TRUE)
table2.isw1[5,3]<-sd(b0.isw1.27,na.rm=TRUE)
table2.isw1[5,4]<-mean(cover.isw1.27[,4],na.rm=TRUE)
table2.isw1[5,5]<-mean(b1.isw1.27,na.rm=TRUE)
table2.isw1[5,6]<-mean(b1.isw1.sd.27,na.rm=TRUE)
table2.isw1[5,7]<-sd(b1.isw1.27,na.rm=TRUE)
table2.isw1[5,8]<-mean(cover.isw1.27[,8],na.rm=TRUE)

#### censoring point at t_0=3 ####
c.0<-5000000
c.1<-61.17
c.3<-22.53
c.5<-13.43
c.7<-8.5

#### t_0=3 & c=0% ####
b0.isw1.30 <- c()
b0.isw1.sd.30 <- c()
b1.isw1.30 <- c()
b1.isw1.sd.30 <- c()
cover.isw1.30<-matrix(NA,2000,8)
colnames(cover.isw1.30)<-c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

set.seed(1)
for (i in 1:2000){
  a <- data.gen(200,c.0)
  fm <- Surv(Z, censored) ~ X
  tryCatch({
    ismb.fit <- qris(fm, data = a, t0 = 3, Q = 0.5, ne = 200, "smooth", "pmb", c(1.040613, 1.0174711))
    b0.isw1.30[i] <- ismb.fit$coefficient[1]
    b0.isw1.sd.30[i] <- ismb.fit$stderr[1]
    b1.isw1.30[i] <- ismb.fit$coefficient[2]
    b1.isw1.sd.30[i] <- ismb.fit$stderr[2]
    # Coverage
    cover.isw1.30[i,1] <- ismb.fit$coefficient[1]-1.96*ismb.fit$stderr[1]
    cover.isw1.30[i,2] <- ismb.fit$coefficient[1]
    cover.isw1.30[i,3] <- ismb.fit$coefficient[1]+1.96*ismb.fit$stderr[1]
    cover.isw1.30[i,4] <- ind(1.040613, cover.isw1.30[i,1], cover.isw1.30[i,3])
    cover.isw1.30[i,5] <- ismb.fit$coefficient[2]-1.96*ismb.fit$stderr[2]
    cover.isw1.30[i,6] <- ismb.fit$coefficient[2]
    cover.isw1.30[i,7] <- ismb.fit$coefficient[2]+1.96*ismb.fit$stderr[2]
    cover.isw1.30[i,8] <- ind(1.0174711, cover.isw1.30[i,5], cover.isw1.30[i,7])}
    , error=function(e){
      b0.isw1.30[i] <- NA
      b0.isw1.sd.30[i] <- NA
      b1.isw1.30[i] <- NA
      b1.isw1.sd.30[i] <- NA
      # Coverage
      cover.isw1.30[i,1] <- NA
      cover.isw1.30[i,2] <- NA
      cover.isw1.30[i,3] <- NA
      cover.isw1.30[i,4] <- NA
      cover.isw1.30[i,5] <- NA
      cover.isw1.30[i,6] <- NA
      cover.isw1.30[i,7] <- NA
      cover.isw1.30[i,8] <- NA
    })
}

# IS beta table
table3.isw1[1,1]<-mean(b0.isw1.30,na.rm=TRUE)
table3.isw1[1,2]<-mean(b0.isw1.sd.30,na.rm=TRUE)
table3.isw1[1,3]<-sd(b0.isw1.30,na.rm=TRUE)
table3.isw1[1,4]<-mean(cover.isw1.30[,4],na.rm=TRUE)
table3.isw1[1,5]<-mean(b1.isw1.30,na.rm=TRUE)
table3.isw1[1,6]<-mean(b1.isw1.sd.30,na.rm=TRUE)
table3.isw1[1,7]<-sd(b1.isw1.30,na.rm=TRUE)
table3.isw1[1,8]<-mean(cover.isw1.30[,8],na.rm=TRUE)

#### t_0=3 & c=10% ####
b0.isw1.31 <- c()
b0.isw1.sd.31 <- c()
b1.isw1.31 <- c()
b1.isw1.sd.31 <- c()
cover.isw1.31<-matrix(NA,2000,8)
colnames(cover.isw1.31)<-c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

set.seed(1)
for (i in 1:2000){
  a<-data.gen(200,c.1)
  fm <- Surv(Z, censored) ~ X
  tryCatch({
    ismb.fit <- qris(fm, data = a, t0 = 3, Q = 0.5, ne = 200, "smooth", "pmb", c(1.040613, 1.0174711))
    b0.isw1.31[i] <- ismb.fit$coefficient[1]
    b0.isw1.sd.31[i] <- ismb.fit$stderr[1]
    b1.isw1.31[i] <- ismb.fit$coefficient[2]
    b1.isw1.sd.31[i] <- ismb.fit$stderr[2]
    # Coverage
    cover.isw1.31[i,1] <- ismb.fit$coefficient[1]-1.96*ismb.fit$stderr[1]
    cover.isw1.31[i,2] <- ismb.fit$coefficient[1]
    cover.isw1.31[i,3] <- ismb.fit$coefficient[1]+1.96*ismb.fit$stderr[1]
    cover.isw1.31[i,4] <- ind(1.040613, cover.isw1.31[i,1], cover.isw1.31[i,3])
    cover.isw1.31[i,5] <- ismb.fit$coefficient[2]-1.96*ismb.fit$stderr[2]
    cover.isw1.31[i,6] <- ismb.fit$coefficient[2]
    cover.isw1.31[i,7] <- ismb.fit$coefficient[2]+1.96*ismb.fit$stderr[2]
    cover.isw1.31[i,8] <- ind(1.0174711, cover.isw1.31[i,5], cover.isw1.31[i,7])}
    , error=function(e){
      b0.isw1.31[i] <- NA
      b0.isw1.sd.31[i] <- NA
      b1.isw1.31[i] <- NA
      b1.isw1.sd.31[i] <- NA
      # Coverage
      cover.isw1.31[i,1] <- NA
      cover.isw1.31[i,2] <- NA
      cover.isw1.31[i,3] <- NA
      cover.isw1.31[i,4] <- NA
      cover.isw1.31[i,5] <- NA
      cover.isw1.31[i,6] <- NA
      cover.isw1.31[i,7] <- NA
      cover.isw1.31[i,8] <- NA
    })
}

# IS beta table
table3.isw1[2,1]<-mean(b0.isw1.31,na.rm=TRUE)
table3.isw1[2,2]<-mean(b0.isw1.sd.31,na.rm=TRUE)
table3.isw1[2,3]<-sd(b0.isw1.31,na.rm=TRUE)
table3.isw1[2,4]<-mean(cover.isw1.31[,4],na.rm=TRUE)
table3.isw1[2,5]<-mean(b1.isw1.31,na.rm=TRUE)
table3.isw1[2,6]<-mean(b1.isw1.sd.31,na.rm=TRUE)
table3.isw1[2,7]<-sd(b1.isw1.31,na.rm=TRUE)
table3.isw1[2,8]<-mean(cover.isw1.31[,8],na.rm=TRUE)

#### t_0=3 & c=30% ####
b0.isw1.33 <- c()
b0.isw1.sd.33 <- c()
b1.isw1.33 <- c()
b1.isw1.sd.33 <- c()
cover.isw1.33<-matrix(NA,2000,8)
colnames(cover.isw1.33)<-c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

set.seed(1)
for (i in 1:2000){
  a<-data.gen(200,c.3)
  fm <- Surv(Z, censored) ~ X
  tryCatch({
    ismb.fit <- qris(fm, data = a, t0 = 3, Q = 0.5, ne = 200, "smooth", "pmb", c(1.040613, 1.0174711))
    b0.isw1.33[i] <- ismb.fit$coefficient[1]
    b0.isw1.sd.33[i] <- ismb.fit$stderr[1]
    b1.isw1.33[i] <- ismb.fit$coefficient[2]
    b1.isw1.sd.33[i] <- ismb.fit$stderr[2]
    # Coverage
    cover.isw1.33[i,1] <- ismb.fit$coefficient[1]-1.96*ismb.fit$stderr[1]
    cover.isw1.33[i,2] <- ismb.fit$coefficient[1]
    cover.isw1.33[i,3] <- ismb.fit$coefficient[1]+1.96*ismb.fit$stderr[1]
    cover.isw1.33[i,4] <- ind(1.040613, cover.isw1.33[i,1], cover.isw1.33[i,3])
    cover.isw1.33[i,5] <- ismb.fit$coefficient[2]-1.96*ismb.fit$stderr[2]
    cover.isw1.33[i,6] <- ismb.fit$coefficient[2]
    cover.isw1.33[i,7] <- ismb.fit$coefficient[2]+1.96*ismb.fit$stderr[2]
    cover.isw1.33[i,8] <- ind(1.0174711, cover.isw1.33[i,5], cover.isw1.33[i,7])}
    , error=function(e){
      b0.isw1.33[i] <- NA
      b0.isw1.sd.33[i] <- NA
      b1.isw1.33[i] <- NA
      b1.isw1.sd.33[i] <- NA
      # Coverage
      cover.isw1.33[i,1] <- NA
      cover.isw1.33[i,2] <- NA
      cover.isw1.33[i,3] <- NA
      cover.isw1.33[i,4] <- NA
      cover.isw1.33[i,5] <- NA
      cover.isw1.33[i,6] <- NA
      cover.isw1.33[i,7] <- NA
      cover.isw1.33[i,8] <- NA
    })
}

# IS beta table
table3.isw1[3,1]<-mean(b0.isw1.33,na.rm=TRUE)
table3.isw1[3,2]<-mean(b0.isw1.sd.33,na.rm=TRUE)
table3.isw1[3,3]<-sd(b0.isw1.33,na.rm=TRUE)
table3.isw1[3,4]<-mean(cover.isw1.33[,4],na.rm=TRUE)
table3.isw1[3,5]<-mean(b1.isw1.33,na.rm=TRUE)
table3.isw1[3,6]<-mean(b1.isw1.sd.33,na.rm=TRUE)
table3.isw1[3,7]<-sd(b1.isw1.33,na.rm=TRUE)
table3.isw1[3,8]<-mean(cover.isw1.33[,8],na.rm=TRUE)

#### t_0=3 & c=50% ####
b0.isw1.35 <- c()
b0.isw1.sd.35 <- c()
b1.isw1.35 <- c()
b1.isw1.sd.35 <- c()
cover.isw1.35<-matrix(NA,2000,8)
colnames(cover.isw1.35)<-c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")
set.seed(1)
for (i in 1:2000){
  a<-data.gen(200,c.5)
  fm <- Surv(Z, censored) ~ X
  tryCatch({
    ismb.fit <- qris(fm, data = a, t0 = 3, Q = 0.5, ne = 200, "smooth", "pmb", c(1.040613, 1.0174711))
    b0.isw1.35[i] <- ismb.fit$coefficient[1]
    b0.isw1.sd.35[i] <- ismb.fit$stderr[1]
    b1.isw1.35[i] <- ismb.fit$coefficient[2]
    b1.isw1.sd.35[i] <- ismb.fit$stderr[2]
    # Coverage
    cover.isw1.35[i,1] <- ismb.fit$coefficient[1]-1.96*ismb.fit$stderr[1]
    cover.isw1.35[i,2] <- ismb.fit$coefficient[1]
    cover.isw1.35[i,3] <- ismb.fit$coefficient[1]+1.96*ismb.fit$stderr[1]
    cover.isw1.35[i,4] <- ind(1.040613, cover.isw1.35[i,1], cover.isw1.35[i,3])
    cover.isw1.35[i,5] <- ismb.fit$coefficient[2]-1.96*ismb.fit$stderr[2]
    cover.isw1.35[i,6] <- ismb.fit$coefficient[2]
    cover.isw1.35[i,7] <- ismb.fit$coefficient[2]+1.96*ismb.fit$stderr[2]
    cover.isw1.35[i,8] <- ind(1.0174711, cover.isw1.35[i,5], cover.isw1.35[i,7])}
    , error=function(e){
      b0.isw1.35[i] <- NA
      b0.isw1.sd.35[i] <- NA
      b1.isw1.35[i] <- NA
      b1.isw1.sd.35[i] <- NA
      # Coverage
      cover.isw1.35[i,1] <- NA
      cover.isw1.35[i,2] <- NA
      cover.isw1.35[i,3] <- NA
      cover.isw1.35[i,4] <- NA
      cover.isw1.35[i,5] <- NA
      cover.isw1.35[i,6] <- NA
      cover.isw1.35[i,7] <- NA
      cover.isw1.35[i,8] <- NA
    })
}

# IS beta table
table3.isw1[4,1]<-mean(b0.isw1.35,na.rm=TRUE)
table3.isw1[4,2]<-mean(b0.isw1.sd.35,na.rm=TRUE)
table3.isw1[4,3]<-sd(b0.isw1.35,na.rm=TRUE)
table3.isw1[4,4]<-mean(cover.isw1.35[,4],na.rm=TRUE)
table3.isw1[4,5]<-mean(b1.isw1.35,na.rm=TRUE)
table3.isw1[4,6]<-mean(b1.isw1.sd.35,na.rm=TRUE)
table3.isw1[4,7]<-sd(b1.isw1.35,na.rm=TRUE)
table3.isw1[4,8]<-mean(cover.isw1.35[,8],na.rm=TRUE)

#### t_0=3 & c=70% ####
b0.isw1.37 <- c()
b0.isw1.sd.37 <- c()
b1.isw1.37 <- c()
b1.isw1.sd.37 <- c()
cover.isw1.37<-matrix(NA,2000,8)
colnames(cover.isw1.37)<-c("lower.b0","b0","upper.b0","I.b0","lower.b1","b1","upper.b1","I.b1")

set.seed(1)
for (i in 1:2000){
  a<-data.gen(200,c.7)
  fm <- Surv(Z, censored) ~ X
  tryCatch({
    ismb.fit <- qris(fm, data = a, t0 = 3, Q = 0.5, ne = 200, "smooth", "pmb", c(1.040613, 1.0174711))
    b0.isw1.37[i] <- ismb.fit$coefficient[1]
    b0.isw1.sd.37[i] <- ismb.fit$stderr[1]
    b1.isw1.37[i] <- ismb.fit$coefficient[2]
    b1.isw1.sd.37[i] <- ismb.fit$stderr[2]
    # Coverage
    cover.isw1.37[i,1] <- ismb.fit$coefficient[1]-1.96*ismb.fit$stderr[1]
    cover.isw1.37[i,2] <- ismb.fit$coefficient[1]
    cover.isw1.37[i,3] <- ismb.fit$coefficient[1]+1.96*ismb.fit$stderr[1]
    cover.isw1.37[i,4] <- ind(1.040613, cover.isw1.37[i,1], cover.isw1.37[i,3])
    cover.isw1.37[i,5] <- ismb.fit$coefficient[2]-1.96*ismb.fit$stderr[2]
    cover.isw1.37[i,6] <- ismb.fit$coefficient[2]
    cover.isw1.37[i,7] <- ismb.fit$coefficient[2]+1.96*ismb.fit$stderr[2]
    cover.isw1.37[i,8] <- ind(1.0174711, cover.isw1.37[i,5], cover.isw1.37[i,7])}
    , error=function(e){
      b0.isw1.37[i] <- NA
      b0.isw1.sd.37[i] <- NA
      b1.isw1.37[i] <- NA
      b1.isw1.sd.37[i] <- NA
      # Coverage
      cover.isw1.37[i,1] <- NA
      cover.isw1.37[i,2] <- NA
      cover.isw1.37[i,3] <- NA
      cover.isw1.37[i,4] <- NA
      cover.isw1.37[i,5] <- NA
      cover.isw1.37[i,6] <- NA
      cover.isw1.37[i,7] <- NA
      cover.isw1.37[i,8] <- NA
    })
}

# IS beta table
table3.isw1[5,1]<-mean(b0.isw1.37,na.rm=TRUE)
table3.isw1[5,2]<-mean(b0.isw1.sd.37,na.rm=TRUE)
table3.isw1[5,3]<-sd(b0.isw1.37,na.rm=TRUE)
table3.isw1[5,4]<-mean(cover.isw1.37[,4],na.rm=TRUE)
table3.isw1[5,5]<-mean(b1.isw1.37,na.rm=TRUE)
table3.isw1[5,6]<-mean(b1.isw1.sd.37,na.rm=TRUE)
table3.isw1[5,7]<-sd(b1.isw1.37,na.rm=TRUE)
table3.isw1[5,8]<-mean(cover.isw1.37[,8],na.rm=TRUE)
