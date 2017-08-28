###
###1. Store True values for X0, X1, Z0, Z1, Z2, P
###2. Store proportion censored pre & post
###3. Simmulate n of 1000, do this 1000 times
###4. Estimate cox, weibull, store all relevant coefficient estimates (exponentiate p's where applicable)
###5. For each value in 4, calculate CPs and RMSEs, store.
###6. Estimate cure exp and zombie weibull-->store all relevant coefficient estimates (exponentiate p's where applicable).
###7. For each value in 6, calculate CPs and RMSEs, store.

##############
####Set Up####
##############
  
  #clear memory
  rm( list=ls() )
  
  #load necessary libraries 						                                 
  library(foreign)
  library(car)
  library(MASS)
  library(VGAM)
  library(survival)
  library(msm)
  library(verification)
  library(corpcor)
  library(Design)
  library(SpatialCure)
  library(MCMCpack)
  #set working directory
  setwd("/Users/bomin8319/Desktop/SpatialCureSurv/simulation")
  
  ##########################################################################
  ##########################################################################
  ############################Monte Carlo###################################
  ##########################################################################
  
  #set seed
  set.seed(123)   
  
  #set the number of observations
  n<-1000
  
  #set the number of simulations, and create matrices to store the results
  nsims<-3
  
  #history matrix for true estimates
  tru.est<-matrix(NA,nrow=nsims,ncol=8)
  #history matrix for true Spatial estimates
  tru.est2<-matrix(NA,nrow=nsims,ncol=6)
  
  #history matrix for exp estimates
  exp.est<-matrix(NA,nrow=nsims,ncol=36)
  #history matrix for weibull estimates
  weib.est<-matrix(NA,nrow=nsims,ncol=42)
  #history matrix for exp RMSE
  exp.rmse<-matrix(NA,nrow=nsims,ncol=18)
  #history matrix for exp RMSE
  weib.rmse<-matrix(NA,nrow=nsims,ncol=21)
  #history matrix for exp CP
  exp.cp<-matrix(NA,nrow=nsims,ncol=18)
  #history matrix for exp CP
  weib.cp<-matrix(NA,nrow=nsims,ncol=21)
  
  #create covariates
  x<-runif(n, min=-2.5, max=12)
  z<-log(runif(n, min=1, max=100))
  s<- sample(1:5, n, replace = TRUE)
  a<- matrix(0, 5, 5)
  a[1,2]=a[2,1] = 1
  a[2,3]=a[3,2] = 1 
  a[5,3]=a[3,5] = 1 
  a[1,4]=a[4,1] = 1 
  w <- rep(0,5)
  w <- vapply(1:5, function(i){rnorm(1, mean(w[which(a[i,]==1)]), sqrt(1/(1*rowSums(a)[i])))}, c(1))
  w <- w - mean(w)
  #create a dependent variable, begin the simmulations
  for(i in 1:nsims){
  
  #Assign parameter values
  tru.est[i,1]<-1
  tru.est[i,2]<-3.5
  tru.est[i,3]<--2
  tru.est[i,4]<--2
  tru.est[i,5]<-3
  tru.est[i,6]<-1
  
  tru.est2[i,1:5] <- w
  tru.est2[i,6]<-1   #lambda
  
  W <- w[s]
  myrates <- exp(tru.est[i,1]+(tru.est[i,2]*x + W)) 
  y <- rexp(n, rate = myrates) # generates the r.v.
  cen <- rexp(n, rate = 1 )
  ycen <- pmin(y, cen)
  di <- as.numeric(y <= cen)
  tru.est[i,7]<-table(di)[1]
  
  
  #create parameters for ZG
  phi<-1/(1+exp(-(tru.est[i,3]+tru.est[i,4]*z+tru.est[i,5]*x)))
  yzero<-matrix(0,n,1)
  error<--1*rlogis(n)
  flag<-error>qlogis(phi)
  yzero[flag]<-error[flag]
  flag<-yzero==0
  ycen[flag]<-ifelse(di[flag]==1,cen[flag],ycen[flag])
  di[flag]<-ifelse(di[flag]==1,yzero[flag],di[flag])
  tru.est[i,8]<-table(di)[1]
  data<-cbind(ycen,di,x,z, s)
  
  
  #############################################################################
  ########################Simple Exponential Model#############################
  #############################################################################
  Y<-ycen
  C<-di
  X<-cbind(1,x)
  #X <- matrix(x, ncol = 1)
  Exponential <- mcmcSurv(Y, C, X, 550, 50, 1, form = "Exponential")
  output.Exponential = summary(mcmc(Exponential$betas))
  
  #store betas and ses
  exp.est[i,1]<-output.Exponential[[1]][1,1]
  exp.est[i,2]<-output.Exponential[[1]][1,4]
  exp.est[i,3]<-output.Exponential[[1]][2,1]
  exp.est[i,4]<-output.Exponential[[1]][2,4]
  
  #store rmse
  exp.rmse[i,1]<-sqrt((tru.est[i,1]-exp.est[i,1])^2)
  exp.rmse[i,2]<-sqrt((tru.est[i,2]-exp.est[i,3])^2)
  
  #calculate upper and lower 95% CI's
  b0.lower<-output.Exponential[[2]][1,1]
  b0.upper<-output.Exponential[[2]][1,5]
  b1.lower<-output.Exponential[[2]][2,1]
  b1.upper<-output.Exponential[[2]][2,5]
  
  
  #store coverage parameters
  exp.cp[i,1]<-ifelse(tru.est[i,1]>b0.lower & tru.est[i,1]<b0.upper, 1,0)
  exp.cp[i,2]<-ifelse(tru.est[i,2]>b1.lower & tru.est[i,2]<b1.upper, 1,0)
  
  
  ################################################################################
  #########################Simple Weibull Model ##################################
  ################################################################################
  Y<-ycen
  C<-di
  X<-cbind(1,x)
  #X <- matrix(x, ncol = 1)
  Weibull <- mcmcSurv(Y, C, X, 550, 50, 1, form = "Weibull")
  output.Weibull = summary(mcmc(Weibull$betas))
  output.Weibull2 = summary(mcmc(Weibull$rho))
  
  #store betas and ses
  weib.est[i,1]<-output.Weibull[[1]][1,1]
  weib.est[i,2]<-output.Weibull[[1]][1,4]
  weib.est[i,3]<-output.Weibull[[1]][2,1]
  weib.est[i,4]<-output.Weibull[[1]][2,4]
  weib.est[i,5]<-output.Weibull2[[1]][1]
  weib.est[i,6]<-output.Weibull2[[1]][4]
  
  #store rmse
  weib.rmse[i,1]<-sqrt((tru.est[i,1]-weib.est[i,1])^2)
  weib.rmse[i,2]<-sqrt((tru.est[i,2]-weib.est[i,3])^2)
  weib.rmse[i,3]<-sqrt((tru.est[i,6]-weib.est[i,5])^2)
  
  #calculate upper and lower 95% CI's
  b0.lower<-output.Weibull[[2]][1,1]
  b0.upper<-output.Weibull[[2]][1,5]
  b1.lower<-output.Weibull[[2]][2,1]
  b1.upper<-output.Weibull[[2]][2,5]
  p.lower<-output.Weibull2[[2]][1]
  p.upper<-output.Weibull2[[2]][5]
  
  #store coverage parameters
  weib.cp[i,1]<-ifelse(tru.est[i,1]>b0.lower & tru.est[i,1]<b0.upper, 1,0)
  weib.cp[i,2]<-ifelse(tru.est[i,2]>b1.lower & tru.est[i,2]<b1.upper, 1,0)
  weib.cp[i,3]<-ifelse(tru.est[i,6]>p.lower & tru.est[i,6]<p.upper, 1,0)
  
  ###############################################################################
  ############################Cure Exponential Model#############################
  ###############################################################################
  Y<-ycen
  C<-di
  X<-cbind(1,x)
  #X <- matrix(x, ncol = 1)
  Z<-cbind(1,z,x)
  
  CExponential <- mcmcCure(Y, C, X, Z, 3000, 500, 5, form = "Exponential")
  output.CExponential = summary(mcmc(CExponential$betas))
  output.CExponential2 = summary(mcmc(CExponential$gammas))
  
  #store betas and ses
  exp.est[i,5]<-output.CExponential[[1]][1,1]
  exp.est[i,6]<-output.CExponential[[1]][1,4]
  exp.est[i,7]<-output.CExponential[[1]][2,1]
  exp.est[i,8]<-output.CExponential[[1]][2,4]
  exp.est[i,9]<-output.CExponential2[[1]][1,1]
  exp.est[i,10]<-output.CExponential2[[1]][1,4]
  exp.est[i,11]<-output.CExponential2[[1]][2,1]
  exp.est[i,12]<-output.CExponential2[[1]][2,4]
  exp.est[i,13]<-output.CExponential2[[1]][3,1]
  exp.est[i,14]<-output.CExponential2[[1]][3,4]
  
  #store rmse
  exp.rmse[i,3]<-sqrt((tru.est[i,1]-exp.est[i,5])^2)
  exp.rmse[i,4]<-sqrt((tru.est[i,2]-exp.est[i,7])^2)
  exp.rmse[i,5]<-sqrt((tru.est[i,3]-exp.est[i,9])^2)
  exp.rmse[i,6]<-sqrt((tru.est[i,4]-exp.est[i,11])^2)
  exp.rmse[i,7]<-sqrt((tru.est[i,5]-exp.est[i,13])^2)
  
  #calculate upper and lower 95% CI's
  b0.lower<-output.CExponential[[2]][1,1]
  b0.upper<-output.CExponential[[2]][1,5]
  b1.lower<-output.CExponential[[2]][2,1]
  b1.upper<-output.CExponential[[2]][2,5]
  g0.lower<-output.CExponential2[[2]][1,1]
  g0.upper<-output.CExponential2[[2]][1,5]
  g1.lower<-output.CExponential2[[2]][2,1]
  g1.upper<-output.CExponential2[[2]][2,5]
  g2.lower<-output.CExponential2[[2]][3,1]
  g2.upper<-output.CExponential2[[2]][3,5]
  
  #store coverage parameters
  exp.cp[i,3]<-ifelse(tru.est[i,1]>b0.lower & tru.est[i,1]<b0.upper, 1,0)
  exp.cp[i,4]<-ifelse(tru.est[i,2]>b1.lower & tru.est[i,2]<b1.upper, 1,0)
  exp.cp[i,5]<-ifelse(tru.est[i,3]>g0.lower & tru.est[i,3]<g0.upper, 1,0)
  exp.cp[i,6]<-ifelse(tru.est[i,4]>g1.lower & tru.est[i,4]<g1.upper, 1,0)
  exp.cp[i,7]<-ifelse(tru.est[i,5]>g2.lower & tru.est[i,5]<g2.upper, 1,0)
  
  
  #####################################################################################
  ############################Cure Weibull Model ######################################
  #####################################################################################
  Y<-ycen
  C<-di
  X<-cbind(1,x)
  #X <- matrix(x, ncol = 1)
  Z<-cbind(1,z,x)
  
  CWeibull <- mcmcCure(Y, C, X, Z, 3000, 500, 5, form = "Weibull")
  output.CWeibull = summary(mcmc(CWeibull$betas))
  output.CWeibull2 = summary(mcmc(CWeibull$gammas))
  output.CWeibull3 = summary(mcmc(CWeibull$rho))
  
  #store betas and ses
  weib.est[i,7]<-output.CWeibull[[1]][1,1]
  weib.est[i,8]<-output.CWeibull[[1]][1,4]
  weib.est[i,9]<-output.CWeibull[[1]][2,1]
  weib.est[i,10]<-output.CWeibull[[1]][2,4]
  weib.est[i,11]<-output.CWeibull2[[1]][1,1]
  weib.est[i,12]<-output.CWeibull2[[1]][1,4]
  weib.est[i,13]<-output.CWeibull2[[1]][2,1]
  weib.est[i,14]<-output.CWeibull2[[1]][2,4]
  weib.est[i,15]<-output.CWeibull2[[1]][3,1]
  weib.est[i,16]<-output.CWeibull2[[1]][3,4]
  weib.est[i,17]<-output.CWeibull3[[1]][1]
  weib.est[i,18]<-output.CWeibull3[[1]][4]
   
  #store rmse
  weib.rmse[i,4]<-sqrt((tru.est[i,1]-weib.est[i,7])^2)
  weib.rmse[i,5]<-sqrt((tru.est[i,2]-weib.est[i,9])^2)
  weib.rmse[i,6]<-sqrt((tru.est[i,3]-weib.est[i,11])^2)
  weib.rmse[i,7]<-sqrt((tru.est[i,4]-weib.est[i,13])^2)
  weib.rmse[i,8]<-sqrt((tru.est[i,5]-weib.est[i,15])^2)
  weib.rmse[i,9]<-sqrt((tru.est[i,6]-weib.est[i,17])^2)
  
  #calculate upper and lower 95% CI's
  b0.lower<-output.CWeibull[[2]][1,1]
  b0.upper<-output.CWeibull[[2]][1,5]
  b1.lower<-output.CWeibull[[2]][2,1]
  b1.upper<-output.CWeibull[[2]][2,5]
  g0.lower<-output.CWeibull2[[2]][1,1]
  g0.upper<-output.CWeibull2[[2]][1,5]
  g1.lower<-output.CWeibull2[[2]][2,1]
  g1.upper<-output.CWeibull2[[2]][2,5]
  g2.lower<-output.CWeibull2[[2]][3,1]
  g2.upper<-output.CWeibull2[[2]][3,5]
  p.lower<-output.CWeibull3[[2]][1]
  p.upper<-output.CWeibull3[[2]][5]
  
  #store coverage parameters
  weib.cp[i,4]<-ifelse(tru.est[i,1]>b0.lower & tru.est[i,1]<b0.upper, 1,0)
  weib.cp[i,5]<-ifelse(tru.est[i,2]>b1.lower & tru.est[i,2]<b1.upper, 1,0)
  weib.cp[i,6]<-ifelse(tru.est[i,3]>g0.lower & tru.est[i,3]<g0.upper, 1,0)
  weib.cp[i,7]<-ifelse(tru.est[i,4]>g1.lower & tru.est[i,4]<g1.upper, 1,0)
  weib.cp[i,8]<-ifelse(tru.est[i,5]>g2.lower & tru.est[i,5]<g2.upper, 1,0)
  weib.cp[i,9]<-ifelse(tru.est[i,6]>p.lower & tru.est[i,6]<p.upper, 1,0)
  
  
  ###############################################################################
  ########################Spatial Cure Exponential Model#########################
  ###############################################################################
  Y<-ycen
  C<-di
  X<-cbind(1,x)
  #X <- matrix(x, ncol = 1)
  Z<-cbind(1, z,x)
  S<-s
  A<-a
  
  SCExponential <- mcmcSpatialCure2(Y, C, X, Z, S, A, 3000, 500, 5, form = "Exponential", prop.var = 0.05)
  output.SCExponential = summary(mcmc(SCExponential$betas))
  output.SCExponential2 = summary(mcmc(SCExponential$gammas))
  output.SCExponential3 = summary(mcmc(SCExponential$W))
  output.SCExponential4 = summary(mcmc(SCExponential$lambda))
  
  #store betas and ses
  exp.est[i,15]<-output.SCExponential[[1]][1,1]
  exp.est[i,16]<-output.SCExponential[[1]][1,4]
  exp.est[i,17]<-output.SCExponential[[1]][2,1]
  exp.est[i,18]<-output.SCExponential[[1]][2,4]
  exp.est[i,19]<-output.SCExponential2[[1]][1,1]
  exp.est[i,20]<-output.SCExponential2[[1]][1,4]
  exp.est[i,21]<-output.SCExponential2[[1]][2,1]
  exp.est[i,22]<-output.SCExponential2[[1]][2,4]
  exp.est[i,23]<-output.SCExponential2[[1]][3,1]
  exp.est[i,24]<-output.SCExponential2[[1]][3,4]
  
  #W's and V's
  exp.est[i,25]<-output.SCExponential3[[1]][1,1]
  exp.est[i,26]<-output.SCExponential3[[1]][1,4]
  exp.est[i,27]<-output.SCExponential3[[1]][2,1]
  exp.est[i,28]<-output.SCExponential3[[1]][2,4]
  exp.est[i,29]<-output.SCExponential3[[1]][3,1]
  exp.est[i,30]<-output.SCExponential3[[1]][3,4]
  exp.est[i,31]<-output.SCExponential3[[1]][4,1]
  exp.est[i,32]<-output.SCExponential3[[1]][4,4]
  exp.est[i,33]<-output.SCExponential3[[1]][5,1]
  exp.est[i,34]<-output.SCExponential3[[1]][5,4]

  exp.est[i,35]<-output.SCExponential4[[1]][1]
  exp.est[i,36]<-output.SCExponential4[[1]][4]
  
  
  #store rmse
  exp.rmse[i,8]<-sqrt((tru.est[i,1]-exp.est[i,15])^2)
  exp.rmse[i,9]<-sqrt((tru.est[i,2]-exp.est[i,17])^2)
  exp.rmse[i,10]<-sqrt((tru.est[i,3]-exp.est[i,19])^2)
  exp.rmse[i,11]<-sqrt((tru.est[i,4]-exp.est[i,21])^2)
  exp.rmse[i,12]<-sqrt((tru.est[i,5]-exp.est[i,23])^2)
  exp.rmse[i,13]<-sqrt((tru.est2[i,1]-exp.est[i,25])^2)
  exp.rmse[i,14]<-sqrt((tru.est2[i,2]-exp.est[i,27])^2)
  exp.rmse[i,15]<-sqrt((tru.est2[i,3]-exp.est[i,29])^2)
  exp.rmse[i,16]<-sqrt((tru.est2[i,4]-exp.est[i,31])^2)
  exp.rmse[i,17]<-sqrt((tru.est2[i,5]-exp.est[i,33])^2)
  exp.rmse[i,18]<-sqrt((tru.est2[i,6]-exp.est[i,45])^2)
  
  #calculate upper and lower 95% CI's
  b0.lower<-output.SCExponential[[2]][1,1]
  b0.upper<-output.SCExponential[[2]][1,5]
  b1.lower<-output.SCExponential[[2]][2,1]
  b1.upper<-output.SCExponential[[2]][2,5]
  g0.lower<-output.SCExponential2[[2]][1,1]
  g0.upper<-output.SCExponential2[[2]][1,5]
  g1.lower<-output.SCExponential2[[2]][2,1]
  g1.upper<-output.SCExponential2[[2]][2,5]
  g2.lower<-output.SCExponential2[[2]][3,1]
  g2.upper<-output.SCExponential2[[2]][3,5]
  
  w1.lower<-output.SCExponential3[[2]][1,1]
  w1.upper<-output.SCExponential3[[2]][1,5]
  w2.lower<-output.SCExponential3[[2]][2,1]
  w2.upper<-output.SCExponential3[[2]][2,5]
  w3.lower<-output.SCExponential3[[2]][3,1]
  w3.upper<-output.SCExponential3[[2]][3,5]
  w4.lower<-output.SCExponential3[[2]][4,1]
  w4.upper<-output.SCExponential3[[2]][4,5]
  w5.lower<-output.SCExponential3[[2]][5,1]
  w5.upper<-output.SCExponential3[[2]][5,5]
  l1.lower<-output.SCExponential5[[2]][1]
  l1.upper<-output.SCExponential5[[2]][5]
  #store coverage parameters
  exp.cp[i,8]<-ifelse(tru.est[i,1]>b0.lower & tru.est[i,1]<b0.upper, 1,0)
  exp.cp[i,9]<-ifelse(tru.est[i,2]>b1.lower & tru.est[i,2]<b1.upper, 1,0)
  exp.cp[i,10]<-ifelse(tru.est[i,3]>g0.lower & tru.est[i,3]<g0.upper, 1,0)
  exp.cp[i,11]<-ifelse(tru.est[i,4]>g1.lower & tru.est[i,4]<g1.upper, 1,0)
  exp.cp[i,12]<-ifelse(tru.est[i,5]>g2.lower & tru.est[i,5]<g2.upper, 1,0)
  
  exp.cp[i,13]<-ifelse(tru.est2[i,1]>w1.lower & tru.est2[i,1]<w1.upper, 1,0)
  exp.cp[i,14]<-ifelse(tru.est2[i,2]>w2.lower & tru.est2[i,2]<w2.upper, 1,0)
  exp.cp[i,15]<-ifelse(tru.est2[i,3]>w3.lower & tru.est2[i,3]<w3.upper, 1,0)
  exp.cp[i,16]<-ifelse(tru.est2[i,4]>w4.lower & tru.est2[i,4]<w4.upper, 1,0)
  exp.cp[i,17]<-ifelse(tru.est2[i,5]>w5.lower & tru.est2[i,5]<w5.upper, 1,0)
  exp.cp[i,18]<-ifelse(tru.est2[i,6]>l1.lower & tru.est2[i,6]<l1.upper, 1,0)
  
  
  #####################################################################################
  #######################Spatial Cure Weibull Model ###################################
  #####################################################################################
  Y<-ycen
  C<-di
  X<-cbind(1,x)
  #X <- matrix(x, ncol = 1)
  Z<-cbind(1,z,x)
  S<-s
  A<-a
  
  SCWeibull <- mcmcSpatialCure2(Y, C, X, Z, S, A, 3000, 500, 5, form = "Weibull", prop.var = 0.05)
  output.SCWeibull = summary(mcmc(SCWeibull$betas))
  output.SCWeibull2 = summary(mcmc(SCWeibull$gammas))
  output.SCWeibull3 = summary(mcmc(SCWeibull$rho))
  output.SCWeibull4 = summary(mcmc(SCWeibull$W))
  output.SCWeibull5 = summary(mcmc(SCWeibull$lambda))
  
  #store betas and ses
  weib.est[i,19]<-output.SCWeibull[[1]][1,1]
  weib.est[i,20]<-output.SCWeibull[[1]][1,4]
  weib.est[i,21]<-output.SCWeibull[[1]][2,1]
  weib.est[i,22]<-output.SCWeibull[[1]][2,4]
  weib.est[i,23]<-output.SCWeibull2[[1]][1,1]
  weib.est[i,24]<-output.SCWeibull2[[1]][1,4]
  weib.est[i,25]<-output.SCWeibull2[[1]][2,1]
  weib.est[i,26]<-output.SCWeibull2[[1]][2,4]
  weib.est[i,27]<-output.SCWeibull2[[1]][3,1]
  weib.est[i,28]<-output.SCWeibull2[[1]][3,4]
  weib.est[i,29]<-output.SCWeibull3[[1]][1]
  weib.est[i,30]<-output.SCWeibull3[[1]][4]
  
  weib.est[i,31]<-output.SCWeibull4[[1]][1,1]
  weib.est[i,32]<-output.SCWeibull4[[1]][1,4]
  weib.est[i,33]<-output.SCWeibull4[[1]][2,1]
  weib.est[i,34]<-output.SCWeibull4[[1]][2,4]
  weib.est[i,35]<-output.SCWeibull4[[1]][3,1]
  weib.est[i,36]<-output.SCWeibull4[[1]][3,4]
  weib.est[i,37]<-output.SCWeibull4[[1]][4,1]
  weib.est[i,38]<-output.SCWeibull4[[1]][4,4]
  weib.est[i,39]<-output.SCWeibull4[[1]][5,1]
  weib.est[i,40]<-output.SCWeibull4[[1]][5,4]

  weib.est[i,41]<-output.SCWeibull5[[1]][1]
  weib.est[i,42]<-output.SCWeibull5[[1]][4]
  
  #store rmse
  weib.rmse[i,10]<-sqrt((tru.est[i,1]-weib.est[i,19])^2)
  weib.rmse[i,11]<-sqrt((tru.est[i,2]-weib.est[i,21])^2)
  weib.rmse[i,12]<-sqrt((tru.est[i,3]-weib.est[i,23])^2)
  weib.rmse[i,13]<-sqrt((tru.est[i,4]-weib.est[i,25])^2)
  weib.rmse[i,14]<-sqrt((tru.est[i,5]-weib.est[i,27])^2)
  weib.rmse[i,15]<-sqrt((tru.est[i,6]-weib.est[i,29])^2)
  weib.rmse[i,16]<-sqrt((tru.est2[i,1]-weib.est[i,31])^2)
  weib.rmse[i,17]<-sqrt((tru.est2[i,2]-weib.est[i,33])^2)
  weib.rmse[i,18]<-sqrt((tru.est2[i,3]-weib.est[i,35])^2)
  weib.rmse[i,19]<-sqrt((tru.est2[i,4]-weib.est[i,37])^2)
  weib.rmse[i,20]<-sqrt((tru.est2[i,5]-weib.est[i,39])^2)
  weib.rmse[i,21]<-sqrt((tru.est2[i,11]-weib.est[i,51])^2)
  
  #calculate upper and lower 95% CI's
  b0.lower<-output.SCWeibull[[2]][1,1]
  b0.upper<-output.SCWeibull[[2]][1,5]
  b1.lower<-output.SCWeibull[[2]][2,1]
  b1.upper<-output.SCWeibull[[2]][2,5]
  g0.lower<-output.SCWeibull2[[2]][1,1]
  g0.upper<-output.SCWeibull2[[2]][1,5]
  g1.lower<-output.SCWeibull2[[2]][2,1]
  g1.upper<-output.SCWeibull2[[2]][2,5]
  g2.lower<-output.SCWeibull2[[2]][3,1]
  g2.upper<-output.SCWeibull2[[2]][3,5]
  p.lower<-output.SCWeibull3[[2]][1]
  p.upper<-output.SCWeibull3[[2]][5]
  
  w1.lower<-output.SCWeibull4[[2]][1,1]
  w1.upper<-output.SCWeibull4[[2]][1,5]
  w2.lower<-output.SCWeibull4[[2]][2,1]
  w2.upper<-output.SCWeibull4[[2]][2,5]
  w3.lower<-output.SCWeibull4[[2]][3,1]
  w3.upper<-output.SCWeibull4[[2]][3,5]
  w4.lower<-output.SCWeibull4[[2]][4,1]
  w4.upper<-output.SCWeibull4[[2]][4,5]
  w5.lower<-output.SCWeibull4[[2]][5,1]
  w5.upper<-output.SCWeibull4[[2]][5,5]
  l1.lower<-output.SCWeibull5[[2]][1]
  l1.upper<-output.SCWeibull5[[2]][5]
  
  #store coverage parameters
  weib.cp[i,10]<-ifelse(tru.est[i,1]>b0.lower & tru.est[i,1]<b0.upper, 1,0)
  weib.cp[i,11]<-ifelse(tru.est[i,2]>b1.lower & tru.est[i,2]<b1.upper, 1,0)
  weib.cp[i,12]<-ifelse(tru.est[i,3]>g0.lower & tru.est[i,3]<g0.upper, 1,0)
  weib.cp[i,13]<-ifelse(tru.est[i,4]>g1.lower & tru.est[i,4]<g1.upper, 1,0)
  weib.cp[i,14]<-ifelse(tru.est[i,5]>g2.lower & tru.est[i,5]<g2.upper, 1,0)
  weib.cp[i,15]<-ifelse(tru.est[i,6]>p.lower & tru.est[i,6]<p.upper, 1,0)
  
  weib.cp[i,16]<-ifelse(tru.est2[i,1]>w1.lower & tru.est2[i,1]<w1.upper, 1,0)
  weib.cp[i,17]<-ifelse(tru.est2[i,2]>w2.lower & tru.est2[i,2]<w2.upper, 1,0)
  weib.cp[i,18]<-ifelse(tru.est2[i,3]>w3.lower & tru.est2[i,3]<w3.upper, 1,0)
  weib.cp[i,19]<-ifelse(tru.est2[i,4]>w4.lower & tru.est2[i,4]<w4.upper, 1,0)
  weib.cp[i,20]<-ifelse(tru.est2[i,5]>w5.lower & tru.est2[i,5]<w5.upper, 1,0)
  weib.cp[i,21]<-ifelse(tru.est2[i,6]>l1.lower & tru.est2[i,6]<l1.upper, 1,0)
  
  }
  #combine matrices and label variables
  main.data<-cbind(tru.est,exp.est,weib.est,exp.rmse,weib.rmse,exp.cp,weib.cp)
  colnames(main.data)<-c("true.x0","true.x1","true.z0","true.z1","true.z2","true.p","cen.lat","cen.obs",
  	"exp.x0","exp.x0.se","exp.x1","exp.x1.se",
  	"cexp.x0","cexp.x0.se","cexp.x1","cexp.x1.se","cexp.z0","cexp.z0.se","cexp.z1","cexp.z1.se","cexp.z2","cexp.z2.se",
  	"scexp.x0","scexp.x0.se","scexp.x1","scexp.x1.se","scexp.z0","scexp.z0.se","scexp.z1","scexp.z1.se","scexp.z2","scexp.z2.se",
  	"scexp.w1","scexp.w1.se","scexp.w2","scexp.w2.se","scexp.w3","scexp.w3.se","scexp.w4","scexp.w4.se","scexp.w5","scexp.w5.se",
  	"scexp.l","scexp.l.se",
  	"wei.x0","wei.x0.se","wei.x1","wei.x1.se","wei.p","wei.p.se",
  	"cwei.x0","cwei.x0.se","cwei.x1","cwei.x1.se","cwei.z0","cwei.z0.se","cwei.z1","cwei.z1.se","cwei.z2","cwei.z2.se","cwei.p","cwei.p.se",
  	"scwei.x0","scwei.x0.se","scwei.x1","scwei.x1.se","scwei.z0","scwei.z0.se","scwei.z1","scwei.z1.se","scwei.z2","scwei.z2.se","scwei.p","scwei.p.se",
  	"scwei.w1","scwei.w1.se","scwei.w2","scwei.w2.se","scwei.w3","scwei.w3.se","scwei.w4","scwei.w4.se","scwei.w5","scwei.w5.se",
  	"scwei.l","scwei.l.se",
  	"exp.x0.rmse","exp.x1.rmse","cexp.x0.rmse","cexp.x1.rmse","cexp.z0.rmse","cexp.z1.rmse","cexp.z2.rmse",
  	"scexp.x0.rmse","scexp.x1.rmse","scexp.z0.rmse","scexp.z1.rmse","scexp.z2.rmse",
  	"scexp.w1.rmse","scexp.w2.rmse","scexp.w3.rmse","scexp.w4.rmse","scexp.w5.rmse",
  	"scexp.v1.rmse","scexp.v2.rmse","scexp.v3.rmse","scexp.v4.rmse","scexp.v5.rmse",
  	"scexp.l.rmse",
  	"wei.x0.rmse","wei.x1.rmse","wei.p.rmse","cwei.x0.rmse","cwei.x1.rmse","cwei.z0.rmse","cwei.z1.rmse","cwei.z2.rmse","cwei.p.rmse",
  	"scwei.x0.rmse","scwei.x1.rmse","scwei.z0.rmse","scwei.z1.rmse","scwei.z2.rmse","scwei.p.rmse",
  	"scwei.w1.rmse","scwei.w2.rmse","scwei.w3.rmse","scwei.w4.rmse","scwei.w5.rmse",
  	"scwei.l.rmse",	
  	"exp.x0.cp","exp.x1.cp","cexp.x0.cp","cexp.x1.cp","cexp.z0.cp","cexp.z1.cp","cexp.z2.cp",
  	"wei.x0.cp","wei.x1.cp","wei.p.cp", "cwei.x0.cp","cwei.x1.cp","cwei.z0.cp","cwei.z1.cp","cwei.z2.cp","cwei.p.cp",
  	"scwei.x0.cp","scwei.x1.cp","scwei.z0.cp","scwei.z1.cp","scwei.z2.cp","scwei.p.cp",
  	"scwei.w1.cp","scwei.w2.cp","scwei.w3.cp","scwei.w4.cp","scwei.w5.cp",
  	"scwei.l.cp"
  	)
  
  #save dataset
  main.data<-as.data.frame(main.data)
  write.dta(main.data,"main.data.dta")
  
  #the end
