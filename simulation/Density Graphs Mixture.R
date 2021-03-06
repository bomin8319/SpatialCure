#graph B2

##Ben Bagozzi
##
##June 8th 2011
##
##Combine MC results
##
##
##
##

#clear memory  
rm(list=ls())                                           
library(car)
library(Hmisc)                                                     
library(mvtnorm)
library(foreign)
library(graphics)
library(MASS)
library(lattice)
library(tseries)
library(Matrix)
library(Design)
library(msm)
library(corpcor)
library(Zelig)


Results<-read.dta("/Users/bomin8319/Desktop/SpatialCure/simulation/main.data_U.dta")
Results<-as.matrix(Results, )
resize.win <- function(Width=6, Height=6)
{
        # works for windows
    dev.off(); # dev.new(width=6, height=6)
    windows(record=TRUE, width=Width, height=Height)
}
resize.win(7,4)

par(mfrow=c(2,3))
par(cex.lab=1)
par(cex.axis=1)
par(cex.main=1)

par(mar=c(5.1,4.1,2.1,2.1))
#sets the bottom, left, top and right 

#B0
local.xlim<-c(0,2.5)
local.ylim<-c(0,4)
plot(density(Results[,15],na.rm=TRUE), main = "",  ylab = "", xlab = "", xlim=local.xlim, ylim=local.ylim,col="blue", xaxt='n', yaxt='n')
abline(v=1,lty=3)
par(new=TRUE)
plot(density(Results[,19],na.rm=TRUE), main = "",  ylab = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="red", xaxt='n', yaxt='n')
par(new=TRUE)
plot(density(Results[,29],na.rm=TRUE), main = "",  ylab = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="pink", xaxt='n', yaxt='n')
par(new=TRUE)
plot(density(Results[,51],na.rm=TRUE), main = "",  ylab = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="forestgreen", xaxt='n', yaxt='n')
par(new=TRUE)
plot(density(Results[,57],na.rm=TRUE), main = "",  ylab = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="orange", xaxt='n', yaxt='n')
par(new=TRUE)
plot(density(Results[,69],na.rm=TRUE), main = "", xlab = "Beta 0",xlim=local.xlim, ylim=local.ylim,col="purple")
text (x=.6, y =3.2, "Cure.Exp",col="red")
text (x=1.33, y =3.4, "Wei",col="forestgreen")
text (x=1.33, y =3.2, "Exp",col="blue")
text (x=0.6, y =3.4, "Spat.Cure.Wei",col="purple")
text (x=0.6, y =3.6, "Spat.Cure.Exp",col="pink")
text (x=0.6, y =3.0, "Cure.Wei",col="orange")
par(new=FALSE)
 

#B1
local.xlim<-c(3,3.9)
local.ylim<-c(0,7.5)
plot(density(Results[,17],na.rm=TRUE), main = "",  ylab = "", xlab = "", xlim=local.xlim, ylim=local.ylim,col="blue", xaxt='n', yaxt='n')
abline(v=3.5,lty=3)
par(new=TRUE)
plot(density(Results[,21],na.rm=TRUE), main = "",  ylab = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="red", xaxt='n', yaxt='n')
par(new=TRUE)
plot(density(Results[,31],na.rm=TRUE), main = "",  ylab = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="pink", xaxt='n', yaxt='n')
par(new=TRUE)
plot(density(Results[,53],na.rm=TRUE), main = "",  ylab = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="forestgreen", xaxt='n', yaxt='n')
par(new=TRUE)
plot(density(Results[,59],na.rm=TRUE), main = "",  ylab = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="orange", xaxt='n', yaxt='n')
par(new=TRUE)
plot(density(Results[,71],na.rm=TRUE), main = "", xlab = "Beta 1",xlim=local.xlim, ylim=local.ylim,col="purple")
text (x=3.75, y =5, "Cure.Exp",col="red")
text (x=3.15, y =5, "Wei",col="forestgreen")
text (x=3.15, y =4.75, "Exp",col="blue")
text (x=3.75, y =4.75, "Spat.Cure.Wei",col="purple")
text (x=3.75, y =4.5, "Spat.Cure.Exp",col="pink")
text (x=3.75, y =4.25, "Cure.Wei",col="orange")
par(new=FALSE)

#B1
local.ylim<-c(0,13)
local.xlim<-c(0.6,1.2)
plot(density(Results[,55],na.rm=TRUE), main = "",  ylab = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="forestgreen", xaxt='n', yaxt='n')
abline(v=1,lty=3)
par(new=TRUE)
plot(density(Results[,67],na.rm=TRUE), main = "",  ylab = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="orange", xaxt='n', yaxt='n')
par(new=TRUE)
plot(density(Results[,79],na.rm=TRUE), main =" ", xlab = "P",xlim=local.xlim, ylim=local.ylim,col="purple")
text (x=0.7, y =11.6, "Wei",col="forestgreen")
text (x=0.7, y =11, "Cure.Wei",col="orange")
text (x=1.1, y =6, "Spat.Cure.Wei",col="purple")
par(new=FALSE)

#B1
local.ylim<-c(0,.5)
local.xlim<-c(-5.5,3)
plot(density(Results[,23],na.rm=TRUE), main = "",  ylab = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="red", xaxt='n', yaxt='n')
abline(v=-2,lty=3)
par(new=TRUE)
plot(density(Results[,33],na.rm=TRUE), main = "",  ylab = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="pink", xaxt='n', yaxt='n')
par(new=TRUE)
plot(density(Results[,61],na.rm=TRUE), main = "",  ylab = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="orange", xaxt='n', yaxt='n')
par(new=TRUE)
plot(density(Results[,73],na.rm=TRUE), main = "", xlab = "Gamma 0",xlim=local.xlim, ylim=local.ylim,col="purple")
text (x=0.5, y =.42, "Cure.Exp",col="red")
text (x=0.5, y =.4, "Spat.Cure.Exp",col="pink")
text (x=0.5, y =.38, "Cure.Wei",col="orange")
text (x=0.5, y =0.36, "Spat.Cure.Wei",col="purple")
par(new=FALSE)


#B1
local.ylim<-c(0,1.2)
local.xlim<-c(-4,1)
plot(density(Results[,25],na.rm=TRUE), main = "",  ylab = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="red", xaxt='n', yaxt='n')
abline(v=-2,lty=3)
par(new=TRUE)
plot(density(Results[,35],na.rm=TRUE), main = "",  ylab = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="pink", xaxt='n', yaxt='n')
par(new=TRUE)
plot(density(Results[,63],na.rm=TRUE), main = "",  ylab = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="orange", xaxt='n', yaxt='n')
par(new=TRUE)
plot(density(Results[,75],na.rm=TRUE), main = "", xlab = "Gamma 1",xlim=local.xlim, ylim=local.ylim,col="purple")
text (x=-3, y =1.0, "Cure.Exp",col="red")
text (x=-3, y =0.95, "Spat.Cure.Exp",col="pink")
text (x=-3, y =0.9, "Cure.Wei",col="orange")
text (x=-3, y =0.85, "Spat.Cure.Wei",col="purple")
par(new=FALSE)

#B1
local.ylim<-c(0,1.1)
local.xlim<-c(1,5)
plot(density(Results[,27],na.rm=TRUE), main = "",  ylab = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="red", xaxt='n', yaxt='n')
abline(v=3,lty=3)
par(new=TRUE)
plot(density(Results[,37],na.rm=TRUE), main = "",  ylab = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="pink", xaxt='n', yaxt='n')
par(new=TRUE)
plot(density(Results[,65],na.rm=TRUE), main = "",  ylab = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="orange", xaxt='n', yaxt='n')
par(new=TRUE)
plot(density(Results[,77],na.rm=TRUE), main = "", xlab = "Gamma 2",xlim=local.xlim, ylim=local.ylim,col="purple")
text (x=2, y =1.0, "Cure.Exp",col="red")
text (x=2, y =0.95, "Spat.Cure.Exp",col="pink")
text (x=2, y =0.9, "Cure.Wei",col="orange")
text (x=2, y =0.85, "Spat.Cure.Wei",col="purple")
par(new=FALSE)


#B1
local.ylim<-c(0,1.1)
local.xlim<-c(0,3.5)
plot(density(Results[,49],na.rm=TRUE), main = "",  ylab = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="pink", xaxt='n', yaxt='n')
abline(v=1,lty=3)
par(new=TRUE)
plot(density(Results[,91],na.rm=TRUE), main = "", xlab = "lambda",xlim=local.xlim, ylim=local.ylim,col="purple")
text (x=2.5, y =1.0, "Spat.Cure.Exp",col="pink")
text (x=2.5, y =0.95, "Spat.Cure.Wei",col="purple")
par(new=FALSE)


#B1
local.ylim<-c(0,2.8)
local.xlim<-c(-1.5,1)
plot(density(Results[,39],na.rm=TRUE), main = "",  ylab = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="pink", xaxt='n', yaxt='n')
abline(v=-0.6428763,lty=3)
par(new=TRUE)
plot(density(Results[,81],na.rm=TRUE), main = "", xlab = "w1",xlim=local.xlim, ylim=local.ylim,col="purple")
text (x=0.2, y =2, "Spat.Cure.Exp",col="pink")
text (x=0.2, y =1.8, "Spat.Cure.Wei",col="purple")
par(new=FALSE)

#B1
local.ylim<-c(0,2.8)
local.xlim<-c(-1.3,1)
plot(density(Results[,41],na.rm=TRUE), main = "",  ylab = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="pink", xaxt='n', yaxt='n')
abline(v=-0.2796147 ,lty=3)
par(new=TRUE)
plot(density(Results[,83],na.rm=TRUE), main = "", xlab = "w2",xlim=local.xlim, ylim=local.ylim,col="purple")
text (x=0.3, y =2, "Spat.Cure.Exp",col="pink")
text (x=0.3, y =1.8, "Spat.Cure.Wei",col="purple")
par(new=FALSE)

#w3
local.ylim<-c(0,2.8)
local.xlim<-c(-1.3,0.8)
plot(density(Results[,43],na.rm=TRUE), main = "",  ylab = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="pink", xaxt='n', yaxt='n')
abline(v=-0.7002307 ,lty=3)
par(new=TRUE)
plot(density(Results[,85],na.rm=TRUE), main = "", xlab = "w3",xlim=local.xlim, ylim=local.ylim,col="purple")
text (x=0, y =2, "Spat.Cure.Exp",col="pink")
text (x=0, y =1.8, "Spat.Cure.Wei",col="purple")
par(new=FALSE)


#w3
local.ylim<-c(0,3.2)
local.xlim<-c(-0.5,1.2)
plot(density(Results[,45],na.rm=TRUE), main = "",  ylab = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="pink", xaxt='n', yaxt='n')
abline(v=0.5647177 ,lty=3)
par(new=TRUE)
plot(density(Results[,87],na.rm=TRUE), main = "", xlab = "w4",xlim=local.xlim, ylim=local.ylim,col="purple")
text (x=0, y =2, "Spat.Cure.Exp",col="pink")
text (x=0, y =1.8, "Spat.Cure.Wei",col="purple")
par(new=FALSE)


#w3
local.ylim<-c(0,3.7)
local.xlim<-c(-0.7,1.6)
plot(density(Results[,47],na.rm=TRUE), main = "",  ylab = "", xlab = "",xlim=local.xlim, ylim=local.ylim,col="pink", xaxt='n', yaxt='n')
abline(v=1.058004 ,lty=3)
par(new=TRUE)
plot(density(Results[,89],na.rm=TRUE), main = "", xlab = "w5",xlim=local.xlim, ylim=local.ylim,col="purple")
text (x=0, y =2, "Spat.Cure.Exp",col="pink")
text (x=0, y =1.8, "Spat.Cure.Wei",col="purple")
par(new=FALSE)

