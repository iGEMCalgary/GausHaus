#####################
# Gaus house Proof of concept scripts
# Andrew Symes
# iGEM 2020
#####################
set.seed(16)
library(GPfit)
library(animation)
library(plot3D)


#ssn is a dataframe of the euclidean coordinates
#ssn <- read.csv("~/Desktop/GAUSHOUSE/SEGI8Dynamic20.csv")
ssn <- read.csv("~/Desktop/EG1Dynamic20.csv")

index <- c(1:length(ssn$resi))
atomsperchain <- min(ifelse(ssn$resi[100:10000] == 1,index,10000)) + 98


#Generate dataframes for every axis of all points
for(ii in 1:atomsperchain){
  tempxvec <- c()
  tempyvec <- c()
  tempzvec <- c()
  
  for(i in 0:20){
    tempx <- ssn$x[atomsperchain*i + ii]
    tempxvec <- c(tempxvec,tempx)
    
    tempy <- ssn$y[atomsperchain*i + ii]
    tempyvec <- c(tempyvec,tempy)
    
    tempz <- ssn$z[atomsperchain*i + ii]
    tempzvec <- c(tempzvec,tempz)
    
  }
  
  if(ii == 1){
    x.df <- data.frame(tempxvec)
    y.df <- data.frame(tempyvec)
    z.df <- data.frame(tempzvec)
    
  }
  
  if(ii != 1){
    x.df <- data.frame(x.df,tempxvec)
    y.df <- data.frame(y.df,tempyvec)
    z.df <- data.frame(z.df,tempzvec)
    print(ii)
    
  }
}

#Name them after your atom names
names(x.df) <- ssn$name[1:21]
names(y.df) <- ssn$name[1:21]
names(z.df) <- ssn$name[1:21]


#1 if it is an alpha carbon 0 if it is not
alphacarb <- ifelse(ssn$name== "CA",1,0)

#Generate dataframes necessary for getting x y and z for each point in time.

aCarbx.df <- data.frame((x.df[,10]))
aCarby.df <- data.frame(y.df[,10])
aCarbz.df <- data.frame(z.df[,10])

for (i in 10:4114){
  if(alphacarb[i] == 1){
    aCarbx.df <- data.frame(aCarbx.df,x.df[,i])
    aCarby.df <- data.frame(aCarby.df,y.df[,i])
    aCarbz.df <- data.frame(aCarbz.df,z.df[,i])
    
  }
  
}

###X measurement
xmat <- as.matrix(aCarbx.df[1:10,])
xmat <- t(xmat) /100

xfit <- as.matrix((aCarbx.df[11,]))
xfit <- t(xfit) /100

xGP1 <- GP_fit(xmat,xfit,trace = T,maxit=10)

#y measure
ymat <- as.matrix(aCarby.df[1:10,])
ymat <- t(ymat) /100

yfit <- as.matrix((aCarby.df[11,]))
yfit <- t(yfit) /100

yGP1 <- GP_fit(ymat,yfit,trace = T,maxit=10)

#z measure
zmat <- as.matrix(aCarbz.df[1:10,])
zmat <- t(zmat) /100

zfit <- as.matrix((aCarbz.df[11,]))
zfit <- t(zfit) /100

zGP1 <- GP_fit(zmat,zfit,trace = T,maxit=10)

#Print out the readings
print(xGP)
print(yGP)
print(zGP)

print(xGP1)
print(yGP1)
print(zGP1)

print(data.frame(c('x','y','z'),c(xGP$sig2,yGP$sig2,zGP$sig2)))
print(data.frame(c('x','y','z'),c(xGP1$sig2,yGP1$sig2,zGP1$sig2)))
