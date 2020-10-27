#####################
# Gaus house Proof of concept scripts
# Andrew Symes
# iGEM 2020
#####################
library(GPfit)
library(animation)
library(plot3D)
ssn1 <- read.csv("~/Desktop/GAUSHOUSE/1SSNDATA.csv")

#Generate dataframes for every axis of all points
for(ii in 1:1408){
  tempxvec <- c()
  tempyvec <- c()
  tempzvec <- c()
  
  for(i in 0:101){
    tempx <- ssn$x[1408*i + ii]
    tempxvec <- c(tempxvec,tempx)
    
    tempy <- ssn$y[1408*i + ii]
    tempyvec <- c(tempyvec,tempy)
    
    tempz <- ssn$z[1408*i + ii]
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
    
  }
}

#Name them after your atom names
names(x.df) <- ssn$name[1:102]
names(y.df) <- ssn$name[1:102]
names(z.df) <- ssn$name[1:102]


#1 if it is an alpha carbon 0 if it is not
alphacarb <- ifelse(ssn$name== "C",1,0)

#Generate dataframes necessary for getting x y and z for each point in time.

aCarbx.df <- data.frame((x.df[,10]))
aCarby.df <- data.frame(y.df[,10])
aCarbz.df <- data.frame(z.df[,10])

for (i in 10:1408){
  if(alphacarb[i] == 1){
    aCarbx.df <- data.frame(aCarbx.df,x.df[,i])
    aCarby.df <- data.frame(aCarby.df,y.df[,i])
    aCarbz.df <- data.frame(aCarbz.df,z.df[,i])
    
  }
  
}

#IMAGE GENERATING FUNCTION
# Plots each alphacarbon then connects them with a line.
# It does this for every frame
makeIMG <- function(){
for(i in 1:102){
  tempxx <- c()
  tempyy <- c()
  tempzz <- c()
  for(j in 1:136){
    tempxx <- c(tempxx,as.numeric(c(aCarbx.df[i,])[j]))
    tempyy <- c(tempyy,as.numeric(c(aCarby.df[i,])[j]))
    tempzz <- c(tempzz,as.numeric(c(aCarbz.df[i,])[j]))

  }
  scatter3D(tempxx,tempyy,tempzz,type="l")
  
}
}

#Makes the final gif using the IMG function seen above
$saveGIF(makeIMG(),"~/Desktop/testGif.gif")


##################
###Generalizable graphing function
GmakeIMG <- function(aCarx.df,aCary.df,aCarz.df){
  for(i in 1:length(aCarx.df[,1])){
    tempxx <- c()
    tempyy <- c()
    tempzz <- c()
    for(j in 1:length(aCarx.df[1,])){
      tempxx <- c(tempxx,as.numeric(c(aCarx.df[i,])[j]))
      tempyy <- c(tempyy,as.numeric(c(aCary.df[i,])[j]))
      tempzz <- c(tempzz,as.numeric(c(aCarz.df[i,])[j]))
      
    }
    scatter3D(tempxx,tempyy,tempzz,type="l")
    
  }
}

#Makes the final gif using the IMG function seen above
saveGIF(GmakeIMG(),"~/Desktop/testGif.gif")

############


####

#GPfit Time

####
c9x <- as.vector(x.df[,9])
c9y <- as.vector(y.df[,9])
c9z <- as.vector(z.df[,9])
c9Out <- c(c9x[102],c9y[102],c9z[102]) /100
c9outmat <- matrix(c9Out,3,1)

c9x <- as.vector(x.df[,9])[1:12] /100
c9y <- as.vector(y.df[,9])[1:12] /100
c9z <- as.vector(z.df[,9])[1:12] /100
c9Mat <- matrix(c(c9x,c9y,c9z),3,12,byrow = T)



xProGP <- GP_fit(c9Mat,c9Out,trace = T,maxit=5)


c9newx <- as.vector(x.df[,9])[90:101] /100
c9newy <- as.vector(y.df[,9])[90:101] /100
c9newz <- as.vector(z.df[,9])[90:101] /100
c9newMat <- matrix(c(c9newx,c9newy,c9newz),3,12,byrow = T)

plot.GP(xProGP,1)


preds <- predict.GP(xProGP,c9newMat)

#### Function for predicting more

Crabbe <- function(Gaus,xvec,yvec,zvec,leny){
  #Select out the first 12 for testing
  # leny <- 12
  # xvec <- TM11[1:12]
  # yvec <- TM11[13:24]
  # zvec <- TM11[25:36]
  # Gaus <- TGP1
  lenyy <- leny -11
  cx <- xvec[lenyy:leny]
  cy <- yvec[lenyy:leny] 
  cz <- zvec[lenyy:leny] 
  cMat <- matrix(c(cx,cy,cz),3,12,byrow=T)
  
  temppred <- predict.GP(Gaus,cMat)
  
  return(temppred$Y_hat)
  
}

Lucius <- function(Gaus,xvec,xvec2,xvec3,yvec,yvec2,yvec3,zvec,zvec2,zvec3,leny){
  #Select out the first 12 for testing
  # leny <- 12
  # xvec <- TM11[1:12]
  # yvec <- TM11[13:24]
  # zvec <- TM11[25:36]
  # Gaus <- TGP1
  lenyy <- leny -11
  cx <- xvec[lenyy:leny]
  cy <- yvec[lenyy:leny] 
  cz <- zvec[lenyy:leny] 
  
  cx2 <- xvec2[lenyy:leny]
  cy2 <- yvec2[lenyy:leny] 
  cz2 <- zvec2[lenyy:leny]
  
  cx3 <- xvec3[lenyy:leny]
  cy3 <- yvec3[lenyy:leny] 
  cz3 <- zvec3[lenyy:leny] 
  
  cMat <- matrix(c(cx,cy,cz),3,12,byrow=T)
  
  temppred <- predict.GP(Gaus,cMat)
  
  return(temppred$Y_hat)
  
}

Crabbe(xProGP,as.vector(x.df[,9]),as.vector(y.df[,9]),as.vector(z.df[,9]),102)


Goyle <- function(Gaus,xvec,yvec,zvec,leny,KG){
  
  initial <- Crabbe(Gaus,xvec,yvec,zvec,leny)
  xvec <- c(xvec,initial[1])
  yvec <- c(yvec,initial[2])
  zvec <- c(zvec,initial[3])
  Saved <- data.frame(initial)
  
  for(i in 1:KG){
    temprun <- Crabbe(Gaus,xvec,yvec,zvec,leny+ i)
    xvec <- c(xvec,temprun[1])
    yvec <- c(yvec,temprun[2])
    zvec <- c(zvec,temprun[3])
    Saved <- data.frame(Saved,temprun)
  }
  
  return(Saved)
  
}

Snape <- function(Gaus,xvec,xvec2,xvec3,yvec,yvec2,yvec3,zvec,zvec2,zvec3,leny,KG){
  
  initial <- Lucius(Gaus,xvec,xvec2,xvec3,yvec,yvec2,yvec3,zvec,zvec2,zvec3,leny)
  xvec <- c(xvec,initial[1])
  yvec <- c(yvec,initial[2])
  zvec <- c(zvec,initial[3])
  Saved <- data.frame(initial)
  
  for(i in 1:KG){
    temprun <- Lucius(Gaus,xvec,xvec2,xvec3,yvec,yvec2,yvec3,zvec,zvec2,zvec3,leny+ i)
    xvec <- c(xvec,temprun[1])
    yvec <- c(yvec,temprun[2])
    zvec <- c(zvec,temprun[3])
    Saved <- data.frame(Saved,temprun)
  }
  
  return(Saved)
  
}

# ####TESTING SECTION
# iii <- Goyle(xProGP,as.vector(x.df[,9]),as.vector(y.df[,9]),as.vector(z.df[,9]),102,5)
# 
# 
# #Test design matrix
# TM1 <- c(as.numeric(aCarbx.df[1,1:12]),as.numeric(aCarby.df[1,1:12]),as.numeric(aCarbz.df[1,1:12])) /100
# TM2 <- c(as.numeric(aCarbx.df[1:12,2]),as.numeric(aCarby.df[1:12,2]),as.numeric(aCarbz.df[1:12,2])) /100
# TM3 <-  c(as.numeric(aCarbx.df[3,1:12]),as.numeric(aCarby.df[3,1:12]),as.numeric(aCarbz.df[3,1:12])) /100
# TM4 <-  c(as.numeric(aCarbx.df[4,1:12]),as.numeric(aCarby.df[4,1:12]),as.numeric(aCarbz.df[4,1:12])) /100
# TM5 <-  c(as.numeric(aCarbx.df[5,1:12]),as.numeric(aCarby.df[5,1:12]),as.numeric(aCarbz.df[5,1:12])) /100
# 
# Mat1 <- matrix(TM1,3,12,byrow = T)
# Mat2 <- matrix(TM2,3,12,byrow = T)
# Mat3 <- matrix(TM3,3,12,byrow = T)
# Mat4 <- matrix(TM4,3,12,byrow = T)
# Mat5 <- matrix(TM5,3,12,byrow = T)
# 
# #Test fitting outputs
# Out1 <- matrix(c(as.numeric(aCarbx.df[1,13]),as.numeric(aCarby.df[1,13]),as.numeric(aCarbz.df[1,13])),3,1,byrow = T) /100
# Out2 <- matrix(c(as.numeric(aCarbx.df[13,2]),as.numeric(aCarby.df[13,2]),as.numeric(aCarbz.df[13,2])),3,1,byrow = T) /100
# Out3 <- matrix(c(as.numeric(aCarbx.df[3,13]),as.numeric(aCarby.df[3,13]),as.numeric(aCarbz.df[3,13])),3,1,byrow = T) /100
# Out4 <- matrix(c(as.numeric(aCarbx.df[4,13]),as.numeric(aCarby.df[4,13]),as.numeric(aCarbz.df[4,13])),3,1,byrow = T) /100
# Out5 <- matrix(c(as.numeric(aCarbx.df[5,13]),as.numeric(aCarby.df[5,13]),as.numeric(aCarbz.df[5,13])),3,1,byrow = T) /100
# 
# 
# 
# TGP1 <- GP_fit(Mat1,Out1,trace = T,maxit=5)
# TGP2 <- GP_fit(Mat2,Out2,trace = T,maxit=5)
# TGP3 <- GP_fit(Mat3,Out3,trace = T,maxit=5)
# TGP4 <- GP_fit(Mat4,Out4,trace = T,maxit=5)
# TGP5 <- GP_fit(Mat5,Out5,trace = T,maxit=5)
# 
# 
# TM11 <- c(as.numeric(aCarbx.df[1,2:13]),as.numeric(aCarby.df[1,2:13]),as.numeric(aCarbz.df[1,2:13])) /100
# Mat11 <- Mat1 <- matrix(TM11,3,12,byrow = T)
# 
# TM22 <- c(as.numeric(aCarbx.df[2,2:13]),as.numeric(aCarby.df[2,2:13]),as.numeric(aCarbz.df[2,2:13])) /100
# 
# 
# PredT1 <- Goyle(TGP1,TM1[1:12],TM1[13:24],TM1[25:36],12,10)
# predict(TGP1,Mat11)
# Crabbe(TGP1,TM11[1:12],TM11[13:24],TM11[25:36],12)
# Goyle(TGP1,TM11[1:12],TM11[13:24],TM11[25:36],12,10)

####


Draco <- function(Gaus,xdf,ydf,zdf,leny,KG,joints){
  lenyy <- leny + 1
  xvec <- as.vector(xdf[,1]) / 100
  yvec <- as.vector(xdf[,1]) / 100
  zvec <- as.vector(xdf[,1]) / 100
  initialize <- Goyle(Gaus,xvec,yvec,zvec,leny,KG)
  
  xdfNew <- data.frame(initialize[1,])
  ydfNew <- data.frame(initialize[2,])
  zdfNew <- data.frame(initialize[3,])
  
  for(i in 2:joints){
    xvec <- as.vector(xdf[,i]) /100
    yvec <- as.vector(ydf[,i]) /100
    zvec <- as.vector(zdf[,i]) /100
    tempMat <- matrix(c(xvec[1:leny],yvec[1:leny],zvec[1:leny]),3,leny,byrow=T)
    tempout <- matrix(c(xvec[lenyy],yvec[lenyy],zvec[lenyy]),3,1)
    print(tempMat)
    print(tempout)
    tempGP <- GP_fit(tempMat,tempout,trace = T,maxit=5)
    tempests <- Goyle(tempGP,xvec,yvec,zvec,leny,KG)
    
    xdfNew <- data.frame(xdfNew,tempests[1,])
    ydfNew <- data.frame(ydfNew,tempests[2,])
    zdfNew <- data.frame(zdfNew,tempests[3,])
    
    print(paste("run #",  i+1))
  }
  
  return(data.frame(xdfNew,ydfNew,zdfNew))
  
}


####WIP
Bellatrix <- function(xdf,ydf,zdf,leny,KG,joints){
  xvec <- as.vector(xdf[,1]) / 100
  yvec <- as.vector(ydf[,1]) / 100
  zvec <- as.vector(zdf[,1]) / 100
  xvec2 <- as.vector(xdf[,2]) / 100
  yvec2 <- as.vector(ydf[,2]) / 100
  zvec2 <- as.vector(zdf[,2]) / 100
  xvec3 <- as.vector(xdf[,3]) / 100
  yvec3 <- as.vector(ydf[,3]) / 100
  zvec3 <- as.vector(zdf[,3]) / 100
  
  iniMat <- matrix(c(xvec[1:12],yvec[1:12],zvec[1:12],xvec2[1:12],yvec2[1:12],zvec2[1:12],xvec3[1:12],yvec3[1:12],zvec3[1:12]),9,12,byrow=T)
  iniO <- matrix(c(xvec[13],yvec[13],zvec[13],xvec2[13],yvec2[13],zvec2[13],xvec3[13],yvec3[13],zvec3[13]),9,1)
  print(iniMat)
  print(iniO)
  Gaus <- GP_fit(iniMat,iniO,trace = T,maxit = 5)
  initialize <- Snape(Gaus,xvec,xvec2,xvec3,yvec,yvec2,yvec3,zvec,zvec2,zvec3,leny,KG)
  
  xdfNew <- data.frame(initialize[1,])
  ydfNew <- data.frame(initialize[2,])
  zdfNew <- data.frame(initialize[3,])
  
  for(i in 2:joints){
    xvec <- as.vector(xdf[,i]) /100
    yvec <- as.vector(ydf[,i]) /100
    zvec <- as.vector(zdf[,i]) /100
    xvec2 <- as.vector(xdf[,i+1]) /100
    yvec2 <- as.vector(ydf[,i+1]) /100
    zvec2 <- as.vector(zdf[,i+1]) /100
    xvec3 <- as.vector(xdf[,i+2]) /100
    yvec3 <- as.vector(ydf[,i+2]) /100
    zvec3 <- as.vector(zdf[,i+2]) /100
    
    tempMat <- matrix(c(xvec[1:12],yvec[1:12],zvec[1:12],
                        xvec2[1:12],yvec2[1:12],zvec2[1:12],
                        xvec3[1:12],yvec3[1:12],zvec3[1:12]),9,12,byrow=T)
    tempout <- matrix(c(xvec[13],yvec[13],zvec[13],xvec2[13],yvec2[13],zvec2[13],xvec3[13],yvec3[13],zvec3[13]),9,1)
    print(tempMat)
    print(tempout)
    tempGP <- GP_fit(tempMat,tempout,trace = T,maxit=5)
    tempests <- Snape(tempGP,xvec,xvec2,xvec3,yvec,yvec2,yvec3,zvec,zvec2,zvec3,leny,KG)
    
    xdfNew <- data.frame(xdfNew,tempests[1,])
    ydfNew <- data.frame(ydfNew,tempests[2,])
    zdfNew <- data.frame(zdfNew,tempests[3,])
    
    print(paste("run #",  i+1))
  }
  
  return(data.frame(xdfNew,ydfNew,zdfNew))
  
}




aCarbx.df[1:12,1:5][,1]


#yay <- Draco(xProGP,aCarbx.df[1:13,],aCarby.df[1:13,],aCarbz.df[1:13,],12,59,136)
#yay <- Draco(xProGP,aCarbx.df[1:51,],aCarby.df[1:51,],aCarbz.df[1:51,],50,59,20)
#yay3 <- Bellatrix(aCarbx.df[1:13,],aCarby.df[1:51,],aCarbz.df[1:13,],12,59,20)
yay102frames <- Draco(xProGP,aCarbx.df[1:101,],aCarby.df[1:101,],aCarbz.df[1:25,],101,59,20)

Dobby <- function(yay,acids,frames,NAME){
  i <- length(yay[1,]) /3
  j <- 2*i
  k <- 3*i
  
  Gx.v <- as.numeric(yay[,1:i])
  Gy.v <- as.numeric(yay[,i+1:j])
  Gz.v <- as.numeric(yay[,(j+1):k])

  Gx.df <- data.frame(t(matrix(Gx.v,acids,frames,byrow = T)))
  Gy.df <- data.frame(t(matrix(Gy.v,acids,frames,byrow = T)))
  Gz.df <- data.frame(t(matrix(Gz.v,acids,frames,byrow = T)))
  saveGIF(GmakeIMG(Gx.df,Gy.df,Gz.df),NAME)
}

Dobby(yay24frames,20,60,"~/Desktop/testGifGP24frames.gif")

#Rename eachtimewar
saveGIF(GmakeIMG(Gx.df,Gy.df,Gz.df),"~/Desktop/testGifGP50train.gif")


##CAS thesis script
# disti <- c()
# for(i in 0:101){
#   x1 <- ssn$x[1408*i + 1]
#   y1 <- ssn$y[1408*i + 1]
#   z1 <- ssn$z[1408*i + 1]
#   
#   x2 <- ssn$x[1408*(i + 1)]
#   y2 <- ssn$y[1408*(i + 1)]
#   z2 <- ssn$z[1408*(i + 1)]
#   
#   temp <- c()
#   disti <- c(disti,sqrt((x1-x2)^2 + (y1 -y2)^2 + (z1 - z2)^2))
# }
# 
# mean(disti)
# #Angstrums apart
# 143616/1408
# 
# 
# 1408 , 1
# 
# 102 times
#

# 136 different amino acids


