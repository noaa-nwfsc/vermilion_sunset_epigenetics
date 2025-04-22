# age prediction methylation
setwd("~/Documents/Anita_Wray/methylation_analysis/age_prediction")
meth <- read.csv("methylation_file_all_loci_all_samples_with_multiplex.csv",header=TRUE)
depth <- read.csv("coverage_file_all_loci_all_samples_with_multiplex.csv",header=TRUE)
age <- read.table("ages.txt",header=TRUE)
indDepth <- rowMeans(depth[,3:ncol(depth)],na.rm=TRUE)

minDepth <- 10 # minimum read depth to keep methylation score
maxMiss <- 0.5 # maximum locus missingness tolerated

# fix the locus names
library(stringr)
ids <- colnames(meth[3:ncol(meth)])
ids2 <- rep(NA,length(ids))
for(i in 1:length(ids)){
  splitID <- strsplit(ids[i],split=".",fixed=TRUE)
  id2start <- paste(splitID[[1]][1:(length(splitID[[1]])-1)],"-",collapse="",sep="")
  ids2[i] <- paste(id2start,splitID[[1]][length(splitID[[1]])],sep="")
}

colnames(meth)[3:ncol(meth)] <- ids2
colnames(meth)[1] <- "locus_num"
colnames(depth)[3:ncol(meth)] <- ids2
colnames(depth)[1] <- "locus_num"

# reorder age so to be consistent with the order of the methylation data
newAge <- NULL
for(i in 1:length(ids2)){
  newAge <- rbind(newAge,age[which(age$Specimen.Num == ids2[i]),])
}

#############################################################
# filter methylation scores based on read depth
#############################################################
for(i in 3:ncol(meth)){
  thisDepth <- depth[,i]
  if(sum(thisDepth < minDepth,na.rm=TRUE) > 0){
    print(sum(thisDepth < minDepth,na.rm=TRUE)/length(thisDepth))
    meth[which(thisDepth < minDepth),i] <- NA
  }
}

#############################################
# format glmnet input data and load glmnet
#############################################
ageVec <- newAge$Ages
predMat <- as.data.frame(t(meth[,3:ncol(meth)]))
miss <- colSums(is.na(predMat) == TRUE)/nrow(predMat)
rem <- which(round(miss,digits=2) > maxMiss)
predMat2 <- predMat[,-rem]

pVec <- NULL
# plot the data 1
par(mfrow=c(3,3),mar=c(4,6,1,1))
for(i in 1:ncol(predMat2)){
  plot(predMat2[,i],ageVec,main=rownames(predMat2)[i],xlab="methylation",ylab="age")
  mod <- lm(ageVec~predMat2[,i])
  abline(mod)
  r2 <- as.character(round(summary(mod)$r.squared,digits=2))
  p <- as.character(round(summary(mod)$coefficients[2,4],digits=3))
  pVec <- c(pVec,p)
  legend(x=min(predMat2[,i],na.rm=TRUE),y=60,legend=c(paste("r2 = ",r2,sep=""),
      paste("P = ",p,sep="")),xjust=FALSE,yjust=FALSE,bty="n")
}
  
# # plot the data 2
# for(i in 10:ncol(predMat2)){
#   plot(predMat2[,i],ageVec,main=rownames(predMat2)[i],xlab="methylation",ylab="age")
#   mod <- lm(ageVec~predMat2[,i])
#   abline(mod)
#   r2 <- as.character(round(summary(mod)$r.squared,digits=2))
#   p <- as.character(round(summary(mod)$coefficients[2,4],digits=3))
#   pVec <- c(pVec,p)
#   legend(x=min(predMat2[,i],na.rm=TRUE),y=60,legend=c(paste("r2 = ",r2,sep=""),
#        paste("P = ",p,sep="")),xjust=FALSE,yjust=FALSE,bty="n")
# }




################################################################################################
# keep only loci that seem to be correlated with age
################################################################################################
library(glmnet)

keep <- sort(which(as.numeric(pVec) <= 0.1))
methFilt <- predMat2[,keep]      # filtered methylation data

meth_imput <- makeX(methFilt,na.impute=TRUE)  # imputed methylation data

####################################################################################################
# K-fold cross validation with cv.glmnet. Imputation is necessary as missing values are not allowed
####################################################################################################
# make the model
mod <- cv.glmnet(x=meth_imput,y=ageVec,nfolds=20)

# predicted values for all individuals
pred <- predict(mod, newx = meth_imput, s = "lambda.min")

# plot the results
library(scales)
plot(ageVec,pred,xlab="Real age",ylab="Predicted age",main="glmnet k-fold cross validation",col=alpha("blue",alpha=0.4),pch=16,cex=1.5)

lines(c(0,100),c(0,100))
####################################################################################################
# K-fold cross validation with cv.glmnet. Imputation is necessary as missing values are not allowed
####################################################################################################


keep <- sort(which(as.numeric(pVec) <= 0.4))
methFilt <- predMat2[,keep]      # filtered methylation data

meth_imput <- makeX(methFilt,na.impute=TRUE)  # imputed methylation data
library(glmnet)


# make the model
mod <- cv.glmnet(x=meth_imput,y=ageVec,nfolds=20)

# predicted values for all individuals
pred <- predict(mod, newx = meth_imput, s = "lambda.min")

# plot the results
library(scales)
plot(ageVec,pred,xlab="Real age",ylab="Predicted age",main="glmnet k-fold cross validation",col=alpha("blue",alpha=0.4),pch=16,cex=1.5)


##############################################
# multiple linear regression with imputation
##############################################
mod <- lm(ageVec~meth_imput[,1] + meth_imput[,2] + meth_imput[,3] + meth_imput[,4] + meth_imput[,5] + meth_imput[,6])
hist(mod$residuals)
# leave-one-out using only loci that are individually statistically significantly (P <= 0.05) correlated with age

preds <- rep(NA,length(ageVec))
for(i in 1:length(preds)){
  # leave one out
  thisAge <- ageVec[-i]
  thisMeth <- meth_imput[-i,]
  
  # make the model
  mod <- lm(thisAge~thisMeth[,1] + thisMeth[,2] + thisMeth[,3] + thisMeth[,4] + thisMeth[,5] + thisMeth[,6])
  intercept <- as.numeric(mod$coefficients[1])
  coefs <- as.numeric(mod$coefficients[2:7])
  
  # predict the age of the one left out
  newVals <- as.numeric(meth_imput[i,])
  preds[i] <- intercept + sum(coefs*newVals)
}

plot(ageVec,preds,xlab="True age",ylab="Predicted age",col=alpha("blue",alpha=0.4),pch=16,cex=1.5,main="Leave-one-out with lm models")



##############################################
# multiple linear regression without imputation
##############################################

keep <- sort(which(as.numeric(pVec) <= 0.1))
methFilt <- predMat2[,keep]      # filtered methylation data
meth_imput <- methFilt          # no imputation this time


mod <- lm(ageVec~meth_imput[,1] + meth_imput[,2] + meth_imput[,3] + meth_imput[,4] + meth_imput[,5] + meth_imput[,6])
hist(mod$residuals)
# leave-one-out using only loci that are individually statistically significantly (P <= 0.05) correlated with age

preds <- rep(NA,length(ageVec))
for(i in 1:length(preds)){
  # leave one out
  thisAge <- ageVec[-i]
  thisMeth <- meth_imput[-i,]
  
  # make the model
  mod <- lm(thisAge~thisMeth[,1] + thisMeth[,2] + thisMeth[,3] + thisMeth[,4] + thisMeth[,5] + thisMeth[,6])
  intercept <- as.numeric(mod$coefficients[1])
  coefs <- as.numeric(mod$coefficients[2:7])
  
  # predict the age of the one left out
  newVals <- as.numeric(meth_imput[i,])
  preds[i] <- intercept + sum(coefs*newVals)
}

plot(ageVec,preds,xlab="True age",ylab="Predicted age",col=alpha("blue",alpha=0.4),pch=16,cex=1.5,main="Leave-one-out with lm models")






##############################################
# multiple linear regression without imputationm using only loci with p-value 0.05 or less
##############################################

keep <- sort(which(as.numeric(pVec) <= 0.05))
methFilt <- predMat2[,keep]      # filtered methylation data
meth_imput <- methFilt          # no imputation this time


mod <- lm(ageVec~meth_imput[,1] + meth_imput[,2] + meth_imput[,3] + meth_imput[,4] + meth_imput[,5] + meth_imput[,6])
hist(mod$residuals)
# leave-one-out using only loci that are individually statistically significantly (P <= 0.05) correlated with age

preds <- rep(NA,length(ageVec))
for(i in 1:length(preds)){
  # leave one out
  thisAge <- ageVec[-i]
  thisMeth <- meth_imput[-i,]
  
  # make the model
  mod <- lm(thisAge~thisMeth[,1] + thisMeth[,2] + thisMeth[,3] + thisMeth[,4] + thisMeth[,5] + thisMeth[,6])
  intercept <- as.numeric(mod$coefficients[1])
  coefs <- as.numeric(mod$coefficients[2:6])
  
  # predict the age of the one left out
  newVals <- as.numeric(meth_imput[i,])
  preds[i] <- intercept + sum(coefs*newVals)
}

plot(ageVec,preds,xlab="True age",ylab="Predicted age",col=alpha("blue",alpha=0.4),pch=16,cex=1.5,main="Leave-one-out with lm models")



#####################################################
# multiple linear regression with log transformation
#####################################################

multRegMod <- lm(ageVec~meth_imput[,1] + meth_imput[,2] + meth_imput[,3] + meth_imput[,4] + meth_imput[,5] + meth_imput[,6] + meth_imput[,7])


# leave-one-out using only loci that are individually statistically significantly (P <= 0.05) correlated with age

preds <- rep(NA,length(lmAge))
for(i in 1:length(preds)){
  # leave one out
  thisAge <- log(ageVec[-i]+0.001)
  thisMeth <- meth_imput[-i,]
  
  # make the model
  mod <- lm(thisAge~thisMeth[,1] + thisMeth[,2] + thisMeth[,3] + thisMeth[,4] + thisMeth[,5] + thisMeth[,6] + thisMeth[,7])
  intercept <- as.numeric(mod$coefficients[1])
  coefs <- as.numeric(mod$coefficients[2:7])
  # predict the age of the one left out
  newVals <- as.numeric(lmDat[i,])
  preds[i] <- exp(intercept + sum(coefs*newVals))
}

plot(lmAge,preds,xlab="True age",ylab="Predicted age",col=alpha("blue",alpha=0.4),pch=16,cex=1.5,main="Leave-one-out with lm models")
