
load("/Volumes/LATEST/killaPack/Genes41019.RData")
load("/Volumes/LATEST/killaPack/repbalClust.RData")
rm(twStool,twStool4,twStool5,twTax)
rm(ttwStool,ttwStool2,ttwStool3,ttwStool4)
putspp<-rep_spplevel2R1
head(colnames(putspp),3)
colnames(putspp)<-substr(colnames(putspp), 1,nchar(colnames(putspp))-9)
head(colnames(putspp),3)
cmon4<-intersect(rownames(combo2),colnames(putspp) )
head(rownames(combo2),3)
colnames(putspp)<-paste0 ("11624.",colnames(putspp),sep="")
cmon4<-intersect(rownames(combo2),colnames(putspp) )
which(!rownames(combo2)%in%cmon4)
rownames(combo2)[73]
which(!colnames(putspp)%in%cmon4)
colnames(putspp)[76]
colnames(putspp)[76]<-"11624.J40.363.78"
cmon4<-intersect(rownames(combo2),colnames(putspp) )
putAlltog<-combo2[cmon4,]
putspp<-putspp[,cmon4]
putspp<-t(putspp)
putspp<-as.data.frame(putspp)
cmon4<-intersect(rownames(putAlltog),rownames(putspp) )
putAlltog<-putAlltog[cmon4,]
putspp<-putspp[cmon4,]
head(colnames(putAlltog),3)
colnames(putAlltog)<- gsub("X","EC",colnames(putAlltog))
colnames(putAlltog)<- gsub("..","",colnames(putAlltog),fixed=T)

putspp2 <- putspp
putspp2<-apply(putspp,1,function(x) x/sum(x))
dim(putspp2)
putspp2<-t(putspp2)
sum(putspp2[1,])
sum(putspp2[50,])
library(compositions)
putspp2<-clr(putspp2)
putspp2<-as.data.frame.rmult(putspp2)
detach("package:compositions",unload=TRUE)
putAlltog<-cbind.data.frame(putAlltog,putspp2)
library(biomformat)
putr1clust<-read_biom("/Volumes/LATEST/killaPack/clean/replic_clust/R1/feature-table.biom")
putr1clust<-as.matrix(biom_data(putr1clust))
putr1clust<-as.data.frame(putr1clust)
colnames(putr1clust)<-substr(colnames(putr1clust), 1,nchar(colnames(putr1clust))-9)
colnames(putr1clust)<-paste0 ("11624.",colnames(putr1clust),sep="")
colnames(putr1clust)[76]
colnames(putr1clust)[76]<-"11624.J40.363.78"
cmon4<-intersect(rownames(putAlltog),rownames(putr1clust))
putr1clust<-t(putr1clust)
putr1clust<-as.data.frame(putr1clust)
cmon4<-intersect(rownames(putAlltog),rownames(putr1clust))
putAlltog<-putAlltog[cmon4,]
putr1clust<-putr1clust[cmon4,]
colnames(putr1clust)<-gsub("y","Cb",colnames(putr1clust))
putAlltog<-cbind.data.frame(putAlltog,putr1clust)

metad<-putAlltog[,c(212:227)]



#################
library(caret)
metad<-putAlltog[,c(212:227)]
dat<-putAlltog[,-c(212:227)]
dat$Class<-metad$Age
set.seed(159)
inTrain3 <- createDataPartition(y=dat$Class,p=.80,list=F)
training3 <-dat[ inTrain3,]
testing3 <-dat[ -inTrain3,]
ctrl <- trainControl(method="repeatedcv", repeats = 3,classProbs = T)
plsFit <- train (Class~., data=training3, method="ridge", tuneLength=15, trControl=ctrl,preProc=c("center","scale"), importance=T) 
plsFit
plsClasses <- predict(plsFit, newdata = testing3)
confusionMatrix(data = plsClasses, testing3$Class)
plot(plsFit)
plot(varImp(plsFit))
combo <- data.frame(testing3$Class, plsClasses)
colnames(combo)<- c("Class","predClass")
combo$absDiff <- abs(combo$Class - combo$predClass)
combo$Diff <- (combo$Class - combo$predClass)
mean(combo$absDiff)
sd(combo$Class)
sqrt(mean(combo$Diff*combo$Diff))
mean(combo$absDiff)


####pick varying data
asv2<-matrix(0,nrow=1,ncol=396)
datvar<- dat
for (i in datvar[,1:227]) {
  asv<-var(i)
  asv2[,i]<-asv
  
}
range(asv2[1,])
mean(asv2[1,])
asv3<-asv2[,which(colMeans(asv2)>1)]
asv4<-asv2[,which(colMeans(asv2)<(-1))]
asv5<- append(colnames(asv3), colnames(asv4))
datvar2<-dat[,which(colnames(dat)%in%asv5)]
###
datvar2<-datvar[,which(apply(datvar,2,function(x) var(x))>0.0001)]
identical(rownames(datvar2),rownames(metad))
all.equal(rownames(datvar2),rownames(metad))
datvar2$Class<-metad$Age
set.seed(159)
inTrain3 <- createDataPartition(y=datvar2$Class,p=.80,list=F)
training3 <-datvar2[ inTrain3,]
testing3 <-datvar2[ -inTrain3,]
ctrl <- trainControl(method="repeatedcv", repeats = 3)
plsFit <- train (Class~., data=training3, method="rf", 
                 trControl=ctrl, importance=T,preProcess=c("center","scale")) 
plsFit
plsClasses <- predict(plsFit, newdata = testing3)
confusionMatrix(data = plsClasses, testing3$Class)
plot(plsFit)
plot(varImp(plsFit))
combo <- data.frame(testing3$Class, plsClasses)
colnames(combo)<- c("Class","predClass")
combo$absDiff <- abs(combo$Class - combo$predClass)
combo$Diff <- (combo$Class - combo$predClass)
mean(combo$absDiff)
sd(combo$Class)
sqrt(mean(combo$Diff*combo$Diff))
mean(combo$absDiff)

