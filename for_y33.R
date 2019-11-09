library(ape)
library(caret)
cognClust<-read.csv("/Volumes/LATEST/cognClust.csv",row.names = 1,check.names = F)
y33<- extract.clade(imp20tr,"y33")
y33<-y33$tip.label

cmon2 <- intersect(colnames(totuImp), y33)
y33CLRut <- totuImp[,cmon2]
cmon <- intersect(rownames(y33CLRut), rownames(cognClust))
y33CLRut <- y33CLRut[cmon,]
cognClust<-cognClust[cmon,]
y33CLRut$Class <- cognClust$`cognPam1$clustering`
df<-colnames(y33CLRut)
y33CLRut <- y33CLRut[complete.cases(y33CLRut$Class),]
df
colnames(y33CLRut)<-c("Negativicoccus_sp", "Sutterella_sp","Bacteroides_sp","Prevotella_sp","Actinomycetaceae","Actinomyces_sp","Alloscardovia","Pepton.coxii","Class")
y33CLRut$Class<-gsub("1","Group1", y33CLRut$Class)
y33CLRut$Class<-gsub("2","Group2", y33CLRut$Class)
y33CLRut$Class<-gsub("3","Group3", y33CLRut$Class)

set.seed(159)
inTrain3 <- createDataPartition(y=y33CLRut$Class,p=.75,list=F)
training3 <-y33CLRut[ inTrain3,]
testing3 <-y33CLRut[ -inTrain3,]
set.seed(159)
ctrl <- trainControl(method="repeatedcv",repeats = 10,classProbs = T)
plsFit <- train(Class ~ ., data=training3,method="rda", tuneLength=15,trControl=ctrl,importance=T)
#add importance=T to plsFit if model needs it
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

unique(testing3$Class<
