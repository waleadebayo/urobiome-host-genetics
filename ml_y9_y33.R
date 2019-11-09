library(ape,rJava,bartMachine,randomForest,caret)

load("abund.RData")
 library(ape)
 y9<- extract.clade(imp20tr,"y9")
y9<-y9$tip.label
 
 cmon2 <- intersect(colnames(totuImp), y9)
 y9CLRut <- totuImp[,cmon2]
 cmon <- intersect(rownames(y9CLRut), rownames(UTI_subclass))
 y9CLRut <- y9CLRut[cmon,]
 UTI_subclass<-UTI_subclass[cmon,]
 y9CLRut$Class <- UTI_subclass$UTI_subclass
 y9CLRut <- y9CLRut[complete.cases(y9CLRut$Class),]
 y9CLRut$Class <- as.factor(y9CLRut$Class)
#colnames(y9CLRut)<-c("Dialister sp","Peptoniphilus sp","Campylobacter sp","Corynbacterium sp","Class")
set.seed(159)
inTrain3 <- createDataPartition(y=y9CLRut$Class,p=.75,list=F)
training3 <-y9CLRut[ inTrain3,]
testing3 <-y9CLRut[ -inTrain3,]
set.seed(159)
ctrl <- trainControl(method="repeatedcv",repeats = 10,classProbs=T,summaryFunction = prSummary, sampling = "up")
plsFit <- train(Class ~ ., data=training3, method="pls", tuneLength=15,metric="ROC",preProc=c("center","scale"),trControl=ctrl,importance=T)
plsFit
plsClasses <- predict(plsFit, newdata = testing3)
confusionMatrix(data = plsClasses, testing3$Class)

#####
#ncomp  AUC        Precision  Recall     F        
#1      0.6483809  0.6399174  0.5571494  0.5800237
#2      0.6592804  0.6472130  0.5977241  0.5920469
#3      0.6567592  0.6501921  0.6239080  0.6106583

#AUC was used to select the optimal model using the
#largest value.
#The final value used for the model was ncomp = 2.
##Confusion Matrix and Statistics

#Reference
#Prediction  none recurrent
#none        78        38
#recurrent   20        12

#Accuracy : 0.6081          
#95% CI : (0.5246, 0.6872)
#No Information Rate : 0.6622          
#P-Value [Acc > NIR] : 0.9289          

#Kappa : 0.0394          

#Mcnemar's Test P-Value : 0.0256          

#Sensitivity : 0.7959          
#Specificity : 0.2400          
#Pos Pred Value : 0.6724          
#Neg Pred Value : 0.3750          
#Prevalence : 0.6622          
#Detection Rate : 0.5270          
#Detection Prevalence : 0.7838          
#Balanced Accuracy : 0.5180          
#'Positive' Class : none
plot(plsFit)
plot(varImp(plsFit))