##Using random Forest machine learning lagorithm in caret package within R 

# For python scripts, please check, https://pycaret.readthedocs.io/en/latest/tutorials.html#classification

####The data below was on urinary incontinence, with a sample size of ~1800 with microbial data as features
# For the current project , prepare your data input 
#For help on reading in and preparing your data input, rf check the Intro to R in the emails sent to you

#ensure clinical phenotype data are aligned to omics data and combine: 
#add to the omics data with column of the phenotype  named Class, e.g.treatment/control, disease/nodisease
metsdata=read.delim('metsdata.txt', row.names=1)
phenodata=read.delim('phenodata.txt',row.names = 1)
identical(rownames(metsdata),rownames(phenodata))#should be TRUE
qq=metsdata ; qq$Class=phenodata$GroupsoOfInterest
dim(qq)
table(qq$Class)
# In this case 'Class' includes two groups of samples and what we have here is a binary classification problem i.e. predicting just 2 groups

#The data is relatively 'big' and but not too big, so we split into test set and validation set
#by seprating out 25% of the data randomly but with data coverage
insub <- createDataPartition(y=qq$Class,p=.75,list=F)
datasub <-qq[ insub,]
valsub <-qq[ -insub,]
dim(datasub)
dim(valsub)

#we set a seed for reproducibility in computer-generated random seeds
set.seed(1659)

# create models, model parameters, 
#Here, its a five-fold cross validated model with 10 repeats model that prints what its doing
ctrl <- trainControl(method="repeatedcv", repeats = 10,number=5,verboseIter=T)
#optional #classProbs = T
#Here we use random forest algorithm specifically from randomForest package with data centered around mean and quite granular
rfFit <- train(Class ~., data=datasub,method="rf", 
               tuneLength=32,trControl=ctrl,preProcess=c("center"),importance=T)

##see output and accuracy
rfFit
###plot feature Importance
vi=varImp(rfFit)
plot(varImp(rfFit))

###plot rocs
roc_imp <- filterVarImp(x = datasub[, -ncol(datasub)], y = datasub$Class)

##apply model on the validation set i.e. the 25% of data that was separated initially
rfClasses <- predict(rfFit, newdata = valsub)
#print and plot accuracy metrics on the validation set
confusionMatrix(data = rfClasses, factor(valsub$Class))
fourfoldplot(cc$table)

##
timestamp()
savehistory('test_machlearn_rf.R')
