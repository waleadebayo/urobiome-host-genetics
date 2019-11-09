
library(factoextra)
library(ggplot2)
library(vegan)
pca<-prcomp(dist((geneSig)))
fviz_eig(pca,addlabels=TRUE) + th
#plot2: pca. What can we conclude from the plot?
dataGG = data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2])
percentVar <- round(100*pca$sdev^2/sum(pca$sdev^2),1)
ggplot(dataGG,aes(x=dataGG$PC1,y=dataGG$PC2)) +  geom_point(size=6,shape=1) + labs(x = paste0("PC1, VarExp:", round(percentVar[1],4)), y = paste0("PC2, VarExp:", round(percentVar[2],4))) + th
head(rownames(dataGG),3)
rownames(dataGG)<-substr(rownames(dataGG), 1,nchar(rownames(dataGG))-9)
rownames(dataGG)<-paste0 ("11624.",rownames(dataGG),sep="")
rownames(dataGG)[73]
rownames(dataGG)[73]<-"11624.J40.363.78"
cmon3<-intersect(rownames(dataGG), rownames(putspp2))
dataGG<-dataGG[cmon3,]
dataSPP<-putspp2[cmon3,]
####
clrInd <- dataSPP
clrInd$PC1<-dataGG$PC1
cmon4<-intersect(rownames(clrInd), rownames(sigClustMg))
clrInd<-clrInd[cmon4,] ; sigClustMg<-sigClustMg[cmon4,]
clrInd$scUTI<-sigClustMg$scUTI
clrInd1 <- clrInd
modInd <- lm(PC1 ~ ., data=clrInd1, na.action = na.omit)
modInd2<-summary(modInd)$coefficients
modInd3<-as.data.frame(modInd2)
modInd3<-modInd3[which(modInd3$`Pr(>|t|)`<0.05),]

################## many in geneSig are highly skewed, so normalise compositionally and re-run
gps<-geneSig
rownames(gps)<-substr(rownames(gps), 1,nchar(rownames(gps))-9)
rownames(gps)<-paste0 ("11624.",rownames(gps),sep="")
rownames(gps)[73]
rownames(gps)[73]<-"11624.J40.363.78"
cmon5<-intersect(rownames(gps),rownames(clrInd))
cmon5<-intersect(rownames(clrInd),rownames(gps))
gps<-gps[cmon5,]
library(compositions)
gps<-clr(gps)
gps<-as.data.frame.rmult(gps)
detach("package:compositions",unload=TRUE)
pca<-prcomp(gps)
fviz_eig(pca,addlabels=TRUE) + th
##########continue as before


#signifi
modInd <- lm(PC1 ~ ., data=clrInd1, na.action = na.omit)
modInd2<-summary(modInd)$coefficients
modInd3<-as.data.frame(modInd2)
modInd3<-modInd3[which(modInd3$`Pr(>|t|)`<0.1),]
#rsquared
 sumb<-clrInd1[,1:94]
 sumb2<-sumb
 for (i in 1:ncol(sumb)){ 
   b2<-sumb[,i]
   r2vals<-summary(lm(clrInd1[,"PC1"]~b2))
   r2v=r2vals$r.squared
   sumb2[,i]=r2v
   }

 
sumb2<-sumb2[1,]
sumb2[2,]<-colnames(sumb2) 
sumb2<-t(sumb2)
sumb2<-as.data.frame(sumb2)
colnames(sumb2)<-c("Rsq","Species")
View(sumb2)
rownames(sumb2)<-sumb2$Species
cmon2<-intersect(rownames(sumb2), rownames(modInd3))
sumb3<-sumb2[cmon2,]
modInd3<-modInd3[cmon2,]
#sumb3$Rsq<-as.numeric(sumb3$Rsq)
#sumb3$names<-reorder(sumb3$Species,-sumb3$Rsq)
sumb3$betaCoef<-modInd3$Estimate
sumb3$betaCoef2<-abs(sumb3$betaCoef)
library(lattice)
barchart(betaCoef2~names, data=sumb3,
         scales=list(x=list(rot=35)),
         ylab="Effect size (on age-related genes)",
         xlab="Species")
#########
###for Rsq
barchart(betaCoef2~names, data=sumb3,
         scales=list(x=list(rot=90)),ylab="Impact on age-related microbial genes",
         xlab="Effective Species")
#also for beta but on other PCs
barchart(betaCoef2~names, data=sumb3,
         ylab="Effect size (on age-related genes)[PC4]",
         xlab="Species")
########
#sumb4$Species<-reorder(sumb4$Species,-sumb4$betaCoef2)
sumb4<-sumb3[order(sumb3$betaCoef2,decreasing=T),]
sumb4<-sumb4[1:20,]
sumb4$names<-reorder(sumb4$Species,-sumb4$betaCoef2)
barchart(betaCoef2~names, data=sumb4,
         scales=list(x=list(rot=45,cex=0.6)),
         ylab="Effect size (on age-related genes)",
         xlab="Species")
####adj for the choice of samples,discordance, family,zygosity
forwgs2<-read.csv("/Users/macuser/Downloads/for_wgs_2.csv",row.names = 1)

  rownames(forwgs2)<-paste0 ("11624.",rownames(forwgs2),sep="")

  rownames(forwgs2)[105] <-"11624.J40.363.78"
cmonn<-intersect(rownames(forwgs2),rownames(clrInd))
forwgs2<-forwgs2[cmonn,]
um<-clrInd[cmonn,]
um$discord<-forwgs2$Discordance
um$fam_id<-forwgs2$fam_id

  um$Zygosity<-forwgs2$Zygosity
clrInd1 <- um
modInd <- lm(PC1 ~ ., data=clrInd1, na.action = na.omit)
modInd2<-summary(modInd)$coefficients
modInd3<-as.data.frame(modInd2)
modInd3<-modInd3[which(modInd3$`Pr(>|t|)`<0.05),]

  View(modInd3)
#rsquared
   sumb<-clrInd1[,1:94]
sumb2<-sumb
for (i in 1:ncol(sumb)){ 
      b2<-sumb[,i]
      r2vals<-summary(lm(clrInd1[,"PC1"]~b2))
      r2v=r2vals$r.squared
     sumb2[,i]=r2v
      }
sumb2<-sumb2[1,]
sumb2[2,]<-colnames(sumb2)
sumb2<-t(sumb2)
sumb2<-as.data.frame(sumb2)
colnames(sumb2)<-c("Rsq","Species")
rownames(sumb2)<-sumb2$Species
cmon2<-intersect(rownames(sumb2), rownames(modInd3))
sumb3<-sumb2[cmon2,]
modInd3<-modInd3[cmon2,]
#sumb3$Rsq<-as.numeric(sumb3$Rsq)
  #sumb3$names<-reorder(sumb3$Species,-sumb3$Rsq)
  sumb3$betaCoef<-modInd3$Estimate
sumb3$betaCoef2<-abs(sumb3$betaCoef)
########
#sumb4$Species<-reorder(sumb4$Species,-sumb4$betaCoef2)
  sumb4<-sumb3[order(sumb3$betaCoef2,decreasing=T),]
sumb4$names<-reorder(sumb4$Species,-sumb4$betaCoef2)
barchart(betaCoef2~names, data=sumb4,
                     scales=list(x=list(rot=45)),
                     ylab="Effect size (on age-related genes)",
                     xlab="Species")
