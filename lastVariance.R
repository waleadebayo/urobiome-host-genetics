#####use everyone who has all data
####import,arrange,fill in
kID2<- read.csv("/Volumes/LATEST/kinship/allKin3.csv", row.names=1)
load("/Users/macuser/Documents/uromg/Variance/hoci2k.RData")
hoci2k$rname<-rownames(hoci2k)
rownames(hoci2k)<-hoci2k$study_id
pckin<-read.csv("/Volumes/LATEST/kinship/PCdistKinAdd.csv",row.names = 1)
sdev<-read.csv("/Volumes/LATEST/kinship/PCdistKinAddSdev.csv")
head(round(100*sdev$x^2/sum(sdev$x^2),1))
##sort them out
cmon<-intersect(rownames(kID2),rownames(pckin))
kID2<-kID2[cmon,]
pckin<-pckin[cmon,]
pckin2<-pckin
rownames(pckin2)<-kID2$study_id2
kid3<-kID2
rownames(kid3)<-kid3$study_id2
cmon2<-intersect(rownames(kid3),rownames(hoci2k))
hoci2k<-hoci2k[cmon2,]
kid3<-kid3[cmon2,]

cmon3<-intersect(rownames(hoci2k), rownames(pckin2))
hoci2k<-hoci2k[cmon3,]
pckin3<-pckin2[cmon3,]
hoci2k$genetickinPC1<-pckin3$PC1
hoci2k$genetickinPC2<-pckin3$PC2


##
compAll<- hoci2k[complete.cases(hoci2k$meno_stat),]
compAll<- compAll[complete.cases(compAll$utihis),]
compAll<- compAll[complete.cases(compAll$scqFI),]
compAll<- compAll[complete.cases(compAll$HEI),]
compAll<- compAll[complete.cases(compAll$antib),]
compAll<- compAll[complete.cases(compAll$nBirths),]
###
load("/Users/macuser/Documents/uromg/Variance/mgd2t_f2RR.RData")
##
library(vegan)
rownames(compAll)<-compAll$rname
bdivUH<-prune_samples(rownames(compAll),mgd2t_f2)
bdivUHd <-distance(bdivUH,method="bray")
hoci2kP<-sample_data(bdivUH)
hoci2kP<-as.data.frame(hoci2kP)
cmon<-intersect(rownames(hoci2kP), rownames(compAll))
compAll<-compAll[cmon,]
set.seed(1659);adonis2(bdivUHd~genetickinPC1+genetickinPC2+scAge+scUTI+scale(nBirths)+antib_use+scale(HEI)+scqFI+runSeq+reads_filt+extractionkit_lot+mastermix_lot+processing_robot+extraction_robot,data=compAll)

set.seed(1659);pem1<-adonis2(bdivUHd~genetickinPC1+genetickinPC2+scAge+scUTI+scale(nBirths)+antib_use+scale(HEI)+scqFI+runSeq+reads_filt+extractionkit_lot+mastermix_lot+processing_robot+extraction_robot,data=compAll)

###reorder
set.seed(1659);pem1<-adonis2(bdivUHd~genetickinPC1+genetickinPC2+scAge+scUTI+scale(nBirths)+antib_use+scale(HEI)+meno_stat+scqFI+runSeq+reads_filt+primer_plate+extractionkit_lot+mastermix_lot+processing_robot+extraction_robot,data=compAll)
set.seed(1659);pem2<-adonis2(bdivUHd~genetickinPC2+genetickinPC1+scAge+scUTI+scale(nBirths)+antib_use+scale(HEI)+meno_stat+scqFI+runSeq+reads_filt+primer_plate+extractionkit_lot+mastermix_lot+processing_robot+extraction_robot,data=compAll)
set.seed(1659);pem3<-adonis2(bdivUHd~genetickinPC2+scAge+genetickinPC1+scUTI+scale(nBirths)+antib_use+scale(HEI)+meno_stat+scqFI+primer_plate+extractionkit_lot+mastermix_lot+processing_robot+extraction_robot+reads_filt+runSeq,data=compAll)
set.seed(1659);pem4<-adonis2(bdivUHd~genetickinPC2+scAge+scUTI+genetickinPC1+scale(nBirths)+antib_use+scale(HEI)+meno_stat+scqFI+primer_plate+extractionkit_lot+mastermix_lot+processing_robot+extraction_robot+reads_filt+runSeq,data=compAll)
set.seed(1659);pem6<-adonis2(bdivUHd~scAge+scUTI+scale(nBirths)+antib_use+scale(HEI)+meno_stat+scqFI+genetickinPC1+genetickinPC2+extractionkit_lot+extraction_robot+primer_plate+mastermix_lot+processing_robot+reads_filt+runSeq,data=compAll)
set.seed(1659);pem5<-adonis2(bdivUHd~genetickinPC2+scAge+scUTI+scale(nBirths)+antib_use+scale(HEI)+meno_stat+scqFI+genetickinPC1+primer_plate+extractionkit_lot+mastermix_lot+processing_robot+extraction_robot+reads_filt+runSeq,data=compAll)
set.seed(1659);pem7<-adonis2(bdivUHd~scUTI+scAge+meno_stat+antib_use+scale(nBirths)+scale(HEI)+scqFI+genetickinPC1+genetickinPC2+extractionkit_lot+extraction_robot+primer_plate+mastermix_lot+processing_robot+reads_filt+runSeq,data=compAll)
set.seed(1659);pem8<-adonis2(bdivUHd~scAge+scUTI+meno_stat+scale(nBirths)+scale(HEI)+scqFI+genetickinPC1+genetickinPC2+antib_use+extractionkit_lot+extraction_robot+primer_plate+mastermix_lot+processing_robot+reads_filt+runSeq,data=compAll)
 set.seed(1659);pem9<-adonis2(bdivUHd~scAge+scUTI+meno_stat+scale(nBirths)+genetickinPC1+scale(HEI)+scqFI+genetickinPC2+antib_use+extractionkit_lot+extraction_robot+primer_plate+mastermix_lot+processing_robot+reads_filt+runSeq,data=compAll)
 set.seed(1659);pem9<-adonis2(bdivUHd~meno_stat+scale(nBirths)+genetickinPC1+scUTI+scale(HEI)+scqFI+genetickinPC2+antib_use+scAge+extractionkit_lot+extraction_robot+primer_plate+mastermix_lot+processing_robot+reads_filt+runSeq,data=compAll)
 set.seed(1659);pem10<-adonis2(bdivUHd~meno_stat+scale(nBirths)+genetickinPC1+scUTI+scale(HEI)+scqFI+genetickinPC2+antib_use+scAge+extractionkit_lot+extraction_robot+primer_plate+mastermix_lot+processing_robot+reads_filt+runSeq,data=compAll)
 set.seed(1659);pem11<-adonis2(bdivUHd~meno_stat+genetickinPC1+scale(nBirths)+scUTI+scale(HEI)+scqFI+genetickinPC2+antib_use+scAge+extractionkit_lot+extraction_robot+primer_plate+mastermix_lot+processing_robot+reads_filt+runSeq,data=compAll)
 set.seed(1659);pem12<-adonis2(bdivUHd~meno_stat+scqFI+scale(nBirths)+genetickinPC1+scUTI+scale(HEI)+genetickinPC2+antib_use+scAge+extractionkit_lot+extraction_robot+primer_plate+mastermix_lot+processing_robot+reads_filt+runSeq,data=compAll)
 set.seed(1659);pem13<-adonis2(bdivUHd~meno_stat+scqFI+scale(nBirths)+scUTI+scale(HEI)+genetickinPC2+antib_use+scAge+extractionkit_lot+extraction_robot+runSeq+genetickinPC1+primer_plate+mastermix_lot+processing_robot+reads_filt,data=compAll)
 set.seed(1659);pem14<-adonis2(bdivUHd~meno_stat+scqFI+scale(nBirths)+scUTI+scale(HEI)+antib_use+scAge+extractionkit_lot+extraction_robot+runSeq+genetickinPC2+genetickinPC1+primer_plate+mastermix_lot+processing_robot+reads_filt,data=compAll)
 set.seed(1659);pem15<-adonis2(bdivUHd~scAge+scqFI+scUTI+antib_use+meno_stat+scale(nBirths)+scale(HEI)+genetickinPC1+extractionkit_lot+extraction_robot+runSeq+genetickinPC2+primer_plate+mastermix_lot+processing_robot+reads_filt,data=compAll)
 set.seed(1659);pem16<-adonis2(bdivUHd~scAge+genetickinPC1+scqFI+scUTI+genetickinPC2+antib_use+meno_stat+scale(nBirths)+scale(HEI)+extractionkit_lot+extraction_robot+runSeq+primer_plate+mastermix_lot+processing_robot+reads_filt,data=compAll)
 set.seed(1659);pem17<-adonis2(bdivUHd~scAge+scqFI+scUTI+antib_use+meno_stat+scale(nBirths)+scale(HEI)+genetickinPC1+extractionkit_lot+extraction_robot+genetickinPC2+primer_plate+mastermix_lot+processing_robot+reads_filt+runSeq,data=compAll)
 set.seed(1659);pem18<-adonis2(bdivUHd~scale(nBirths)+scale(HEI)+scAge+scqFI+antib_use+meno_stat+scUTI+genetickinPC1+extractionkit_lot+extraction_robot+genetickinPC2+primer_plate+mastermix_lot+processing_robot+reads_filt+runSeq,data=compAll)
set.seed(1659);pem19<-adonis2(bdivUHd~scale(nBirths)+scale(HEI)+scAge+scqFI+antib_use+meno_stat+genetickinPC1+scUTI+extractionkit_lot+extraction_robot+genetickinPC2+primer_plate+mastermix_lot+processing_robot+reads_filt+runSeq,data=compAll)
 set.seed(1659);pem20<-adonis2(bdivUHd~antib_use+scale(nBirths)+genetickinPC1+scale(HEI)+scAge+scqFI+meno_stat+scUTI+extraction_robot+genetickinPC2+extractionkit_lot+primer_plate+processing_robot+reads_filt+mastermix_lot+runSeq,data=compAll)
set.seed(1659);pem21<-adonis2(bdivUHd~genetickinPC1+scale(nBirths)+scale(HEI)+scqFI+scUTI+scAge+antib_use+meno_stat+extraction_robot+genetickinPC2+extractionkit_lot+primer_plate+processing_robot+reads_filt+mastermix_lot+runSeq,data=compAll)
 pemAll <- rbind(pem1,pem2,pem3,pem4,pem5,pem6,pem7,pem8,pem9,pem10,pem11,pem12,pem13,pem14,pem15,pem16,pem17,pem18,pem19,pem20,pem21)
 save(bdivUHd,file="pemBray533D.RData")
 save(compAll,file="pemBray533S.RData")

 ####capture all data and plot
 kID2<- read.delim("/Volumes/LATEST/kinship/hociKin.txt", sep= " ")
 pckin<-read.csv("/Volumes/LATEST/kinship/PCdistKin.csv",row.names = 1)
 sdev<-read.csv("/Volumes/LATEST/kinship/PCdistKinSdev.csv")
 head(round(100*sdev$x^2/sum(sdev$x^2),1))
 cmon<-intersect(rownames(kID2),rownames(pckin))
 kID2<-kID2[cmon,]
 pckin<-pckin[cmon,]
 pckin2<-pckin
 head(rownames(pckin2),3)
 head(rownames(kID2),3)
 rownames(pckin2)<-kID2$study_id
 kid3<-kID2
 rownames(kid3)<-kid3$study_id
 hoci2k<-kid3
 hoci2k$genetickinPC1<-pckin2$PC1
 hoci2k$genetickinPC2<-pckin2$PC2
  compAll2<- hoci2k[complete.cases(hoci2k$meno_stat),]
  compAll2<- compAll2[complete.cases(compAll2$utihis),]
  compAll2<- compAll2[complete.cases(compAll2$scqFI),]
  compAll2<- compAll2[complete.cases(compAll2$HEI),]
  compAll2<- compAll2[complete.cases(compAll2$antib),]
 compAll2<- compAll2[complete.cases(compAll2$nBirths),]
   rownames(compAll2)<-compAll2$ucsd_Name
 bdivUH<-prune_samples(rownames(compAll2),mgd2t_f2)
  bdivUHd <-distance(bdivUH,method="bray")
  hoci2kP<-sample_data(bdivUH)
  hoci2kP<-as.data.frame(hoci2kP)
  cmon<-intersect(rownames(hoci2kP), rownames(compAll2))
  compAll2<-compAll2[cmon,]
  set.seed(1659);adonis2(bdivUHd~genetickinPC1+runSeq+reads_filt+extractionkit_lot+mastermix_lot+processing_robot+extraction_robot,data=compAll2) 
 set.seed(1659);adonis2(bdivUHd~genetickinPC1+scAge+scUTI+meno_stat+scqFI+antib_use+scale(nBirths)+runSeq+reads_filt+extractionkit_lot+mastermix_lot+processing_robot+extraction_robot,data=compAll2)
 save(bdivUHd,file="pemBray545D.RData")
 save(compAll2,file="pemBray545S.RData")
 
 hociHEI<-read.csv("/Users/macuser/Desktop/divers/hociHEI.csv",row.names = 1)
head(rownames(hociHEI),4)
compAll3<-compAll2
rownames(compAll3)<-compAll3$study_id
cmon2<-intersect(rownames(compAll3), rownames(hociHEI))
hociHEI<-hociHEI[cmon2,]
compAll3<-compAll3[cmon2,]
compAll3$HEI<-hociHEI$HEI
compAll3$surg<-hociHEI$surgical_pertubation
compAll3$surgCaes<-hociHEI$surg_hist_withCaes
compAll3<- compAll3[complete.cases(compAll3$HEI),] 
compAll3<- compAll3[complete.cases(compAll3$surg),] 

#########use original 835 but many do not have data
kID2<- read.delim("/Volumes/LATEST/kinship/hociKin.txt", sep= " ")
pckin<-read.csv("/Volumes/LATEST/kinship/PCdistKin.csv",row.names = 1)
sdev<-read.csv("/Volumes/LATEST/kinship/PCdistKinSdev.csv")
head(round(100*sdev$x^2/sum(sdev$x^2),1))
cmon<-intersect(rownames(kID2),rownames(pckin))
kID2<-kID2[cmon,]
pckin<-pckin[cmon,]
pckin2<-pckin
head(rownames(pckin2),3)
head(rownames(kID2),3)
rownames(pckin2)<-kID2$study_id
kid3<-kID2
rownames(kid3)<-kid3$study_id
hoci2k<-kid3
hoci2k$genetickinPC1<-pckin2$PC1
hoci2k$genetickinPC2<-pckin2$PC2
compAll<- hoci2k[complete.cases(hoci2k$age),]
rownames(compAll)<-compAll$ucsd_Name
bdivUH<-prune_samples(rownames(compAll),mgd2t_f2)
bdivUH
bdivUHd <-distance(bdivUH,method="bray")
bdivUHd <-distance(bdivUH,method="wunifrac")
hoci2kP<-sample_data(bdivUH)
hoci2kP<-as.data.frame(hoci2kP)
cmon<-intersect(rownames(hoci2kP), rownames(compAll))
compAll<-compAll[cmon,]
set.seed(1659);adonis2(bdivUHd~genetickinPC1+runSeq+reads_filt+extractionkit_lot+mastermix_lot+processing_robot+extraction_robot,data=compAll)

set.seed(1659);adonis2(bdivUHd~genetickinPC1+genetickinPC2+scAge+scUTI+scale(nBirths)+antib_use+scale(HEI)+scqFI+runSeq+reads_filt+extractionkit_lot+mastermix_lot+processing_robot+extraction_robot,data=compAll)

#######atchinsonComposition
load("/Users/macuser/Desktop/package/rds/mgd2t.RData")
bdivAt<-prune_samples(rownames(compAll),mgd2t)
bdivAt<-prune_taxa(taxa_names(bdivUH),bdivAt)
otuSub2 <- otu_table(bdivAt) 
dim(otuSub2) ; 
otuSub2<-t(otuSub2)
dim(otuSub2)
otuSub2<-as.data.frame(otuSub2)

d1<-compositions::acomp(otuSub2,MNAR = 0)
dim(d1)
d1<-dist((d1))

hoci2kP<-sample_data(bdivAt)
hoci2kP<-as.data.frame(hoci2kP)
dim(hoci2kP)
cmon<-intersect(rownames(hoci2kP), rownames(compAll))
compAll<-compAll[cmon,]

set.seed(1659);adonis2(d1~genetickinPC1+genetickinPC2+scAge+antib_use+as.factor(meno_stat)+scUTI+scale(nBirths)+antib_use+scale(HEI)+scqFI+runSeq+reads_filt+extractionkit_lot+mastermix_lot+processing_robot+extraction_robot,data=compAll)

###

set.seed(1659);pemMWU178 <- adonis2(bdivUHd~geneticPC1+geneticPC2+geneticPC3+geneticPC4+geneticPC5+geneticPC6+geneticPC7+geneticPC8+geneticPC9+geneticPC10+scAge+antib_use+scUTI+scale(nBirths)+antib_use+scale(HEI)+scqFI+runSeq+reads_filt+extractionkit_lot+mastermix_lot+processing_robot+extraction_robot,data=compAll)
save(pemMWU178,file="pemMWU178.RData")
save(pemMWU302,file="pemMWU302.RData")
save(pemMWU302D,file="pemMWU302D.RData")
save(pemMWU302S,file="pemMWU302S.RData")
load("pemMWU302.RData")
View(pemMWU302)
pemMWU178<-pemMWU302[c(1:13),]
View(pemMWU178)
pemMWU178<-pemMWU302[c(1:12),]
View(pemMWU178)
sum(pemMWU178$R2)
pemMWU178$percentage_contribution<-(pemMWU178$R2 / 0.06575365)
rownames(pemMWU178)<- c("GeneticPC1","GeneticPC2","GeneticPC3","GeneticPC4","GeneticPC5","GeneticPC6","GeneticPC7","GeneticPC8","GeneticPC9","GeneticPC10","Age","Antibiotics usage")
pemMWU178$names<-rownames(pemMWU178)
pemMWU178$names<-reorder(pemMWU178$names, pemMWU178$R2)
barchart(pemMWU178$percentage_contribution~pemMWU178$names, data=pemMWU178,ylab="Contribution to variance(R2=0.11)",xlab="urinary microbiome factors",scales=list(x=list(rot=35)),main="microbiome variation(n=302)")
barchart(pemMWU178$percentage_contribution~pemMWU178$names, data=pemMWU178,ylab="Contribution to variance(R2=0.11)",xlab="urinary microbiome factors",scales=list(x=list(rot=35)))
##remove2
pemM<-read.csv("pemM.csv",row.names = 1)
pemM<-pemM[-2,]
library(lattice) ; library(ggplot2)
pemM$names<-reorder(pemM$names, pemM$R2)
##sum(pemM$R2);pemM$average / 0.03844487
pemM$percentage_contribution <- (pemM$percentage_contribution*100)
barchart(pemM$percentage_contribution~pemM$names, data=pemM,ylab="Contribution to variance (R2=10.3%)",xlab="urinary microbiome factors",scales=list(x=list(rot=35))