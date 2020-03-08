
#distance of pair
load("/Users/adewaleadebayo/Desktop/postproc/re_metad/hoci_paired.RData")
sam<- read.csv("/Users/adewaleadebayo/Desktop/postproc/RthGd/hociHEI.csv",row.names = 1)
#load("/Users/adewaleadebayo/Desktop/postproc/compareUrSt/mgd2t_f2NV.RData")
#mgd2t_f2 <- phyloseq::prune_samples(rownames(sampdat),mgd2t_f2)
ord<-phyloseq::ordinate(mgd2t_f2,method = "PCoA",distance = "bray",formula=~Family_ID)
#phyloseq::plot_ordination(mgd2t_f2,ordination = ord)
ord2 <- phyloseq::distance(mgd2t_f2,method = "bray")
##
ss<-rownames(as.matrix(ord2))
 ss <- as.data.frame(ss)
 rownames(ss)<-ss$ss
 ss$names<-ss$ss
twins <- sam
twinwide <- fast.reshape(twins, id="fam_id",varying= "study_id")
twinwide <- twinwide[complete.cases(twinwide$study_id1, twinwide$study_id2),]
x <- twinwide$fam_id
sampdat <- twins[which(twins$fam_id%in%x),]
 cmon <-intersect(rownames(ss),rownames(sampdat))
 ss$Zygosity <- sampdat$Zygosity
 ss$fam_id <- sampdat$fam_id
 ss$study_id <- sampdat$study_id
 DispZyg2 <- vegan::betadisper(ord2, ss$fam_id)
  ppD <- as.data.frame(DispZyg2$distances)
  ppD$Distance.to.PC.median <- ppD$`DispZyg2$distances`
  cmon<-intersect(rownames(ppD),rownames(ss))
 ppD<-ppD[cmon,]
  ss <- ss[cmon,]
  boxplot(ppD$Distance.to.PC.median~ss$Zygosity)
  boxplot(ppD$Distance.to.PC.median~ss$Zygosity,outline=F,notch=T,boxwex=0.5,xlab="Twin type",ylab="Dispersal from pair")

###### within family segregation of urobiome
