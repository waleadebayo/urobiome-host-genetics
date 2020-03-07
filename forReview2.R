
library(mets)
twins <- sam
twinwide <- fast.reshape(twins, id="fam_id",varying= "study_id")

twinwide <- twinwide[complete.cases(twinwide$study_id1, twinwide$study_id2),]
x <- twinwide$fam_id
sampdat <- twins[which(twins$fam_id%in%x),]
table(duplicated(sampdat$Family_No))

load("/Users/adewaleadebayo/Desktop/postproc/re_metad/hoci_paired.RData")
load("~/Desktop/postproc/bdivBr.RData")
load("~/Desktop/postproc/bdivBrd.RData")
sam<- read.csv("/Users/adewaleadebayo/Desktop/postproc/RthGd/hociHEI.csv",row.names = 1)
 bb<-rownames(as.matrix(bdivBrd))
bb <- as.data.frame(bb)
rownames(bb)<-bb$bb
bb$names<-bb$bb
twins <- sam
rownames(twins) <- twins$ucsd_Name
cmon<-intersect(rownames(twins), rownames(bb))
twins <- twins[cmon,]
bb<-bb[cmon,]

 twinwide <- fast.reshape(twins, id="fam_id",varying= "study_id")
 twinwide <- twinwide[complete.cases(twinwide$study_id1, twinwide$study_id2),]
 x <- twinwide$fam_id
 sampdat <- twins[which(twins$fam_id%in%x),]
sampdat<-sampdat[order(sampdat$fam_id),]

bb<- ape::pcoa(bdivBrd,correction = "lingoes")
bb<-as.data.frame(bb$vectors)
cmon<-intersect(rownames(sampdat), rownames(bb))
sampdat <- sampdat[cmon,]
bb<-bb[cmon,]
asv<-bb$Axis.2
mod <- lm(asv~runSeq + reads_filt+extractionkit_lot + 
            mastermix_lot + extraction_robot + 
            processing_robot, data=sampdat)
resids=summary(mod)$residuals
mod1 <- twinlm(resids ~ 1, data=sampdat, DZ="DZ", zyg="Zygosity", id="fam_id")
summary(mod1)
#load("/Users/adewaleadebayo/Desktop/postproc/compareUrSt/mgd2t_f2NV.RData")
#mgd2t_f2 <- phyloseq::prune_samples(rownames(sampdat),mgd2t_f2)
#ord<-phyloseq::ordinate(mgd2t_f2,method = "PCoA",distance = "bray")
#phyloseq::plot_ordination(mgd2t_f2,ordination = ord)
ord <- phyloseq::distance(mgd2t_f2,method = "bray")
total = phyloseq::sample_sums(mgd2t_f2)
total = median(total)
standf = function(x, t=total) round(t * (x / sum(x)))
gps = phyloseq::transform_sample_counts(mgd2t_f2, standf)
ord2 <- phyloseq::distance(gps,method = "bray")
ss<-rownames(as.matrix(ord2))
ss <- as.data.frame(ss)
rownames(ss)<-ss$ss
ss$names<-ss$ss
cmon <-intersect(rownames(ss),rownames(sampdat))
ss$Zygosity <- sampdat$Zygosity
ss$fam_id <- sampdat$fam_id
ss$study_id <- sampdat$study_id
DispZyg <- vegan::betadisper(ord2, ss$fam_id)
permutest(DispZyg,permutations = 1000)
plot(TukeyHSD(DispZyg))
#DispZyg2 <- vegan::betadisper(ord, ss$,bias.adjust = T)
pp <- as.data.frame(DispZyg$centroids)
tt<- fast.reshape(ss, id="fam_id",varying= "study_id")
rownames(tt)<-tt$fam_id
cmon<-intersect(rownames(pp),rownames(tt))
pp<-pp[cmon,]
tt <- tt[cmon,]
boxplot(pp$PCoA2~tt$Zygosity)

ppD <- as.data.frame(DispZyg2$distances)
ppD$Distance.to.PC.median <-ppD$`DispZyg2$distances`
cmon<-intersect(rownames(ppD),rownames(ss))
ppD<-ppD[cmon,]
ss <- ss[cmon,]
boxplot(ppD$Distance.to.PC.median~ss$Zygosity)
kruskal.test(ppD$Distance.to.PC.median~ss$Zygosity)
#ord<- phyloseq::
##
ss<-rownames(as.matrix(ord2))
 ss <- as.data.frame(ss)
 rownames(ss)<-ss$ss
 ss$names<-ss$ss
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

######plotting from desktop at work

 load("/Users/macuser/Documents/uromg/Variance/hoci2k.RData")
 dispZyg<-read.csv("/Users/macuser/Downloads/reheritscri-2/DisppairedTwinZyg.csv",row.names=1)
library(mets)
 twinwide <- fast.reshape(hoci2k,varying="study_id",id="fam_id")
 twinwide<-twinwide[complete.cases(twinwide$study_id1,twinwide$study_id2),]
pair1<-hoci2k[which(hoci2k$study_id%in%twinwide$study_id1),]
 pair2<-hoci2k[which(hoci2k$study_id%in%twinwide$study_id2),]
 pair1<- pair1[order(pair1$study_id),]
pair2<- pair2[order(pair2$study_id),]
cmon<- intersect(rownames(pair1),rownames(dispZyg))
pair1<-pair1[cmon,]
dispZyg1<-dispZyg[cmon,]
cmon<- intersect(rownames(pair2),rownames(dispZyg))
pair2<-pair2[cmon,]
dispZyg2<-dispZyg[cmon,]
 pair1$Dispersal <- dispZyg1
 pair2$Dispersal <- dispZyg2
 DispersalDiscord <- abs(pair1$Dispersal - pair2$Dispersal)
 pair1$DispersalDiscord <- DispersalDiscord
 pair2$DispersalDiscord <- DispersalDiscord
pairJoin <- rbind(pair1,pair2)
 boxplot(pairJoin$DispersalDiscord ~ pairJoin$Zygosity)
 boxplot(pairJoin$DispersalDiscord ~ pairJoin$Zygosity,xlab="Twin Type",ylab="Dispersal of microbiome within pair",boxwex=0.6)
 boxplot(pairJoin$DispersalDiscord ~ pairJoin$Zygosity,xlab="Twin Type",ylab="Dispersal of microbiome within pair",boxwex=0.5,notch=T)
 save.image("/Users/macuser/Desktop/forReview/DispersalPair.RData")