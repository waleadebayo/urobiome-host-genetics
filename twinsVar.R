
###
load("/Users/macuser/Documents/uromg/Variance/mgd2t_f2RR.RData")
load("/Users/macuser/Documents/uromg/Variance/hoci2k.RData")
rem <- read.csv(file="/Users/macuser/Documents/uromg/Variance/twins_all7factors_2k.csv",row.names = 1)
length(unique(rem$fam_id))
bdivUH<-prune_samples(rownames(rem),mgd2t_f2)
bdivUHd <-ordinate(bdivUH,method="PCoA", distance="bray")
bdivUHd<-as.data.frame(bdivUHd$vectors)
cmon<-intersect(rownames(bdivUHd), rownames(rem))
rem<-rem[cmon,]
bdivUHd<-bdivUHd[cmon,]
TwLME<-lme4::lmer(bdivUHd$Axis.1~scale(ageyear_urine)+scale(utihis)+scale(nBirths)+antib_use+scale(HEI)+scale(scqFI)+meno_stat+runSeq+scale(reads_filt)+extractionkit_lot+mastermix_lot+processing_robot+extraction_robot+(1|fam_id),REML=F,data=rem)
summary(TwLME)
TwLME2<-lme4::lmer(bdivUHd$Axis.1~scAge+scUTI+scale(nBirths)+antib_use+scale(HEI)+scqFI+meno_stat+runSeq+reads_filt+extractionkit_lot+mastermix_lot+processing_robot+extraction_robot,REML=F,data=rem)
######res
TwLME<-lme4::lmer(bdivUHd$Axis.2~scale(ageyear_urine)+scale(utihis)+scale(nBirths)+antib_use+scale(HEI)+scale(scqFI)+meno_stat+runSeq+scale(reads_filt)+extractionkit_lot+mastermix_lot+processing_robot+extraction_robot+(1|fam_id),REML=F,data=rem)
summary(TwLME)
TwLME2<-lm(bdivUHd$Axis.2~scale(ageyear_urine)+scale(utihis)+scale(nBirths)+antib_use+scale(HEI)+scale(scqFI)+meno_stat+runSeq+scale(reads_filt)+extractionkit_lot+mastermix_lot+processing_robot+extraction_robot,data=rem)
anova(TwLME,TwLME2)