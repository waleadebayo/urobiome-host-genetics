library(broom);
library(ggplot2)
; library(compositions)
cmon2<-intersect(rownames(putr1clust),rownames(sam))
sigClustMg<-putr1clust[cmon2,]
sam<-sam[cmon2,]
sigClustMg$fam_id<-sam$fam_id
sigClustMg$study_id<-sam$study_id
sigClustMg$age<-sam$ageyear_urine
sigClustMg$utihis<-sam$utihis
sigClustMg$meno_stat<-sam$meno_stat
sigClustMg$Zygosity<-sam$Zygosity
sigClustMg<-cbind.data.frame(sigClustMg,sam[,90:161])
# sumSam<-rowSums(putspp)
# sumSam<-as.data.frame(sumSam)
# cmon3<-intersect(rownames(sigClustMg), rownames(sumSam))
# sumSam<-sumSam[cmon3,]
# sigClustMg$sumSamMG20perc<-sumSam
#cmonn<-intersect(rownames(forwgs2),rownames(sigClustMg))
#forwgs2<-forwgs2[cmonn,]
#sigClustMg<-sigClustMg[cmonn,]
#sigClustMg$discord<-forwgs2$Discordance
sigClustMg$scAge<- scale(sigClustMg$age)
sigClustMg$scUTI<- scale(sigClustMg$utihis)

###
totuImp <- sigClustMg[,1:90]
clrInd <- totuImp
clrInd1 <- clrInd
clrInd2 <- clrInd1
clrInd3 <- clrInd1
###
for(i in 1:ncol(clrInd1)) {
       sp <- clrInd1[,i]
       modInd <- lm(sp ~ scAge+discord+scUTI+meno_stat+antib_use2+sumSamMG20perc, data=sigClustMg, na.action = na.omit)
       summ <- tidy(modInd)
       summ <- summ[2,5]
       clrInd1[,i]= summ}
for(i in 1:ncol(clrInd2)) {
       sp <- clrInd2[,i]
       modInd <- lm(sp ~ scAge+discord+scUTI+meno_stat+antib_use2+sumSamMG20perc, data=sigClustMg, na.action = na.omit)
       summ <- tidy(modInd)
       summ <- summ[2,2]
      clrInd2[,i]= summ}
for(i in 1:ncol(clrInd3)) {
       sp <- clrInd3[,i]
       modInd <- lm(sp ~ scAge+discord+scUTI+meno_stat+antib_use2+sumSamMG20perc, data=sigClustMg, na.action = na.omit)
       summ <- tidy(modInd)
       summ <- summ[2,3]
       clrInd3[,i]= summ}
clrInd1 <- clrInd1[1,]
clrInd2 <- clrInd2[1,]
clrInd3 <- clrInd3[1,]
###
clrInd1 <- t(clrInd1)
clrInd2 <- t(clrInd2)
clrInd3 <- t(clrInd3)
####
colnames(clrInd1) <- "age_pval"
colnames(clrInd2) <- "age_estimate"
colnames(clrInd3) <- "age_st.error"
clrInd1 <- data.frame(clrInd1)
clrInd2 <- data.frame(clrInd2)
clrInd3 <- data.frame(clrInd3)
####
clrInd1$estimate <- clrInd2$age_estimate
clrInd1$st.error <- clrInd3$age_st.error
clrInd1$fdr<-p.adjust(clrInd1$age_pval,method = c("fdr"))
quantile(clrInd1$fdr,probs=c(0,0.25,0.5,1))


ageASV3 <- clrInd1[clrInd1$age_pval<0.05,]
ageASV3$min <- (ageASV3$estimate - ageASV3$st.error)
ageASV3$max <- (ageASV3$estimate + ageASV3$st.error)

library(ggplot2)
c<- ggplot(ageASV3) + geom_bar(aes(x=rownames(ageASV3), y=estimate, 
                                    fill=1/fdr), stat="identity") + 
  xlab("Species")+ylab("Effect on average proportions")+ 
  scale_fill_gradient(aes(), low="aquamarine", high = "aquamarine4", 
                      guide = FALSE) + 
  geom_errorbar(aes(x=rownames(ageASV3),ymin=min, ymax=max), width=0.1, size=0.3, alpha=0.9) + 
  theme_minimal() + coord_flip()
 
 c+ annotate("segment",size=0.6,alpha=0.6, x=14.05,y=0.35,yend=0.58,xend=14.05,arrow=arrow())+
   annotate("text",x=14.40,y=0.40,label="Age",size=3.5)
 
 #cmon2<-intersect(rownames(putspp2),rownames(sigClustMg))
 # putspp2<-putspp2[cmon2,]
 # sigClustMg<-sigClustMg[cmon2,]