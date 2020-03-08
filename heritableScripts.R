##########heritability scripts
####import ASV table
totuImp <-read.csv("totuImp.csv", row.names = 1, check.names = F)
totuImp <- t(totuImp)
dim(totuImp)
totuImp <- data.frame(totuImp, check.names = F)
####match  and sort microbiome data ID, sort with paired twins in the dataset (the twin of some individual were nt available)
cmon <- intersect(rownames(twins),rownames(totuImp))
twins<-twins[cmon,]
totuImp<-totuImp[cmon,]
identical(rownames(totuImp),rownames(totuImp))
 cmon2 <- intersect(colnames(totuImp), rownames(tax_nv))
 totuImp <- totuImp[,which(colnames(totuImp)%in%cmon2)]
 #####transform data
 ASV20clr = (totuImp + 1)
 library(compositions)
 ASV20clr = clr(ASV20clr)
 ASV20clr = as.data.frame.rmult(ASV20clr, check.names=F)
 detach("package:compositions", unload = TRUE)
 detach("package:tensorA", unload = TRUE)
 library(car): library(mets)
 clrb <- ASV20clr
 clrbresids <- ASV20clr
 #####loop over models and run analyses
 for(i in 1:ncol(clrb)){
     asv <- clrb[,i]
     #use clr only later u can use bcnpower with clr
     mod <- lm(asv~runSeq + reads_filt+extractionkit_lot + mastermix_lot + extraction_robot + processing_robot, data=twins)
     resids=summary(mod)$residuals
     #overwrite the column with residuals
     clrbresids[,i]=resids
 }

asv.table2 <- data.frame(clrbresids, check.names=F)
ace <- apply(asv.table2, 2, function(x) {mod1 <- twinlm(x ~ 1, data=twins, DZ="DZ", zyg="Zygosity", id="fam_id")
 coefs <- summary(mod1)$coef
 coefs2 <- summary(mod1)[[1]][14:16]
 coefs3 <- cbind(coefs, coefs2)
 return(coefs3)
 })

dft <- data.frame(ace, row.names = c("A_Estimate", "C_Estimate", "E_Estimate","A_2.5%", "C_2.5%", "E_2.5%","A_97.5%", "C_97.5%", "E_97.5%", "A_model.P.val","C_model.P.val","E_model.P.val"), check.names = F)
 dft <- t(dft)
dft <- data.frame(dft,check.names=F)
cmon<- intersect(rownames(dft), rownames(tax_5perc))
dft <- dft[cmon,]
tax_herit <- tax_5perc[cmon,]
tax_herit$A_Estimate <- dft$A_Estimate
tax_herit$A_low <- dft$A_2.5.
tax_herit$A_upp <- dft$A_97.5.
tax_herit$A_pval <- dft$A_model.P.val
tax_herit$C_Estimate <- dft$C_Estimate
tax_herit$E_Estimate <- dft$E_Estimate
 View(tax_herit)
write.csv(tax_herit, file="/Users/adewaleadebayo/Desktop/new_herit.csv")

###########repeat for clusters of core urobiome
residImp20 <-read.csv("/Users/adewaleadebayo/Desktop/exported_demux/Imp/imp20/Paird/core/Paird_residuals.csv",row.names=1, check.names=F)
residImp20 <- t(residImp20)
residImp20 <- data.frame(residImp20, check.names = F)
cmon2 <- intersect(rownames(twins), rownames(residImp20))
twins2 <- twins[cmon2,] 
twinwide <- fast.reshape(twins2, id="fam_id",varying="study_id")
twinwide <- twinwide[complete.cases(twinwide$study_id1, twinwide$study_id2),]
y <- twinwide$fam_id
twins2 <- twins2[which(twins2$fam_id%in%y),]
cmon2 <- intersect(rownames(twins2), rownames(residImp20))
residImp20 <- residImp20[cmon2,]
twins2 <- twins2[cmon2,]
identical(rownames(twins2), rownames(residImp20))
all.equal(rownames(twins2), rownames(residImp20))
head(twins2$fam_id, 10)
head(twins2$study_id, 10)
asv.table2 <- data.frame(residImp20, check.names=F)
ace <- apply(asv.table2, 2, function(x) {mod1 <- twinlm(x ~ 1, data=twins2, DZ="DZ", zyg="Zygosity", id="fam_id")
 coefs <- summary(mod1)$coef
 coefs2 <- summary(mod1)[[1]][14:16]
 coefs3 <- cbind(coefs, coefs2)
 return(coefs3)
 })
dft <- data.frame(ace, row.names = c("A_Estimate", "C_Estimate", "E_Estimate","A_2.5%", "C_2.5%", "E_2.5%","A_97.5%", "C_97.5%", "E_97.5%","A_model.P.val","C_model.P.val","E_model.P.val"), check.names = F)
 dft <- t(dft)
 dft <- data.frame(dft,check.names=F)
 View(dft)
