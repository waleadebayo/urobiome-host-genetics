

urob <- read.delim("/Users/macuser/Desktop/GWAS/urobHLI2.txt", header=F)
colnames(urob) <- c("IID","FID","COV1","COV2","COV3")
rawBal <- read.csv("/Volumes/LATEST/gwas2/imp20residuals.csv", row.names = 1,check.names = F)
rawBal<-t(rawBal)
rawBal<-as.data.frame(rawBal)
samBal <- read.csv( "/Users/macuser/Desktop/package/earlier/hociHEI.csv", row.names = 1)
rownames(samBal)<-samBal$rname
cmon<-intersect(rownames(samBal), rownames(rawBal))
samBal<-samBal[cmon,]
rawBal<-rawBal[cmon,]
rownames(samBal)<-samBal$study_id
rownames(rawBal)<-rownames(samBal)
x<-urob$IID
samBal<-samBal[which(rownames(samBal)%in%x),]
rawBal<-rawBal[which(rownames(rawBal)%in%x),]
cmon <- intersect(rownames(rawBal), rownames(samBal))
rawBal <- rawBal[cmon,]
samBal <- samBal[cmon,]
########if it is not residuals
resBal <- rawBal
for (i in 1:ncol(rawBal)){
  asv <- rawBal[,i]
  mod <- lm(asv ~ runSeq+extractionkit_lot+mastermix_lot+processing_robot+extraction_robot+primer_plate, data=samBal)
  resids=summary(mod)$residuals
  resBal[,i]=resids
}
######################
resBal<-rawBal
cmon <- intersect(rownames(resBal), rownames(samBal))
resBal<-resBal[cmon,]
samBal<-samBal[cmon,]
resBal$scAge <- scale(samBal$ageyear_urine)
resBal$scVFL <- scale(samBal$VFL_score)
resBal$scMM <- scale(samBal$total_mmse)
resBal$scFI <- scale(samBal$newBCQFI)
resBal$scUTI <- scale(samBal$utihis)
#if needed rownames(resBal)<-samBal$study_id
UROB<-urob
rownames(UROB) <- UROB$IID
cmon <- intersect(rownames(UROB), rownames(resBal))
UROB <- UROB[cmon,]
resBal<- resBal[cmon,]
UROB <- cbind(UROB, resBal)
####
rownames(UROB)<-UROB$FID
UROB$SID<-UROB$IID
UROB$IID<-NULL
###this works
UROB2<-UROB
UROB2<-UROB2[,c(2,1,3:ncol(UROB2))]
write.table(UROB2,file="/Users/macuser/Desktop/UROB7", sep= " " ,row.names = F)

###########terminal2Rosalind
#scp  /Users/macuser/Desktop/UROB6 k1787722@login.rosalind.kcl.ac.uk:brc_scratch/HLI/urob.unique2
##remove quotes and replace y with Cb
#######inrosalind
#cd brc_scratch
#mkdir HLI
#
#plink --file ~/ehlai/HLI/new/newww  --keep urob.unique2 --make-bed  --out urobhli2 --memory 56000
#plink --bfile urobhli2 --pca --out PCAurobhli --memory 56000
#plink --bfile urobhli2 --pheno urob.unique2   --all-pheno  --linear --covar PCAurobhli.eigenvec  --pfilter 1e-5  --out Cball --memory 56000
#put in sh file
#scp  /Users/macuser/Desktop/hello2.sh k1787722@login.rosalind.kcl.ac.uk:brc_scratch/HLI/hello2.sh
#scp  /Users/macuser/Desktop/hello2_y33.sh k1787722@login.rosalind.kcl.ac.uk:brc_scratch/HLI/hello2_y33.sh

#qsub -V -N urogw2  -o logurob -cwd -pe smp 6 -l h_vmem=32G hello2.sh
#qsub -V -N urogw  -o logurob -cwd -l h_vmem=32G hello2_y33.sh
#plink --bfile urobhli2 --pheno urob.unique2   --pheno-name y33  --linear --covar PCAurobhli.eigenvec  --pfilter 1e-2  --out Cball_y33 --memory 56000
#use 302 set because 407 has samples with zero microbiome
#do the pca, send to shell and queue it on hpc,
#shell taking time do the guys u were interested in,replace output names with names now being used in new tree
#plink --bfile urobhli2 --pheno ~/ehlai/HLI/new/subset2k/pheno   --pheno-name Cb52,Cb7,Cb10  --linear --adjust --covar PCAurobhli.eigenvec  --pfilter 1e-4  --out cb33,cb9,cb13_302 --memory 56000 

#in rosalindlinuxR
#module add general/R/3.5.2
#cd ~/ehlai/HLI/new/subset2k
#R
#module add general/R/3.5.0
#results_as<-read.table("~/brc_scratch/output/cb13_302.assoc.linear.adjusted" ,head=T)
# dim(results_as)
# results_as2 <- results_as[order(results_as$FDR_BH, decreasing=F),]
# results_as2[1:10,]
#plot qq and mh
#results was filtered at p=1e-4 so plots will not have full xsome cannot
#results_as<-read.table("~/brc_scratch/output/cb13_302.assoc.linear" ,head=T)
#qq(results_as$P)
#manhattan(results_as,chr=“CHR”,bp=“BP”,p=“P”,snp=“SNP”)
 



################clinical measure
pheno <- read.delim("/Users/macuser/Desktop/GWAS/nGWASwtCovar/pheno",sep=" ", header=T)
rownames(pheno)<-pheno$IID
cmon <- intersect(rownames(pheno), rownames(samBal))
pheno<-pheno[cmon,]
samBal<-samBal[cmon,]
pheno$scAge <- scale(samBal$ageyear_urine)
pheno$scVFL <- scale(samBal$VFL_score)
pheno$scMM <- scale(samBal$total_mmse)
pheno$scFI <- scale(sqrt(samBal$newBCQFI))
pheno$scUTI <- scale(samBal$utihis)
pheno$scqVFL <- scale(samBal$scqVFL)
pheno$scqMM <- scale(samBal$scqMMSE)
pheno$scqFI <- scale((samBal$scqFI))
UROB2<-pheno
UROB2<-UROB2[,c(2,1,3:ncol(UROB2))]
write.table(UROB2,file="/Users/macuser/Desktop/UROB9", sep= " " ,row.names = F)
###add
#add --prune to remove missing phenotypes#


