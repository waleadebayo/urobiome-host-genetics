
load("/Volumes/LATEST/ethnic.RData")
 load("/Volumes/LATEST/Ancest2.RData")
load("/Volumes/LATEST/hoci2kAnc.RData")
 load("/Users/macuser/Desktop/divers/bdivBrd.RData")
pc <- cmdscale(bdivBrd)
pc <- as.data.frame(pc)
cmon<-intersect(rownames(pc), rownames(hoci2kAnc))
s1 <- hoci2kAnc[cmon,]
s2 <- pc[cmon,]
kruskal.test(s2$V1 ~ as.factor(s1$ancestry))
kruskal.test(s2$V2 ~ as.factor(s1$ancestry))
s1$PC1 <- s2$V1
 s1$PC2 <- s2$V2
 s11<-s1[which(!(s1$ancestry)=="British"),]
 s22<-s1[which((s1$ancestry)=="British"),]
set.seed(1659)
set1 <- caret::createFolds(rownames(s22),k=20)
sample11 <- s22[ set1$Fold11 , , drop=F]
sample12 <- s22[ set1$Fold12 , , drop=F]
sample13 <- s22[ set1$Fold13 , , drop=F]
sample14 <- s22[ set1$Fold14 , , drop=F]
sample15 <- s22[ set1$Fold15 , , drop=F]
sample16 <- s22[ set1$Fold16 , , drop=F]
sample17 <- s22[ set1$Fold17 , , drop=F]
sample18 <- s22[ set1$Fold18 , , drop=F]
sample19 <- s22[ set1$Fold19 , , drop=F]
sample20 <- s22[ set1$Fold20 , , drop=F]

dim(sample1)

sample11 <- rbind(sample11,s11)
sample12 <- rbind(sample12,s11)
sample13 <- rbind(sample13,s11)
sample14 <- rbind(sample14,s11)
sample15 <- rbind(sample15,s11)
sample16 <- rbind(sample16,s11)
sample17 <- rbind(sample17,s11)
sample18 <- rbind(sample18,s11)
sample19 <- rbind(sample19,s11)
sample20 <- rbind(sample20,s11)

indices=rownames(sample1)
foo <- function(data, indices){
 dt<-data[indices,]
 c(kruskal.test(dt[,"PC2"] ~ as.factor(dt[,"ancestry"]))$p.value
 )
 }
 library(boot)
 
set.seed(1659);
 myBootsrap <- boot(s1,foo,R=100)
myBootsrap$t
myBootsrap$t0
plot(myBootsrap)
boot.ci(myBootsrap,type="basic")
tableOfIndices<-boot.array(myBootstrap, indices=T)

set.seed(1659); myBootsrap <- boot(sample1,foo,R=10)




set.seed(1659);boot.ci(boot(sample11,foo,R=1000),type="perc") 
set.seed(1659); boot.ci(boot(sample12,foo,R=1000),type="perc")
set.seed(1659); boot.ci(boot(sample13,foo,R=1000),type="perc")
set.seed(1659); boot.ci(boot(sample14,foo,R=1000),type="perc")
set.seed(1659); boot.ci(boot(sample15,foo,R=1000),type="perc")
set.seed(1659); boot.ci(boot(sample16,foo,R=1000),type="perc")
set.seed(1659); boot.ci(boot(sample17,foo,R=1000),type="perc")
set.seed(1659); boot.ci(boot(sample18,foo,R=1000),type="perc")
set.seed(1659); boot.ci(boot(sample19,foo,R=1000),type="perc")
set.seed(1659); boot.ci(boot(sample20,foo,R=1000),type="perc")

set.seed(1659); mean(boot(sample11,foo,R=1000)$t)
set.seed(1659); mean(boot(sample12,foo,R=1000)$t)
set.seed(1659); mean(boot(sample13,foo,R=1000)$t)
set.seed(1659); mean(boot(sample14,foo,R=1000)$t)
set.seed(1659); mean(boot(sample15,foo,R=1000)$t)
set.seed(1659); mean(boot(sample16,foo,R=1000)$t)
set.seed(1659); mean(boot(sample17,foo,R=1000)$t)
set.seed(1659); mean(boot(sample18,foo,R=1000)$t)
set.seed(1659); mean(boot(sample19,foo,R=1000)$t)
set.seed(1659); mean(boot(sample20,foo,R=1000)$t)