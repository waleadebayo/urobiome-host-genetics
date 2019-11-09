comboData<-read.csv("/Users/macuser/Downloads/comboData4.csv", row.names=1)
forwgs<-read.csv("/Users/macuser/Downloads/for_wgs_2.csv", row.names=1)
load("/Users/macuser/Documents/uromg/Variance/hoci2k.RData")
head(rownames(forwgs),3)
head(rownames(comboData),3)
hoci2k$rname<-rownames(hoci2k)
rownames(hoci2k)<-hoci2k$study_id
forwgs$rname<-rownames(forwgs)
rownames(forwgs)<-forwgs$study_id
cmon<-intersect(rownames(forwgs),rownames(hoci2k))
forwgs<-forwgs[cmon,]
hociF<-hoci2k[cmon,]
forwgs$rrname<-hociF$rname
rownames(forwgs)<-forwgs$rrname
cmon<-intersect(rownames(forwgs), rownames(comboData))
head(forwgs$fam_id,3)
forwgs<-forwgs[cmon,]
comboData<-comboData[cmon,]
comboData$Div_Discordance <- forwgs$Discordance
comboData$Discordance <- forwgs$category
comboData$Zygosity <- forwgs$Zygosity
forwgs$genID<-paste("S",seq(1:nrow(forwgs)),sep="")
##heatmap(as.matrix(comboData[,-c(39:41)]), scale = "row",labRow = forwgs$genID,col=cm.colors(256),Colv=NA,Rowv=NA)
##reorderfun = function(d,w)reorder(Rowv,forwgs$Zygosity)
##heatmap(as.matrix(comboData[,-c(39:41)]), scale = "none",)
mat1<-comboData[,-c(19:21)]
rownames(mat1)<-forwgs$genID


mat2<-as.data.frame(mat2)
mat3<-mat2
for (i in 1:ncol(mat3)){
  asv<-mat3[,i]
  mod1<-scale(asv)
  mat2[,i]=mod1}


mat2<-t(mat1)
mat2<-as.data.frame(mat2)
mat3<-mat2
for (i in 1:ncol(mat3)){
  asv<-mat3[,i]
  mod1<-scale(asv)
  mat2[,i]=mod1}
mat2<-as.matrix(mat2)
mat2<-t(mat2)

Heatmap(mat2,name="proportion", row_order = rownames(mat2),column_order = colnames(mat2),row_names_gp = gpar(fontsize = 7))
Heatmap(mat2,split = forwgs$Zygosity,row_order = rownames(mat2),column_order = colnames(mat2),row_names_gp = gpar(fontsize = 7))
Heatmap(mat2,row_order = rownames(mat2),column_order = colnames(mat2),row_names_gp = gpar(fontsize = 7))

Heatmap(mat2,split = data.frame(type=forwgs$Zygosity,fam=forwgs$fam_id),row_order = rownames(mat2),column_order = colnames(mat2),row_names_gp = gpar(fontsize = 7))

row_names_gp = gpar(fontsize = 5),row_dend_gp = gpar(fontsize = 5),column_names_gp = gpar(fontsize = 5)