load("/Users/macuser/Downloads/comparedata.RData")
library(phyloseq)
Family.sum = tapply(taxa_sums(pearceStudy), tax_table(pearceStudy)[, "Genus"], sum, na.rm=TRUE)
 top15phyla = names(sort(Family.sum, TRUE))[1:15] ; top15phyla
 Family.sum = tapply(taxa_sums(whiteStudy), tax_table(whiteStudy)[, "Genus"], sum, na.rm=TRUE)
 top15phyla = append(top15phyla ,names(sort(Family.sum, TRUE))[1:15]) ; top15phyla
 Family.sum = tapply(taxa_sums(TwU_Study), tax_table(TwU_Study)[, "Genus"], sum, na.rm=TRUE)
 top15phyla = append(top15phyla ,names(sort(Family.sum, TRUE))[1:15]) ; top15phyla
 topphyla<-unique(top15phyla)
 topphyla <- topphyla[nchar(topphyla)>1]
 taxcombo2
 ordFam15 = subset_taxa(taxcombo2, Genus==topphyla)
 ordFam15 <- transform_sample_counts(ordFam15, function(x) x/sum(x))
 ##ssT10 <- subset_samples(taxcombo2, Study%in%c("Pearce_Ur","T_White_Ur", "Tw_Uro100"))
 ##ordFam15 <- transform_sample_counts(ssT10, function(x) x/sum(x))
 ##ordFam15 = subset_taxa(ssT10, Genus==topphyla)
 colrCount <- length(unique(get_taxa_unique(ordFam15, "Genus")))
 data10 <- psmelt(ordFam15)
 #data101 <- data10
 data101<-data10[which((data10$Study)%in%c("Pearce_Ur","T_White_Ur", "Tw_Uro100")),]
 data101$Study<-gsub("Pearce_Ur", "Urine1",data101$Study)
 data101$Study<-gsub("T_White_Ur", "Urine2",data101$Study)
 data101$Study<-gsub("Tw_Uro100", "Urine3",data101$Study)
 tiff("Taxcombo_rawRawTop10each.tiff", width=4.6, height=5, res=300, units="in") 
 ggplot(data101, aes(x=Study, y=Abundance, fill=Genus)) 
 +geom_bar(aes(), stat="identity") + scale_fill_manual(values=getPalette(colrCount)) 
 + theme(panel.background = element_rect(fill="white"), axis.text.y = element_blank(),
         axis.ticks.y.left = element_blank(), legend.position="right") + guides(fill=guide_legend(ncol=1)) +ylab("Abundance") ; dev.off()
 tiff("Taxcombo_relRawTop10each.tiff", width=4.6, height=5, res=300, units="in") 
 
 tiff("Taxcombo_rawFillTop10each.tiff", width=4.6, height=5, res=300, units="in") 
 ggplot(data101, aes(x=Study, y=Abundance, fill=Genus)) +geom_bar(aes(),position="fill", stat="identity") + scale_fill_manual(values=getPalette(colrCount)) + theme(panel.background = element_rect(fill="white"), axis.text.y = element_blank(),axis.ticks.y.left = element_blank(), legend.position="right") + guides(fill=guide_legend(ncol=1)) +ylab("Relative abundance") ;dev.off()
 tiff("Taxcombo_relFillTop10each.tiff", width=4.6, height=5, res=300, units="in") 
 
 
 
 
 
 ###jumboTax <- rbind(TwU_biom1501,white_tax,pearce_tax,prico_tax,TwSt_tax1)
 ##comboSet2wtTaxa <- merge_phyloseq(comboSet2, jumboTax)
 
 ###ade4 plotting
 unif1kRR <- distance(comboSet2wtTaxa, method = "unifrac")
 ssT<-sample_data(comboSet2wtTaxa)[,"Study"]
 ssT<-data.frame(ssT)
 unif<-cmdscale(unif1kRR)
 pdf("uni_study2Set2.pdf", height=5.4, width=7)
 ade4::s.class(unif,as.factor(ssT$Study),col = c("darkgoldenrod1","dodgerblue2","darkorange3","forestgreen","firebrick3"),label=c("Urine1","Vagina","Urine2","Gut","Urine3"),cpoint=0.4);dev.off()
 pdf("uni_study2.pdf", height=5.4, width=7)
 ade4::s.class(unif,ssT$Study,col = c("darkgoldenrod1","dodgerblue2","darkorange3","forestgreen","firebrick3"),
               label=c("Urine1","Vagina","Urine2","Gut","Urine3"),cpoint=0.4) ; dev.off()
 
 ###
 ssT10 <- subset_samples(comboSet2wtTaxa, Study%in%c("Pearce_Ur","T_White2_Ur", "Tw2_Uro100"))
 ssT10<- transform_sample_counts(ssT10, function(x) x/sum(x))
 ssT10 = subset_taxa(ssT10, Genus%in%topphyla)
 data10 <- psmelt(ssT10)
  data101 <- data10
  data101$Study<-gsub("Pearce_Ur", "Urine1",data101$Study)
  data101$Study<-gsub("T_White2_Ur", "Urine2",data101$Study)
  data101$Study<-gsub("Tw2_Uro100", "Urine3",data101$Study)
 tiff("comboSet2.tiff", width=4.6, height=5.2, res=300, units="in")
 ggplot(data101, aes(x=Study, y=Abundance, fill=Genus)) +geom_bar(aes(),position="fill", stat="identity") + scale_fill_manual(values=getPalette(colrCount)) + theme(panel.background = element_rect(fill="white"), axis.text.y = element_blank(),axis.ticks.y.left = element_blank(), legend.position="right") + guides(fill=guide_legend(ncol=1)) +ylab("Relative abundance") ; dev.off()
 
 ###take unrarefied raw data from the three urine
 ordFam15 <- transform_sample_counts(whiteStudy, function(x) x/sum(x))
 ordFam15 = subset_taxa(ordFam15, Genus%in%topphyla)
 ordFam152 <- transform_sample_counts(pearceStudy, function(x) x/sum(x))
 ordFam152 = subset_taxa(ordFam152, Genus%in%topphyla)
 ordFam153 <- transform_sample_counts(mgd2t_f2, function(x) x/sum(x))
 ordFam153 = subset_taxa(ordFam153, Genus%in%topphyla)
 colnames(tax_table(ordFam153)) <- c("Domain" , "Phylum",  "Class" ,  "Order" ,  "Family" , "Genus",   "Species")
 urines_unrarefO<-merge_phyloseq(otu_table(ordFam15),otu_table(ordFam152),otu_table(ordFam153))
 dim(urines_unrarefO)
 urines_unrarefT<-merge_phyloseq(tax_table(ordFam15),tax_table(ordFam152),tax_table(ordFam153))
 dim(urines_unrarefT)
 oo1<-sample_data(ordFam15)
 oo1$StudyType<-"Urine2"
 oo2<-sample_data(ordFam152)
 oo2$StudyType<-"Urine1"
 oo3<-sample_data(ordFam153)
 oo3$StudyType<-"Urine3"
 oo1<-oo1[,c("Study","StudyType")] 
 oo2<-oo2[,c("Study","StudyType")] 
 oo3<-oo3[,c("Study","StudyType")] 
 oo1<-sample_data(oo1)
 oo2<-sample_data(oo2)
 oo3<-sample_data(oo3)
 urines_unrarefS<-merge_phyloseq(oo1,oo2,oo3)
 dim(urines_unrarefS)
urines_unraref <- phyloseq(urines_unrarefO,urines_unrarefT,urines_unrarefS) 
urines_unraref
data10 <- psmelt(urines_unraref)
tiff("comboRAW_Urines.tiff", width=4.6, height=5.2, res=300, units="in")
ggplot(data10, aes(x=Study, y=Abundance, fill=Genus)) +geom_bar(aes(),position="fill", stat="identity") + scale_fill_manual(values=getPalette(colrCount)) + theme(panel.background = element_rect(fill="white"), axis.text.y = element_blank(),axis.ticks.y.left = element_blank(), legend.position="right") + guides(fill=guide_legend(ncol=1)) +ylab("Relative abundance") 
; dev.off()

########people who have never had UTI
ss<-read.csv("/Users/macuser/Desktop/package/earlier/hociHEI.csv",row.names = 1)
rownames(ss)<-ss$ucsd_Name
cmon<-intersect(rownames(oo3),rownames(ss))
oo3<-oo3[cmon,]
ss <-ss[cmon,]
oo3$Prior_UTI <- ss$utihis
oo3<-sample_data(oo3)
ss <- merge_phyloseq(mgd2t_f2, oo3)
ss <- subset_samples(ss, Prior_UTI=="0" )
ordFam154 <- transform_sample_counts(ss, function(x) x/sum(x))
ordFam154 = subset_taxa(ordFam154, Genus%in%topphyla)
tt<-otu_table(ordFam154)
oo4<-sample_data(ordFam154)
colnames(tt)<-paste("NH",colnames(tt), sep="_")
rownames(oo4)<-paste("NH",rownames(oo4), sep="_")
oo4$StudyType<-"Urine3_NoHistUTI"
oo4<-oo4[,c("Study","StudyType")]
urines_unrarefS<-merge_phyloseq(oo1,oo2,oo3,oo4)
urines_unrarefO<-merge_phyloseq(otu_table(ordFam15),otu_table(ordFam152),otu_table(ordFam153),tt)
urines_unraref <- phyloseq(urines_unrarefO,urines_unrarefT,urines_unrarefS) 
data101 <- psmelt(urines_unraref)
tiff("comboNoRaref_Urines_4type.tiff", width=4.6, height=5.2, res=300, units="in")
ggplot(data101, aes(x=StudyType, y=Abundance, fill=Genus)) +geom_bar(aes(),position="fill", stat="identity") + scale_fill_manual(values=getPalette(colrCount)) + theme(panel.background = element_rect(fill="white"), axis.text.y = element_blank(),axis.ticks.y.left = element_blank(), legend.position="right") + guides(fill=guide_legend(ncol=1)) +ylab("Relative abundance") 
; dev.off()

###############
#####Equalize number of samples with rarefied subsamples data with the noutihis group added#
set.seed(165387) ;  rr <- sample(rownames(oo4), size=100, replace = F)
dim(tt)
tt <-prune_samples(rr, tt)
dim(tt)
oo4 <-prune_samples(rr, oo4)
ssT10 <- subset_samples(comboSet2wtTaxa, Study%in%c("Pearce_Ur","T_White2_Ur", "Tw2_Uro100"))
ssT10<- transform_sample_counts(ssT10, function(x) x/sum(x))
ssT10 = subset_taxa(ssT10, Genus%in%topphyla)
oo4$Study<-oo4$StudyType
rr <- merge_phyloseq(otu_table(ssT10),tt)
axx<-merge_phyloseq(tax_table(ssT10),urines_unrarefT)
ass<-merge_phyloseq(sample_data(ssT10), oo4)
urines_raref<-phyloseq(rr,axx,ass)
data10 <- psmelt(urines_raref)
range(sample_sums(urines_raref))
data101 <- data10
unique(sample_data(urines_raref)[,"Study"])
data101$Study<-gsub("Pearce_Ur", "Urine1",data101$Study)
data101$Study<-gsub("T_White2_Ur", "Urine2",data101$Study)
data101$Study<-gsub("Tw2_Uro100", "Urine3",data101$Study)

tiff("comboSet2_Urines_4type.tiff", width=4.6, height=5.5, res=300, units="in")
ggplot(data101, aes(x=Study, y=Abundance, fill=Genus)) +geom_bar(aes(),position="fill", stat="identity") + 
  scale_fill_manual(values=getPalette(colrCount)) + theme(panel.background = element_rect(fill="white"),
                                                          axis.text.x = element_text(angle = 20), axis.text.y = element_blank(),axis.ticks.y.left = element_blank(), legend.position="right") + 
  guides(fill=guide_legend(ncol=1)) +ylab("Relative abundance");dev.off()