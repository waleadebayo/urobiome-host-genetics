#!/bin/bash

conda install -c biobakery graphlan

###prepare annotation files to ASVs at different levels of annotations
annot1,annot2, annot3
#####to ensure proper encoding re-write annotations as tab-delimited from R
#tt<-read.csv("/Users/adewaleadebayo/Desktop/GraphlAN/final/annot3.csv" ,head=F)
#write.table(tt,file="/Users/adewaleadebayo/Desktop/GraphlAN/final/annot3.txt", quote=F, row.names = F,col.names =F,sep="\t")
 #####
graphlan_annotate.py --annot annot_0.txt origM.txt tre12x.xml
graphlan.py tre12x.xml test2.pdf --dpi 300 --size 8
graphlan_annotate.py --annot annot1.txt tre12x.xml tre12x1.xml
graphlan_annotate.py --annot annot2.txt tre12x1.xml tre12x2.xml
graphlan_annotate.py --annot annot3.txt tre12x2.xml tre12x2.xml

graphlan.py tre12x2.xml test31.pdf --dpi 300 --size 5
graphlan_annotate.py --annot annot1.txt tre12x3.xml tre12x3.xml
graphlan.py tre12x3.xml test33.pdf --dpi 300 --size 3.5