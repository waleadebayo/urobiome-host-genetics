
##import newick format of core urobiome associations with ageing
install.packages("Phinch")
load("/Volumes/LATEST/tre13.RData")
?tidytree::as.treedata
?tidytree::data_frame
treD<-tidytree::as_tibble(tre13)
View(treD)
write.csv(treD,file="treD.csv")

############
library(devtools)
install_github("Displayr/flipPlots")
Fig1sk<-read.csv("/Volumes/LATEST/Fig1sankey.csv",row.names = 1)
library(flipPlots)
Fig2sk<-read.csv("/Volumes/LATEST/Fig1sankey2.csv",row.names = 1)

SankeyDiagram(Fig1sk[, -c(7,8)],link.color = "Source",weights = Fig1sk$logInvFDR)
snk_node<-SankeyDiagram(Fig1sk[, -c(7,8)],link.color = "Source",weights = Fig1sk$logInvFDR,label.show.varname = F)
SankeyDiagram(Fig1sk[, -c(7,8)],link.color = "Source",weights = Fig1sk$logInvFDR,label.show.varname = F)

#snk2<-SankeyDiagram(Fig2sk[, -c(7,8)],link.color = "Source",weights = Fig2sk$logInvFDR,label.show.varname = F,font.size=10,max.categories = 50)

library(htmlwidgets)
install.packages("webshot")
library(webshot)
webshot::install_phantomjs()

saveWidget(snk_node,'temp.html')
saveWidget(snk2,'temp2.html')
webshot('temp.html',file='trySankey_colNode.pdf',vwidth = 1210,vheight = 920)
webshot('temp2.html',file='trySankey_colNode2.pdf',vwidth = 1000,vheight = 1250)

##################final

library(devtools)
install_github("Displayr/flipPlots")
library(flipPlots)
library(htmlwidgets)
install.packages("webshot")
library(webshot)
webshot::install_phantomjs()

Fig3sk<-read.csv("/Volumes/LATEST/Fig1sankey.csv",row.names = 1)
snk3<-SankeyDiagram(Fig3sk[, -c(7,8)],link.color = "Source",weights = Fig3sk$logInvFDR,label.show.varname = F,font.size=15,max.categories = 50)
saveWidget(snk3,'temp3.html')
webshot('temp3.html',file='Fig1Sankey_uro.pdf',vwidth = 1210,vheight = 920)

Fig4sk<-read.csv("/Volumes/LATEST/Fig1sankey4.csv",row.names = 1); snk4<-SankeyDiagram(Fig4sk[, -c(6,7)],link.color = "Source",weights = Fig4sk$logInvFDR,label.show.varname = F,font.size=15,max.categories = 50)
saveWidget(snk4,'temp4.html')
webshot('temp4.html',file='trySankey5.pdf',vwidth = 750,vheight = 500)
