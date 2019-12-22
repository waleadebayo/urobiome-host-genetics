
library(dplyr)
library(tidyr)
z=sepMed
m<-c("retard","conti","maleate","fumarate","nanogram","hydrochloride",
     "solo","sodium","potassium","pill","implant")
strpAttach<- function(j,z,c){z %>% mutate(c = gsub("\\j\\d+$", "", c, ignore.case=TRUE))
  return(z)}
strpAttach<- function(j,z,c){z %>% mutate(c = gsub("j+$", "", c,ignore.case=TRUE))
  return(z)}
res=z
for (i in 1:ncol(z)){
  c=z[,i]
  for(j in m){
    j=m[j]
    res=strpAttach(j,z,c)
  }
}

strpAttach<- function(j,z,c){z %>% mutate(c = gsub("j+$", "", c,ignore.case=TRUE))
  return(z)}
##########new part
m <- as.data.frame(m)
fr=m
 for(i in 1:ncol(m)){
      r=m[,i]
      fr1=paste0(r,"+$")
      fr[,i]=fr1
    }
strpAttach<- function(j,z,c){z %>% mutate(c = gsub(j, "", c,ignore.case=TRUE))
  return(z)}
strpAttach<- function(j,z,c){z %>% mutate(c = gsub(j, "", c,ignore.case=TRUE))
  return(c)}

res=z
res1=res
for (i in 1:ncol(z)){
  c=z[,i]
  for(j in fr1){
    j=m[j]
    res=strpAttach(j,z,c)
    res1[,i]=res  }
}

for(j in fr1){
  res=strpAttach(j,z,c)
   }
for(j in fr1){
  res1=strpAttach2(j,res)
  y2[,i]=res1}

res=z
zz=z
strpAttach4<- function(fr1){for(j in fr1){res %>% mutate(c = gsub(j, "", c,ignore.case=TRUE))}}

for( j in fr1){
   res1=strpAttach3(j)
    }
res1=strpAttach3("conti+$")
res1=strpAttach3("ID10+$")
res1
#####################################################
###########fresh starts
##########new part
library(tibble)
library(dplyr)
library(fuzzyjoin)
library(tidytext)
library(tidyr)
library(stringdist)
sepMed0<-bcqMed[,3:12]
sepMed<-sepMed0
z=sepMed
m<-c("retard","conti","maleate","fumarate","nanogram","hydrochloride",
     "solo","sodium","potassium","pill","implant","fumarate","nocte",......)
mm <- as.data.frame(m)
fr=mm
for(i in 1:ncol(mm)){
  r=mm[,i]
  fr1=paste0(r,"+$")
  fr[,i]=fr1
}
#z
#dput(z)
#solution1
fr2 <- paste0(fr1,collapse = "|")
df4 <- z %>%
  mutate(c = stringr::str_replace_all(c,fr2,""))
df4
#solution 2
df6<-z
for(pat in fr1) df6$c <- sub(pat, "", df6$c)
#########
bnf2<-read.csv("/Users/macuser/Downloads/BNF_mod32.csv")
zif1<-data.frame(study_id=bcqMed$ParticipantID,value=zif$Y4,stringsAsFactors = F)
zif1[is.na(zif1)] <- ""
back2 <- as_tibble(zif1)

backdict2 <- data.frame(correct=bnf2$BNF.Product_Chemical,Class=bnf2$BNF.Subparagraph, stringsAsFactors = F)
backdict2 <- as_tibble(backdict2)
string_annotated <- back2 %>% 
  stringdist_left_join(y = backdict2, by = c("value"="correct"), max_dist = 1.5, ignore_case = TRUE,method="osa",weight=c(0.7,0.7,0.9,0.1)) %>% 
  mutate(match = !is.na(correct))

string_annotated_col <- string_annotated %>% 
  group_by(study_id) %>% 
  summarise(
    match = sum(match),
    keyword = toString(unique(na.omit(correct))),
Class = toString(unique(na.omit(Class))))
string_annotated_col3<-string_annotated_col
#Description = paste(word, collapse = " ")
write.csv(string_annotated_col,file="bcq_annotatedMatchesD1.5.csv")
newmed<-string_annotated_col
rownames(newmed)<-newmed$study_id
combo<-bcqMed
combo<-combo[which(!duplicated(combo$ParticipantID)),]
rownames(combo)<-combo$ParticipantID
cmon<-intersect(rownames(newmed),rownames(combo))
combo<-combo[cmon,]
newmed<-newmed[cmon,]
combo <- cbind(combo,newmed)
write.csv(combo,file="comboD1.5_ExpandedDictsa.csv")


##50zeros300threeormores
#stringdist_left_join(y = backdict2, by = c("value"="correct"), max_dist = 0.12,p=0.1,bt=0.1, ignore_case = TRUE,method="jw",weight=c(1,1,1)) %>% 
#..zeros...threeormores
# stringdist_left_join(y = backdict2, by = c("value"="correct"), max_dist = 0.15,p=0.1,bt=0.1, ignore_case = TRUE,method="jw",weight=c(0.7,1,0.8)) %>% 
#380zeros0threeormores
# stringdist_left_join(y = backdict2, by = c("value"="correct"), max_dist = 0.2, ignore_case = TRUE,method="osa",weight=c(1,1,0.5,0.5)) %>% 
#120zeros100threeormores
#stringdist_left_join(y = backdict2, by = c("value"="correct"), max_dist = 1.2, ignore_case = TRUE,method="osa",weight=c(1,1,0.5,0.5)) %>% 
#160zeros20threeormores
# stringdist_left_join(y = backdict2, by = c("value"="correct"), max_dist = 1.5, ignore_case = TRUE,method="osa",weight=c(1,1,1,0.1)) %>% 
##120zeros50threemores
#stringdist_left_join(y = backdict2, by = c("value"="correct"), max_dist = 1.7, ignore_case = TRUE,method="osa",weight=c(0.7,1,1,0.7)) %>% 
##130zeros50threemores
#max_dist = 1.7, ignore_case = TRUE,method="osa",weight=c(0.7,0.9,1,0.7))
#160zeros20threeormores
#max_dist = 1.5, ignore_case = TRUE,method="osa",weight=c(0.7,0.9,0.9,0.7)
##140zeros30threemores
#max_dist = 1.5, ignore_case = TRUE,method="osa",weight=c(0.7,0.7,0.9,0.1)
##108zeros100threeormores
#max_dist = 1.5, ignore_case = TRUE,method="osa",weight=c(0.7,0.5,0.9,0.1)


########################
##############now loop through all columns
#####
backdict2 <- data.frame(correct=bnf2$BNF.Product_Chemical,Class=bnf2$BNF.Subparagraph, stringsAsFactors = F)
backdict2 <- as_tibble(backdict2)
cll<-list()
for (j in 1:ncol(zif)){
  ll<-zif[,j] 
  zif2<-data.frame(study_id=bcqMed$ParticipantID,value=ll,stringsAsFactors = F)
  zif2[is.na(zif2)] <- ""
  back2 <- as_tibble(zif2)
  

  string_annotated <- back2 %>% 
    stringdist_left_join(y = backdict2, by = c("value"="correct"), max_dist = 1.5, ignore_case = TRUE,method="osa",weight=c(0.7,0.7,0.9,0.1)) %>% 
    mutate(match = !is.na(correct))
  
  string_annotated$i<-i
  cll[[i]]<-string_annotated
}  
  string_annotated_col <- string_annotated %>% 
    group_by(study_id) %>% 
    summarise(
      match = sum(match),
      keyword = toString(unique(na.omit(correct))),
      Class = toString(unique(na.omit(Class))))


string_annotated_col$i<-i
cll[[i]]<-string_annotated_col
