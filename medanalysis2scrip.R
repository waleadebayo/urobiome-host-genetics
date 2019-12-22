

library(tibble)
library(dplyr)
library(fuzzyjoin)
library(tidytext)
library(tidyr)
z=sepMed
m<-c("retard","conti","maleate","fumarate","nanogram","hydrochloride",
     "solo","sodium","potassium","pill","implant")
#strpAttach<- function(j,z,c){z %>% mutate(c = gsub("\\j\\d+$", "", c))
 # return(z)}
#res=z
#for (i in 1:ncol(z)){
#  c=z[,i]
#  for(j in m){
#    j=m[j]
#    res=strpAttach(j,z,c)
#  }
#}

y2<-z
#fret<-sapply(m,function(e) {y2 %>% mutate(c = gsub("\\ed+$", "", c))})
#strpAttach<- function(j,y,c){y %>% mutate(c = gsub("j+$", "", c))
#  return(y)}
#y=z
#res2=z
#res=z
#for (i in 1:ncol(y)){
 #   for(j in m){
  #  c=y[,"c"]
   # res=strpAttach(j,y,c)
    #res2[,i]=res
#  }
#}
###################go on
cat(m,sep="|")
y<-z
#y2<-y %>% mutate(c = gsub("(retard|conti|maleate|fumarate|nanogram|hydrochloride|solo|sodium|potassium|pill|implant)+$", "", c))
sepMed0<-bcqMed[,3:12]
sepMed<-sepMed0
#sepMed2<-sepMed0
#for (i in 1:ncol(sepMed)) {
#  c<-sepMed[,i]
#  zi<- sepMed %>% mutate(c = gsub("(retard|conti|maleate|fumarate|nanogram|hydrochloride|solo|sodium|potassium|pill|implant)+$", "", c,ignore.case = T))
#}
zif<-sepMed
for (i in 1:ncol(sepMed)) {
zi<-as.data.frame(sepMed[,i])
 colnames(zi)<-"c"
zi<- zi %>% mutate(c= gsub("(retard|conti|maleate|nocte|fumarate|nanogram|hydrochloride|solo|sodium|potassium|pill|implant)+$", "",c,ignore.case = TRUE))
zif[,i]=zi$c
}
zif1<-data.frame(study_id=bcqMed$ParticipantID,value=zif$Y1,stringsAsFactors = F)
zif1[is.na(zif1)] <- ""
back2 <- as_tibble(zif1)

backdict2 <- data.frame(correct=bnf$BNF.Product,stringsAsFactors = F)
backdict2 <- as_tibble(backdict2)
string_annotated <- back2 %>% 
  stringdist_left_join(y = backdict2, by = c("value"="correct"), max_dist = 2, ignore_case = TRUE,method="lv") %>% 
  mutate(match = !is.na(correct))

string_annotated_col <- string_annotated %>% 
  group_by(study_id) %>% 
  summarise(
            match = sum(match),
            keyword = toString(unique(na.omit(correct))))
#Event_Name = toString(unique(na.omit(Event_Name)))
#Description = paste(word, collapse = " ")