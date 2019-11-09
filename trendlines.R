 
odkegg<-read.csv("/Users/macuser/Downloads/peptide_metabolism/AgeR2.order2.csv",row.names = 1)
mat1<-t(odkegg)
 mat2<-(mat1[2:19, ] + 1 )
 mat2<-apply(mat2,2,function(x) x/sum(x))

 mat3<-rbind(mat2,mat1[1,])
 mat3<-t(mat3)
 mat3<-as.data.frame(mat3)

 library(ggplot2)
ggplot(data=mat3, aes(x=V19, y=X2.3.1..Transferring.groups.other.than.aminoacyl.groups)) +
  geom_smooth(method="lm",se=T)+theme(axis.text.x = element_text(size=10))+xlab("Age")

ggplot(data=mat3, aes(x=V19, y=X2.3.2..Aminoacyltransferases)) +
  geom_smooth(method="lm",se=T)+theme(axis.text.x = element_text(size=10))+xlab("Age")


ggplot(data=mat3, aes(x=V19, y=X3.4.11..Aminopeptidases)) +
  geom_smooth(method="lm",se=T)+theme(axis.text.x = element_text(size=10))+xlab("Age")


ggplot(data=mat3, aes(x=V19, y=X3.4.13..Dipeptidases)) +
  geom_smooth(method="lm",se=T)+theme(axis.text.x = element_text(size=10))+xlab("Age")



ggplot(data=mat3, aes(x=V19, y=X3.4.14..Dipeptidyl.peptidases.and.tripeptidyl.peptidases)) +
  geom_smooth(method="lm",se=T)+theme(axis.text.x = element_text(size=10))+xlab("Age")



ggplot(data=mat3, aes(x=V19, y=X3.4.15..Peptidyl.dipeptidases)) +
  geom_smooth(method="lm",se=T)+theme(axis.text.x = element_text(size=10))+xlab("Age")


ggplot(data=mat3, aes(x=V19, y=X3.4.16..Serine.type.carboxypeptidases)) +
  geom_smooth(method="lm",se=T)+theme(axis.text.x = element_text(size=10))+xlab("Age")


ggplot(data=mat3, aes(x=V19, y=X3.4.17..Metallocarboxypeptidases)) +
  geom_smooth(method="lm",se=T)+theme(axis.text.x = element_text(size=10))+xlab("Age")


ggplot(data=mat3, aes(x=V19, y=X3.4.19..Omega.peptidases)) +
  geom_smooth(method="lm",se=T)+theme(axis.text.x = element_text(size=10))+xlab("Age")



ggplot(data=mat3, aes(x=V19, y=X3.4.21..Serine.endopeptidases)) +
  geom_smooth(method="lm",se=T)+theme(axis.text.x = element_text(size=10))+xlab("Age")



ggplot(data=mat3, aes(x=V19, y=X3.4.23..Aspartic.endopeptidases)) +
  geom_smooth(method="lm",se=T)+theme(axis.text.x = element_text(size=10))+xlab("Age")




ggplot(data=mat3, aes(x=V19, y=X3.4.22..Cysteine.endopeptidases)) +
  geom_smooth(method="lm",se=T)+theme(axis.text.x = element_text(size=10))+xlab("Age")



ggplot(data=mat3, aes(x=V19, y=X3.4.22..Cysteine.endopeptidases)) +
  geom_smooth(method="lm",se=T)+theme(axis.text.x = element_text(size=10))+xlab("Age")



ggplot(data=mat3, aes(x=V19, y=X3.4.24..Metalloendopeptidases)) +
  geom_smooth(method="lm",se=T)+theme(axis.text.x = element_text(size=10))+xlab("Age")


ggplot(data=mat3, aes(x=V19, y=X3.4.25..Threonine.endopeptidases)) +
  geom_smooth(method="lm",se=T)+theme(axis.text.x = element_text(size=10))+xlab("Age")



ggplot(data=mat3, aes(x=V19, y=X3.4....Acting.on.peptide.bonds..peptidases.)) +
  geom_smooth(method="lm",se=T)+theme(axis.text.x = element_text(size=10))+xlab("Age")


ggplot(data=mat3, aes(x=V19, y=X3.5....Acting.on.carbon.nitrogen.bonds..other.than.peptide.bonds)) +
  geom_smooth(method="lm",se=T)+theme(axis.text.x = element_text(size=10))+xlab("Age")


ggplot(data=mat3, aes(x=V19, y=X6.3.2..Acid.D.amino.acid.ligases..peptide.synthases.)) +
  geom_smooth(method="lm",se=T)+theme(axis.text.x = element_text(size=10))+xlab("Age")


ggplot(data=mat3, aes(x=V19, y=X5.1.1..Acting.on.amino.acids.and.derivatives)) +
  geom_smooth(method="lm",se=T)+theme(axis.text.x = element_text(size=10))+xlab("Age")

###
plots.dir.path <- list.files(tempdir(), pattern="rs-graphics", full.names = TRUE); 
plots.png.paths <- list.files(plots.dir.path, pattern=".png", full.names = TRUE)
file.copy(from=plots.png.paths, to="peptide2")

plots.png.detials <- file.info(plots.png.paths)
plots.png.detials <- plots.png.detials[order(plots.png.detials$mtime),]
sorted.png.names <- gsub(plots.dir.path, "peptide2", row.names(plots.png.detials), fixed=TRUE)
numbered.png.names <- paste0("peptide2/", 1:length(sorted.png.names), ".png")

# Rename all the .png files as: 1.png, 2.png, 3.png, and so on.
file.rename(from=sorted.png.names, to=numbered.png.names)

###tryheatmap # heatmap(as.matrix(mat3[,1:18]), scale = "column",labRow = mat3[,19],Colv=NA,Rowv=NA)
################supplementary
summary(lm(data=mat3, V19~X3.4.13..Dipeptidases)) 
#for R1 significant p<0.1
#Coefficients:
#  Estimate
#(Intercept)                                                        67.0246
#X3.5....Acting.on.carbon.nitrogen.bonds..other.than.peptide.bonds  53.8549
#Std. Error
#(Intercept)                                                           0.6563
#X3.5....Acting.on.carbon.nitrogen.bonds..other.than.peptide.bonds    39.1012
#t value
#(Intercept)                                                       102.126
#X3.5....Acting.on.carbon.nitrogen.bonds..other.than.peptide.bonds   1.377
#Pr(>|t|)
#(Intercept)                                                         <2e-16
#X3.5....Acting.on.carbon.nitrogen.bonds..other.than.peptide.bonds     0.17

#(Intercept)                                                       ***
#  X3.5....Acting.on.carbon.nitrogen.bonds..other.than.peptide.bonds    
#---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 7.893 on 175 degrees of freedom
#Multiple R-squared:  0.01072,	Adjusted R-squared:  0.005071 
#F-statistic: 1.897 on 1 and 175 DF,  p-value: 0.1702


#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)                     66.5125     0.8622  77.146   <2e-16 ***
# X3.4.15..Peptidyl.dipeptidases  61.7542    43.0123   1.436    0.153    
#---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 7.889 on 175 degrees of freedom
#Multiple R-squared:  0.01164,	Adjusted R-squared:  0.005994 
#F-statistic: 2.061 on 1 and 175 DF,  p-value: 0.1529

#Coefficients:
#  Estimate
#(Intercept)                                                66.6617
#X3.4.14..Dipeptidyl.peptidases.and.tripeptidyl.peptidases  49.0238
#Std. Error
#(Intercept)                                                   0.8275
#X3.4.14..Dipeptidyl.peptidases.and.tripeptidyl.peptidases    37.7106
#t value
#(Intercept)                                                 80.56
#X3.4.14..Dipeptidyl.peptidases.and.tripeptidyl.peptidases    1.30
#Pr(>|t|)    
#(Intercept)                                                 <2e-16 ***
#  X3.4.14..Dipeptidyl.peptidases.and.tripeptidyl.peptidases    0.195    
#---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 7.897 on 175 degrees of freedom
#Multiple R-squared:  0.009565,	Adjusted R-squared:  0.003905 

#Coefficients:
# Estimate Std. Error t value Pr(>|t|)    
#(Intercept)             65.880      1.067  61.760   <2e-16 ***
# X3.4.13..Dipeptidases   62.552     36.262   1.725   0.0863 .  
#---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 7.869 on 175 degrees of freedom
#Multiple R-squared:  0.01672,	Adjusted R-squared:  0.0111 
#F-statistic: 2.976 on 1 and 175 DF,  p-value: 0.08629

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)                71.652      1.549  46.269  < 2e-16 ***
# X3.4.11..Aminopeptidases  -45.274     15.321  -2.955  0.00356 ** 
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 7.745 on 175 degrees of freedom
#Multiple R-squared:  0.04753,	Adjusted R-squared:  0.04208 
#F-statistic: 8.732 on 1 and 175 DF,  p-value: 0.003557
#######for R2 significant p=<0.1
#mat3$X3.4.11..Aminopeptidases  -46.928      0.00307
#mat3$X3.4.13..Dipeptidases    65.34       0.0668 .
#mat3$X3.4.14..Dipeptidyl.peptidases.and.tripeptidyl.peptidases  61.1623
#mat3$X3.4.14..Dipeptidyl.peptidases.and.tripeptidyl.peptidases   0.0988

#mat3$X3.4.15..Peptidyl.dipeptidases  84.7261
#mat3$X3.4.15..Peptidyl.dipeptidases   0.0455 * 
#mat3$X3.4.17..Metallocarboxypeptidases  57.4795    37.7816   1.521     0.13 
#mat3$X3.4.21..Serine.endopeptidases  -51.165
#mat3$X3.4.21..Serine.endopeptidases   0.0269 * 
#mat3$X3.4.25..Threonine.endopeptidases  53.6078
#mat3$X3.4.25..Threonine.endopeptidases   0.0964 .  
#mat3$X3.5....Acting.on.carbon.nitrogen.bonds..other.than.peptide.bonds  74.9713 beta
#mat3$X3.5....Acting.on.carbon.nitrogen.bonds..other.than.peptide.bonds   0.0491 (pval)



#######
colnames(mat3)
[1] "X2.3.1..Transferring.groups.other.than.aminoacyl.groups"          
[2] "X2.3.2..Aminoacyltransferases"                                    
[3] "X3.4.11..Aminopeptidases"                                         
[4] "X3.4.13..Dipeptidases"                                            
[5] "X3.4.14..Dipeptidyl.peptidases.and.tripeptidyl.peptidases"        
[6] "X3.4.15..Peptidyl.dipeptidases"                                   
[7] "X3.4.16..Serine.type.carboxypeptidases"                           
[8] "X3.4.17..Metallocarboxypeptidases"                                
[9] "X3.4.19..Omega.peptidases"                                        
[10] "X3.4.21..Serine.endopeptidases"                                   
[11] "X3.4.23..Aspartic.endopeptidases"                                 
[12] "X3.4.22..Cysteine.endopeptidases"                                 
[13] "X3.4.24..Metalloendopeptidases"                                   
[14] "X3.4.25..Threonine.endopeptidases"                                
[15] "X3.4....Acting.on.peptide.bonds..peptidases."                     
[16] "X3.5....Acting.on.carbon.nitrogen.bonds..other.than.peptide.bonds"
[17] "X6.3.2..Acid.D.amino.acid.ligases..peptide.synthases."            
[18] "X5.1.1..Acting.on.amino.acids.and.derivatives"                    
[19] "V19"                                                              

#####################combining plots
draw.data <- function(xy, add = FALSE) {
  x.lab <- "concentration [M]"
  y.lab <- "normalised luminescence [%]"
  my.data <- data.frame(xy)
  my_labels <- parse(text = paste("1E", seq(-10, -4, 1), sep = ""))
  y.max <- 1
  y.min <- 0
  diff <- y.max - y.min
  my.data$y <- apply(my.data, 1, function(z) ((z["y"] - y.min)/diff)*100)
  if (!add) {
    ggplot(my.data, aes(x,y)) +
      geom_line() +
      geom_hline(yintercept = c(50, 90), linetype = "dotted") +
      scale_x_continuous(x.lab, breaks = seq(-10, -4, 1), 
                         labels = my_labels) + 
      labs(title = "Graph", x = x.lab, y = y.lab)
  } else {
    geom_line(aes(x, y), data = my.data )
  }
}

###################
###############another way of combining

plot_aligned_series.R https://gist.github.com/tomhopper/faa24797bb44addeba79
#' When plotting multiple data series that share a common x axis but different y axes,
#' we can just plot each graph separately. This suffers from the drawback that the shared axis will typically
#' not align across graphs due to different plot margins.
#' One easy solution is to reshape2::melt() the data and use ggplot2's facet_grid() mapping. However, there is
#' no way to label individual y axes.
#' facet_grid() and facet_wrap() were designed to plot small multiples, where both x- and y-axis ranges are
#' shared acros all plots in the facetting. While the facet_ calls allow us to use different scales with
#' the \code{scales = "free"} argument, they should not be used this way.
#' A more robust approach is to the grid package grid.draw(), rbind() and ggplotGrob() to create a grid of 
#' individual plots where the plot axes are properly aligned within the grid.
#' Thanks to https://rpubs.com/MarkusLoew/13295 for the grid.arrange() idea.

library(ggplot2)
library(grid)
library(dplyr)

#' Create some data to play with. Two time series with the same timestamp.
df <- data.frame(DateTime = ymd("2010-07-01") + c(0:8760) * hours(2), series1 = rnorm(8761), series2 = rnorm(8761, 100))

#' Create the two plots.
plot1 <- df %>%
  select(DateTime, series1) %>%
  na.omit() %>%
  ggplot() +
  geom_point(aes(x = DateTime, y = series1), size = 0.5, alpha = 0.75) +
  ylab("Red dots / m") +
  theme_minimal() +
  theme(axis.title.x = element_blank())

plot2 <- df %>%
  select(DateTime, series2) %>%
  na.omit() %>%
  ggplot() +
  geom_point(aes(x = DateTime, y = series2), size = 0.5, alpha = 0.75) +
  ylab("Blue drops / L") +
  theme_minimal() +
  theme(axis.title.x = element_blank())

grid.newpage()
grid.draw(rbind(ggplotGrob(plot1), ggplotGrob(plot2), size = "last"))

###############
