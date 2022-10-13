#Set the library
library(readxl)
library(ade4)
library(factoextra)
library(vegan)
library(corrplot)
library(gridExtra)
library(dplyr)
library(magrittr)
library(psych)
library(tidyverse)
library(ggpubr)
library(Hmisc)
library(viridis)
library(ggcorrplot)
library(cowplot)
library(lemon)
library(colorspace)
library(egg)


## Set the working directory and read in the data
dir <- "C:/Users/nicol/Documents/experimental_analysis"
setwd(dir)
df<-read.table("GlobalAnalysis.csv", header = TRUE, sep = ",")

# Convert data to data frame
df<-data.frame(df)


# subset the data by depth
Global_0.5<- subset(df, df$Depth == 0.5)
Global_1 <- subset(df, df$Depth == 1)
Global_1.5 <- subset(df, df$Depth == 1.5)
Global_2 <- subset(df, df$Depth == 2)

#Histograms of each environmental variable
hist(Global_0.5$CC, xlab = "Cloud cover (*100%)", ylim=range(0,120))
hist(Global_0.5$WH,  xlab = "Wave height (cm)", ylim=range(0,120))
hist(Global_0.5$WS, xlab = "Wind speed (km/hr)", ylim=range(0,120))
hist(Global_0.5$AVG_SA, xlab = "Sun alitude (degrees)", ylim=range(0,120), xlim=range(0,60))
hist(Global_0.5$Secchi_Adj,xlab = "Secchi depth (m)", ylim=range(0,120))
hist(Global_0.5$T_H, xlab = "Tidal height (m)", ylim=range(0,120))

# Set row names
rownames(Global_0.5)<-Global_0.5$ID
rownames(Global_1)<-Global_1$ID
rownames(Global_1.5)<-Global_1.5$ID
rownames(Global_2)<-Global_2$ID

#Summary statistics
mean(Global_0.5$SPA_P)
sd(Global_0.5$SPA_P)
mean(Global_0.5$SPB_P)
sd(Global_0.5$SPB_P)
mean(Global_0.5$SPC_P)
sd(Global_0.5$SPC_P)
mean(Global_0.5$SPD_P)
sd(Global_0.5$SPD_P)

mean(Global_1$SPA_P)
sd(Global_1$SPA_P)
mean(Global_1$SPB_P)
sd(Global_1$SPB_P)
mean(Global_1$SPC_P)
sd(Global_1$SPC_P)
mean(Global_1$SPD_P)
sd(Global_1$SPD_P)

mean(Global_1.5$SPA_P)
sd(Global_1.5$SPA_P)
mean(Global_1.5$SPB_P)
sd(Global_1.5$SPB_P)
mean(Global_1.5$SPC_P)
sd(Global_1.5$SPC_P)
mean(Global_1.5$SPD_P)
sd(Global_1.5$SPD_P)

mean(Global_2$SPA_P)
sd(Global_2$SPA_P)
mean(Global_2$SPB_P)
sd(Global_2$SPB_P)
mean(Global_2$SPC_P)
sd(Global_2$SPC_P)
mean(Global_2$SPD_P)
sd(Global_2$SPD_P)

mean(Global_0.5$Global_P)
sd(Global_0.5$Global_P)
mean(Global_1$Global_P)
sd(Global_1$Global_P)
mean(Global_1.5$Global_P)
sd(Global_1.5$Global_P)
mean(Global_2$Global_P)
sd(Global_2$Global_P)

##SUMMARY STATS FOR ENVIRONMENTAL VARIABLES
median(Global_0.5$AVG_SA)
sd(Global_0.5$AVG_SA)
range(Global_0.5$AVG_SA)
mean(Global_0.5$AVG_SA)

median(Global_0.5$CC)
sd(Global_0.5$CC)
range(Global_0.5$CC)
mean(Global_0.5$CC)

median(Global_0.5$WS)
sd(Global_0.5$WS)
range(Global_0.5$WS)
mean(Global_0.5$WS)

median(Global_0.5$WH)
sd(Global_0.5$WH)
range(Global_0.5$WH)
mean(Global_0.5$WH)

median(Global_0.5$Secchi_Adj)
sd(Global_0.5$Secchi_Adj)
range(Global_0.5$Secchi_Adj)
mean(Global_0.5$Secchi_Adj)

View(Global_0.5)
# Subset each data set so it only has specific columns (Chosen based on iterations of PCA)
Global_0.5<-select(Global_0.5, c( Secchi_Adj, CC, WH,WS, AVG_SA, LureA))
Global_1<-select(Global_1, c(Secchi_Adj, CC, WH, WS, AVG_SA, LureA))
Global_1.5<-select(Global_1.5, c(Secchi_Adj, CC, WH, WS, AVG_SA, LureA, LureB, LureC, LureD))
Global_2<-select(Global_2, c(Secchi_Adj,CC, WH,WS, AVG_SA, LureB, LureC, LureD))


#Global PCA for each depth 
# columns 1-5 include the environmental variables and column 6 has the global % 
PCA_0.5<-dudi.pca(Global_0.5 [,],scale = TRUE,center = TRUE,scannf = FALSE,nf=5)
PCA_1<-dudi.pca(Global_1 [,],scale = TRUE,center = TRUE,scannf = FALSE,nf=5)
PCA_1.5<-dudi.pca(Global_1.5 [,],scale = TRUE,center = TRUE,scannf = FALSE,nf=5)
PCA_2<-dudi.pca(Global_2 [,],scale = TRUE,center = TRUE,scannf = FALSE,nf=5)

#This checks the % of variance explained by each PCs
# If the label shows an eigen % of >10%, The PC accounts for more variance than one of the original (scaled) variables
# PC's with >10% should be kept 
fviz_eig(PCA_0.5, addlabels = TRUE, ylim=c(0,60), main= "PCA 0.5m")
fviz_eig(PCA_1, addlabels = TRUE, ylim=c(0,60),main= "PCA 1m")
fviz_eig(PCA_1.5, addlabels = TRUE, ylim=c(0,60),main= "PCA 1.5m" )
fviz_eig(PCA_2, addlabels = TRUE, ylim=c(0,60),main= "PCA 2m" )

## Extract results
var_0.5<- get_pca_var(PCA_0.5)
var_1<- get_pca_var(PCA_1)
var_1.5<- get_pca_var(PCA_1.5)
var_2<- get_pca_var(PCA_2)

#correlation circle
fviz_pca_var(PCA_0.5, col.var="cos2", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")) 
fviz_pca_var(PCA_1, col.var="cos2", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07") )
fviz_pca_var(PCA_1.5, col.var="cos2", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"))
fviz_pca_var(PCA_2,  col.var="cos2", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"))

# Quality of representation
var_0.5$cos2<-data.matrix(var_0.5$cos2)
corrplot(var_0.5$cos2, is.corr=FALSE)

var_1$cos2<-data.matrix(var_1$cos2)
corrplot(var_1$cos2, is.corr=FALSE)

var_1.5$cos2<-data.matrix(var_1.5$cos2)
corrplot(var_1.5$cos2, is.corr=FALSE)

var_2$cos2<-data.matrix(var_2$cos2)
corrplot(var_2$cos2, is.corr=FALSE)

# COntributions of vairables to each PC
# Contributions of variables in accounting for the variablilty in a given PC (%)
#Corrplot to visualize the most contributing variables to each dimension/PC.
#The larger the value, the more a variable contributes to the PC
corrplot(var_0.5$contrib, is.corr = FALSE)
corrplot(var_1$contrib, is.corr = FALSE)
corrplot(var_1.5$contrib, is.corr = FALSE)
corrplot(var_2$contrib, is.corr = FALSE)

# A way to view the most important contributing variable to each PC
# add "axes=2" to show the most important for dimension/PC 2
# the redline shows variables that contribute more than the expected average (1/6=17%) 
# this can be considered a cut off point for analysis - if contribute more, than it is important contributing factorto the PC
fviz_contrib(PCA_0.5, choice = "var", axes=1:3)
fviz_contrib(PCA_1, choice = "var", axes = 1:3)
fviz_contrib(PCA_1.5, choice = "var", axes = 1:3)
fviz_contrib(PCA_2,choice = "var", axes = 1:2)

# Another way to visualize contribution to PC 1 and 2 
fviz_pca_var(PCA_0.5, title= "PCA 0.5m", col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"))
fviz_pca_var(PCA_1, title="PCA 1m", col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"))
fviz_pca_var(PCA_1.5, title="PCA 1.5m",col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"))
fviz_pca_var(PCA_2, title="PCA 2m", col.var = "contrib",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"))


#Plot Diatom
fviz_pca_biplot(PCA_0.5,axes=c(1,2), alpha.var ="contrib",col.var = "contrib",
                gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                pointsize = 4,label="var", labelsize = 6,
                addEllipses = FALSE,ellipse.level=0.6,col="black",
                palette = c("red","navy"))+theme_gray()+xlim(-6,6)+ylim(-5,5)+
  theme(axis.text.x = element_text(angle = 90, size = 28, colour = "black", vjust = 0.5, hjust = 1, face= "bold"),
        axis.title.x = element_text(size = 28, face = "bold"),
        axis.text.y = element_text(colour = "black", size = 28, face = "bold"),
        legend.text = element_text(size = 28, face = "bold", colour = "black"),
        axis.title.y = element_text(size = 28, face = "bold"), legend.title = element_text(size = 28, face = "bold"))+
  #GGtitle
  ggtitle("PCA 0.5m") +
  theme(plot.title = element_text(size=40,face="bold",hjust = 0.5))


dev.off()

fviz_pca_biplot(PCA_1,alpha.var ="contrib",col.var = "contrib",
                gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                pointsize = 4,label="var", labelsize = 6,
                addEllipses = FALSE,ellipse.level=0.6,col="black",
                palette = c("red","navy"))+theme_gray()+xlim(-6,6)+ylim(-5,5)+
  theme(axis.text.x = element_text(angle = 90, size = 28, colour = "black", vjust = 0.5, hjust = 1, face= "bold"),
        axis.title.x = element_text(size = 28, face = "bold"),
        axis.text.y = element_text(colour = "black", size = 28, face = "bold"),
        legend.text = element_text(size = 28, face = "bold", colour = "black"),
        axis.title.y = element_text(size = 28, face = "bold"), legend.title = element_text(size = 28, face = "bold"))+
  #GGtitle
  ggtitle("PCA 1m") +
  theme(plot.title = element_text(size=40,face="bold",hjust = 0.5))


dev.off()

fviz_pca_biplot(PCA_1.5,alpha.var ="contrib",col.var = "contrib",
                gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                pointsize = 4,label="var", labelsize = 6,
                addEllipses = FALSE,ellipse.level=0.6,col="black",
                palette = c("red","navy"))+theme_gray()+xlim(-6,6)+ylim(-5,5)+
  theme(axis.text.x = element_text(angle = 90, size = 28, colour = "black", vjust = 0.5, hjust = 1, face= "bold"),
        axis.title.x = element_text(size = 28, face = "bold"),
        axis.text.y = element_text(colour = "black", size = 28, face = "bold"),
        legend.text = element_text(size = 28, face = "bold", colour = "black"),
        axis.title.y = element_text(size = 28, face = "bold"), legend.title = element_text(size = 28, face = "bold"))+
  #GGtitle
  ggtitle("PCA 1.5m") +
  theme(plot.title = element_text(size=40,face="bold",hjust = 0.5))


dev.off()

fviz_pca_biplot(PCA_2,alpha.var ="contrib",col.var = "contrib",
                gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                pointsize = 4,label="var", labelsize = 6,
                addEllipses = FALSE,ellipse.level=0.6,col="black",
                palette = c("red","navy"))+theme_gray()+xlim(-6,6)+ylim(-5,5)+
  theme(axis.text.x = element_text(angle = 90, size = 28, colour = "black", vjust = 0.5, hjust = 1, face= "bold"),
        
        axis.title.x = element_text(size = 28, face = "bold"),
        axis.text.y = element_text(colour = "black", size = 28, face = "bold"),
        legend.text = element_text(size = 28, face = "bold", colour = "black"),
        axis.title.y = element_text(size = 28, face = "bold"), legend.title = element_text(size = 28, face = "bold"))+
  #GGtitle
  ggtitle("PCA 2m") +
  theme(plot.title = element_text(size=40,face="bold",hjust = 0.5))


dev.off()

##### CORRELATION #####
## Plotting fish visibility vs each environmental variable to determine if the relationships are linear and have good spread of data
# EXPLORING RELATIONSHIPS THAT WERE SIGNIFICANT IN PEARSONS REGRESSION
ggplot(Global_0.5, aes(WH, WS)) +geom_point(position=position_jitter(h=0.1, w=0.1),shape = 21, alpha = 0.5, size = 3, color="#01587A", fill="#077DAA")+geom_smooth(method = lm)
ggplot(Global_0.5, aes(LureA, CC)) +geom_point(position=position_jitter(h=0.1, w=0.1),shape = 21, alpha = 0.5, size = 3, color="#01587A", fill="#077DAA")+geom_smooth(method = lm)
ggplot(Global_0.5, aes(LureA, Secchi_Adj)) +geom_point(position=position_jitter(h=0.1, w=0.1), shape = 21, alpha = 0.5, size = 3, color="#01587A", fill="#077DAA")+geom_smooth(method = lm)
ggplot(Global_1, aes(LureA, CC)) +geom_point(position=position_jitter(h=0.1, w=0.1),shape = 21, alpha = 0.5, size = 3, color="#01587A", fill="#077DAA")+geom_smooth(method = lm)
ggplot(Global_1.5, aes(LureA, CC)) +geom_point(position=position_jitter(h=0.1, w=0.1), shape = 21, alpha = 0.5, size = 3, color="#01587A", fill="#077DAA")+geom_smooth(method = lm)


ggplot(Global_1.5, aes(LureB, WH)) +geom_point(position=position_jitter(h=0.1, w=0.1), shape = 21, alpha = 0.5, size = 3, color="#01587A", fill="#077DAA")+geom_smooth(method = lm)
ggplot(Global_1.5, aes(LureC, WH)) +geom_point(position=position_jitter(h=0.1, w=0.1), shape = 21, alpha = 0.5, size = 3, color="#01587A", fill="#077DAA")+geom_smooth(method = lm)
ggplot(Global_1.5, aes(LureD, WH)) +geom_point(position=position_jitter(h=0.1, w=0.1), shape = 21, alpha = 0.5, size = 3, color="#01587A", fill="#077DAA")+geom_smooth(method = lm)

ggplot(Global_1.5, aes(LureB, WS)) +geom_point(position=position_jitter(h=0.1, w=0.1), shape = 21, alpha = 0.5, size = 3, color="#01587A", fill="#077DAA")+geom_smooth(method = lm)
ggplot(Global_1.5, aes(LureC, WS)) +geom_point(position=position_jitter(h=0.1, w=0.1), shape = 21, alpha = 0.5, size = 3, color="#01587A", fill="#077DAA")+geom_smooth(method = lm)
ggplot(Global_1.5, aes(LureD, WS)) +geom_point(position=position_jitter(h=0.1, w=0.1), shape = 21, alpha = 0.5, size = 3, color="#01587A", fill="#077DAA")+geom_smooth(method = lm)

ggplot(Global_1.5, aes(LureB, AVG_SA)) +geom_point(position=position_jitter(h=0.1, w=0.1), shape = 21, alpha = 0.5, size = 3, color="#01587A", fill="#077DAA")+geom_smooth(method = lm)
ggplot(Global_1.5, aes(LureC, AVG_SA)) +geom_point(position=position_jitter(h=0.1, w=0.1), shape = 21, alpha = 0.5, size = 3, color="#01587A", fill="#077DAA")+geom_smooth(method = lm)
ggplot(Global_1.5, aes(LureD, AVG_SA)) +geom_point(position=position_jitter(h=0.1, w=0.1), shape = 21, alpha = 0.5, size = 3, color="#01587A", fill="#077DAA")+geom_smooth(method = lm)

ggplot(Global_2, aes(LureB, WH)) +geom_point(position=position_jitter(h=0.1, w=0.1), shape = 21, alpha = 0.5, size = 3, color="#01587A", fill="#077DAA")+geom_smooth(method = lm)
ggplot(Global_2, aes(LureC, WH)) +geom_point(position=position_jitter(h=0.1, w=0.1), shape = 21, alpha = 0.5, size = 3, color="#01587A", fill="#077DAA")+geom_smooth(method = lm)
ggplot(Global_2, aes(LureD, WH)) +geom_point(position=position_jitter(h=0.1, w=0.1), shape = 21, alpha = 0.5, size = 3, color="#01587A", fill="#077DAA")+geom_smooth(method = lm)

ggplot(Global_2, aes(LureB, WS)) +geom_point(position=position_jitter(h=0.1, w=0.1), shape = 21, alpha = 0.5, size = 3, color="#01587A", fill="#077DAA")+geom_smooth(method = lm)
ggplot(Global_2, aes(LureC, WS)) +geom_point(position=position_jitter(h=0.1, w=0.1), shape = 21, alpha = 0.5, size = 3, color="#01587A", fill="#077DAA")+geom_smooth(method = lm)
ggplot(Global_2, aes(LureD, WS)) +geom_point(position=position_jitter(h=0.1, w=0.1), shape = 21, alpha = 0.5, size = 3, color="#01587A", fill="#077DAA")+geom_smooth(method = lm)

ggplot(Global_2, aes(LureB, AVG_SA)) +geom_point(position=position_jitter(h=0.1, w=0.1), shape = 21, alpha = 0.5, size = 3, color="#01587A", fill="#077DAA")+geom_smooth(method = lm)
ggplot(Global_2, aes(LureC, AVG_SA)) +geom_point(position=position_jitter(h=0.1, w=0.1), shape = 21, alpha = 0.5, size = 3, color="#01587A", fill="#077DAA")+geom_smooth(method = lm)
ggplot(Global_2, aes(LureD, AVG_SA)) +geom_point(position=position_jitter(h=0.1, w=0.1), shape = 21, alpha = 0.5, size = 3, color="#01587A", fill="#077DAA")+geom_smooth(method = lm)
ggplot(Global_2, aes(LureD, Secchi_Adj)) +geom_point(position=position_jitter(h=0.1, w=0.1), shape = 21, alpha = 0.5, size = 3, color="#01587A", fill="#077DAA")+geom_smooth(method = lm)


pairs(~LureA+CC+Secchi_Adj,data=Global_0.5, main="0.5m depth")
pairs(~LureA+CC,data=Global_1, main="1.0m depth")
pairs(~LureA+LureB+LureC+LureD+CC+WS+WH+AVG_SA+Secchi_Adj,data=Global_1.5, main="1.5m depth")
pairs(~LureB+LureC+LureD+WS+WH+AVG_SA+Secchi_Adj,data=Global_2, main="2.0 depth")


# PEARSONS correlation matrix for each depth 
#Matrix of correlation and p-value
cor_0.5<-rcorr(type="pearson", as.matrix(Global_0.5))
#Calculating correlation
M_0.5<-cor_0.5$r
#calculating P-value
p_mat_0.5<-cor_0.5$P
ggcorrplot(M_0.5, method="square", type="upper",tl.srt = 50,lab=TRUE,lab_size = 4,p.mat = p_mat_0.5, sig.level = 0.05,insig = "blank")+theme_update()+
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))+
  theme(axis.text.x = element_text(angle = 90, size = 10, colour = "black", vjust = 0.5, hjust = 1, face= "bold"),
        plot.margin = margin(0.5,0.5,0.5,0.5,"cm"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 10, face = "bold", colour = "black"),
        #legend.text = element_blank(),legend.position = "none",
        legend.text = element_text(size = 10, face = "bold", colour = "black"),
        legend.position = "bottom",legend.spacing = unit(0.5,"cm"),legend.key.width = unit(1, 'cm'),
        legend.key.height = unit(0.5, 'cm'),legend.spacing.x = unit(0.5,"cm"),legend.spacing.y = unit(0.5,"cm"),
        axis.text.y = element_text(colour = "black", size = 10, face = "bold"))+
  ggtitle("Pearson's R 0.5m")+ 
  theme(plot.title = element_text(size=15,face="bold",hjust = 0.5))

dev.off()


cor_1<-rcorr(as.matrix(Global_1))
#Calculating correlation
M_1<-cor_1$r
#calculating P-value
p_mat_1<-cor_1$P
ggcorrplot(M_1, method="square", type="upper",tl.srt = 50,lab=TRUE,lab_size = 4,p.mat = p_mat_1, sig.level = 0.05,insig = "blank")+theme_update()+
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))+
  theme(axis.text.x = element_text(angle = 90, size = 10, colour = "black", vjust = 0.5, hjust = 1, face= "bold"),
        plot.margin = margin(0.5,0.5,0.5,0.5,"cm"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 10, face = "bold", colour = "black"),
        #legend.text = element_blank(),legend.position = "none",
        legend.text = element_text(size = 10, face = "bold", colour = "black"),
        legend.position = "bottom",legend.spacing = unit(0.5,"cm"),legend.key.width = unit(1, 'cm'),
        legend.key.height = unit(0.5, 'cm'),legend.spacing.x = unit(0.5,"cm"),legend.spacing.y = unit(0.5,"cm"),
        axis.text.y = element_text(colour = "black", size = 10, face = "bold"))+
  ggtitle("Pearson's R 1m")+ 
  theme(plot.title = element_text(size=15,face="bold",hjust = 0.5))

dev.off()


#Matrix of correlation and p-value
cor_1.5<-rcorr(as.matrix(Global_1.5))
#Calculating correlation
M_1.5<-cor_1.5$r
#calculating P-value
p_mat_1.5<-cor_1.5$P
ggcorrplot(M_1.5, method="square", type="upper",tl.srt = 50,lab=TRUE,lab_size = 3,p.mat = p_mat_1.5, sig.level = 0.05,insig = "blank")+theme_update()+
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))+
  theme(axis.text.x = element_text(angle = 90, size = 10, colour = "black", vjust = 0.5, hjust = 1, face= "bold"),
        plot.margin = margin(0.5,0.5,0.5,0.5,"cm"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 10, face = "bold", colour = "black"),
        #legend.text = element_blank(),legend.position = "none",
        legend.text = element_text(size = 10, face = "bold", colour = "black"),
        legend.position = "bottom",legend.spacing = unit(0.5,"cm"),legend.key.width = unit(1, 'cm'),
        legend.key.height = unit(0.5, 'cm'),legend.spacing.x = unit(0.5,"cm"),legend.spacing.y = unit(0.5,"cm"),
        axis.text.y = element_text(colour = "black", size = 10, face = "bold"))+
  ggtitle("Pearson's R 1.5m")+ 
  theme(plot.title = element_text(size=15,face="bold",hjust = 0.5))

dev.off()


#Matrix of correlation and p-value
cor_2<-rcorr(as.matrix(Global_2))
#Calculating correlation
M_2<-cor_2$r
#calculating P-value
p_mat_2<-cor_2$P
ggcorrplot(M_2, method="square", type="upper",tl.srt = 50,lab=TRUE,lab_size = 4,p.mat = p_mat_2, sig.level = 0.05,insig = "blank")+theme_update()+
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))+
  theme(axis.text.x = element_text(angle = 90, size = 10, colour = "black", vjust = 0.5, hjust = 1, face= "bold"),
        plot.margin = margin(0.5,0.5,0.5,0.5,"cm"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 10, face = "bold", colour = "black"),
        #legend.text = element_blank(),legend.position = "none",
        legend.text = element_text(size = 10, face = "bold", colour = "black"),
        legend.position = "bottom",legend.spacing = unit(0.5,"cm"),legend.key.width = unit(1, 'cm'),
        legend.key.height = unit(0.5, 'cm'),legend.spacing.x = unit(0.5,"cm"),legend.spacing.y = unit(0.5,"cm"),
        axis.text.y = element_text(colour = "black", size = 10, face = "bold"))+
  ggtitle("Pearson's R 2m")+ 
  theme(plot.title = element_text(size=15,face="bold",hjust = 0.5))

dev.off()


### SPEARMANS correlation and p value

cor_0.5s<-rcorr(type="spearman", as.matrix(Global_0.5))
#Calculating correlation
M_0.5s<-cor_0.5s$r
#calculating P-value
p_mat_0.5s<-cor_0.5s$P
ggcorrplot(M_0.5s, method="square", type="upper",tl.srt = 50,lab=TRUE,lab_size = 6,p.mat = p_mat_0.5s, sig.level = 0.05,insig = "blank")+theme_update()+
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))+
  theme(axis.text.x = element_text(angle = 90, size = 10, colour = "black", vjust = 0.5, hjust = 1, face= "bold"),
        plot.margin = margin(0.5,0.5,0.5,0.5,"cm"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 10, face = "bold", colour = "black"),
        #legend.text = element_blank(),legend.position = "none",
        legend.text = element_text(size = 10, face = "bold", colour = "black"),
        legend.position = "bottom",legend.spacing = unit(0.5,"cm"),legend.key.width = unit(2, 'cm'),
        legend.key.height = unit(0.5, 'cm'),legend.spacing.x = unit(0.5,"cm"),legend.spacing.y = unit(0.5,"cm"),
        axis.text.y = element_text(colour = "black", size = 10, face = "bold"))+
  ggtitle("Spearman's R 0.5m")+ 
  theme(plot.title = element_text(size=20,face="bold",hjust = 0.5))

dev.off()

### 1m
cor_1s<-rcorr(type= "spearman", as.matrix(Global_1))
#Calculating correlation
M_1s<-cor_1s$r
#calculating P-value
p_mat_1s<-cor_1s$P
ggcorrplot(M_1s, method="square", type="upper",tl.srt = 50,lab=TRUE,lab_size = 6,p.mat = p_mat_1s, sig.level = 0.05,insig = "blank")+theme_update()+
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))+
  theme(axis.text.x = element_text(angle = 90, size = 10, colour = "black", vjust = 0.5, hjust = 1, face= "bold"),
        plot.margin = margin(0.5,0.5,0.5,0.5,"cm"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 10, face = "bold", colour = "black"),
        #legend.text = element_blank(),legend.position = "none",
        legend.text = element_text(size = 10, face = "bold", colour = "black"),
        legend.position = "bottom",legend.spacing = unit(0.5,"cm"),legend.key.width = unit(2, 'cm'),
        legend.key.height = unit(0.5, 'cm'),legend.spacing.x = unit(0.5,"cm"),legend.spacing.y = unit(0.5,"cm"),
        axis.text.y = element_text(colour = "black", size = 10, face = "bold"))+
  ggtitle("Spearman's R 1m")+ 
  theme(plot.title = element_text(size=20,face="bold",hjust = 0.5))

dev.off()


#Matrix of correlation and p-value
cor_1.5s<-rcorr(type="spearman",as.matrix(Global_1.5))
#Calculating correlation
M_1.5s<-cor_1.5s$r
#calculating P-value
p_mat_1.5s<-cor_1.5s$P
ggcorrplot(M_1.5s, method="square", type="upper",tl.srt = 50,lab=TRUE,lab_size = 4,p.mat = p_mat_1.5s, sig.level = 0.05,insig = "blank")+theme_update()+
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))+
  theme(axis.text.x = element_text(angle = 90, size = 10, colour = "black", vjust = 0.5, hjust = 1, face= "bold"),
        plot.margin = margin(0.5,0.5,0.5,0.5,"cm"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 10, face = "bold", colour = "black"),
        #legend.text = element_blank(),legend.position = "none",
        legend.text = element_text(size = 10, face = "bold", colour = "black"),
        legend.position = "bottom",legend.spacing = unit(0.5,"cm"),legend.key.width = unit(2, 'cm'),
        legend.key.height = unit(0.5, 'cm'),legend.spacing.x = unit(0.5,"cm"),legend.spacing.y = unit(0.5,"cm"),
        axis.text.y = element_text(colour = "black", size = 10, face = "bold"))+
  ggtitle("Spearman's R 1.5m")+ 
  theme(plot.title = element_text(size=20,face="bold",hjust = 0.5))

dev.off()


#Matrix of correlation and p-value
cor_2s<-rcorr(type= "spearman", as.matrix(Global_2))
#Calculating correlation
M_2s<-cor_2s$r
#calculating P-value
p_mat_2s<-cor_2s$P
ggcorrplot(M_2s, method="square", type="upper",tl.srt = 50,lab=TRUE,lab_size = 4,p.mat = p_mat_2s, sig.level = 0.05,insig = "blank")+theme_update()+
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"))+
  theme(axis.text.x = element_text(angle = 90, size = 10, colour = "black", vjust = 0.5, hjust = 1, face= "bold"),
        plot.margin = margin(0.5,0.5,0.5,0.5,"cm"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 10, face = "bold", colour = "black"),
        #legend.text = element_blank(),legend.position = "none",
        legend.text = element_text(size = 10, face = "bold", colour = "black"),
        legend.position = "bottom",legend.spacing = unit(0.5,"cm"),legend.key.width = unit(2, 'cm'),
        legend.key.height = unit(0.5, 'cm'),legend.spacing.x = unit(0.5,"cm"),legend.spacing.y = unit(0.5,"cm"),
        axis.text.y = element_text(colour = "black", size = 10, face = "bold"))+
  ggtitle("Spearmans's R 2m")+ 
  theme(plot.title = element_text(size=20,face="bold",hjust = 0.5))

dev.off()


