## Libraries 
library(ggplot2)
library(scales)
library(animation)
library(tidyr)
library(ggpubr)
library(devtools)
library(dplyr)
library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) # clustering algorithms & visualization
library(partykit)
library(sf)
library(tmap)    # for static and interactive maps
library(leaflet) # for interactive maps
library(FactoClass)
library(factoextra)

## Creating working directory and adding the data in
dir <- "C:/Users/nicol/Documents/experimental_analysis"
setwd(dir)
df<-read.table("GlobalAnalysis.csv", header = TRUE, sep = ",")

## removing all data with empty cells (Ale said this was necessary)
df <- (na.omit(df))

## scaling all the variables 
rescale_df <- df %>%
  mutate(WS_scal = scale(WS),
         WH_scal = scale(WH),
         CC_scal = scale(CC),
         SA_scal = scale(AVG_SA),
         SD_scal = scale(Secchi_Adj))

# final check to make sure that you have no na values
rescale_df <- (na.omit(rescale_df))

#goal in this section of code (above) is to only take the rescaled columns
rescale_df<- select(rescale_df, c(WS_scal, WH_scal, CC_scal, SA_scal, SD_scal))

##Choosing how many clusters to have
set.seed(123)
fviz_nbclust(rescale_df, kmeans, method = "wss")
fviz_nbclust(rescale_df, kmeans, method = "silhouette") 
fviz_nbclust(rescale_df, kmeans, method = "gap_stat")

# Conclusion: I think 5 classes makes sense

## Cluster analysis
#Goal is to minimize SS within groups and maximize between group variation. Groups should have similar sample sizes. The plot is good to visualize thier overlap (minimal overlap is best)
#test with 5 classes
set.seed(123)
final <- kmeans(rescale_df, 5, nstart = 25)
print(final)

#final plot
fviz_cluster(final, data = rescale_df, lablesize=0.1)

rescale_df %>%
  mutate(Cluster = final$cluster) %>%
  group_by(Cluster) %>%
  summarise_all("mean")


## Put cluster as a new column back into both data frames so you can see the non-scaled data in thier clusters
df$cluster <- final$cluster
df$cluster =as.factor(df$cluster)
rescale_df$cluster <- final$cluster
rescale_df$cluster =as.factor(rescale_df$cluster)

## Create box and wisker plots for each cluster for each variable 
p <- ggplot(df, aes(x=cluster, y=CC, fill= cluster)) + geom_boxplot()
p
q <- ggplot(df, aes(x=cluster, y=WS, fill= cluster)) +  geom_boxplot()
q
r <- ggplot(df, aes(x=cluster, y=WH, fill= cluster)) + geom_boxplot()
r
s <- ggplot(df, aes(x=cluster, y=AVG_SA, fill= cluster)) +  geom_boxplot()
s
t <- ggplot(df, aes(x=cluster, y=Secchi_Adj, fill= cluster)) + geom_boxplot()
t
# Box and wisker plots for each cluster (this includes all depths, so probably not the best measure)
u <- ggplot(df, aes(x=cluster, y=Global_P, fill= cluster)) + geom_boxplot()
u

# Subseting each depth
Global_0.5<- subset(df, df$Depth == 0.5)
Global_1 <- subset(df, df$Depth == 1)
Global_1.5 <- subset(df, df$Depth == 1.5)
Global_2 <- subset(df, df$Depth == 2)

# box and wisker plot for each depth
v <- ggplot(Global_0.5, aes(x=cluster, y=Global_P, fill= cluster)) + geom_boxplot()
v
w <- ggplot(Global_1, aes(x=cluster, y=Global_P, fill= cluster)) + geom_boxplot()
w
x <- ggplot(Global_1.5, aes(x=cluster, y=Global_P, fill= cluster)) + geom_boxplot()
x
y <- ggplot(Global_2, aes(x=cluster, y=Global_P, fill= cluster)) + geom_boxplot()
y

### SPECIES SPECIFIC
## SP A
# ALL counts, not specific to depth
A <- ggplot(df, aes(x=cluster, y=SPA_P, fill= cluster)) + geom_boxplot()

# for each depth
B <- ggplot(Global_0.5, aes(x=cluster, y=LureA, fill= cluster)) + geom_boxplot()+ylim(0, 100)
B + labs( y="0.5m" )
C <- ggplot(Global_1, aes(x=cluster, y=LureA, fill= cluster)) + geom_boxplot()+ylim(0, 100)
C + labs( y="1m")
D <- ggplot(Global_1.5, aes(x=cluster, y=LureA, fill= cluster)) + geom_boxplot()+ylim(0, 100)
D + labs( y="1.5m")
E <- ggplot(Global_2, aes(x=cluster, y=LureA, fill= cluster)) + geom_boxplot()+ylim(0, 100)
E+ labs(y="2m")

## SP B
# ALL counts, not specific to depth
G <- ggplot(df, aes(x=cluster, y=SPB_P, fill= cluster)) + geom_boxplot()
G
# for each depth
H <- ggplot(Global_0.5, aes(x=cluster, y=LureB, fill= cluster)) + geom_boxplot()+ylim(0, 100)
H
I <- ggplot(Global_1, aes(x=cluster, y=LureB, fill= cluster)) + geom_boxplot()+ylim(0, 100)
I
J <- ggplot(Global_1.5, aes(x=cluster, y=LureB, fill= cluster)) + geom_boxplot()+ylim(0, 100)
J
K <- ggplot(Global_2, aes(x=cluster, y=LureB, fill= cluster)) + geom_boxplot()+ylim(0, 100)
K

## SP C
# ALL counts, not specific to depth
L <- ggplot(df, aes(x=cluster, y=SPC_P, fill= cluster)) + geom_boxplot()
L
# for each depth
M <- ggplot(Global_0.5, aes(x=cluster, y=LureC, fill= cluster)) + geom_boxplot()+ylim(0, 100)
M+ labs(x="Cluster", y="Lure visibility (%)")
N <- ggplot(Global_1, aes(x=cluster, y=LureC, fill= cluster)) + geom_boxplot()+ylim(0, 100)
N+ labs(x="Cluster", y="Lure visibility (%)")
O <- ggplot(Global_1.5, aes(x=cluster, y=LureC, fill= cluster)) + geom_boxplot()+ylim(0, 100)
O+ labs(x="Cluster", y="Lure visibility (%)")
P <- ggplot(Global_2, aes(x=cluster, y=LureC, fill= cluster)) + geom_boxplot()+ylim(0, 100)
P+ labs(x="Cluster", y="Lure visibility (%)")

## SP D
# ALL counts, not specific to depth
Q <- ggplot(df, aes(x=cluster, y=SPD_P, fill= cluster)) + geom_boxplot()
Q
# for each depth
R <- ggplot(Global_0.5, aes(x=cluster, y=LureD, fill= cluster)) + geom_boxplot()+ylim(0, 100)
R+ labs(x="Cluster", y="Lure visibility (%)")
S <- ggplot(Global_1, aes(x=cluster, y=LureD, fill= cluster)) + geom_boxplot()+ylim(0, 100)
S+ labs(x="Cluster", y="Lure visibility (%)")
U <- ggplot(Global_1.5, aes(x=cluster, y=LureD, fill= cluster)) + geom_boxplot()+ylim(0, 100)
U+ labs(x="Cluster", y="Lure visibility (%)")
V <- ggplot(Global_2, aes(x=cluster, y=LureD, fill= cluster)) + geom_boxplot()+ylim(0, 100)
V+ labs(x="Cluster", y="Lure visibility (%)")

