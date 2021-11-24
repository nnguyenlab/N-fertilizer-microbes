#This R code is for comparing beta diversity of prokaryotic samples in the radish experiment. The code includes:
#i Transforming abundance data based on feature table data of OTU frequencies by sample
#ii finding distances/similarity metrics
#iii performing statistics- PERMANOVAs and Pairwise Comparisons of those distances grouped by treatment and time
#iv plotting PCA ordinations in ggplot
#combining plots using grids

library(ggplot2) #graphing
library(grid) #graphing
library(gridBase) #graphing
library(gridExtra) #graphing
library(vegan) #vegdist, Adonis2
library(RVAideMemoire) #pairwise.perm.manova
library(stats) #prcomp
library(ggfortify) #for autoplot

#-------plot BACTERIA data------
#read csv files
otus.b.1=read.table(file="unrarefied-feature-table-b-1.tsv",sep="\t", header=TRUE, check.names=FALSE, row.names=1)
t.otus.b.1<-t(otus.b.1)
otus.b.1<-data.frame(t.otus.b.1)
metadata.b.1=read.csv(file='metadata-radish-16S_1.csv', header=TRUE, check.names=FALSE, row.names=1)

#this manually organizes the treatments for the legend
metadata.b.1$Treatment <- factor(metadata.b.1$Treatment, c("Compost","Urea","Control"))

#color pallete vector
pallete1 <- c("#a62900","#ffcb04","#036900")#kabir
pallete2<- c("#E69F00","#56B4E9","#000000")#Steven's pallete
pallete3 <- c("#de663e","#fecb00","#3e99d2")#(brown,yellow,blue)
pallete4 <- c("#de663e","#fecb00","#4470cf")#(brown,yellow,purple)
pallete7 <- c("#f45c85","#f2e127","#0a4d8c")#(pink,yellow,blue)
custom_labels <- theme(axis.text=element_text(size=8),axis.title=element_text(size=12),plot.title = element_text(hjust=1))
pallete <- scale_colour_manual(values=pallete7)

#apply clr standardization
stan.bac.1 <- decostand(otus.b.1, "hell")
#grab euclidean distances
dis.bac.1 <-vegdist(stan.bac.1, method= "euclidean")

#------------first ordination-----
#Creating the first ordinations. PCA chosen.
pca.bac.1<- autoplot(prcomp(dis.bac.1), data=metadata.b.1, show.legend=FALSE)

p.b.1<- pca.bac.1 + ggtitle("C1") +
  geom_point(aes(shape=Treatment, fill="black"), show.legend=FALSE, size=4) +
  geom_point(aes(shape=Treatment, color=metadata.b.1$Treatment), show.legend=FALSE, size=3) +
  geom_polygon(aes(color=metadata.b.1$Treatment), fill=NA,show.legend=FALSE,size=1) +
  theme_classic() + custom_labels + pallete
g7 <- ggplotGrob(p.b.1)#Create grob vector for future combining of graphs.

#Permutational Analysis of Variance. 
#first unnested
perm= how(nperm=999)
adonis2(dis.bac.1 ~ Treatment*TimePoint, data=metadata.b.1, permutations= perm)

#now nested
setBlocks(perm) <- with(metadata.b.1, TimePoint)
adonis2(dis.bac.1 ~ Treatment, data=metadata.b.1, permutations= perm)

#Now pairwise comparisons
pairwise.perm.manova(dis.bac.1, metadata.b.1$Treatment, test="Wilks", nperm=999, p.method="fdr")

#------------1Thru2-----------
#read csv files
otus.b.1thru2=read.table(file="unrarefied-feature-table-b-1thru2.tsv",sep="\t", header=TRUE, check.names=FALSE, row.names=1)
t.otus.b.1thru2<-t(otus.b.1thru2)
otus.b.1thru2<-data.frame(t.otus.b.1thru2)
metadata.b.1thru2=read.csv(file='metadata-radish-16S_1thru2.csv', header=TRUE, check.names=FALSE, row.names=1)
metadata.b.1thru2$Treatment <- factor(metadata.b.1$Treatment, c("Compost","Urea","Control"))

#apply clr standardization
stan.bac.1thru2 <- decostand(otus.b.1thru2, "hell")
#grab euclidean distances
dis.bac.1thru2 <-vegdist(stan.bac.1thru2, method= "euclidean")
pca.bac.1thru2<- autoplot(prcomp(dis.bac.1thru2), data=metadata.b.1thru2)

p.b.1thru2<- pca.bac.1thru2 + ggtitle("C1-2*") +
  geom_point(aes(shape=Treatment, fill="black"), show.legend=FALSE, size=4) +
  geom_point(aes(shape=Treatment, color=metadata.b.1thru2$Treatment), show.legend=FALSE, size=3) + 
  stat_ellipse(geom="polygon",aes(colour=metadata.b.1thru2$Treatment),fill=NA, linetype= 1, size= 1, show.legend=FALSE) + 
  theme_classic() + custom_labels + pallete
g8 <- ggplotGrob(p.b.1thru2)

#first none nested
perm= how(nperm=999)
adonis2(dis.bac.1thru2 ~ Treatment*TimePoint, data=metadata.b.1thru2, permutations= perm)
setBlocks(perm) <- with(metadata.b.1thru2, TimePoint)
adonis2(dis.bac.1thru2 ~ Treatment, data=metadata.b.1thru2, permutations= perm)
pairwise.perm.manova(dis.bac.1thru2, metadata.b.1thru2$Treatment, test="Wilks", nperm=999, p.method="fdr")
pairwise.perm.manova(dis.bac.1thru2, metadata.b.1thru2$TimePoint, test="Wilks", nperm=999, p.method="fdr")

#------------1Thru3-----------
#read csv files
otus.b.1thru3=read.table(file="unrarefied-feature-table-b-1thru3.tsv",sep="\t", header=TRUE, check.names=FALSE, row.names=1)
t.otus.b.1thru3<-t(otus.b.1thru3)
otus.b.1thru3<-data.frame(t.otus.b.1thru3)
metadata.b.1thru3=read.csv(file='metadata-radish-16S_1thru3.csv', header=TRUE, check.names=FALSE, row.names=1)
metadata.b.1thru3$Treatment <- factor(metadata.b.1thru3$Treatment, c("Compost","Urea","Control"))
stan.bac.1thru3 <- decostand(otus.b.1thru3, "hell")
#grab euclidean distances
dis.bac.1thru3 <-vegdist(stan.bac.1thru3, method= "euclidean")
pca.bac.1thru3<- autoplot(prcomp(dis.bac.1thru3), data=metadata.b.1thru3)

p.b.1thru3<- pca.bac.1thru3 + ggtitle("C1-3*") +
  geom_point(aes(shape=Treatment, fill="black"), show.legend=FALSE,size=4) +
  geom_point(aes(shape=Treatment, color=metadata.b.1thru3$Treatment), show.legend=FALSE,size=3) + 
  stat_ellipse(geom="polygon",aes(colour=metadata.b.1thru3$Treatment),fill=NA, linetype= 1, size= 1, show.legend=FALSE) + 
  theme_classic() + custom_labels + pallete
g9 <- ggplotGrob(p.b.1thru3)

#first none nested
perm= how(nperm=999)
adonis2(dis.bac.1thru3 ~ Treatment*TimePoint, data=metadata.b.1thru3, permutations= perm)

setBlocks(perm) <- with(metadata.b.1thru3, TimePoint)
adonis2(dis.bac.1thru3 ~ Treatment, data=metadata.b.1thru3, permutations= perm)
pairwise.perm.manova(dis.bac.1thru3, metadata.b.1thru3$Treatment, test="Wilks", nperm=999, p.method="fdr")
pairwise.perm.manova(dis.bac.1thru3, metadata.b.1thru3$TimePoint, test="Wilks", nperm=999, p.method="fdr")

#------------1Thru4-----------
#read csv files
otus.b.1thru4=read.table(file="unrarefied-feature-table-b-1thru4.tsv",sep="\t", header=TRUE, check.names=FALSE, row.names=1)
t.otus.b.1thru4<-t(otus.b.1thru4)
otus.b.1thru4<-data.frame(t.otus.b.1thru4)
metadata.b.1thru4=read.csv(file='metadata-radish-16S_1thru4.csv', header=TRUE, check.names=FALSE, row.names=1)
metadata.b.1thru4$Treatment <- factor(metadata.b.1thru4$Treatment, c("Compost","Urea","Control"))

#apply standardization
stan.bac.1thru4 <- decostand(otus.b.1thru4, "hell")
#grab euclidean distances
dis.bac.1thru4 <-vegdist(stan.bac.1thru4, method="euclidean")
pca.bac.1thru4<- autoplot(prcomp(dis.bac.1thru4), data=metadata.b.1thru4)

p.b.1thru4<- pca.bac.1thru4 + ggtitle("C1-4*") +
  geom_point(aes(shape=Treatment, fill="black"), show.legend=FALSE,size=4) +
  geom_point(aes(shape=Treatment, color=metadata.b.1thru4$Treatment), show.legend=FALSE, size=3) + 
  stat_ellipse(geom="polygon",aes(colour=metadata.b.1thru4$Treatment),fill=NA, linetype= 1, size= 1, show.legend=FALSE) + 
  theme_classic() + custom_labels + pallete
g10 <- ggplotGrob(p.b.1thru4)

#first none nested
perm= how(nperm=999)
adonis2(dis.bac.1thru4 ~ Treatment*TimePoint, data=metadata.b.1thru4, permutations= perm)

setBlocks(perm) <- with(metadata.b.1thru4, TimePoint)
adonis2(dis.bac.1thru4 ~ Treatment, data=metadata.b.1thru4, permutations= perm)
pairwise.perm.manova(dis.bac.1thru4, metadata.b.1thru4$Treatment, test="Wilks", nperm=999, p.method="fdr")
pairwise.perm.manova(dis.bac.1thru4, metadata.b.1thru4$TimePoint, test="Wilks", nperm=999, p.method="fdr")

#------------1Thru5-----------
#read csv files
otus.b.1thru5=read.table(file="unrarefied-feature-table-b-1thru5.tsv",sep="\t", header=TRUE, check.names=FALSE, row.names=1)
t.otus.b.1thru5<-t(otus.b.1thru5)
otus.b.1thru5<-data.frame(t.otus.b.1thru5)
metadata.b.1thru5=read.csv(file='metadata-radish-16S_1thru5.csv', header=TRUE, check.names=FALSE, row.names=1)
metadata.b.1thru5$Treatment <- factor(metadata.b.1thru5$Treatment, c("Compost","Urea","Control"))

#apply clr standardization
stan.bac.1thru5 <- decostand(otus.b.1thru5, "hell")
#grab euclidean distances
dis.bac.1thru5 <-vegdist(stan.bac.1thru5, method= "euclidean")
pca.bac.1thru5<- autoplot(prcomp(dis.bac.1thru5), data=metadata.b.1thru5)

p.b.1thru5<- pca.bac.1thru5 + ggtitle("C1-5*") +
  geom_point(aes(shape=Treatment, fill="black"), show.legend=FALSE, size=4) +
  geom_point(aes(shape=Treatment, color=metadata.b.1thru5$Treatment), show.legend=FALSE, size=3) + 
  stat_ellipse(geom="polygon", aes(colour=metadata.b.1thru5$Treatment), fill=NA, linetype= 1, size= 1, show.legend=FALSE) + 
  theme_classic() + custom_labels + pallete
g11 <- ggplotGrob(p.b.1thru5)

#first none nested
perm= how(nperm=999)
adonis2(dis.bac.1thru5 ~ Treatment*TimePoint, data=metadata.b.1thru5, permutations= perm)
setBlocks(perm) <- with(metadata.b.1thru5, TimePoint)
adonis2(dis.bac.1thru5 ~ Treatment, data=metadata.b.1thru5, permutations= perm)
pairwise.perm.manova(dis.bac.1thru5, metadata.b.1thru5$Treatment, test="Wilks", nperm=999, p.method="fdr")
pairwise.perm.manova(dis.bac.1thru5, metadata.b.1thru5$TimePoint, test="Wilks", nperm=999, p.method="fdr")

#------------1Thru6-----------
otus.b.1thru6=read.table(file="unrarefied-feature-table-b-1thru6.tsv",sep="\t", header=TRUE, check.names=FALSE, row.names=1)
t.otus.b.1thru6<-t(otus.b.1thru6)
otus.b.1thru6<-data.frame(t.otus.b.1thru6)
metadata.b.1thru6=read.csv(file='metadata-radish-16S_1thru6.csv', header=TRUE, check.names=FALSE, row.names=1)
metadata.b.1thru6$Treatment <- factor(metadata.b.1thru6$Treatment, c("Compost","Urea","Control"))

#apply clr standardization
stan.bac.1thru6 <- decostand(otus.b.1thru6, "hell")
#grab euclidean distances
dis.bac.1thru6 <-vegdist(stan.bac.1thru6, method= "euclidean")
pca.bac.1thru6<- autoplot(prcomp(dis.bac.1thru6), data=metadata.b.1thru6)

p.b.1thru6<- pca.bac.1thru6 + ggtitle("C1-6*") +
  geom_point(aes(shape=Treatment, fill="black"), show.legend=FALSE, size=4) +
  geom_point(aes(shape=Treatment, color=metadata.b.1thru6$Treatment), show.legend=FALSE, size=3) + 
  stat_ellipse(geom="polygon",aes(colour=metadata.b.1thru6$Treatment),fill=NA, linetype= 1, size= 1, show.legend=FALSE) + 
  theme_classic() + custom_labels + pallete
g12 <- ggplotGrob(p.b.1thru6)

#first none nested
perm= how(nperm=999)
adonis2(dis.bac.1thru6 ~ Treatment*TimePoint, data=metadata.b.1thru6, permutations= perm)

setBlocks(perm) <- with(metadata.b.1thru6, TimePoint)
adonis2(dis.bac.1thru6 ~ Treatment, data=metadata.b.1thru6, permutations= perm)
pairwise.perm.manova(dis.bac.1thru6, metadata.b.1thru6$Treatment, test="Wilks", nperm=999, p.method="fdr")
pairwise.perm.manova(dis.bac.1thru6, metadata.b.1thru6$TimePoint, test="Wilks", nperm=999, p.method="fdr")

#===========Plotting FUNGAL data===========
#read csv files
otus.f.1=read.csv(file="feature-table-f-1.csv", header=TRUE, check.names=FALSE, row.names=1)
metadata.f.1=read.csv(file='metadata-radish-ITS_1.csv', header=TRUE, check.names=FALSE, row.names=1)
metadata.f.1$Treatment <- factor(metadata.f.1$Treatment, c("Compost","Urea","Control"))

#------------first ordination-----
#apply clr standardization
stan.fun.1 <- decostand(otus.f.1, "hell")
#grab euclidean distances
dis.fun.1 <-vegdist(stan.fun.1, method="euclidean")

#ordination plot. chose PCA 
pca.fun.1<- autoplot(prcomp(dis.fun.1), data=metadata.f.1, show.legend=FALSE)

p.f.1<- pca.fun.1 + ggtitle("C1") +
  geom_point(aes(shape=Treatment), show.legend=FALSE, size=4) +
  geom_point(aes(shape=Treatment, color=metadata.f.1$Treatment), show.legend=FALSE, size=3) +
  geom_polygon(aes(color=metadata.f.1$Treatment), fill=NA, show.legend=FALSE, size=1) + 
  theme_classic() + custom_labels + pallete
g1 <- ggplotGrob(p.f.1)#creating a graph vector for later combining

#------------1Thru2-----------
#read csv files
otus.f.1thru2=read.csv(file="feature-table-f-1thru2.csv", header=TRUE, check.names=FALSE, row.names=1)
metadata.f.1thru2=read.csv(file='metadata-radish-ITS_1thru2.csv', header=TRUE, check.names=FALSE, row.names=1)
metadata.f.1thru2$Treatment <- factor(metadata.f.1thru2$Treatment, c("Compost","Urea","Control"))

#apply clr standardization
stan.fun.1thru2 <- decostand(otus.f.1thru2, "hell")
#grab euclidean distances
dis.fun.1thru2 <-vegdist(stan.fun.1thru2, method= "euclidean")

#first PCA Graph
pca.fun.1thru2 <- autoplot(prcomp(dis.fun.1thru2), data=metadata.f.1thru2)

p.f.1thru2 <- pca.fun.1thru2 + ggtitle("C1-2*") +
  geom_point(aes(shape=Treatment), show.legend=FALSE, size=4) +
  geom_point(aes(shape=Treatment, color=metadata.f.1thru2$Treatment), show.legend=FALSE,size=3) + 
  stat_ellipse(geom="polygon",aes(colour=metadata.f.1thru2$Treatment), fill=NA, linetype=1, size=1, show.legend=FALSE) + 
  theme_classic() + custom_labels + pallete
g2 <- ggplotGrob(p.f.1thru2)


#first none nested
perm= how(nperm=999)
adonis2(dis.fun.1thru2 ~ Treatment*TimePoint, data=metadata.f.1thru2, permutations= perm)
setBlocks(perm) <- with(metadata.f.1thru2, TimePoint)
adonis2(dis.fun.1thru2 ~ Treatment, data=metadata.f.1thru2, permutations= perm)
pairwise.perm.manova(dis.fun.1thru2, metadata.f.1thru2$Treatment, test="Wilks", nperm=999, p.method="fdr")
pairwise.perm.manova(dis.fun.1thru2, metadata.f.1thru2$TimePoint, test="Wilks", nperm=999, p.method="fdr")

#------------1thru3---------
#read csv files
otus.f.1thru3=read.csv(file="feature-table-f-1thru3.csv", header=TRUE, check.names=FALSE, row.names=1)
metadata.f.1thru3=read.csv(file='metadata-radish-ITS_1thru3.csv', header=TRUE, check.names=FALSE, row.names=1)
metadata.f.1thru3$Treatment <- factor(metadata.f.3$Treatment, c("Compost","Urea","Control"))

#apply clr standardization
stan.fun.1thru3 <- decostand(otus.f.1thru3, "hell")
#grab euclidean distances
dis.fun.1thru3 <-vegdist(stan.fun.1thru3, method="euclidean")
pca.fun.1thru3<- autoplot(prcomp(dis.fun.1thru3), data=metadata.f.1thru3)

p.f.1thru3<- pca.fun.1thru3 + ggtitle("C1-3*") +
  geom_point(aes(shape=Treatment, fill="black"), show.legend=FALSE, size=4) +
  geom_point(aes(shape=Treatment, color=metadata.f.1thru3$Treatment), show.legend=FALSE, size=3) + 
  stat_ellipse(geom="polygon",aes(colour=metadata.f.1thru3$Treatment),fill=NA, linetype=1, size=1, show.legend=FALSE) + 
  theme_classic() + custom_labels + pallete
g3 <- ggplotGrob(p.f.1thru3)

pairwise.perm.manova(dis.fun.1thru3, metadata.f.1thru3$Treatment, test="Wilks", nperm=999, p.method="fdr")
#first none nested
perm= how(nperm=999)
adonis2(dis.fun.1thru3 ~ Treatment*TimePoint, data=metadata.f.1thru3, permutations= perm)
setBlocks(perm) <- with(metadata.f.1thru3, TimePoint)
adonis2(dis.fun.1thru3 ~ Treatment, data=metadata.f.1thru3, permutations= perm)

#------------1thru4-----------------
#read csv files
otus.f.1thru4=read.csv(file="feature-table-f-1thru4.csv", header=TRUE, check.names=FALSE, row.names=1)
metadata.f.1thru4=read.csv(file='metadata-radish-ITS_1thru4.csv', header=TRUE, check.names=FALSE, row.names=1)
metadata.f.1thru4$Treatment <- factor(metadata.f.4$Treatment, c("Compost","Urea","Control"))

#apply clr standardization
stan.fun.1thru4 <- decostand(otus.f.1thru4, "norm")
#grab euclidean distances
dis.fun.1thru4 <-vegdist(stan.fun.1thru4, method= "euclidean")

pca.fun.1thru4<- autoplot(prcomp(dis.fun.1thru4), data=metadata.f.1thru4)

p.f.1thru4<- pca.fun.1thru4 + ggtitle("C1-4*") +
  geom_point(aes(shape=Treatment, fill="black"), show.legend=FALSE, size=4) +
  geom_point(aes(shape=Treatment, color=metadata.f.1thru4$Treatment), show.legend=FALSE, size=3) + 
  stat_ellipse(geom="polygon", aes(colour=metadata.f.1thru4$Treatment), fill=NA, linetype=1, size=1, show.legend=FALSE) + 
  theme_classic() + custom_labels + pallete
g4 <- ggplotGrob(p.f.1thru4)


#first none nested
perm= how(nperm=999)
adonis2(dis.fun.1thru4 ~ Treatment*TimePoint, data=metadata.f.1thru4, permutations= perm)

setBlocks(perm) <- with(metadata.f.1thru4, TimePoint)
adonis2(dis.fun.1thru4 ~ Treatment, data=metadata.f.1thru4, permutations= perm)
pairwise.perm.manova(dis.fun.1thru4, metadata.f.1thru4$Treatment, test="Wilks", nperm=999, p.method="fdr")
pairwise.perm.manova(dis.fun.1thru4, metadata.f.1thru4$TimePoint, test="Wilks", nperm=999, p.method="fdr")

#------------1thru5-----------------
#read csv files
otus.f.1thru5=read.csv(file="feature-table-f-1thru5.csv", header=TRUE, check.names=FALSE, row.names=1)
metadata.f.1thru5=read.csv(file='metadata-radish-ITS_1thru5.csv', header=TRUE, check.names=FALSE, row.names=1)
metadata.f.1thru5$Treatment <- factor(metadata.f.5$Treatment, c("Compost","Urea","Control"))

#apply clr standardization
stan.fun.1thru5 <- decostand(otus.f.1thru5, "hell")
#grab euclidean distances
dis.fun.1thru5 <-vegdist(stan.fun.1thru5, method= "euclidean")

pca.fun.1thru5<- autoplot(prcomp(dis.fun.1thru5), data=metadata.f.1thru5)

p.f.1thru5<- pca.fun.1thru5 + ggtitle("C1-5*") +
  geom_point(aes(shape=Treatment, fill="black"), show.legend=FALSE, size=4) +
  geom_point(aes(shape=Treatment, color=metadata.f.1thru5$Treatment), show.legend=FALSE, size=3) + 
  stat_ellipse(geom="polygon",aes(colour=metadata.f.1thru5$Treatment), fill=NA, linetype=1, size=1, show.legend=FALSE) + 
  theme_classic() + custom_labels + pallete
g5 <- ggplotGrob(p.f.1thru5)


#first none nested
perm= how(nperm=999)
adonis2(dis.fun.1thru5 ~ Treatment*TimePoint, data=metadata.f.1thru5, permutations= perm)

setBlocks(perm) <- with(metadata.f.1thru5, TimePoint)
adonis2(dis.fun.1thru5 ~ Treatment, data=metadata.f.1thru5, permutations= perm)

pairwise.perm.manova(dis.fun.1thru5, metadata.f.1thru5$Treatment, test="Wilks", nperm=999, p.method="fdr")
#------------1thru6-----
#read csv files
otus.f.1thru6=read.csv(file="feature-table-f-1thru6.csv", header=TRUE, check.names=FALSE, row.names=1)
metadata.f.1thru6=read.csv(file='metadata-radish-ITS_1thru6.csv', header=TRUE, check.names=FALSE, row.names=1)
metadata.f.1thru6$Treatment <- factor(metadata.f.6$Treatment, c("Compost","Urea","Control"))

#apply clr standardization
stan.fun.1thru6 <- decostand(otus.f.1thru6, "hell")
#grab euclidean distances
dis.fun.1thru6 <-vegdist(stan.fun.1thru6, method= "euclidean")
pca.fun.1thru6<- autoplot(prcomp(dis.fun.1thru6), data=metadata.f.1thru6)

p.f.1thru6<- pca.fun.1thru6 + ggtitle("C1-6*") +
  geom_point(aes(shape=Treatment, fill="black"), show.legend=FALSE, size=4) +
  geom_point(aes(shape=Treatment, color=metadata.f.1thru6$Treatment), show.legend=FALSE, size=3) + 
  stat_ellipse(geom="polygon", aes(colour=metadata.f.1thru6$Treatment), fill=NA, linetype=1, size=1, show.legend=FALSE) + 
  theme_classic() + custom_labels + pallete
g6 <- ggplotGrob(p.f.1thru6)

#first none nested
perm= how(nperm=999)
adonis2(dis.fun.1thru6 ~ Treatment*TimePoint, data=metadata.f.1thru6, permutations= perm)

setBlocks(perm) <- with(metadata.f.1thru6, TimePoint)
adonis2(dis.fun.1thru6 ~ Treatment, data=metadata.f.1thru6, permutations= perm)

pairwise.perm.manova(dis.fun.1thru6, metadata.f.1thru6$Treatment, test="Wilks", nperm=999, p.method="fdr")



#------------MAKING FIGURE 3 ----------------
pca.fun.leg<- autoplot(prcomp(dis.fun.1thru6), data=metadata.f.1thru6)
p.f.1leg<- pca.fun.1thru6 + theme_classic()+ scale_fill_viridis_d() + geom_point(aes(shape=Treatment)) + stat_ellipse(geom="polygon",aes(colour=Treatment),fill=NA, linetype= 1, size= 1) + scale_colour_manual(values=pallete2)

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

new.legend<-g_legend(p.f.1leg)

#standardizing width, probably not needed for this example
#g2$heights<-unit.pmax(g3$heights)  
#g3$heights<-unit.pmax(g3$heights) 
#g4$heights<-unit.pmax(g3$heights) 
#g5$heights<-unit.pmax(g3$heights) 
#g6$heights<-unit.pmax(g3$heights) 
#g2$widths<-unit.pmax(g3$widths)  
#g3$widths<-unit.pmax(g3$widths)   
#g4$widths<-unit.pmax(g3$widths)  
#g5$widths<-unit.pmax(g3$widths)  
#g6$widths<-unit.pmax(g3$widths)  

#plot all
grid.arrange(g1,g2,g3,g4,g5,g6,nrow=2,ncol=3)


#Making figure for significant predicted values
#with legend to grab
pca.bac.leg<- autoplot(prcomp(dis.bac.1thru6), data=metadata.b.1thru6)
p.bac.leg<- pca.bac.leg + 
  theme_classic() + 
  geom_point(aes(shape=Treatment, fill="black"),size=3,show.legend=FALSE) +
  geom_point(aes(shape=Treatment, color=Treatment), size=1.5,show.legend=FALSE) + 
  stat_ellipse(geom="polygon",aes(color=Treatment),fill=NA, linetype= 1, size= 1) + 
  scale_colour_manual(values= pallete2)

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

new.legend<-g_legend(p.bac.leg)
plot(new.legend)



#standardizing width, probably not needed for this example
#g7$heights<-unit.pmax(g8$heights)  
#g8$heights<-unit.pmax(g8$heights) 
#g9$heights<-unit.pmax(g8$heights) 
#g10$heights<-unit.pmax(g8$heights) 
#g11$heights<-unit.pmax(g8$heights) 
#g12$heights<-unit.pmax(g8$heights) 
 
#g7$widths<-unit.pmax(g8$widths)  
#g8$widths<-unit.pmax(g8$widths)   
#g9$widths<-unit.pmax(g8$widths)  
#g10$widths<-unit.pmax(g8$widths)  
#g11$widths<-unit.pmax(g8$widths) 
#g12$widths<-unit.pmax(g8$widths) 

#plot
grid.arrange(g7,g8,g9,g10,g11,g12,nrow=2,ncol=3)

ggtitle()
#Bacteria and Fungi together in one plot
grid.arrange(g7,g8,g9,g10,g11,g12,g1,g2,g3,g4,g5,g6,nrow=4,ncol=3)
