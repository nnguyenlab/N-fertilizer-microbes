#amended by NN 7/25/21

#Set and check the working directory.
setwd("/Volumes/GoogleDrive/Shared drives/Radish amendment microbes/Manuscript/Prelim Data and Figures")
getwd()

library(Rmisc)
library(ggplot2)
library(utils) #sum
library(stats) #AOV
library(agricolae) #LSD test
library(gmodels) #contrasts
library(dplyr) #ggplot
library(lsmeans)
library(nlme)
library(rcompanion)
library(grid)
library(gridBase)
library(gridExtra) #for grid.arrange

#import excel sheets
library(readxl)
Diversity_Metrics_bacteria <- read_excel("Diversity_Metrics_bacteria.xlsx")
Diversity_Metrics_fungi <- read_excel("Diversity_Metrics_fungi.xlsx")

data_b <-Diversity_Metrics_bacteria
data_f <- Diversity_Metrics_fungi

pallete1 <- c("#a62900","#ffcb04","#036900")#kabir
pallete2<- c("#E69F00","#56B4E9","#000000")#Steven's pallete
pallete3 <- c("#de663e","#fecb00","#3e99d2")#(brown,yellow,blue)
pallete4 <- c("#de663e","#fecb00","#4470cf")#(brown,yellow,purple)
pallete7 <- c("#f45c85","#f2e127","#0a4d8c")#(pink,yellow,blue)
custom_labels <- theme(axis.text=element_text(size=11),axis.title=element_text(size=14), plot.title=element_text(size=20))

#------------------------------------------------------------
#summarize each attribute
#for bacteria
sum.richness <- summarySE(data_b, measurevar="richness", groupvars=c("Treatment","TimePoint"))
sum.taxa <- summarySE(data_b, measurevar="shannon", groupvars=c("Treatment","TimePoint"))
sum.phylo <- summarySE(data_b, measurevar="faith_pd", groupvars=c("Treatment","TimePoint"))

#for fungi
sum.f.richness <- summarySE(data_f, measurevar="richness", groupvars=c("Treatment","TimePoint"))
sum.f.taxa <- summarySE(data_f, measurevar="shannon", groupvars=c("Treatment","TimePoint"))

#-=================GRAPHING================
#Graph each attribute for bacteria
#pd will be a position dodge for the points so that they don't overlap
pd <- position_dodge(0.65)

#--------next plot richness---------
#First we are turning the character column into a vector so we can arrange and color the legend as we see fit
#Turn your 'treatment' column into a character vector
sum.richness$Treatment <- as.character(sum.richness$Treatment)
#Then turn it back into a factor with the levels in the correct order
sum.richness$Treatment <- factor(sum.richness$Treatment, levels=c("Compost", "Urea", "Control"))

p1 <- ggplot(sum.richness, aes(x=TimePoint, y=richness, group=Treatment, shape=Treatment, color=Treatment)) +
  geom_errorbar(aes(ymin=richness-se, ymax=richness+se), colour="#696268", width=.1, position=pd) +
  geom_point(color="black", size=5, show.legend=FALSE, position=pd) +
  geom_point(size=4, show.legend=FALSE, position=pd) +
  scale_colour_manual(values=pallete7) +
  labs(title="A", y="Observed OTU richness", x=NULL) +
  theme_classic() + custom_labels

#--------now taxonomic diversity------
#Turn your 'treatment' column into a character vector
sum.taxa$Treatment <- as.character(sum.taxa$Treatment)
#Then turn it back into a factor with the levels in the correct order
sum.taxa$Treatment <- factor(sum.taxa$Treatment, levels=c("Compost", "Urea", "Control"))

p2 <- ggplot(sum.taxa, aes(x=TimePoint, y=shannon, group=Treatment, shape=Treatment, color=Treatment)) +
  geom_errorbar(aes(ymin=shannon-se, ymax=shannon+se), colour="#696268", width=.1, position=pd) +
  geom_point(color="black", size=5, show.legend=FALSE, position=pd) +
  geom_point(size=4, show.legend=FALSE, position=pd) +
  scale_colour_manual(values=pallete7) +
  labs(title=NULL, y="Shannon Diversity (H')", x="Crop Cycle") +
  theme_classic() + custom_labels

#--------now phylogenetic diversity------
#Turn your 'treatment' column into a character vector
sum.phylo$Treatment <- as.character(sum.phylo$Treatment)
#Then turn it back into a factor with the levels in the correct order
sum.phylo$Treatment <- factor(sum.phylo$Treatment, levels=c("Compost", "Urea", "Control"))

p3 <- ggplot(sum.phylo, aes(x=TimePoint, y=faith_pd, group=Treatment, shape=Treatment, color=Treatment)) +
  geom_errorbar(aes(ymin=faith_pd-se, ymax=faith_pd+se), colour="#696268", width=.1, position=pd) +
  geom_point(color="black", size=5, show.legend=FALSE, position=pd) +
  geom_point(size=4, show.legend=FALSE, position=pd) +
  scale_colour_manual(values=pallete7) +
  labs(title=NULL, y="Faith's PD", x=NULL) +
  theme_classic() + custom_labels


#==========Graphing the Fungi============
#Graph each attribute
#pd will be a position dodge for the points
pd <- position_dodge(0.65)

#-----------plot richness----------------
#First we are turning the character column into a vector so we can arrange and color the legend as we see fit
#Turn your 'treatment' column into a character vector
sum.f.richness$Treatment <- as.character(sum.f.richness$Treatment)
#Then turn it back into a factor with the levels in the correct order
sum.f.richness$Treatment <- factor(sum.f.richness$Treatment, levels=c("Compost", "Urea", "Control"))

p4 <- ggplot(sum.f.richness, aes(x=TimePoint, y=richness, group=Treatment, shape=Treatment, color=Treatment)) +
  geom_errorbar(aes(ymin=richness-se, ymax=richness+se), colour="#696268", width=.1, position=pd) +
  geom_point(color="black", size=5, show.legend=FALSE, position=pd) +
  geom_point(size=4, show.legend=FALSE, position=pd) +
  scale_colour_manual(values=pallete7) +
  labs(title="B", y="Observed OTU richness", x=NULL) +
  theme_classic() + custom_labels

#-----------now taxonomic diversity----------------
#Turn your 'treatment' column into a character vector
sum.f.taxa$Treatment <- as.character(sum.f.taxa$Treatment)
#Then turn it back into a factor with the levels in the correct order
sum.f.taxa$Treatment <- factor(sum.f.taxa$Treatment, levels=c("Compost", "Urea", "Control"))

p5 <- ggplot(sum.f.taxa, aes(x=TimePoint, y=shannon, group=Treatment, shape=Treatment, color=Treatment)) +
  geom_errorbar(aes(ymin=shannon-se, ymax=shannon+se), colour="#696268", width=.1, position=pd) +
  geom_point(color="black", size=5, show.legend=FALSE, position=pd) +
  geom_point(size=4, show.legend=FALSE, position=pd) +
  scale_colour_manual(values=pallete7) +
  labs(title=NULL, y="Shannon Diversity (H')", x="Crop Cycle") +
  theme_classic() + custom_labels

#=============================================
#Arranging Everything into a single Figure
#Making figure legend
pleg <- richness.scat +scale_colour_manual(values=pallete2) + theme_classic() +labs(title="Bacteria", y="Richness", x= NULL) 

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(pleg)


g1 <- ggplotGrob(p1)
g2 <- ggplotGrob(p2)
g3 <- ggplotGrob(p3)
g4 <- ggplotGrob(p4)
g5 <- ggplotGrob(p5)

g2$widths<-unit.pmax(g1$widths)
g3$widths<-unit.pmax(g1$widths)
g4$widths<-unit.pmax(g1$widths)
g5$widths<-unit.pmax(g1$widths)

#grid.arrange(g1,g4,g2,g5,g3,mylegend,ncol=2, nrow=3)old plot with lots of graphs
grid.arrange(g1,g4,g2,g5,ncol=2, nrow=2)#more simplified version

#DIFFERENCES==============================================
#I decided on fitting a linear model using generalized least squares
#this is then followed by a means separation and a pairwise comparison with Tukey correction.
#--------------------------------------------

#Bacteria
#Richness
b.fit.richness <- gls(richness~ TimePoint+Treatment+TimePoint*Treatment, data = data_b, corr = corAR1(, form= ~ 1 |Pot))
summary(b.fit.richness)
anova(b.fit.richness)
plot(b.fit.richness)
#No significance
library(sjPlot)
plot_model(b.fit.richness)
plot_model(b.fit.richness, type="slope")
#matches visually
#--------
#Taxa
b.fit.taxa <- gls(shannon~ TimePoint+Treatment+TimePoint*Treatment, data = data_b, corr = corAR1(, form= ~ 1 |Pot))
summary(b.fit.taxa)
anova(b.fit.taxa)
plot(b.fit.taxa)
#No significance
library(sjPlot)
plot_model(b.fit.taxa)
plot_model(b.fit.taxa, type="slope")
#-----
#Phylo
b.fit.phylo <- gls(faith_pd~ TimePoint+Treatment+timePoint*Treatment, data = data_b, corr = corAR1(, form= ~ 1 |Pot))
summary(b.fit.phylo)
anova(b.fit.phylo)
plot(b.fit.phylo)
#No significance
plot_model(b.fit.phylo)
plot_model(b.fit.phylo, type="slope")
plot_model(b.fit.phylo, type="int")
plot_model(b.fit.phylo, type="pred", terms = c("Treatment", "TimePoint"), colors= c("#eff3ff","#c6dbef","#9ecae1","#6baed6","#4292c6","#2171b5","#084594"))+ theme_dark()

#---------------Fungi
#Richness
f.fit.richness <- gls(richness~ TimePoint+Treatment+TimePoint*Treatment, data = data_f, corr = corAR1(, form= ~ 1 |TimePoint))
summary(f.fit.richness)
anova(f.fit.richness)
plot(f.fit.richness)
plot(f.fit.richness, richness ~ fitted(.) | TimePoint, abline = c(0,1))
#time significance
plot_model(f.fit.richness)
plot_model(f.fit.richness, type="slope")
plot_model(f.fit.richness, type="pred", terms = c("Treatment", "TimePoint"), colors="gs") +theme_classic() +labs(title="Prediction of Fungal Richness by Treatment", y="Number of OTUs", x= "Treatment") 
plot(f.fit.richness$fitted)
marginal <- lsmeans(f.fit.richness, ~ TimePoint:Treatment) 
plot(marginal)
#optional marginal
marginal2<-lsmeans(f.fit.richness, ~TimePoint)
marginal3<-lsmeans(f.fit.richness, ~Treatment) 
#create the ls means comparisons your intersted in 
fit.f.rich.cld<- cld(marginal, alpha   = 0.05, Letters = letters,
         adjust  = "tukey")
fit.f.rich.time.cld<- cld(marginal2, alpha   = 0.05, Letters = letters,
                     adjust  = "tukey")
plot(fit.f.rich.time.cld)
#-------------
#taxa
f.fit.taxa <- gls(shannon~ TimePoint * Treatment, data = data_f, corr = corAR1(, form= ~ 1 |Pot))
summary(f.fit.taxa)
anova(f.fit.taxa)
plot(f.fit.taxa)
#time significance
plot_model(f.fit.taxa)
plot_model(f.fit.taxa, type="slope")
plot_model(f.fit.taxa, type="pred", terms = c("Treatment", "TimePoint"), colors=c("#eff3ff","#c6dbef","#9ecae1","#6baed6","#4292c6","#2171b5","#084594")) +theme_dark()

marginal <- lsmeans(f.fit.taxa, ~ TimePoint:Treatment) 
plot(marginal)
#optional marginal
marginal2<-lsmeans(f.fit.taxa, ~TimePoint)
marginal3<-lsmeans(f.fit.taxa, ~Treatment) 
#create the ls means comparisons your intersted in 
fit.f.taxa.cld<- cld(marginal, alpha   = 0.05, Letters = letters,
                     adjust  = "tukey")
fit.f.taxa.time.cld<- cld(marginal2, alpha   = 0.05, Letters = letters,
                          adjust  = "tukey")
plot(fit.f.taxa.time.cld)

#=================
#Making figure for significan predicted values

#with legend to grab
f.p.fit.rich.leg<-plot_model(f.fit.richness, type="pred", terms = c("Treatment", "TimePoint"), colors="gs") +theme_classic() 

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

new.legend<-g_legend(f.p.fit.rich.leg)
#-------

f.rich.p1 <-plot_model(f.fit.richness, type="pred", terms = c("Treatment", "TimePoint"), colors="gs") +theme_classic() + theme(legend.position = "none") +labs(title="Prediction", y="Number of OTUs", x="Treatment") 

f.taxa.p2<-plot_model(f.fit.taxa, type="pred", terms = c("Treatment", "TimePoint"), colors=c("gs")) +theme_classic()+ theme(legend.position = "none") +labs(title="Prediction", y="Shannon Index", x="Treatment")

#making grobs
g6 <- ggplotGrob(f.rich.p1)
g7 <- ggplotGrob(f.taxa.p2)


#standardizing width, probably not needed for this example
g7$widths<-unit.pmax(g6$widths)  

#plot
grid.arrange(g6,g7,new.legend,ncol=3)
