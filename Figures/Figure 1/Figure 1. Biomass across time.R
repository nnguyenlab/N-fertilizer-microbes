library(readxl)
library(Rmisc)
library(ggplot2)
library(gridExtra) #for grid.arrange

#import excel sheets
biomass <- read_excel("Biomass_Data.xlsx")
abv <- summarySE(biomass, measurevar="ABV_Dry", groupvars=c("Treatment","TimePoint"))
blw <- summarySE(biomass, measurevar="BLW_Dry", groupvars=c("Treatment","TimePoint"))

#Graph each attribute for bacteria
#pd will be a position dodge for the points
pd <- position_dodge(0.65)
#First we are turning the character column into a vector so we can arrange and color the legend as we see fit
#Turn your 'treatment' column into a character vector
abv$Treatment <- as.character(abv$Treatment)
#Then turn it back into a factor with the levels in the correct order
abv$Treatment <- factor(abv$Treatment, levels=c("Compost", "Urea", "Control"))
blw$Treatment <- as.character(blw$Treatment)
#Then turn it back into a factor with the levels in the correct order
blw$Treatment <- factor(blw$Treatment, levels=c("Compost", "Urea", "Control"))

pallete1 <- c("#a62900","#ffcb04","#036900")#kabir
pallete2<- c("#E69F00","#56B4E9","#000000")#Steven's pallete
pallete3 <- c("#de663e","#fecb00","#3e99d2")
pallete4 <- c("#de663e","#fecb00","#4470cf")
pallete5 <- c("#7fc97f","#beaed4","#ffff99")
pallete6 <- c("#035c8c","#00a587","#f2c849")
pallete7 <- c("#f45c85","#f2e127","#0a4d8c")#(pink,yellow,blue)
custom_labels <- theme(axis.text=element_text(size=12),axis.title=element_text(size=14))

#next create the basic ggplot
#starting with richness
p1 <- ggplot(abv, aes(x=TimePoint, y=ABV_Dry, group=Treatment, shape=Treatment, color=Treatment)) +
      geom_errorbar(aes(ymin=ABV_Dry-se, ymax=ABV_Dry+se), colour="#696268", width=.1) +
      geom_point(color="black", size=5, show.legend=FALSE) +
      geom_point(size=4, show.legend=FALSE) +
      scale_colour_manual(values=pallete7) +
      labs(title="A", y="Weight (g)", x= NULL) +
      theme_classic() + custom_labels

#now for belowground
p2 <- ggplot(blw, aes(x=TimePoint, y=BLW_Dry, group=Treatment, shape=Treatment, color=Treatment)) +
      geom_errorbar(aes(ymin=BLW_Dry-se, ymax=BLW_Dry+se), colour="#696268", width=.1) +
      geom_point(color="black", size=5, show.legend=FALSE) +
      geom_point(size=4, show.legend=FALSE) +
      scale_colour_manual(values=pallete7) +
      labs(title="B", y="Weight (g)", x= NULL) +
      theme_classic() + custom_labels

#======Arranging Everything into a single Figure======================
#grid the plots
g1 <- ggplotGrob(p1)
g2 <- ggplotGrob(p2)
grid.arrange(g1,g2,ncol=1,nrow=2)

#======Statistical comparison for Figure======================

#boxcox transformation of data
out1<-aov(ABV_Dry~Pot+Treatment*TimePoint, data=biomass)
bc<-boxcox(out1)
bc$x[which.max(bc$y)]


biomass$trans<- (biomass$ABV_Dry)^-0.989899
bartlett.test(trans~interaction(TimePoint, Treatment), biomass)
leveneTest(trans~interaction(TimePoint, Treatment), biomass)


model.01bc = lme(fixed = trans~ Treatment*TimePoint,  
                 random = ~ 1 | Pot,
                 data = biomass,
                 method = 'REML')
summary(model.01bc)
anova(model.01bc)

shapiro.test(model.01bc$residuals)
plot(model.01bcb)
plot_model(model.01bl, type='diag')



#Tukey for by comparing treatment by time level
posthoc<-glht(model.01bc, linfct = K %*% X)
summary(glht(model.01bc, linfct = K %*% X),test=adjusted("bonferroni"))

#Box cox
out2<-aov(BLW_Dry~Pot+Treatment*TimePoint, data=biomass)
bc2<-boxcox(out2)
bc2$x[which.max(bc2$y)]
biomass$transb<- (biomass$BLW_Dry)^-0.1010101

bartlett.test(BLW_Dry~Treatment, biomass)
leveneTest(BLW_Dry~Treatment, biomass)

bartlett.test(transb~Treatment, biomass)
leveneTest(transb~Treatment, biomass)

model.01bcb = lme(fixed = transb~ Treatment*TimePoint,
                  random = ~ 1 | Pot,
                  data = biomass,
                  method = 'REML')
summary(model.01bcb)
anova(model.01bcb)
VarCorr(model.01bcb)


#Tukey for by comparing treatment by time level
posthoc<-glht(model.01bcb, linfct = K %*% X)
summary(glht(model.01bcb, mcp(Treatment = "Tukey")))
summary(glht(model.01bcb, linfct = K %*% X),test=adjusted("bonferroni"))

biomass$Treatment<-as.factor(biomass$Treatment)
biomass$TimePoint<-as.factor(biomass$TimePoint)
