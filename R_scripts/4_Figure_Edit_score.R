##################################################
# Plastid transcript editing across dinoflagellate
# lineages shows lineage-specific application 
# but conserved trends
##################################################

#  Wrote by Lucas Paoli. Contact : lucas.paoli@ens.fr | lucas.paoli@gmail.com

###########
rm(list=ls())
# Directory
setwd("/Users/Lucas/Documents/Publications/2017_Dino/ANALYSIS")
# Variables
table=read_csv(file="../DATA/AA_analysis.csv")
temp = table %>%
  group_by(Gene, Organism, Reference) %>%
  summarize(Num.AA.changes = sum(genomic.AA != transcript.AA),
            mean.edit.score = mean(edit.score.difference),
            sd.edit.score = sd(edit.score.difference))
table = left_join(table,temp, by = c('Gene', 'Organism', 'Reference'))
table$Reference = as.factor(table$Reference)
table = subset(table, Lineage %in% c('Fucoxanthin', 'Peridinin') & Gene.type == 'Photosynthetic')
attach(table)
View(table) # Inspecting the Data
col2=c("#1f78b4","#fb9a99")
col2F=c("#1f78b4","#a6cee3")
col2P=c("#e31a1c","#fb9a99")
col4=c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3")
col5=c("#8dd3c7","#ffffb3","#bebada","#fb8072","#80b1d3")
table$Organism = as.factor(table$Organism)
levels(table$Organism) <- c("K. mikimotoi", "K. veneficum", "P. lunula", "S. minutum")
ylab = "Editing score diff."
ylim1 = ylim(c(-7.5,NA))
theme.set = theme_grey() + 
  theme(text = element_text(size = 10), 
        legend.position = "none", 
        plot.title = element_text(face="bold", color = "black", size=12),
        legend.title = NULL, legend.text = element_text(size = 7), 
        axis.title.x=element_blank(), 
        plot.subtitle = element_text(size=8))
theme.set2 = theme_grey() + theme(legend.position = "top", plot.title = element_text(face="bold", color = "black", size=12))
xlab2 = xlab("Number of Amino Acid changes")
ylab2 = ylab("Per gene mean edit score differences")
###########

###########
# Libraries
###########
library(tidyverse)
library(cowplot)
library(car)
library(lsmeans)
library(nlme)
library(lme4)
library(multcomp)
library(lmerTest)
###########

###########
# Figure 1
###########

# We first investigate the effect of reference on the data.
hist(edit.score.difference)
lm.test = lm(edit.score.difference ~ Reference*(Gene.family+Lineage), data=table)
plot(lm.test)
summary(lm.test)
Anova(lm.test)
TukeyHSD(aov(edit.score.difference ~ Reference, data=table))
# No reference effect!

test.var = var.test(edit.score.difference~Lineage)
test.mean = t.test(edit.score.difference~Lineage)
summary(glht(lmer(edit.score.difference ~ Lineage + (1|Reference), data=table),linfct=mcp(Lineage="Tukey")))
pval = 2e-16

p1 = ggplot(table,aes(x=Lineage,y=edit.score.difference))+geom_boxplot(aes(fill=Lineage),outlier.shape=NA)+theme.set
p1 = p1 + ggtitle("Editing scores for each plastid type", 
                  subtitle = paste0("Mean : fuco = ",format(test.mean$estimate[[1]],digits=4),
                                    " ; pere = ",format(test.mean$estimate[[2]],digits=4),
                                    " ; p-val = ",format(pval,digits = 2),
                                    " (Tukey)\nVariance : fuco/pere = F = ",format(test.var$estimate[[1]],digits=4)," ; p-val = ",format(test.var$p.value,digits = 2)," (F.test)")) +
  ylab(ylab) + ylim1 +
  scale_fill_manual(values=col2) + 
  annotate("text", x = c(1,2), y = rep(10,2), size = 3, fontface = "italic",
           label=c(paste0("n = ",nrow(subset(table,Lineage=="Fucoxanthin"))),
                   paste0("n = ",nrow(subset(table,Lineage=="Peridinin")))))
p1 = p1 + annotate("segment", x = 1, xend = 2, y = 9, yend = 9) + annotate("text", x = 1.5, y = 9.5,label="***")

p1


p2 = ggplot(table,aes(x=Organism,y=edit.score.difference))+geom_boxplot(aes(fill=Organism),outlier.shape = NA)+theme.set
p2 = p2 + ggtitle("Editing scores for each Organism") + ylab(ylab) + ylim1 +
  scale_fill_manual(values=c(col2F,col2P)) + 
  annotate("text", x = c(1,2,3,4), y = rep(18,4), size = 3, fontface = "italic",
           label=c(paste0("n = ",nrow(subset(table,Organism=="K. mikimotoi"))),
                   paste0("n = ",nrow(subset(table,Organism=="K. veneficum"))),
                   paste0("n = ",nrow(subset(table,Organism=="P. lunula"))),
                   paste0("n = ",nrow(subset(table,Organism=="S. minutum")))))
summary(glht(lmer(edit.score.difference ~ Organism + (1|Reference), data=table),linfct=mcp(Organism="Tukey")))
p2 = p2 + 
  annotate("segment", x = 1, xend = 2, y = 10.5, yend = 10.5) + annotate("text", x = 1.5, y = 11,label="***") +
  annotate("segment", x = 1, xend = 3, y = 12, yend = 12) + annotate("text", x = 2, y = 12.5,label="***") +
  annotate("segment", x = 1, xend = 4, y = 13.5, yend = 13.5) + annotate("text", x = 2.5, y = 14,label="***") +
  annotate("segment", x = 2, xend = 3, y = 15, yend = 15) + annotate("text", x = 2.5, y = 15.5,label="***") +
  annotate("segment", x = 2, xend = 4, y = 16.5, yend = 16.5) + annotate("text", x = 3, y = 17,label="***")
p2


p3 = ggplot(table,aes(x=Gene.family,y=edit.score.difference))+geom_boxplot(aes(fill=Gene.family),outlier.shape=NA)+theme.set
p3 = p3 + ggtitle("Editing scores for each gene family") + ylab(ylab) + ylim1 + 
  scale_fill_manual(values=col4) + 
  annotate("text", x = c(1,2,3,4), y = rep(13.5,4), size = 3, fontface = "italic",
           label=c(paste0("n = ",nrow(subset(table,Gene.family=="atp"))),
                   paste0("n = ",nrow(subset(table,Gene.family=="pet"))),
                   paste0("n = ",nrow(subset(table,Gene.family=="psa"))),
                   paste0("n = ",nrow(subset(table,Gene.family=="psb")))))
summary(glht(lmer(edit.score.difference ~ Gene.family + (1|Reference), data=table),linfct=mcp(Gene.family="Tukey")))
p3 = p3 + 
  annotate("segment", x = 1, xend = 2, y = 9, yend = 9) + annotate("text", x = 1.5, y = 9.5,label="***") +
  annotate("segment", x = 1, xend = 3, y = 10.5, yend = 10.5) + annotate("text", x = 2, y = 11,label="***") +
  annotate("segment", x = 2, xend = 4, y = 12, yend = 12) + annotate("text", x = 3, y = 12.5,label="***") +
p3


p4 = ggplot(table,aes(x=Reference,y=edit.score.difference))+geom_boxplot(aes(fill=Reference),outlier.shape = NA)+theme.set
p4 = p4 + ggtitle("Editing scores for each reference") + ylab(ylab) + ylim1 + 
  scale_fill_manual(values=col5) + 
  annotate("text", x = c(1,2,3,4,5), y = rep(10,5), size = 3, fontface = "italic",
           label=c(paste0("n = ",nrow(subset(table,Reference=="Pt"))),
                   paste0("n = ",nrow(subset(table,Reference=="Ct"))),
                   paste0("n = ",nrow(subset(table,Reference=="Eh"))),
                   paste0("n = ",nrow(subset(table,Reference=="Ac"))),
                   paste0("n = ",nrow(subset(table,Reference=="Vb")))))
TukeyHSD(aov(edit.score.difference ~ Reference))
p4

plot_save1 = plot_grid(p1,p2,p3,p4, labels = 'AUTO',nrow=1)
ggsave(filename = "../RESULTS/4_Editing_score_1.pdf",plot_save1, width=15, height=5)
###########

###########
# Figure 2
###########

cor.test(mean.edit.score,Num.AA.changes)

# Y ~ Num.AA.changes*Lineage
p1 = ggplot(table,aes(Num.AA.changes,mean.edit.score)) + geom_point(aes(color=factor(Lineage)),alpha=0.9) + theme.set2
summary(lmer(mean.edit.score ~ Lineage*Num.AA.changes + (1|Reference), data=table))
Anova(lmer(mean.edit.score ~ Lineage*Num.AA.changes + (1|Reference), data=table))
p1 = p1 + xlab2 + ylab2 + ggtitle("Per gene mean edit score differences depending\non the number of Amino Acid changes (***)")
p1 = p1 + scale_color_manual(values=col2,name ="Plastid Type :")
p1

# Y ~ Num.AA.changes*Organism
p2 = ggplot(table,aes(Num.AA.changes,mean.edit.score)) + geom_point(aes(color=factor(Organism)),alpha=0.9) + theme.set2
Anova(lmer(mean.edit.score ~ Organism*Num.AA.changes + (1|Reference), data=table))
p2 = p2 + xlab2 + ylab2 + ggtitle("Per gene mean edit score differences depending\non the number of edits and the lineage (***)")
p2 = p2 + facet_grid( . ~ Lineage) + 
  scale_color_manual(values=c(col2F,col2P),name ="Organism :")
p2

# Y ~ Gene.family*Organism
p4 = ggplot(table,aes(Gene.family,mean.edit.score)) + geom_boxplot(aes(fill=factor(Organism)),alpha=0.8) + theme.set2 + theme(axis.title.x=element_blank())
Anova(lmer(mean.edit.score ~ Organism*Gene.family + (1|Reference), data=table))
p4 = p4 + ylab(ylab) + ggtitle("Per gene mean edit score differences depending\non the gene family and the organism (***)")
p4 = p4 + facet_grid( . ~ Lineage) + 
  scale_fill_manual(values=c(col2F,col2P),name ="Organism :")
p4

plot_save2 = plot_grid(p1, p2, p4, labels = 'AUTO', nrow=1)
ggsave(filename = "../RESULTS/4_Editing_score_2.pdf",plot_save2, width=16, height=5)
###########
