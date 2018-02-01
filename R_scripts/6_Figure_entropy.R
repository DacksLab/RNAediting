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
table = read_csv(file="../DATA/entropy_analysis.csv")
attach(table)
theme.set2 = theme_grey() + theme(legend.position = "top", plot.title = element_text(face="bold", color = "black", size=12))
col2=c("#1f78b4","#fb9a99")
col2F=c("#1f78b4","#a6cee3")
col2P=c("#e31a1c","#fb9a99")
col4=c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3")
colCV=c("#c994c7","#94c996")
table$Organism <- as.factor(table$Organism)
levels(table$Organism) <- c("K. mikimotoi","K. veneficum","P. lunula","S. minutum")
theme.set = theme_grey() + theme(text = element_text(size = 10), legend.position = "none", plot.title = element_text(face="bold", color = "black", size=12),
                                 legend.title = NULL, legend.text = element_text(size = 7), axis.title.x=element_blank(), plot.subtitle = element_text(size=8))
ylab = ylab("Positional Entropy")
xlab2 = xlab("Editing score")
ylab2 = ylab("Positional Entropy")
###########

###########
# Libraries
###########
library(tidyverse)
library(cowplot)
library(car)
###########

###########
# Figure 1 
###########

#Boxplot of the Entropy score depending on the plastid type
p1 = ggplot(table,aes(x=Lineage,y=positional.entropy))+geom_boxplot(aes(fill=Lineage))+theme.set
p1 = p1 + ggtitle("Positional Entropy for each plastid type") + ylab(ylab) + ylim(NA,1.1) +
  scale_fill_manual(values=col2) + 
  annotate("text", x = c(1,2), y = rep(1.1,2), size = 3, fontface = "italic",
           label=c(paste0("n = ",nrow(subset(table,Lineage=="Fucoxanthin"))),
                   paste0("n = ",nrow(subset(table,Lineage=="Peridinin")))))
wilcox.test(positional.entropy~Lineage)
p1 = p1 + annotate("segment", x = 1, xend = 2, y = 1.03, yend = 1.03) + annotate("text", x = 1.5, y = 1.05,label="***")
p1

#Boxplot of the Entropy score depending on the Organism
p2 = ggplot(table,aes(x=Organism,y=positional.entropy))+geom_boxplot(aes(fill=Organism))+theme.set
p2 = p2 + ggtitle("Positional Entropy for each Organism") + ylab(ylab) + ylim(NA,1.1) +
  scale_fill_manual(values=c(col2F,col2P)) + 
  annotate("text", x = c(1,2,3,4), y = rep(1.1,4), size = 3, fontface = "italic",
           label=c(paste0("n = ",nrow(subset(table,Organism=="K. mikimotoi"))),
                   paste0("n = ",nrow(subset(table,Organism=="K. veneficum"))),
                   paste0("n = ",nrow(subset(table,Organism=="P. lunula"))),
                   paste0("n = ",nrow(subset(table,Organism=="S. minutum")))))
pairwise.wilcox.test(positional.entropy,Organism)
p2 = p2 + 
  annotate("segment", x = 2, xend = 3, y = 1.02, yend = 1.02) + annotate("text", x = 2.5, y = 1.03,label="***")+
  annotate("segment", x = 3, xend = 4, y = 1.04, yend = 1.04) + annotate("text", x = 3.5, y = 1.05,label="***")
p2

#Boxplot of the Entropy score depending on whether there is edition
p3 = ggplot(table,aes(x=edited.or.not,y=positional.entropy))+geom_boxplot(aes(fill=edited.or.not))+theme.set
p3 = p3 + ggtitle("Positional Entropy depending on the\nedition status") + ylab(ylab) + ylim(NA,1.1) + 
  scale_fill_manual(values=colCV) + scale_x_discrete(labels=c("Non Edited", "Edited")) +
  annotate("text", x = c(1,2), y = rep(1.1,2), size = 3, fontface = "italic",
           label=c(paste0("n = ",nrow(subset(table,edited.or.not=="no"))),
                   paste0("n = ",nrow(subset(table,edited.or.not=="yes")))))
wilcox.test(positional.entropy~edited.or.not)
p3 = p3 + annotate("segment", x = 1, xend = 2, y = 1.05, yend = 1.05) + annotate("text", x = 1.5, y = 1.07,label="***")
p3

plot_save1 = plot_grid(p1,p2,p3,labels='AUTO',nrow=1)
ggsave(filename= "../RESULTS/6_Entropy_1.pdf",plot_save1, width=12, height=5)
###########

###########
# Figure 2
###########

# Y ~ score*Lineage
p1 = ggplot(table,aes(x=average.score.difference,y=positional.entropy)) + geom_point(aes(color=Lineage),alpha = 0.7) + theme.set2
p1 = p1 + xlab2 + ylab2 + ggtitle("Positional Entropy depending on the Editing score and the lineage (**)")
summary(lm(positional.entropy~average.score.difference*Lineage))
p1 = p1 + scale_color_manual(values=col2, name = "Plastid : ") + geom_density(aes(x=average.score.difference,y=..scaled..,color=Lineage))
p1

# Y ~ score*Organism
p2 = ggplot(table,aes(x=average.score.difference,y=positional.entropy)) + geom_point(aes(color=Organism),alpha = 0.7) + theme.set2
p2 = p2 + xlab2 + ylab2 + ggtitle("Positional Entropy depending on the Editing score and the organism (***)")
summary(lm(positional.entropy~average.score.difference*Organism))
p2 = p2 + facet_grid( . ~ Lineage) + 
  scale_color_manual(values=c(col2F,col2P),name ="Organism : ")  + geom_density(aes(x=average.score.difference, y=..scaled.., color=Organism))
p2

plot_save2 = plot_grid(p1, p2, labels = 'AUTO',nrow=1)
ggsave(filename = "../RESULTS/6_Entropy_2.pdf", plot_save2, width=13, height=6)
###########

###########
# Figure Score
###########
table.score = subset(table, edited.or.not == 'yes')
table.score$conservation.status=ifelse(table.score$positional.entropy>=0.95,'Conserved', 'Not Conserved')
table.score$conservation.lineage = paste0(table.score$conservation.status,'\n',table.score$Lineage)

p.score = ggplot(table.score) + geom_boxplot(aes(x=conservation.lineage,y=average.score.difference, fill=Lineage)) + 
  theme.set + scale_fill_manual(values=col2) + ggtitle("Editing score depending on the conservation and the plastid type") +
  annotate("text", x = c(1,2,3,4), y = rep(17,4), size = 3, fontface = "italic",
           label=c(paste0("n = ",nrow(subset(table.score,conservation.lineage=="Conserved\nFucoxanthin"))),
                   paste0("n = ",nrow(subset(table.score,conservation.lineage=="Conserved\nPeridinin"))),
                   paste0("n = ",nrow(subset(table.score,conservation.lineage=="Not Conserved\nFucoxanthin"))),
                   paste0("n = ",nrow(subset(table.score,conservation.lineage=="Not Conserved\nPeridinin")))))
pairwise.t.test(table.score$average.score.difference, table.score$conservation.lineage)
p.score = p.score +
  annotate("segment", x = 1, xend = 2, y = 12, yend = 12) + annotate("text", x = 1.5, y = 12.3,label="***") +
  annotate("segment", x = 2, xend = 3, y = 13.2, yend = 13.2) + annotate("text", x = 2.5, y = 13.5,label="***") +
  annotate("segment", x = 2, xend = 4, y = 14.4, yend = 14.4) + annotate("text", x = 3, y = 14.7,label="***") +
  annotate("segment", x = 1, xend = 4, y = 15.6, yend = 15.6) + annotate("text", x = 2.5, y = 15.9,label="*")
p.score

ggsave(filename= "../RESULTS/6_Editing_Score_Entropy.pdf",p.score, width=7, height=7)
###########