##################################################
# Statistical Analysis of Dinos ARN editing Data # @Name of the article
##################################################

#  Created by Lucas Paoli. Contact : lucas.paoli@ens.fr | lucas.paoli@gmail.com
#  On OS X 10.9.5 - x86_64, darwin13.4.0
#  R version 3.3.1 (2016-06-21)

# Notation :
# We focus on 4 organisms 
# kv : K. veneficum / Fucoxanthin
# km : K. mikimotoi / Fucoxanthin
# sm : S. minutum / Peridinin
# pl : P. lunula / Peridinin

##########################
# SLIDING WINDOW FIGURES #
##########################

library(lattice)
library(ggplot2)
library(cowplot)

################
# LOADING DATA #
################

table = read.table(file="../DATA/entropy_global_table.csv",sep=",",header=T)

attach(table)


###########
# Figures #
###########

col2=c("#1f78b4","#fb9a99")
col2F=c("#1f78b4","#a6cee3")
col2P=c("#e31a1c","#fb9a99")
col4=c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3")
colCV=c("#c994c7","#94c996")

levels(table$plastid) <- c("Fucoxanthin","Peridinin")
levels(table$organism) <- c("K. mikimotoi",
                            "K. veneficum",
                            "P. lunula",
                            "S. minutum")

theme.set = theme_grey() + theme(text = element_text(size = 10), legend.position = "none", plot.title = element_text(face="bold", color = "black", size=12),
                                 legend.title = NULL, legend.text = element_text(size = 7), axis.title.x=element_blank(), plot.subtitle = element_text(size=8))
ylab = "Positional Entropy"

#Boxplot of the Entropy score depending on the plastid type
p1 = ggplot(table,aes(x=plastid,y=positional.entropy))+geom_boxplot(aes(fill=plastid))+theme.set
p1 = p1 + ggtitle("Positional Entropy for each plastid type") + ylab(ylab) + ylim(NA,1.1) +
  scale_fill_manual(values=col2) + 
  annotate("text", x = c(1,2), y = rep(1.1,2), size = 3, fontface = "italic",
           label=c(paste0("n = ",nrow(subset(table,plastid=="Fucoxanthin"))),
                   paste0("n = ",nrow(subset(table,plastid=="Peridinin")))))
p1

#Boxplot of the Entropy score depending on the organism
p2 = ggplot(table,aes(x=organism,y=positional.entropy))+geom_boxplot(aes(fill=organism))+theme.set
p2 = p2 + ggtitle("Positional Entropy for each organism") + ylab(ylab) + ylim(NA,1.1) +
  scale_fill_manual(values=c(col2F,col2P)) + 
  annotate("text", x = c(1,2,3,4), y = rep(1.1,4), size = 3, fontface = "italic",
           label=c(paste0("n = ",nrow(subset(table,organism=="K. mikimotoi"))),
                   paste0("n = ",nrow(subset(table,organism=="K. veneficum"))),
                   paste0("n = ",nrow(subset(table,organism=="P. lunula"))),
                   paste0("n = ",nrow(subset(table,organism=="S. minutum")))))
p2

#Boxplot of the Entropy score depending on whether there is edition
p3 = ggplot(table,aes(x=edited.or.not,y=positional.entropy))+geom_boxplot(aes(fill=edited.or.not))+theme.set
p3 = p3 + ggtitle("Positional Entropy depending on the\nedition status") + ylab(ylab) + ylim(NA,1.1) + 
  scale_fill_manual(values=colCV) + scale_x_discrete(labels=c("Non Edited", "Edited")) +
  annotate("text", x = c(1,2), y = rep(1.1,2), size = 3, fontface = "italic",
           label=c(paste0("n = ",nrow(subset(table,edited.or.not=="no"))),
                   paste0("n = ",nrow(subset(table,edited.or.not=="yes")))))
p3 = p3 + annotate("segment", x = 1, xend = 2, y = 1.05, yend = 1.05) + annotate("text", x = 1.5, y = 1.07,label="***")
p3

pdf(file= "../RESULTS/6_Entropy_1.pdf", width=12, height=5)
plot_grid(p1,p2,p3,labels='AUTO',nrow=1)
dev.off()

###
###
###

theme.set2 = theme_grey() + theme(legend.position = "top", plot.title = element_text(face="bold", color = "black", size=12))
xlab = "Editing score"
ylab = "Positional Entropy"

# Y ~ score
p1 = ggplot(table,aes(x=score,y=positional.entropy)) + geom_point(aes(color=plastid),alpha = 0.7) + theme.set2
p1 = p1 + xlab(xlab) + ylab(ylab) + ggtitle("Positional Entropy depending on the Editing score")
p1 = p1 + scale_color_manual(values=col2, name = "Plastid : ") + geom_density(aes(x=score,y=..scaled..,color=plastid))
p1

# Y ~ score*organism
p3 = ggplot(table,aes(x=score,y=positional.entropy)) + geom_point(aes(color=organism),alpha = 0.7) + theme.set2
p3 = p3 + xlab(xlab) + ylab(ylab) + ggtitle("Positional Entropy depending on the Editing score and the lineage")
p3 = p3 + facet_grid( . ~ plastid) + 
  scale_color_manual(values=c(col2F,col2P),name ="Organism : ")  + geom_density(aes(x=score, y=..scaled.., color=organism))
p3

pdf(file= "../RESULTS/6_Entropy_2.pdf", width=12, height=6)
plot_grid(p1, p3, labels = 'AUTO',nrow=1)
dev.off()
