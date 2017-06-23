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

#########################
# EDITING SCORE FIGURES #
#########################

library(lattice)
library(ggplot2)
library(cowplot)

table=read.delim(file="../DATA/table_analysis_global_lucas.csv",sep=",")

attach(table)

View(table) # Inspecting the Data

levels(table$plastid) <- c("Fucoxanthin","Peridinin")
levels(table$organism) <- c("K. mikimotoi",
                            "K. veneficum",
                            "P. lunula",
                            "S. minutum")
###########
# Figures #
###########



col2=c("#1f78b4","#fb9a99")
col2F=c("#1f78b4","#a6cee3")
col2P=c("#e31a1c","#fb9a99")
col4=c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3")
col5=c("#8dd3c7","#ffffb3","#bebada","#fb8072","#80b1d3")

theme.set = theme_grey() + theme(text = element_text(size = 10), legend.position = "none", plot.title = element_text(face="bold", color = "black", size=12),
                                 legend.title = NULL, legend.text = element_text(size = 7), axis.title.x=element_blank(), plot.subtitle = element_text(size=8))

ylab = "Average editing score diff."

p1 = ggplot(table,aes(x=plastid,y=average.edit.score.diff))+geom_boxplot(aes(fill=plastid))+theme.set
p1 = p1 + ggtitle("Editing scores for each plastid type", subtitle = "Mean : fuco = 0.626 ; pere = 1.601 ; p-val = 8.4e-09 (t.test)\nVariance : fuco/pere = F = 3.137 ; p-val = 6.2e-10 (F.test)") + ylab(ylab) + ylim(-6,7) +
  scale_fill_manual(values=col2) + 
  annotate("text", x = c(1,2), y = rep(7,2), size = 3, fontface = "italic",
           label=c(paste0("n = ",nrow(subset(table,plastid=="Fucoxanthin"))),
                   paste0("n = ",nrow(subset(table,plastid=="Peridinin")))))
p1 = p1 + annotate("segment", x = 1, xend = 2, y = 6, yend = 6) + annotate("text", x = 1.5, y = 6.5,label="***")
var.test(average.edit.score.diff~plastid)
t.test(average.edit.score.diff~plastid)

# Mean -> fuco : 0.6260843 ; pere : 1.6005357 ; p-val = 8.437e-09 (t.test)
# Variance -> fuco/pere = F = 3.1371 ; p-val = 6.181e-10 (F.test)

p1


p2 = ggplot(table,aes(x=organism,y=average.edit.score.diff))+geom_boxplot(aes(fill=organism))+theme.set
p2 = p2 + ggtitle("Editing scores for each organism") + ylab(ylab) + ylim(-6,9) +
  scale_fill_manual(values=c(col2F,col2P)) + 
  annotate("text", x = c(1,2,3,4), y = rep(9,4), size = 3, fontface = "italic",
           label=c(paste0("n = ",nrow(subset(table,organism=="K. mikimotoi"))),
                   paste0("n = ",nrow(subset(table,organism=="K. veneficum"))),
                   paste0("n = ",nrow(subset(table,organism=="P. lunula"))),
                   paste0("n = ",nrow(subset(table,organism=="S. minutum")))))
p2 = p2 + annotate("segment", x = 2, xend = 3, y = 6, yend = 6) + annotate("text", x = 2.5, y = 6.5,label="***") +
  annotate("segment", x = 2, xend = 4, y = 7, yend = 7) + annotate("text", x = 3, y = 7.5,label="***")
p2


p3 = ggplot(table,aes(x=gene_family,y=average.edit.score.diff))+geom_boxplot(aes(fill=gene_family))+theme.set
p3 = p3 + ggtitle("Editing scores for each gene family") + ylab(ylab) + ylim(-6,9) + 
  scale_fill_manual(values=col4) + 
  annotate("text", x = c(1,2,3,4), y = rep(9,4), size = 3, fontface = "italic",
           label=c(paste0("n = ",nrow(subset(table,gene_family=="atp"))),
                   paste0("n = ",nrow(subset(table,gene_family=="pet"))),
                   paste0("n = ",nrow(subset(table,gene_family=="psa"))),
                   paste0("n = ",nrow(subset(table,gene_family=="psb")))))
p3 = p3 + annotate("segment", x = 1, xend = 2, y = 6, yend = 6) + annotate("text", x = 1.5, y = 6.5,label="***") +
  annotate("segment", x = 1, xend = 4, y = 7, yend = 7) + annotate("text", x = 2.5, y = 7.5,label="***")
p3


p4 = ggplot(table,aes(x=reference,y=average.edit.score.diff))+geom_boxplot(aes(fill=reference))+theme.set
p4 = p4 + ggtitle("Editing scores for each reference") + ylab(ylab) + ylim(-6,9) + 
  scale_fill_manual(values=col5) + 
  annotate("text", x = c(1,2,3,4,5), y = rep(9,5), size = 3, fontface = "italic",
           label=c(paste0("n = ",nrow(subset(table,reference=="Pt"))),
                   paste0("n = ",nrow(subset(table,reference=="Ct"))),
                   paste0("n = ",nrow(subset(table,reference=="Eh"))),
                   paste0("n = ",nrow(subset(table,reference=="Ac"))),
                   paste0("n = ",nrow(subset(table,reference=="Vb")))))
p4

pdf(file= "../RESULTS/4_Editing_score_1.pdf", width=14, height=5)
plot_grid(p1,p2,p3,p4, labels = 'AUTO',nrow=1)
dev.off()

####
####
####

theme.set2 = theme_grey() + theme(legend.position = "top", plot.title = element_text(face="bold", color = "black", size=12))
xlab = "Number of Edits"
ylab = "Average edit score differences"

p1 = ggplot(table,aes(num.AA.changes,average.edit.score.diff)) + geom_point(aes(color=factor(plastid)),alpha=0.9) + theme.set2
p1 = p1 + xlab(xlab) + ylab(ylab) + ggtitle("Average score differences depending on\nthe number of edits")
p1 = p1 + scale_color_manual(values=col2,name ="Plastid Type :")
p1

p2 = ggplot(table,aes(num.AA.changes,average.edit.score.diff)) + geom_point(aes(color=factor(organism)),alpha=0.9) + theme.set2
p2 = p2 + xlab(xlab) + ylab(ylab) + ggtitle("Average score differences depending on the number\nof edits and the lineage (***)")
p2 = p2 + facet_grid( . ~ plastid) + 
  scale_color_manual(values=c(col2F,col2P),name ="Organism :")
p2

# Y ~ Gene_family*organism
p4 = ggplot(table,aes(gene_family,average.edit.score.diff)) + geom_boxplot(aes(fill=factor(organism)),alpha=0.8) + theme.set2 + theme(axis.title.x=element_blank())
p4 = p4 + ylab(ylab) + ggtitle("Average edit score differences depending on the gene\nfamily and the lineage (*)")
p4 = p4 + facet_grid( . ~ plastid) + 
  scale_fill_manual(values=c(col2F,col2P),name ="Organism :")
p4

pdf(file= "../RESULTS/4_Editing_score_2.pdf", width=16, height=5)
plot_grid(p1, p2, p4, labels = 'AUTO', nrow=1)
dev.off()

