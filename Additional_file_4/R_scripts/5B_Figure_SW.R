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
# Loading data #
################

SW=read.delim(file="../DATA/SW_analysis_global_lucas.csv",sep=",")

# Filtering out p-value < 0.05
SW=subset(SW,dt(SW$amino_acid_t_value,df=SW$DF_.N.2.)<=0.05)

# Loading the edit score data : goal is to investigate relationship between SW and edit score
table=read.delim(file="../DATA/table_analysis_global_lucas.csv",sep=",")

#Gene as character to avoid potential issues
SW$gene=as.character(SW$gene)
table$gene=as.character(table$gene)

# Adding the pearson correlation to the table with edit scores
for (i in 1:nrow(table)){
  if (length(subset(SW, gene == table[i,"gene"] &
                    organism == table[i,"organism"] &
                    reference == table[i,"reference"])$amino_acid_pearson_correlation)==0){
    table[i,8]=NA
  }else{
    table[i,8]=subset(SW, gene == table[i,"gene"] &
                        organism == table[i,"organism"] &
                        reference == table[i,"reference"])$amino_acid_pearson_correlation
  }
}

table=na.omit(table)
names(table)=c(names(table)[-8],"SW.score")

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

ylab = "Sliding Window score"


p1 = ggplot(table,aes(x=plastid,y=SW.score))+geom_boxplot(aes(fill=plastid))+theme.set
p1 = p1 + ggtitle("Sliding window correlation for each\nplastid type") + ylab(ylab) + ylim(-1,1) +
  scale_fill_manual(values=col2) + 
  annotate("text", x = c(1,2), y = rep(0.9,2), size = 3, fontface = "italic",
           label=c(paste0("n = ",nrow(subset(table,plastid=="Fucoxanthin"))),
                   paste0("n = ",nrow(subset(table,plastid=="Peridinin")))))
p1 = p1 + annotate("segment", x = 1, xend = 2, y = 0.75, yend = 0.75) + annotate("text", x = 1.5, y = 0.8,label="**")
p1


p2 = ggplot(table,aes(x=organism,y=SW.score))+geom_boxplot(aes(fill=organism))+theme.set
p2 = p2 + ggtitle("Sliding window correlation for each\norganism") + ylab(ylab) + ylim(-1,1.5) +
  scale_fill_manual(values=c(col2F,col2P)) + 
  annotate("text", x = c(1,2,3,4), y = rep(1.4,4), size = 3, fontface = "italic",
           label=c(paste0("n = ",nrow(subset(table,organism=="K. mikimotoi"))),
                   paste0("n = ",nrow(subset(table,organism=="K. veneficum"))),
                   paste0("n = ",nrow(subset(table,organism=="P. lunula"))),
                   paste0("n = ",nrow(subset(table,organism=="S. minutum")))))
p2 = p2 + annotate("segment", x = 1, xend = 2, y = 0.75, yend = 0.75) + annotate("text", x = 1.5, y = 0.8,label="***") +
  annotate("segment", x = 1, xend = 3, y = 1, yend = 1) + annotate("text", x = 2, y = 1.05,label="***") +
  annotate("segment", x = 1, xend = 4, y = 1.25, yend = 1.25) + annotate("text", x = 2.5, y = 1.3,label="*")
p2


p3 = ggplot(table,aes(x=gene_family,y=SW.score))+geom_boxplot(aes(fill=gene_family))+theme.set
p3 = p3 + ggtitle("Sliding window correlation for each\ngene family") + ylab(ylab) + ylim(NA,0.8) + 
  scale_fill_manual(values=col4) + 
  annotate("text", x = c(1,2,3,4), y = rep(0.75,4), size = 3, fontface = "italic",
           label=c(paste0("n = ",nrow(subset(table,gene_family=="atp"))),
                   paste0("n = ",nrow(subset(table,gene_family=="pet"))),
                   paste0("n = ",nrow(subset(table,gene_family=="psa"))),
                   paste0("n = ",nrow(subset(table,gene_family=="psb")))))
p3


p4 = ggplot(table,aes(x=reference,y=SW.score))+geom_boxplot(aes(fill=reference))+theme.set
p4 = p4 + ggtitle("Sliding window correlation for each\nrefence") + ylab(ylab) + ylim(NA,0.8) + 
  scale_fill_manual(values=col5) + 
  annotate("text", x = c(1,2,3,4,5), y = rep(0.75,5), size = 3, fontface = "italic",
           label=c(paste0("n = ",nrow(subset(table,reference=="Pt"))),
                   paste0("n = ",nrow(subset(table,reference=="Ct"))),
                   paste0("n = ",nrow(subset(table,reference=="Eh"))),
                   paste0("n = ",nrow(subset(table,reference=="Ac"))),
                   paste0("n = ",nrow(subset(table,reference=="Vb")))))
p4

# Note that under a model, we can identify significance as follows.
# But the model is quite farfetched and biological relevance is far from established.
#annotate("segment", x = 2, xend = 4, y = 0.75, yend = 0.75) + annotate("text", x = 3, y = 0.75,label="**")
#annotate("segment", x = 3, xend = 4, y = 0.75, yend = 0.75) + annotate("text", x = 3.5, y = 0.75,label="**")
#annotate("segment", x = 1, xend = 4, y = 0.75, yend = 0.75) + annotate("text", x = 2.5, y = 0.75,label="*")

pdf(file= "../RESULTS/5_Sliding_Window_1.pdf", width=16, height=5)
plot_grid(p1,p2,p3,p4, labels = 'AUTO',nrow=1)
dev.off()

###
###
###

theme.set2 = theme_grey() + theme(legend.position = "top", plot.title = element_text(face="bold", color = "black", size=12))
xlab = "Number of Edits"
ylab = "Sliding Window scores"

# Y ~ Number of edits
p1=ggplot(table,aes(num.AA.changes,SW.score)) + geom_point(aes(color=factor(plastid))) + theme.set2 
p1 = p1 + xlab(xlab) + ylab(ylab) + ggtitle("Sliding Window scores depending on the number of edits")
p1 = p1 + scale_color_manual(values=col2,name ="Plastid Type :")
p1

# Y ~ Number of edits*organism
p2 = ggplot(table,aes(num.AA.changes,SW.score)) + geom_point(aes(color=factor(organism))) + theme.set2
p2 = p2 + xlab(xlab) + ylab(ylab) + ggtitle("Sliding Window scores depending on the number of edits\nand the lineage")
p2 = p2 + facet_grid( . ~ plastid) +
  scale_color_manual(values=c(col2F,col2P),name ="Organism :")
p2

# Y ~ Average Edit Score*organism
p3 = ggplot(table,aes(average.edit.score.diff,SW.score)) + geom_point(aes(color=factor(organism))) + theme.set2 
p3 = p3 + xlab("Average edit score differences") + ylab(ylab) + ggtitle("Sliding Window scores depending on the\nAverage edit score differences and the lineage (***)")
p3 = p3 + facet_grid( . ~ plastid) + 
  scale_color_manual(values=c(col2F,col2P),name ="Organism :")
p3

# Y ~ Gene_family*organism
p4 = ggplot(table,aes(gene_family,SW.score)) + geom_boxplot(aes(fill=factor(organism)),alpha=0.8) + theme.set2 + theme(axis.title.x=element_blank())
p4 = p4 + ylab(ylab) + ggtitle("Sliding Window scores depending on the gene family\nand the lineage (***)")
p4 = p4 + facet_grid( . ~ plastid) + 
  scale_fill_manual(values=c(col2F,col2P),name ="Organism :")
p4

pdf(file= "../RESULTS/5_Sliding_Window_2.pdf", width=12, height=10)
plot_grid(p1, p2, p3, p4, labels = 'AUTO')
dev.off()

