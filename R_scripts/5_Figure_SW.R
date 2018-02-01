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
SW=read_csv(file="../DATA/SW_analysis.csv")
table=read_csv(file="../DATA/AA_analysis.csv")
temp = table %>%
  group_by(Gene, Organism, Reference) %>%
  summarize(Num.AA.changes = sum(genomic.AA != transcript.AA),
            mean.edit.score = mean(edit.score.difference),
            sd.edit.score = sd(edit.score.difference))
SW = left_join(SW,temp, by = c('Gene', 'Organism', 'Reference'))
SW = subset(SW, Lineage %in% c('Fucoxanthin', 'Peridinin') & Gene.type == 'Photosynthetic')
# Filtering out p-value < 0.05
sum(dt(SW$nucleotide.t.value,df=SW$DF..N.2.)>0.05) ; nrow(SW)
sum(dt(SW$amino.acid.t.value,df=SW$DF..N.2.)>0.05) ; nrow(SW)
SW=subset(SW,dt(SW$amino.acid.t.value,df=SW$DF..N.2.)<=0.05)
SW=na.omit(SW)
SW$Organism <- as.factor(SW$Organism)
levels(SW$Organism) <- c("K. mikimotoi","K. veneficum","P. lunula","S. minutum")
attach(SW)
col2=c("#1f78b4","#fb9a99")
col2F=c("#1f78b4","#a6cee3")
col2P=c("#e31a1c","#fb9a99")
col4=c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3")
col5=c("#8dd3c7","#ffffb3","#bebada","#fb8072","#80b1d3")
theme.set = theme_grey() + theme(text = element_text(size = 10), legend.position = "none", plot.title = element_text(face="bold", color = "black", size=12),
                                 legend.title = NULL, legend.text = element_text(size = 7), axis.title.x=element_blank(), plot.subtitle = element_text(size=8))
ylab = ylab("Sliding Window score")
theme.set2 = theme_grey() + theme(legend.position = "top", plot.title = element_text(face="bold", color = "black", size=12))
xlab2 = xlab("Number of Edits")
ylab2 = ylab("Sliding Window scores")
###########

###########
# Figure 1 
###########

# We first investigate the effect of reference on the data.
hist(amino.acid.pearson.correlation)
lm.test = lm(amino.acid.pearson.correlation ~ Reference*(Gene.family+Lineage), data=SW)
plot(lm.test)
summary(lm.test)
Anova(lm.test)
TukeyHSD(aov(edit.score.difference ~ Reference, data=table))
# No reference effect!

p1 = ggplot(SW,aes(x=Lineage,y=amino.acid.pearson.correlation))+geom_boxplot(aes(fill=Lineage))+theme.set
p1 = p1 + ggtitle("Sliding window correlation for each\nplastid type") + ylab(ylab) +
  scale_fill_manual(values=col2) + 
  annotate("text", x = c(1,2), y = rep(0.9,2), size = 3, fontface = "italic",
           label=c(paste0("n = ",nrow(subset(SW,Lineage=="Fucoxanthin"))),
                   paste0("n = ",nrow(subset(SW,Lineage=="Peridinin")))))
summary(glht(lmer(amino.acid.pearson.correlation ~ Lineage + (1|Reference), data=SW),linfct=mcp(Lineage="Tukey")))
p1

p2 = ggplot(SW,aes(x=Organism,y=amino.acid.pearson.correlation))+geom_boxplot(aes(fill=Organism))+theme.set
p2 = p2 + ggtitle("Sliding window correlation for each\nOrganism") + ylab(ylab) +
  scale_fill_manual(values=c(col2F,col2P)) + 
  annotate("text", x = c(1,2,3,4), y = rep(1.5,4), size = 3, fontface = "italic",
           label=c(paste0("n = ",nrow(subset(SW,Organism=="K. mikimotoi"))),
                   paste0("n = ",nrow(subset(SW,Organism=="K. veneficum"))),
                   paste0("n = ",nrow(subset(SW,Organism=="P. lunula"))),
                   paste0("n = ",nrow(subset(SW,Organism=="S. minutum")))))
summary(glht(lmer(amino.acid.pearson.correlation ~ Organism + (1|Reference), data=SW),linfct=mcp(Organism="Tukey")))
p2 = p2 + 
  annotate("segment", x = 3, xend = 4, y = 1, yend = 1) + annotate("text", x = 3.5, y = 1.05,label="**") +
  annotate("segment", x = 1, xend = 3, y = 0.75, yend = 0.75) + annotate("text", x = 2, y = 0.8,label="***") +
  annotate("segment", x = 2, xend = 4, y = 1.25, yend = 1.25) + annotate("text", x = 3, y = 1.3,label="***")
p2

p3 = ggplot(SW,aes(x=Gene.family,y=amino.acid.pearson.correlation))+geom_boxplot(aes(fill=Gene.family))+theme.set
p3 = p3 + ggtitle("Sliding window correlation for each\ngene family") + ylab(ylab) + ylim(NA,0.8) + 
  scale_fill_manual(values=col4)+
  annotate("text", x = c(1,2,3,4), y = rep(0.75,4), size = 3, fontface = "italic",
           label=c(paste0("n = ",nrow(subset(SW,Gene.family=="atp"))),
                   paste0("n = ",nrow(subset(SW,Gene.family=="pet"))),
                   paste0("n = ",nrow(subset(SW,Gene.family=="psa"))),
                   paste0("n = ",nrow(subset(SW,Gene.family=="psb")))))
summary(glht(lmer(amino.acid.pearson.correlation ~ Gene.family + (1|Reference), data=SW),linfct=mcp(Gene.family="Tukey")))
p3

p4 = ggplot(SW,aes(x=Reference,y=amino.acid.pearson.correlation))+geom_boxplot(aes(fill=Reference))+theme.set
p4 = p4 + ggtitle("Sliding window correlation for each\nrefence") + ylab(ylab) + ylim(NA,0.8) + 
  scale_fill_manual(values=col5) + 
  annotate("text", x = c(1,2,3,4,5), y = rep(0.75,5), size = 3, fontface = "italic",
           label=c(paste0("n = ",nrow(subset(SW,Reference=="Pt"))),
                   paste0("n = ",nrow(subset(SW,Reference=="Ct"))),
                   paste0("n = ",nrow(subset(SW,Reference=="Eh"))),
                   paste0("n = ",nrow(subset(SW,Reference=="Ac"))),
                   paste0("n = ",nrow(subset(SW,Reference=="Vb")))))
p4

plot_save1 = plot_grid(p1,p2,p3,p4, labels = 'AUTO',nrow=1)
ggsave(filename= "../RESULTS/5_Sliding_Window_1.pdf",plot_save1, width=16, height=5)
###########

###########
# Figure 2 
###########

# Y ~ Number of edits*Lineage
p1=ggplot(SW,aes(Num.AA.changes,amino.acid.pearson.correlation)) + geom_point(aes(color=factor(Lineage))) + theme.set2 
Anova(lmer(amino.acid.pearson.correlation ~ Lineage*Num.AA.changes + (1|Reference), data=SW))
p1 = p1 + xlab2 + ylab2 + ggtitle("Sliding Window scores depending on the number of edits\nand the lineage (*)")
p1 = p1 + scale_color_manual(values=col2,name ="Plastid Type :")
p1

# Y ~ Number of edits*Organism
p2 = ggplot(SW,aes(Num.AA.changes,amino.acid.pearson.correlation)) + geom_point(aes(color=factor(Organism))) + theme.set2
Anova(lmer(amino.acid.pearson.correlation ~ Organism*Num.AA.changes + (1|Reference), data=SW))
p2 = p2 + xlab2 + ylab2 + ggtitle("Sliding Window scores depending on the number of edits\nand the organism")
p2 = p2 + facet_grid( . ~ Lineage) +
  scale_color_manual(values=c(col2F,col2P),name ="Organism :")
p2

# Y ~ Average Edit Score*Organism
p3 = ggplot(SW,aes(mean.edit.score,amino.acid.pearson.correlation)) + geom_point(aes(color=factor(Organism))) + theme.set2 
Anova(lmer(amino.acid.pearson.correlation ~ Organism+mean.edit.score + (1|Reference), data=SW))
p3 = p3 + xlab("Average edit score differences") + ylab2 + ggtitle("Sliding Window scores depending on the\nAverage edit score differences and the organism")
p3 = p3 + facet_grid( . ~ Lineage) + 
  scale_color_manual(values=c(col2F,col2P),name ="Organism :")
p3

# Y ~ Gene.family*Organism
p4 = ggplot(SW,aes(Gene.family,amino.acid.pearson.correlation)) + geom_boxplot(aes(fill=factor(Organism)),alpha=0.8) + theme.set2 + theme(axis.title.x=element_blank())
Anova(lmer(amino.acid.pearson.correlation ~ Organism+Gene.family + (1|Reference), data=SW))
p4 = p4 + ylab(ylab) + ggtitle("Sliding Window scores depending on the gene family\nand the organism")
p4 = p4 + facet_grid( . ~ Lineage) + 
  scale_fill_manual(values=c(col2F,col2P),name ="Organism :")
p4

plot_save2 = plot_grid(p1, p2, p3, p4, labels = 'AUTO')
ggsave(filename= "../RESULTS/5_Sliding_Window_2.pdf",plot_save2, width=12, height=10)
###########
