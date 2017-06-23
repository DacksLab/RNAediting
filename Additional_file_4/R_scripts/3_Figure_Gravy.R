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

########################
# GRAVY SCORE ANALYSIS #
########################

############
# PACKAGES #
############

library(lattice)
library(ggplot2)
library(cowplot)
library(car)

table = read.delim(file="../DATA/gravy_global_table.csv",sep=",")
table.edit = read.delim(file="../DATA/table_analysis_global_lucas.csv",sep=",")

for (i in 1:nrow(table)){
  value = mean(subset(table.edit,as.character(gene)==as.character(table[i,"gene"]) & organism==table[i,"organism"])$average.edit.score.diff)
  table[i,"edit"]=value
}

attach(table)

########################
# STATISTICAL ANALYSIS #
########################

shapiro.test(weight)
hist(weight)
t.test(weight)

#Linear models testing all possible explanatory effects
gravy.lm1=lm(gravy~weight+plastid+gene_status+gene_family+
               weight*plastid+weight*gene_status+weight*gene_family+
               plastid*gene_status+plastid*gene_family+
               gene_status*gene_family)
summary(gravy.lm1)
Anova(gravy.lm1)

# Organism level instead of plastid
gravy.lm2=lm(gravy~weight+organism+gene_status+gene_family+
               weight*organism+weight*gene_status+weight*gene_family+
               organism*gene_status+organism*gene_family+
               gene_status*gene_family)
summary(gravy.lm2)
Anova(gravy.lm2)

anova(gravy.lm1,gravy.lm2)

# Simplifying the model
stats::step(gravy.lm1,direction='both')
stats::step(gravy.lm2,direction='both')

gravy.lm0=lm(gravy~weight)
summary(gravy.lm0)
anova(gravy.lm0,gravy.lm1)
Anova(gravy.lm0)

# Only the intercept is significant : Gravy score is better explained by a
# constant than by any of the models.

# Testing Spearman Rank correlation as done in former articles
test.all=cor.test(weight,gravy,method="spearman") # whole dataset

## FUCO
fuco=subset(table,plastid=="fucoxanthin")
test.fuco=cor.test(fuco$weight,fuco$gravy,method="spearman") # Fucoxanthin
test.fuco.CV=cor.test(subset(fuco,gene_status=="CV")$weight, # Fucoxanthin*Conserved
         subset(fuco,gene_status=="CV")$gravy,
         method="spearman")
test.fuco.NCV=cor.test(subset(fuco,gene_status=="NCV")$weight, # Fucoxanthin*Non-Conserved
         subset(fuco,gene_status=="NCV")$gravy,
         method="spearman")
test.fuco.noATPAB=cor.test(subset(fuco,gene_status=="CV" & gene!="atpA" & gene!="atpB")$weight, # Fucoxanthin*subset of Conserved
                           subset(fuco,gene_status=="CV" & gene!="atpA" & gene!="atpB")$gravy,
                           method="spearman")

## PERE
pere=subset(table,plastid=="peridinin") # Peridinin (=Peridinin*Conserved)
test.pere.noATPAB=cor.test(subset(pere,gene!="atpA" & gene!="atpB")$weight, # Fucoxanthin*subset of Conserved
                           subset(pere,gene!="atpA" & gene!="atpB")$gravy,
                           method="spearman")

test.pere=cor.test(pere$weight,pere$gravy,method="spearman")

# Organism level
kv=subset(table,organism=="kv")
km=subset(table,organism=="km")
pl=subset(table,organism=="pl")
sm=subset(table,organism=="sm")

## KV
test.kv=cor.test(kv$weight,kv$gravy,method="spearman")
test.kv.CV=cor.test(subset(kv,gene_status=="CV")$weight,
                    subset(kv,gene_status=="CV")$gravy,method="spearman")
test.kv.noATPAB=cor.test(subset(kv,gene_status=="CV" & gene!="atpA"  &gene!="atpB")$weight,
                    subset(kv,gene_status=="CV" & gene!="atpA"  &gene!="atpB")$gravy,method="spearman")
test.kv.NCV=cor.test(subset(kv,gene_status=="NCV")$weight,
                     subset(kv,gene_status=="NCV")$gravy,method="spearman")
## KM
test.km=cor.test(km$weight,km$gravy,method="spearman")
test.km.CV=cor.test(subset(km,gene_status=="CV")$weight,
                    subset(km,gene_status=="CV")$gravy,method="spearman")
test.km.noATPAB=cor.test(subset(km,gene_status=="CV" & gene!="atpA"  &gene!="atpB")$weight,
                         subset(km,gene_status=="CV" & gene!="atpA"  &gene!="atpB")$gravy,method="spearman")
test.km.NCV=cor.test(subset(km,gene_status=="NCV")$weight,
                     subset(km,gene_status=="NCV")$gravy,method="spearman")

## PL
test.pl=cor.test(pl$weight,pl$gravy,method="spearman")
test.pl.noATPAB=cor.test(subset(pl, gene!="atpA" & gene!="atpB")$weight,
                         subset(pl, gene!="atpA" & gene!="atpB")$gravy,method="spearman")
## SM
test.sm=cor.test(sm$weight,sm$gravy,method="spearman")
test.sm.noATPAB=cor.test(subset(sm, gene!="atpA" & gene!="atpB")$weight,
                         subset(sm, gene!="atpA" & gene!="atpB")$gravy,method="spearman")


##########
# FIGURE #
##########

col2=c("#1f78b4","#fb9a99")
col2F=c("#1f78b4","#a6cee3")
col2P=c("#e31a1c","#fb9a99")

theme.set = theme_grey() + theme(text = element_text(size = 8), legend.position = "top", plot.title = element_text(face="bold", color = "black", size=8),
                                 legend.title = NULL, legend.text = element_text(size = 7))

ylab = "Gravy score difference"
xlab = "Molecular weight difference"
title = "Gravy score difference in function of the mol. weight difference"

# Plot of the whole dataset
p1 = ggplot(table,aes(weight,gravy)) + geom_point(aes(color=factor(plastid)),shape=19, size=2.5) + theme.set
p1 = p1 + xlab(xlab) + ylab(ylab) + ggtitle(title) + 
  scale_color_manual(values=col2, name =NULL, guide = guide_legend(nrow = 2),
                     labels=c(paste("Fucoxanthin","/ rho :",round(test.fuco$estimate,4)," ; p.value :",round(test.fuco$p.value,4)),
                              paste("Peridinin","/ rho :",round(test.pere$estimate,4)," ; p.value :",round(test.pere$p.value,4))))

# Plot of the conserved genes
CV = subset(table, gene_status=="CV")
p2 = ggplot(CV,aes(CV$weight,CV$gravy)) + geom_point(aes(color=factor(CV$plastid)),shape = 15, size=2.5) + theme.set
p2 = p2 + xlab(xlab) + ylab(ylab) + ggtitle(paste0(title,"\n for conserved genes (photosynthetic genes)")) + 
  scale_color_manual(values=col2,name =NULL,guide = guide_legend(nrow = 2),
                     labels=c(paste("Fucoxanthin","/ rho :",round(test.fuco.CV$estimate,4)," ; p.value :",round(test.fuco.CV$p.value,4)),
                              paste("Peridinin","/ rho :",round(test.pere$estimate,4)," ; p.value :",round(test.pere$p.value,4))))

# Plot of the conserved genes w/o ATP A and B
CV.ATPAB = subset(table, gene_status=="CV" & gene!="atpA" & gene!="atpB")
p7 = ggplot(CV.ATPAB,aes(CV.ATPAB$weight,CV.ATPAB$gravy)) + geom_point(aes(color=factor(CV.ATPAB$plastid)), shape = 18, size=3) + theme.set
p7 = p7 + xlab(xlab) + ylab(ylab) + ggtitle(paste0(title,"\n for conserved genes w/o atp A,B")) + 
  scale_color_manual(values=col2,name =NULL,guide = guide_legend(nrow = 2),
                     labels=c(paste("Fucoxanthin","/ rho :",round(test.fuco.noATPAB$estimate,4)," ; p.value :",round(test.fuco.noATPAB$p.value,4)),
                              paste("Peridinin","/ rho :",round(test.pere.noATPAB$estimate,4)," ; p.value :",round(test.pere.noATPAB$p.value,4))))

# Plot of the peridinin lineage
p3 = ggplot(pere,aes(pere$weight,pere$gravy)) + geom_point(aes(color=factor(pere$organism)),shape=15, size=2.5) + theme.set
p3 = p3 + xlab(xlab) + ylab(ylab) + ggtitle(paste0(title, "\n for the peridinin lineage")) +
  scale_color_manual(values=col2P,name =NULL,guide = guide_legend(nrow = 2),
                     labels=c(paste("P. lunula","/ rho :",round(test.pl$estimate,4)," ; p.value :",round(test.pl$p.value,4)),
                              paste("S. minutum","/ rho :",round(test.sm$estimate,4)," ; p.value :",round(test.sm$p.value,4))))

# Plot of the peridinin lineage w/o atp A,B
pere.noATBAB = subset(pere, gene!="atpA" & gene!="atpB")
p8 = ggplot(pere.noATBAB,aes(pere.noATBAB$weight,pere.noATBAB$gravy)) + geom_point(aes(color=factor(pere.noATBAB$organism)),shape=18, size=3) + theme.set
p8 = p8 + xlab(xlab) + ylab(ylab) + ggtitle(paste0(title, "\n for the peridinin lineage w/o atp A,B")) +
  scale_color_manual(values=col2P,name =NULL,guide = guide_legend(nrow = 2),
                     labels=c(paste("P. lunula","/ rho :",round(test.pl.noATPAB$estimate,4)," ; p.value :",round(test.pl.noATPAB$p.value,4)),
                              paste("S. minutum","/ rho :",round(test.sm.noATPAB$estimate,4)," ; p.value :",round(test.sm.noATPAB$p.value,4))))

# plot of the fucoxanthin lineage
p4 = ggplot(fuco,aes(fuco$weight,fuco$gravy)) + geom_point(aes(color=factor(fuco$organism)),shape=19, size=2.5) + theme.set
p4 = p4 + xlab(xlab) + ylab(ylab) + ggtitle(paste0(title, "\n for the fucoxanthin lineage")) +
  scale_color_manual(values=col2F,name =NULL,guide = guide_legend(nrow = 2),
                     labels=c(paste("K. mikimotoi","/ rho :",round(test.km$estimate,4)," ; p.value :",round(test.km$p.value,4)),
                              paste("K. veneficum","/ rho :",round(test.kv$estimate,4)," ; p.value :",round(test.kv$p.value,4))))

# plot of the fucoxanthin lineage and conserved genes
fuco.CV = subset(fuco,gene_status=="CV")
p5 = ggplot(fuco.CV, aes(fuco.CV$weight, fuco.CV$gravy)) + geom_point(aes(color=factor(fuco.CV$organism)),shape=15, size=2.5) + theme.set
p5 = p5 + xlab(xlab) + ylab(ylab) + ggtitle(paste0(title, "\n for the fucoxanthin lineage and conserved genes")) +
  scale_color_manual(values=col2F,name =NULL,guide = guide_legend(nrow = 2),
                     labels=c(paste("K. mikimotoi","/ rho :",round(test.km.CV$estimate,4)," ; p.value :",round(test.km.CV$p.value,4)),
                              paste("K. veneficum","/ rho :",round(test.kv.CV$estimate,4)," ; p.value :",round(test.kv.CV$p.value,4))))

# plot of the fucoxanthin lineage and conserved genes w/o atp A,B
fuco.CV.noATPAB = subset(fuco,gene_status=="CV", gene!="atpA" & gene!="atpB")
p9 = ggplot(fuco.CV.noATPAB, aes(fuco.CV.noATPAB$weight, fuco.CV$gravy)) + geom_point(aes(color=factor(fuco.CV.noATPAB$organism)),shape=18, size=3) + theme.set
p9 = p9 + xlab(xlab) + ylab(ylab) + ggtitle(paste0(title, "\n for the fucoxanthin lineage and conserved genes w/o atp A,B")) +
  scale_color_manual(values=col2F,name =NULL,guide = guide_legend(nrow = 2),
                     labels=c(paste("K. mikimotoi","/ rho :",round(test.km.noATPAB$estimate,4)," ; p.value :",round(test.km.noATPAB$p.value,4)),
                              paste("K. veneficum","/ rho :",round(test.kv.noATPAB$estimate,4)," ; p.value :",round(test.kv.noATPAB$p.value,4))))


# plot of the fucoxanthin lineage and non conserved genes
fuco.NCV = subset(fuco,gene_status=="NCV")
p6 = ggplot(fuco.NCV,aes(fuco.NCV$weight,fuco.NCV$gravy)) + geom_point(aes(color=factor(fuco.NCV$organism)),shape=17, size=2.5) + theme.set
p6 = p6 + xlab(xlab) + ylab(ylab) + ggtitle(paste0(title, "\nfor the fucoxanthin lineage and non conserved genes")) + 
  scale_color_manual(values=col2F,name =NULL,guide = guide_legend(nrow = 2),
                     labels=c(paste("K. mikimotoi","/ rho :",round(test.km.NCV$estimate,4)," ; p.value :",round(test.km.NCV$p.value,4)),
                              paste("K. veneficum","/ rho :",round(test.kv.NCV$estimate,4)," ; p.value :",round(test.kv.NCV$p.value,4))))


########################
# Gravy and Edit score #
########################

lm.w.1 = lm(weight ~ plastid + organism + edit + edit*plastid + edit*organism)
summary(lm.w.1)
Anova(lm.w.1)

stats::step(lm.w.1, direction="both")
lm.w.0 = lm(weight ~ edit)
summary(lm.w.0)
Anova(lm.w.0)

lm.g.1 = lm(gravy ~ plastid + organism + edit + edit*plastid + edit*organism)
summary(lm.g.1)
Anova(lm.g.1)

stats::step(lm.g.1, direction="both")
lm.g.0 = lm(gravy ~ edit+organism)
summary(lm.g.0)
Anova(lm.g.0)

p.edit.w = ggplot(table, aes(x=edit, y=weight)) + geom_point(aes(color=plastid),size=2.5) + theme.set 
p.edit.w = p.edit.w + xlab("Editing score differences") + ylab("Mol. weight") + ggtitle("Mol. weight as a function of the edit score",subtitle = "No significant pattern") + 
  scale_color_manual(values=col2, name =NULL, guide = guide_legend(nrow = 1),
                     labels=c(paste("Fucoxanthin"),
                              paste("Peridinin")))

p.edit.g = ggplot(table, aes(x=edit, y=gravy)) + geom_point(aes(color=plastid),size=2.5) + theme.set 
p.edit.g = p.edit.g + xlab("Editing score differences") + ylab("Gravy score") + ggtitle("Gravy score as a function of the edit score",subtitle = "No significant pattern") + 
  scale_color_manual(values=col2, name =NULL, guide = guide_legend(nrow = 1),
                     labels=c(paste("Fucoxanthin"),
                              paste("Peridinin")))

pdf(file = "../RESULTS/3_gravy_score.pdf", width=24, height=9)
plot_grid(p1, p2, p7, p3, p8, p4, p5, p9, p6,p.edit.w, p.edit.g, labels = 'AUTO', nrow = 2,hjust=-1)
dev.off()
