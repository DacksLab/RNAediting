##################################################
# Plastid transcript editing across dinoflagellate
# lineages shows lineage-specific application 
# but conserved trends
##################################################

#  Wrote by Lucas Paoli. Contact : lucas.paoli@ens.fr | lucas.paoli@gmail.com

############
rm(list=ls())
# Directory
setwd("/Users/Lucas/Documents/Publications/2017_Dino/ANALYSIS")
# Variables
liste=c("Ch","Ht","Km","Kv","Lp","Pl","Sm")
global.table = read_csv(file="../DATA/pipeline_raw.csv")
table.codon.usage = read_csv(file='../DATA/pipeline_codon_usage.csv')
amino = read_tsv('../DATA/amino_acid_table.csv')
nucl = read_tsv('../DATA/nucleotide_table.csv')
lineage = c('Fucoxanthin','Fucoxanthin','Peridinin','Peridinin', NA, NA, NA)
names(lineage) = c('Kv','Km','Pl','Sm','Lp','Ht','Ch')
theme.set = theme_grey() + theme(legend.position = "top", 
                                 legend.title = element_blank(),
                                 plot.title = element_text(face="bold", color = "black", size=12))
###########

###########
# Libraries
library(tidyverse)
library(cowplot)
###########

##################
# Formatting data
##################
amino.global=NULL
nucl.global=NULL

for (j in liste){
  for (i in 1:nrow(amino)){
    amino[i,'usage']=nrow(subset(global.table,mRNA.amino.acid == amino[[i,'Amino Acid']] & organism == j))/
      nrow(subset(global.table,organism == j))
    amino[i,'expected']=sum(subset(table.codon.usage, `amino acid` == amino[[i,'Amino Acid']] & organism == j)$`genome usage`)/
      sum(subset(table.codon.usage,organism == j)$`genome usage`)
  }
  amino$organism=j
  amino$lineage=lineage[j]
  amino.global=rbind(amino.global,amino)
  
  for (i in 1:nrow(nucl)){
    nucl[i,'usage']=nrow(subset(global.table, mRNA.codon == as.character(nucl[[i,'Nucleotide']]) & organism == j))/
      nrow(subset(global.table,organism == j))
    nucl[i,'expected']=sum(subset(table.codon.usage, codon == as.character(nucl[[i,'Nucleotide']]) & organism == j)$`genome usage`)/
      sum(subset(table.codon.usage,organism == j)$`genome usage`)
  }
  nucl$organism=j
  nucl$lineage=lineage[j]
  nucl.global=rbind(nucl.global,nucl)
}

nucl.global$variation = nucl.global$usage - nucl.global$expected
amino.global$variation = amino.global$usage - amino.global$expected
##################

##################
# Statistics
##################
n.Xt = chisq.test(nucl.global$usage,nucl.global$expected)
n.Xtf = chisq.test(subset(nucl.global, lineage == 'Fucoxanthin')$usage,
           subset(nucl.global, lineage == 'Fucoxanthin')$expected)
n.Xtp = chisq.test(subset(nucl.global, lineage == 'Peridinin')$usage,
           subset(nucl.global, lineage == 'Peridinin')$expected)

a.Xt = chisq.test(amino.global$usage,amino.global$expected)
a.Xtf = chisq.test(subset(amino.global, lineage == 'Fucoxanthin')$usage,
           subset(amino.global, lineage == 'Fucoxanthin')$expected)
a.Xtp = chisq.test(subset(amino.global, lineage == 'Peridinin')$usage,
           subset(amino.global, lineage == 'Peridinin')$expected)

ml1 = lm(variation ~ organism + Degeneracy, data=nucl.global)
summary(ml1)
Anova(ml1)
ml1.1 = lm(variation ~ Degeneracy, data=nucl.global)
summary(ml1.1)
Anova(ml1.1)

ml2 = lm(variation ~ organism + Redundancy, data=amino.global)
summary(ml2)
Anova(ml2)
ml2.1 = lm(variation ~ Redundancy, data=amino.global)
summary(ml2.1)
Anova(ml2.1)
##################

##########
# Figures
##########

# Codon degeneracy
pn1 = ggplot(nucl.global, aes(x=Nucleotide)) + 
  geom_boxplot(aes(y=usage, fill = 'Edited'), alpha = .7) + 
  geom_boxplot(aes(y=expected, fill = 'Genomic'), alpha = .7) + 
  theme.set + theme(axis.text.x = element_text(angle = 90)) +
  ggtitle('Proportions of genomic codons and edited codons',
          subtitle = paste0('Chi-Square test : ', format(n.Xtf$p.value))) +
  xlab('Codons') + ylab('')
pn1

pn2 = ggplot(nucl.global, aes(x=Nucleotide, y=variation, fill='Variation')) + 
  geom_boxplot() + 
  theme.set + theme(axis.text.x = element_text(angle = 90), legend.position = 'none') +
  ggtitle('Variation between the proportion of genomic codons and edited codons') +
  xlab('Codons') + ylab('Difference')
pn2

pn3 = ggplot(nucl.global, aes(x=as.character(Degeneracy))) + 
  geom_boxplot(aes(y=usage, fill='Edited'), alpha = .7) + 
  geom_boxplot(aes(y=expected, fill='Genomic'), alpha = .7) + 
  theme.set +
  ggtitle('Proportion of genomic codons and edited codons based on degeneracy') +
  xlab('n fold degeneracy') + ylab('')
pn3

pn4 = ggplot(nucl.global, aes(x=as.character(Degeneracy), y=variation, fill='Variation')) + 
  geom_boxplot() + 
  theme.set + theme(legend.position = 'none') +
  ggtitle('Difference between the proportion of genomic codons and edited codons\nbased on degeneracy') +
  xlab('n fold degeneracy') + ylab('Difference')
pn4

#Amino Acid redundancy
pa1 = ggplot(amino.global, aes(x=`Amino Acid`)) + 
  geom_boxplot(aes(y=usage, fill='Edited'), alpha = .7) + 
  geom_boxplot(aes(y=expected, fill='Genomic'), alpha = .7) + 
  theme.set + theme(axis.text.x = element_text(angle = 90)) +
  ggtitle('Proportion of genomic amino acids and edited amino acids',
          subtitle = paste0('Chi-Square test : ', format(a.Xtf$p.value))) +
  xlab('Amino Acids') + ylab('')
pa1

pa2 = ggplot(amino.global, aes(x=`Amino Acid`, y=variation, fill='Variation')) + 
  geom_boxplot() + 
  theme.set + theme(axis.text.x = element_text(angle = 90), legend.position='none') +
  ggtitle('Difference between the proportion of genomic amino acids and edited amino acids') +
  xlab('Amino Acids') + ylab('Difference')
pa2

pa3 = ggplot(amino.global, aes(x=as.character(Redundancy))) + 
  geom_boxplot(aes(y=usage, fill='Edited'), alpha = .7) + 
  geom_boxplot(aes(y=expected, fill='Genomic'), alpha = .7) + 
  theme.set + 
  ggtitle('Proportion of genomic amino acids and edited amino acids based on redundancy') +
  xlab('Amino Acids') + ylab('')
pa3

pa4 = ggplot(amino.global, aes(x=as.character(Redundancy), y=variation, fill='Variation')) + 
  geom_boxplot() + 
  theme.set + theme(legend.position='none') +
  ggtitle('Difference between the proportion of genomic amino acids and edited amino acids\nbased on redundancy') +
  xlab('Amino Acids') + ylab('Difference')
pa4

plot = plot_grid(pn1,pn2,pn3,pn4,pa1,pa2,pa3,pa4, labels = 'AUTO', nrow = 4,hjust=-1)
ggsave(filename = "../RESULTS/SX_dgeneracy_redundancy.pdf", plot, width=15, height=18)
##########



