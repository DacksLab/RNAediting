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
liste = c("Km","Kv","Pl","Sm","Ch","Ht","Lp")
liste2 = c("Ch","Ht","Km","Kv","Lp","Pl","Sm")
Species = c("C. horridum",
            "H. triquetra",
            "K. mikimotoi",
            "K. veneficum",
            "L. polyedrum",
            "P. lunula",
            "S. minutum")
table.raw = read.delim(file="../DATA/pipeline_raw.csv", sep = ',')
table.raw$mutation = paste0(table.raw$genome.base,'.to.',table.raw$mRNAbase)
###########

###########
# Libraries
library(tidyverse)
library(cowplot)
###########

table.mutations = as.data.frame(matrix(0,13,7))
names(table.mutations) = liste

for (i in liste){
  table.mutations[1,i]=nrow(subset(table.raw,Organism==i & mutation == 'A.to.T'))
  table.mutations[2,i]=nrow(subset(table.raw,Organism==i & mutation == 'T.to.A'))
  table.mutations[3,i]=nrow(subset(table.raw,Organism==i & mutation == 'A.to.G'))
  table.mutations[4,i]=nrow(subset(table.raw,Organism==i & mutation == 'G.to.A'))
  table.mutations[5,i]=nrow(subset(table.raw,Organism==i & mutation == 'A.to.C'))
  table.mutations[6,i]=nrow(subset(table.raw,Organism==i & mutation == 'C.to.A'))
  table.mutations[7,i]=nrow(subset(table.raw,Organism==i & mutation == 'T.to.C'))
  table.mutations[8,i]=nrow(subset(table.raw,Organism==i & mutation == 'C.to.T'))
  table.mutations[9,i]=nrow(subset(table.raw,Organism==i & mutation == 'G.to.C'))
  table.mutations[10,i]=nrow(subset(table.raw,Organism==i & mutation == 'C.to.G'))
  table.mutations[11,i]=nrow(subset(table.raw,Organism==i & mutation == 'T.to.G'))
  table.mutations[12,i]=nrow(subset(table.raw,Organism==i & mutation == 'G.to.T'))
  table.mutations[13,i]=nrow(subset(table.raw,Organism==i))
}

table.mutations[1:12,]=sweep(table.mutations[1:12,],2,colSums(table.mutations[1:12,]),"/")
write_tsv(table.mutations, '../DATA/Mutation.rates.tsv')

# 1st position
table.mutations = as.data.frame(matrix(0,13,7))
names(table.mutations) = liste2
table.first = subset(table.raw, codon.position == 1)

for (i in liste2){
  table.mutations[1,i]=nrow(subset(table.first,Organism==i & mutation == 'A.to.T'))
  table.mutations[2,i]=nrow(subset(table.first,Organism==i & mutation == 'T.to.A'))
  table.mutations[3,i]=nrow(subset(table.first,Organism==i & mutation == 'A.to.G'))
  table.mutations[4,i]=nrow(subset(table.first,Organism==i & mutation == 'G.to.A'))
  table.mutations[5,i]=nrow(subset(table.first,Organism==i & mutation == 'A.to.C'))
  table.mutations[6,i]=nrow(subset(table.first,Organism==i & mutation == 'C.to.A'))
  table.mutations[7,i]=nrow(subset(table.first,Organism==i & mutation == 'T.to.C'))
  table.mutations[8,i]=nrow(subset(table.first,Organism==i & mutation == 'C.to.T'))
  table.mutations[9,i]=nrow(subset(table.first,Organism==i & mutation == 'G.to.C'))
  table.mutations[10,i]=nrow(subset(table.first,Organism==i & mutation == 'C.to.G'))
  table.mutations[11,i]=nrow(subset(table.first,Organism==i & mutation == 'T.to.G'))
  table.mutations[12,i]=nrow(subset(table.first,Organism==i & mutation == 'G.to.T'))
  table.mutations[13,i]=nrow(subset(table.first,Organism==i))
}

write_tsv(table.mutations, '../DATA/Mutation.rates.first.tsv')

# 2nd position
table.mutations = as.data.frame(matrix(0,13,7))
names(table.mutations) = liste2
table.second = subset(table.raw, codon.position == 2)

for (i in liste2){
  table.mutations[1,i]=nrow(subset(table.second,Organism==i & mutation == 'A.to.T'))
  table.mutations[2,i]=nrow(subset(table.second,Organism==i & mutation == 'T.to.A'))
  table.mutations[3,i]=nrow(subset(table.second,Organism==i & mutation == 'A.to.G'))
  table.mutations[4,i]=nrow(subset(table.second,Organism==i & mutation == 'G.to.A'))
  table.mutations[5,i]=nrow(subset(table.second,Organism==i & mutation == 'A.to.C'))
  table.mutations[6,i]=nrow(subset(table.second,Organism==i & mutation == 'C.to.A'))
  table.mutations[7,i]=nrow(subset(table.second,Organism==i & mutation == 'T.to.C'))
  table.mutations[8,i]=nrow(subset(table.second,Organism==i & mutation == 'C.to.T'))
  table.mutations[9,i]=nrow(subset(table.second,Organism==i & mutation == 'G.to.C'))
  table.mutations[10,i]=nrow(subset(table.second,Organism==i & mutation == 'C.to.G'))
  table.mutations[11,i]=nrow(subset(table.second,Organism==i & mutation == 'T.to.G'))
  table.mutations[12,i]=nrow(subset(table.second,Organism==i & mutation == 'G.to.T'))
  table.mutations[13,i]=nrow(subset(table.second,Organism==i))
}

write_tsv(table.mutations, '../DATA/Mutation.rates.second.tsv')

# Position
table.position = as.data.frame(matrix(0,7,4))
row.names(table.position) = liste

for (i in liste){
  table.position[i,1]=nrow(subset(table.raw,Organism==i & codon.position == 1))
  table.position[i,2]=nrow(subset(table.raw,Organism==i & codon.position == 2))
  table.position[i,3]=nrow(subset(table.raw,Organism==i & codon.position == 3))
  table.position[i,4]=nrow(subset(table.raw,Organism==i))
}
write_tsv(table.position, '../DATA/Positions.rate.tsv')


