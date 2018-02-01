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
liste = c("Ch","Ht","Km","Kv","Lp","Pl","Sm")
Species = c("C. horridum",
            "H. triquetra",
            "K. mikimotoi",
            "K. veneficum",
            "L. polyedrum",
            "P. lunula",
            "S. minutum")
table = read.delim(file="../DATA/pipeline_summary.csv", sep = ',')
table.raw = read.delim(file="../DATA/pipeline_raw.csv", sep = ',')
table.types = read.delim(file="../DATA/pipeline_types.csv", sep = ',')
table.codon.raw = read.delim(file="../DATA/pipeline_codon_usage.csv", sep = ',')
attach(table)
###########

###########
# Libraries
###########
library(tidyverse)
library(cowplot)
###########

###########
# TABLE 1 #
###########

table_1 = as.data.frame(matrix(0,7,10))
row.names(table_1) = liste

for (i in liste){
  table_1[i,1]=nrow(subset(table,Organism==i))
  table_1[i,2]=sum(subset(table,Organism==i)$nucleotide.length)
  table_1[i,3]=sum(subset(table,Organism==i)$amino.acid.length)
  table_1[i,4]=sum(subset(table,Organism==i)$number.edits)
  table_1[i,5]=mean((subset(table,Organism==i)$number.edits)/(subset(table,Organism==i)$nucleotide.length)*100)
  table_1[i,6]=sd((subset(table,Organism==i)$number.edits)/(subset(table,Organism==i)$nucleotide.length)*100)
  table_1[i,7]=mean(subset(table,Organism==i)$percent.non.synonymous.edits)
  table_1[i,8]=sd(subset(table,Organism==i)$percent.non.synonymous.edits)
  table_1[i,9]=mean((subset(table,Organism==i)$number.amino.acid.edits)/(subset(table,Organism==i)$amino.acid.length)*100)
  table_1[i,10]=sd((subset(table,Organism==i)$number.amino.acid.edits)/(subset(table,Organism==i)$amino.acid.length)*100)
}

table_1 = cbind(Species, table_1)
names(table_1)=c("Species", "# Genes","Nuc Length", "AA Length", "# Edits", "Avg % Edits","+-" ,"Avg % Non-Syn","+- (%)" , "Avg % AA Change","+- (%)")

write_tsv(table_1, '../RESULTS/table_1.tsv')

# Stats comparing : Kv, Km and Pl, Sm
###########
t1.kv = subset(table,Organism=='Kv') ; t1.km = subset(table,Organism=='Km')
list.genes = intersect(t1.kv$gene, t1.km$gene)

length(list.genes)

t.test((subset(t1.km, gene %in% list.genes)$number.edits)/(subset(t1.km, gene %in% list.genes)$nucleotide.length)*100,
       (subset(t1.kv, gene %in% list.genes)$number.edits)/(subset(t1.kv, gene %in% list.genes)$nucleotide.length)*100)

t1.pl = subset(table,Organism=='Pl') ; t1.sm = subset(table,Organism=='Sm')
list.genes = intersect(t1.pl$gene, t1.sm$gene)

length(list.genes)

t.test((subset(t1.pl, gene %in% list.genes)$number.edits)/(subset(t1.pl, gene %in% list.genes)$nucleotide.length)*100,
       (subset(t1.sm, gene %in% list.genes)$number.edits)/(subset(t1.sm, gene %in% list.genes)$nucleotide.length)*100)
###########

###########
# TABLE 2 #
###########
table.types_2 = as.data.frame(matrix(0,7,12))
row.names(table.types_2) = liste

for (i in liste){
  table.types_2[i,1]=sum(subset(table.types,Organism==i)$A.to.T)
  table.types_2[i,2]=sum(subset(table.types,Organism==i)$A.to.G)
  table.types_2[i,3]=sum(subset(table.types,Organism==i)$A.to.C)
  table.types_2[i,4]=sum(subset(table.types,Organism==i)$T.to.A)
  table.types_2[i,5]=sum(subset(table.types,Organism==i)$T.to.G)
  table.types_2[i,6]=sum(subset(table.types,Organism==i)$T.to.C)
  table.types_2[i,7]=sum(subset(table.types,Organism==i)$G.to.A)
  table.types_2[i,8]=sum(subset(table.types,Organism==i)$G.to.T)
  table.types_2[i,9]=sum(subset(table.types,Organism==i)$G.to.C)
  table.types_2[i,10]=sum(subset(table.types,Organism==i)$C.to.A)
  table.types_2[i,11]=sum(subset(table.types,Organism==i)$C.to.T)
  table.types_2[i,12]=sum(subset(table.types,Organism==i)$C.to.G)
}

table.types_2=sweep(table.types_2,1,rowSums(table.types_2),"/")*100

table.types_2 = cbind(Species, table.types_2)
names(table.types_2)=c("Species", "A/T","A/G", "A/C", "T/A", "T/G","T/C" ,"G/A","G/T" ,"G/C","C/A","C/T","C/G")

write_tsv(table.types_2, '../RESULTS/table_2.tsv')
###########

###########
# TABLE 3 #
###########

table_3 = as.data.frame(matrix(0,7,5))
row.names(table_3) = liste

for (i in liste){
  table_3[i,1]=mean(subset(table,Organism==i)$GC.before)
  table_3[i,2]=sd(subset(table,Organism==i)$GC.before)
  table_3[i,3]=mean(subset(table,Organism==i)$GC.after)
  table_3[i,4]=sd(subset(table,Organism==i)$GC.after)
  if (i != "Lp"){
    hist(subset(table,Organism==i)$GC.before)
    hist(subset(table,Organism==i)$GC.after)
    table_3[i,5]=t.test(subset(table,Organism==i)$GC.before,subset(table,Organism==i)$GC.after)$p.value}
  else if (i == "Lp"){
    table_3[i,5] = NA
  }
}

table_3 = cbind(Species, table_3)
names(table_3)=c("Species", "Avg %GC Before","+-", "Avg %GC After", "+- (%)", "P value")

write_tsv(table_3, '../RESULTS/table_3.tsv')

# Stats comparing : Fuco and Peri
###########
hist(subset(table, Lineage == 'Fucoxanthin')$GC.after-subset(table, Lineage == 'Fucoxanthin')$GC.before)
hist(subset(table, Lineage == 'Peridinin')$GC.after-subset(table, Lineage == 'Peridinin')$GC.before)

t.test(subset(table, Lineage == 'Fucoxanthin')$GC.after-subset(table, Lineage == 'Fucoxanthin')$GC.before,
       subset(table, Lineage == 'Peridinin')$GC.after-subset(table, Lineage == 'Peridinin')$GC.before)
###########

###########

###########
# TABLE 4 #
###########

table_4=as.data.frame(matrix(0,7,6))
row.names(table_4)=liste

for (i in liste){
  table_4[i,1]=sum(subset(table,Organism==i)$first.position.edits)
  table_4[i,2]=sum(subset(table,Organism==i)$second.position.edits)
  table_4[i,3]=sum(subset(table,Organism==i)$third.position.edits)
  table_4[i,4]=mean(subset(table,Organism==i)$percent.edits.in.first.two.positions)
  table_4[i,5]=sd(subset(table,Organism==i)$percent.edits.in.first.two.positions)
  if (i != "Lp"){
    table_4[i,6]=chisq.test(c(table_4[i,1],table_4[i,2],table_4[i,3]))$p.value}
  else if (i == 'Lp'){
    table_4[i,6]=NA
  }
}

table_4 = cbind(Species, table_4)
names(table_4)=c("Species", "# 1st Pos Edits","# 2nd Pos Edits", "# 3rd Pos Edits","Avg % 1st Two Pos", "+-", "P value")

write_tsv(table_4, '../RESULTS/table_4.tsv')
###########

###########
# TABLE 5 #
###########

table_5=as.data.frame(matrix(0,7,10))
row.names(table_5)=liste

for (i in liste){
  table_5[i,1]=mean(subset(table,Organism==i)$fraction.polyT.before)
  table_5[i,2]=sd(subset(table,Organism==i)$fraction.polyT.before)
  table_5[i,3]=mean(subset(table,Organism==i)$fraction.polyT.after)
  table_5[i,4]=sd(subset(table,Organism==i)$fraction.polyT.after)
  if (i != "Lp"){
    hist(subset(table,Organism==i)$fraction.polyT.before)
    hist(subset(table,Organism==i)$fraction.polyT.after)
    table_5[i,5]=t.test(subset(table,Organism==i)$fraction.polyT.before,
                        subset(table,Organism==i)$fraction.polyT.after)$p.value}
  else if (i == "Lp"){
    table_5[i,5]=NA
  }
  table_5[i,6]=mean(subset(table,Organism==i)$fraction.70.percent.polyT.before)
  table_5[i,7]=sd(subset(table,Organism==i)$fraction.70.percent.polyT.before)
  table_5[i,8]=mean(subset(table,Organism==i)$fraction.70.percent.polyT.after)
  table_5[i,9]=sd(subset(table,Organism==i)$fraction.70.percent.polyT.after)
  if (i != "Lp"){
    hist(subset(table,Organism==i)$fraction.70.percent.polyT.before)
    hist(subset(table,Organism==i)$fraction.70.percent.polyT.after)
    table_5[i,10]=t.test(subset(table,Organism==i)$fraction.70.percent.polyT.before,
                         subset(table,Organism==i)$fraction.70.percent.polyT.after)$p.value}
  else if (i == 'Lp'){
    table_5[i,10]=NA
  }
}

table_5 = cbind(Species, table_5)
names(table_5)=c("Species", "% polyT before","+-","% polyT after","+- (%)", "P value","% 70 polyT Before", "+- (%)","% 70 polyT After", "+- (%)", "P value")

write_tsv(table_5, '../RESULTS/table_5.tsv')
###########

###########
# TABLE S1 #
###########

table.codon_S1=as.data.frame(matrix(0,21,7))
liste.AA=c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","STOP","T","V","W","Y")
row.names(table.codon_S1)=liste.AA
names(table.codon_S1)=liste

table.codon = table.codon.raw %>%
  group_by(amino.acid, codon, Organism) %>%
  summarize(genome.usage = sum(genome.usage),
            mRNA.usage = sum(mRNA.usage))

table.codon = table.codon[-which(table.codon$mRNA.usage == 0),]

for (j in liste){
  for (i in liste.AA){
    g.usage = subset(table.codon,Organism==j & amino.acid==i)$genome.usage
    r.usage = subset(table.codon,Organism==j & amino.acid==i)$mRNA.usage
    if (length(unique(c(g.usage, r.usage))) > 0){
      table.codon_S1[i,j]=chisq.test(cbind(g.usage, r.usage),correct=F)$p.value
    } else {
      table.codon_S1[i,j] = NA
    } 
  }
}

table.codon_S1 = cbind(liste.AA, table.codon_S1)
names(table.codon_S1)=c('Amino Acid', Species)
write_tsv(table.codon_S1, '../RESULTS/table_S1.tsv')
###########

####################
# Table Simulation #
####################
# 
# setwd("/Users/Lucas/Dropbox/Chris/Very Final Stuff/Pipeline_CSV_files_for_Lucas/Simulation_stats")
# table.simulation=read.delim(file="Global_simulation.csv")
# 
# setwd("/Users/Lucas/Dropbox/Chris/Very Final Stuff/Pipeline_CSV_files_for_Lucas/Master_edit_type_outfiles")
# table.edit=read.delim(file="Global_master_edit_editing_types.csv")
# 
# table.edit[,3:length(table.edit)]=sweep(table.edit[,3:length(table.edit)],1,rowSums(table.edit[,3:length(table.edit)]),"/")*100
# is.nan.data.frame <- function(x)
#   do.call(cbind, lapply(x, is.nan))
# table.edit[is.nan(table.edit)] <- 0
# table.edit[,3:length(table.edit)]=round(table.edit[,3:length(table.edit)], digits = 2)
# 
# setwd("/Users/Lucas/Dropbox/Chris/Very Final Stuff/Pipeline_CSV_files_for_Lucas/Master_outfiles")
# table=read.delim(file="global_master_edit.csv")
# 
# table.position=table[,c("Organism","gene","first.position.edits","second.position.edits","third.position.edits")]
# table.position[,3:length(table.position)]=sweep(table.position[,3:length(table.position)],1,rowSums(table.position[,3:length(table.position)]),"/")*100
# table.position[is.nan(table.position)] <- 0
# table.position[,3:length(table.position)]=round(table.position[,3:length(table.position)], digits = 2)
# 
# table.stats=table.simulation[,1:4]
# table.stats[,3]=rep(NA,nrow(table.stats))
# table.stats[,4]=rep(NA,nrow(table.stats))
# names(table.stats)[3:4]=c("edit position","edit type")
# 
# levels(table.simulation$Organism)=unique(c(levels(table.simulation$Organism),levels(table.position$Organism)))
# levels(table.stats$Organism)=unique(c(levels(table.simulation$Organism),levels(table.position$Organism)))
# levels(table.simulation$gene)=unique(c(levels(table.simulation$gene),levels(table.position$gene)))
# levels(table.stats$gene)=unique(c(levels(table.simulation$gene),levels(table.position$gene)))
# 
# 
# for (i in 1:nrow(table.stats)){
#   
#   transit1=cbind(c(subset(table.position,Organism==table.stats[i,"Organism"] & gene==table.stats[i,"gene"])$first.position.edits,
#                    subset(table.position,Organism==table.stats[i,"Organism"] & gene==table.stats[i,"gene"])$second.position.edits,
#                    subset(table.position,Organism==table.stats[i,"Organism"] & gene==table.stats[i,"gene"])$third.position.edits),
#                  
#                  c(subset(table.simulation,Organism==table.stats[i,"Organism"] & gene==table.stats[i,"gene"])$num.1st.pos,
#                    subset(table.simulation,Organism==table.stats[i,"Organism"] & gene==table.stats[i,"gene"])$num.2nd.pos,
#                    subset(table.simulation,Organism==table.stats[i,"Organism"] & gene==table.stats[i,"gene"])$num.3rd.pos))
#   
#   table.stats[i,3]=chisq.test(transit1[rowSums(transit1)>0,],correct=F)$p.value
#   
#   
#   transit2=cbind(c(subset(table.edit,Organism==table.stats[i,"Organism"] & gene==table.stats[i,"gene"])$A.to.T,
#                    subset(table.edit,Organism==table.stats[i,"Organism"] & gene==table.stats[i,"gene"])$A.to.G,
#                    subset(table.edit,Organism==table.stats[i,"Organism"] & gene==table.stats[i,"gene"])$A.to.C,
#                    subset(table.edit,Organism==table.stats[i,"Organism"] & gene==table.stats[i,"gene"])$T.to.A,
#                    subset(table.edit,Organism==table.stats[i,"Organism"] & gene==table.stats[i,"gene"])$T.to.G,
#                    subset(table.edit,Organism==table.stats[i,"Organism"] & gene==table.stats[i,"gene"])$T.to.C,
#                    subset(table.edit,Organism==table.stats[i,"Organism"] & gene==table.stats[i,"gene"])$G.to.A,
#                    subset(table.edit,Organism==table.stats[i,"Organism"] & gene==table.stats[i,"gene"])$G.to.T,
#                    subset(table.edit,Organism==table.stats[i,"Organism"] & gene==table.stats[i,"gene"])$G.to.C,
#                    subset(table.edit,Organism==table.stats[i,"Organism"] & gene==table.stats[i,"gene"])$C.to.A,
#                    subset(table.edit,Organism==table.stats[i,"Organism"] & gene==table.stats[i,"gene"])$C.to.T,
#                    subset(table.edit,Organism==table.stats[i,"Organism"] & gene==table.stats[i,"gene"])$C.to.G),
#                  
#                  c(subset(table.simulation,Organism==table.stats[i,"Organism"] & gene==table.stats[i,"gene"])$A.to.T,
#                    subset(table.simulation,Organism==table.stats[i,"Organism"] & gene==table.stats[i,"gene"])$A.to.G,
#                    subset(table.simulation,Organism==table.stats[i,"Organism"] & gene==table.stats[i,"gene"])$A.to.C,
#                    subset(table.simulation,Organism==table.stats[i,"Organism"] & gene==table.stats[i,"gene"])$T.to.A,
#                    subset(table.simulation,Organism==table.stats[i,"Organism"] & gene==table.stats[i,"gene"])$T.to.G,
#                    subset(table.simulation,Organism==table.stats[i,"Organism"] & gene==table.stats[i,"gene"])$T.to.C,
#                    subset(table.simulation,Organism==table.stats[i,"Organism"] & gene==table.stats[i,"gene"])$G.to.A,
#                    subset(table.simulation,Organism==table.stats[i,"Organism"] & gene==table.stats[i,"gene"])$G.to.T,
#                    subset(table.simulation,Organism==table.stats[i,"Organism"] & gene==table.stats[i,"gene"])$G.to.C,
#                    subset(table.simulation,Organism==table.stats[i,"Organism"] & gene==table.stats[i,"gene"])$C.to.A,
#                    subset(table.simulation,Organism==table.stats[i,"Organism"] & gene==table.stats[i,"gene"])$C.to.T,
#                    subset(table.simulation,Organism==table.stats[i,"Organism"] & gene==table.stats[i,"gene"])$C.to.G))
#   
#   table.stats[i,4]=chisq.test(transit2[rowSums(transit2)>0,],correct=F)$p.value
# }
# 
# setwd("/Users/Lucas/Dropbox/Chris/Very Final Stuff/Pipeline_CSV_files_for_Lucas/Simulation_stats")
# write.table(table.stats,"stats_results.csv",sep=",",row.names=F)
# 
