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
setwd("/Users/Lucas/Documents/Publications/2017_Dino/ANALYSIS")
Short.names = c("Ch","Ht","Km","Kv","Lp","Pl","Sm")
Ext.names = c("C. horridum", "H. triquetra", "K. mikimotoi",
              "K. veneficum", "L. polyedrum", "P. lunula", "S. minutum")
table.codon.usage = read_csv(file="../DATA/pipeline_codon_usage.csv")
summary.edit.table = read_csv(file="../DATA/pipeline_raw.csv")
dnds_input = read_tsv("../DATA/dnds_input.tsv")
dnds_input.sp = read_tsv("../DATA/dnds_input_sp.tsv")
codon.redundancy = read_tsv("../DATA/amino_acid_table.tsv")
nucleotide.degeneracy = read_tsv("../DATA/nucleotide_table.tsv")
###########

###########
# Libraries
library(tidyverse)
library(cowplot)
###########

############
# Functions 
############

get.codon.edit <- function(table.usage, table.edits, loop = c("Ch","Ht","Km","Kv","Lp","Pl","Sm")){
  table.usage$edit.usage = NA
  codons = as.character(unique(table.usage$codon))
  for (i in loop){
    for (j in codons){
      row = which(table.usage$codon == j & table.usage$Organism == i)
      value = nrow(subset(table.edits, genome.codon == j & Organism == i))
      table.usage[row,"edit.usage"]=value
    }
  }
  return(table.usage)
}

add.degeneracy <- function(table, dnds_input){
  for (i in dnds_input$Nucleotide){
    rows = which(table$codon == i)
    table[rows, "synonymous"] = dnds_input[which(dnds_input$Nucleotide == i), "Synonymous"]
    table[rows, "adj.synonymous"] = dnds_input[which(dnds_input$Nucleotide == i), "Adj.Synonymous"]
    table[rows, "very.adj.synonymous"] = dnds_input[which(dnds_input$Nucleotide == i), "Very.Adj.Synonymous"]
    table[rows, "non.synonymous"] = dnds_input[which(dnds_input$Nucleotide == i), "Non.Synonymous"]
    table[rows, "adj.non.synonymous"] = dnds_input[which(dnds_input$Nucleotide == i), "Adj.Non.Synonymous"]
    table[rows, "very.adj.non.synonymous"] = dnds_input[which(dnds_input$Nucleotide == i), "Very.Adj.Non.Synonymous"]
  }
  return(table)
}

add.sp.degeneracy <- function(table, dnds_input){
  for (org in c("Km", "Kv", "Pl", "Sm", "Ch", "Ht")){
    for (i in dnds_input$Nucleotide){
      rows = which(table$codon == i & table$Organism == org)
      table[rows, "sp.synonymous"] = dnds_input[which(dnds_input$Nucleotide == i & dnds_input$Organism == org), "Adj.Synonymous"]
      table[rows, "sp.non.synonymous"] = dnds_input[which(dnds_input$Nucleotide == i & dnds_input$Organism == org), "Adj.Non.Synonymous"]
    }
  }
  return(table)
}

add.redundancy <- function(table, redundancy){
  for (i in redundancy$`Amino Acid`){
    rows = which(table$`amino acid` == i)
    table[rows, "redundancy"] = redundancy[which(redundancy$`Amino Acid` == i), 'Redundancy']
  }
  return(table)
}

compute.dnds <- function(table, edits, loop = c("Ch","Ht","Km","Kv","Lp","Pl","Sm")){
  res = data.frame(Organism = loop, dnds = rep(NA,length(loop)), adj.dnds = rep(NA,length(loop)),
                   very.adj.dnds = rep(NA,length(loop)))
  for (i in loop){
    rows = which(table$Organism == i)
    exp.S = sum(table[rows,'genome usage']*table[rows,'synonymous'])
    exp.N = sum(table[rows,'genome usage']*table[rows,'non.synonymous'])
    adj.exp.S = sum(table[rows,'genome usage']*table[rows,'adj.synonymous'])
    adj.exp.N = sum(table[rows,'genome usage']*table[rows,'adj.non.synonymous'])
    very.adj.exp.S = sum(table[rows,'genome usage']*table[rows,'very.adj.synonymous'])
    very.adj.exp.N = sum(table[rows,'genome usage']*table[rows,'very.adj.non.synonymous'])

    rows = which(edits$Organism == i)
    obs.S = sum(edits[rows,'genome.amino.acid'] == edits[rows,'edit.amino'],na.rm=T)
    obs.N = sum(edits[rows,'genome.amino.acid'] != edits[rows,'edit.amino'],na.rm=T)
    
    dnds = (obs.N/exp.N)/(obs.S/exp.S)
    adj.dnds = (obs.N/adj.exp.N)/(obs.S/adj.exp.S)
    very.adj.dnds = (obs.N/very.adj.exp.N)/(obs.S/very.adj.exp.S)

    res[which(res$Organism == i), 'dnds'] = dnds
    res[which(res$Organism == i), 'adj.dnds'] = adj.dnds
    res[which(res$Organism == i), 'very.adj.dnds'] = very.adj.dnds
  }
  return(res)
}

compute.sp.dnds <- function(table, edits, loop = c("Kv","Km","Pl","Sm","Ht","Ch")){
  res = data.frame(Organism = loop, dnds = rep(NA,length(loop)))
  for (i in loop){
    rows = which(table$Organism == i)
    exp.S = sum(table[rows,'genome usage']*table[rows,'sp.synonymous'])
    exp.N = sum(table[rows,'genome usage']*table[rows,'sp.non.synonymous'])
    
    rows = which(edits$Organism == i)
    obs.S = sum(edits[rows,'genome.amino.acid'] == edits[rows,'edit.amino'],na.rm=T)
    obs.N = sum(edits[rows,'genome.amino.acid'] != edits[rows,'edit.amino'],na.rm=T)
    
    dnds = (obs.N/exp.N)/(obs.S/exp.S)
    
    res[which(res$Organism == i), 'dnds'] = dnds
  }
  return(res)
}
############

############################
# DN/DS = KA/KS Estimation #
############################

# Correcting for multiple edit on the same codon. Taking each one independently
summary.edit.table$edit.codon=NA
summary.edit.table$edit.amino=NA
for (i in 1:nrow(summary.edit.table)){
  x = as.vector(summary.edit.table[[i,"genome.codon"]])
  pos = summary.edit.table[[i,"codon.position"]]
  substr(x,pos,pos) = as.vector(summary.edit.table[[i,"mRNAbase"]])
  summary.edit.table[[i, 'edit.codon']] = x
  amino.acid = as.character(unique(table.codon.usage[which(table.codon.usage$codon == x),]$`amino acid`))
  amino.acid=ifelse(length(amino.acid)==0,NA,amino.acid)
  summary.edit.table[[i, 'edit.amino']] = amino.acid
}

table.codon.usage = get.codon.edit(table.codon.usage, summary.edit.table)
table.codon.usage = add.degeneracy(table.codon.usage, dnds_input)
table.codon.usage = add.sp.degeneracy(table.codon.usage, dnds_input.sp)
table.codon.usage = add.redundancy(table.codon.usage, codon.redundancy)

# Compute dn/ds ratio :
dnds = compute.dnds(table.codon.usage, summary.edit.table)
dnds
write_tsv(dnds, '../DATA/dnds_output.tsv')

sp.dnds = compute.sp.dnds(table.codon.usage, summary.edit.table)
sp.dnds
write_tsv(sp.dnds, '../DATA/dnds_output.sp.tsv')
############################