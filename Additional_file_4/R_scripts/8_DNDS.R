##################################################
# Statistical Analysis of Dinos ARN editing Data #
##################################################

library(lattice)
library(ggplot2)
library(cowplot)

############################
# DN/DS = KA/KS Estimation #
############################

#############
# Functions #
#############

get.synth.table <- function(dir, loop = c("Ch","Ht","Km","Kv","Lp","Pl","Sm"), 
                            suffix = "_edit.csv", sep = ','){
  global.table=NULL
  for (i in loop){
    table = read.delim(file = paste0(dir, i, suffix), sep = sep)
    table$organism = i
    global.table=rbind(global.table,table)
  }
  return(global.table)
}

get.codon.edit <- function(table.usage, table.edits, loop = c("Ch","Ht","Km","Kv","Lp","Pl","Sm")){
  codons = as.character(unique(table.usage$codon))
  for (i in loop){
    for (j in codons){
      row = which(table.usage$codon == j & table.usage$organism == i)
      value = nrow(subset(table.edits, genome.codon == j & organism == i))
      table.usage[row,"edit.usage"]=value
    }
  }
  return(table.usage)
}

add.degeneracy <- function(table, degeneracy){
  for (i in degeneracy$Nucleotide){
    rows = which(table$codon == i)
    table[rows, "synonymous"] = degeneracy[which(degeneracy$Nucleotide == i), "Synonymous"]
    table[rows, "adj.synonymous"] = degeneracy[which(degeneracy$Nucleotide == i), "Adj.Synonymous"]
    table[rows, "very.adj.synonymous"] = degeneracy[which(degeneracy$Nucleotide == i), "Very.Adj.Synonymous"]
    table[rows, "non.synonymous"] = degeneracy[which(degeneracy$Nucleotide == i), "Non.Synonymous"]
    table[rows, "adj.non.synonymous"] = degeneracy[which(degeneracy$Nucleotide == i), "Adj.Non.Synonymous"]
    table[rows, "very.adj.non.synonymous"] = degeneracy[which(degeneracy$Nucleotide == i), "Very.Adj.Non.Synonymous"]
  }
  return(table)
}

add.sp.degeneracy <- function(table, degeneracy){
  for (org in c("Km", "Kv", "Pl", "Sm", "Ch", "Ht")){
    for (i in degeneracy$Nucleotide){
      rows = which(table$codon == i & table$organism == org)
      table[rows, "sp.synonymous"] = degeneracy[which(degeneracy$Nucleotide == i & degeneracy$Organism == org), "Adj.Synonymous"]
      table[rows, "sp.non.synonymous"] = degeneracy[which(degeneracy$Nucleotide == i & degeneracy$Organism == org), "Adj.Non.Synonymous"]
    }
  }
  return(table)
}

add.redundancy <- function(table, redundancy){
  for (i in redundancy$Amino.Acid){
    rows = which(table$amino.acid == i)
    table[rows, "redundancy"] = redundancy[which(redundancy$Amino.Acid == i), 'Degeneracy']
  }
  return(table)
}

compute.dnds <- function(table, edits, loop = c("Ch","Ht","Km","Kv","Lp","Pl","Sm")){
  res = data.frame(organism = loop, dnds = rep(NA,length(loop)), adj.dnds = rep(NA,length(loop)),
                   very.adj.dnds = rep(NA,length(loop)))
  for (i in loop){
    rows = which(table$organism == i)
    exp.S = sum(table[rows,'genome.usage']*table[rows,'synonymous'])
    exp.N = sum(table[rows,'genome.usage']*table[rows,'non.synonymous'])
    adj.exp.S = sum(table[rows,'genome.usage']*table[rows,'adj.synonymous'])
    adj.exp.N = sum(table[rows,'genome.usage']*table[rows,'adj.non.synonymous'])
    very.adj.exp.S = sum(table[rows,'genome.usage']*table[rows,'very.adj.synonymous'])
    very.adj.exp.N = sum(table[rows,'genome.usage']*table[rows,'very.adj.non.synonymous'])

    rows = which(edits$organism == i)
    obs.S = sum(edits[rows,'genome.amino.acid'] == edits[rows,'edit.amino'],na.rm=T)
    obs.N = sum(edits[rows,'genome.amino.acid'] != edits[rows,'edit.amino'],na.rm=T)
    
    dnds = (obs.N/exp.N)/(obs.S/exp.S)
    adj.dnds = (obs.N/adj.exp.N)/(obs.S/adj.exp.S)
    very.adj.dnds = (obs.N/very.adj.exp.N)/(obs.S/very.adj.exp.S)

    res[which(res$organism == i), 'dnds'] = dnds
    res[which(res$organism == i), 'adj.dnds'] = adj.dnds
    res[which(res$organism == i), 'very.adj.dnds'] = very.adj.dnds
  }
  return(res)
}

compute.sp.dnds <- function(table, edits, loop = c("Kv","Km","Pl","Sm","Ht","Ch")){
  res = data.frame(organism = loop, dnds = rep(NA,length(loop)))
  for (i in loop){
    rows = which(table$organism == i)
    exp.S = sum(table[rows,'genome.usage']*table[rows,'sp.synonymous'])
    exp.N = sum(table[rows,'genome.usage']*table[rows,'sp.non.synonymous'])
    
    rows = which(edits$organism == i)
    obs.S = sum(edits[rows,'genome.amino.acid'] == edits[rows,'edit.amino'],na.rm=T)
    obs.N = sum(edits[rows,'genome.amino.acid'] != edits[rows,'edit.amino'],na.rm=T)
    
    dnds = (obs.N/exp.N)/(obs.S/exp.S)
    
    res[which(res$organism == i), 'dnds'] = dnds
  }
  return(res)
}
#############
# Variables #
#############

setwd("/Users/Lucas/Documents/ENS/BIO_M1_STAGE/Paper_Dino/ANALYSIS")

Short.names = c("Ch","Ht","Km","Kv","Lp","Pl","Sm")
Ext.names = c("C. horridum", "H. triquetra", "K. mikimotoi",
              "K. veneficum", "L. polyedrum", "P. lunula", "S. minutum")

input.dir = "/Users/Lucas/Dropbox/Chris/Very Final Stuff/Pipeline_CSV_files_for_Lucas/"

table.codon.usage = read.delim(paste0(input.directory, "Master_codon_pref_outfiles/Global_master_edit_codons.csv"))
levels(table.codon.usage$organism) = Short.names

nucleotide.degeneracy = read.delim("../DATA/nucleotide_table.csv")
nucleotide.degeneracy.sp = read.delim("../DATA/nucleotide_table_sp.csv")

codon.redundancy = read.delim("../DATA/amino_acid_table.csv")

############################
# DN/DS = KA/KS Estimation #
############################

summary.edit.table = get.synth.table(paste0(input.dir,"Gene_specific_edit_files/"))

# Correcting for multiple edit on the same codon. Taking each one independently
summary.edit.table$edit.codon=NA
summary.edit.table$edit.amino=NA
for (i in 1:nrow(summary.edit.table)){
  x = as.vector(summary.edit.table[i,"genome.codon"])
  pos = summary.edit.table[i,"codon.position"]
  substr(x,pos,pos) = as.vector(summary.edit.table[i,"mRNAbase"])
  summary.edit.table[i, 'edit.codon'] = x
  amino.acid = as.character(unique(table.codon.usage[which(table.codon.usage$codon == x), 'amino.acid']))
  amino.acid=ifelse(length(amino.acid)==0,NA,amino.acid)
  summary.edit.table[i, 'edit.amino'] = amino.acid
}

table.codon.usage = get.codon.edit(table.codon.usage, summary.edit.table)

table.codon.usage = add.degeneracy(table.codon.usage, nucleotide.degeneracy)
table.codon.usage = add.sp.degeneracy(table.codon.usage, nucleotide.degeneracy.sp)
table.codon.usage = add.redundancy(table.codon.usage, codon.redundancy)

# Compute dn/ds ratio :
dnds = compute.dnds(table.codon.usage, summary.edit.table)
dnds

sp.dnds = compute.sp.dnds(table.codon.usage, summary.edit.table)
sp.dnds


##############################
# Computing mutation rates : #
##############################

# synonymous all 
edit.1.2 = subset(summary.edit.table, codon.position == 1 | codon.position ==2)
edit.1.2 = subset(summary.edit.table, organism == "Km" | organism == "Kv" | organism == "Sm" | organism == "Pl")
edit.1.2$mutation = paste0(edit.1.2$genome.base," -> ",edit.1.2$mRNAbase)
mutation.rate = table(edit.1.2$mutation)/sum(table(edit.1.2$mutation))
# Ignoring degenerated bases

# Per species
for (i in Short.names){
  edit.1.2 = subset(summary.edit.table, codon.position == 1 | codon.position == 2)
  edit.1.2 = subset(summary.edit.table, organism == i)
  edit.1.2$mutation = paste0(edit.1.2$genome.base," -> ",edit.1.2$mRNAbase)
  mutation.rate = table(edit.1.2$mutation)/sum(table(edit.1.2$mutation))
  assign(paste0('mutation.',i),mutation.rate)
}
View(mutation.Km)
