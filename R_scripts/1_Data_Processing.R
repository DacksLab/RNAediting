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
###########
org = c("Kv", "Km", "Pl", "Sm")
org.all = c("Kv", "Km", "Pl", "Sm", "Ch", "Ht", "Lp")
ref = c("Pt", "Vb", "Ct", "Eh", "Ac")
lineage = c('Fucoxanthin','Fucoxanthin','Peridinin','Peridinin', NA, NA, NA)
names(lineage) = c('Kv','Km','Pl','Sm','Lp','Ht','Ch')
###########

###########
# Libraries
library(tidyverse)
###########

#############
# Functions
format.amino.acid.directory <- function(dir,organism,reference,org,ref){
  for (o in 1:length(organism)){
    for (r in 1:length(reference)){
      fdir = paste0(dir,organism[o],'/',reference[r],'/Spreadsheets/')
      files = list.files(fdir)
      genes = gsub(paste0(org[o],'_|_frame1_aminoacid_changes.csv'),'',files)
      for (g in genes){
        infile = paste0(fdir, org[o], '_', g, '_frame1_aminoacid_changes.csv')
        outfile = paste0(dir, org[o], '_', ref[r], '_', g, '_frame1_aminoacid_changes.csv')
        file.copy(infile, outfile)
      }
    }
  }
}
get.GC <- function(x){
  x = as.character(x)
  GC = sum(x %in% c('G','C'))/length(x)*100
  return(GC)
}

format.pipeline <- function(dir, org, lineage){
  
  table_raw = NULL
  table_codons = NULL
  table_types = NULL
  table_GC_position = NULL
  table_summary = NULL
  
  for (o in 1:length(org)){
    
    #Table Raw
    o_files = list.files(paste0(dir,'/Sheets/'))[grep(paste0(org[o],'.*'),list.files(paste0(dir,'/Sheets/')))]
    genes = gsub(paste0(org[o],'_|_basic_editing.csv'),'',o_files)
    for (g in genes){
      raw = read.delim(file=paste0(dir,'/Sheets/',org[o],'_',g,'_basic_editing.csv'),sep=",",stringsAsFactors=F)
      raw = na.omit(raw)
      if (nrow(raw) > 0){
        raw$Gene = g
        raw$Organism = org[o]
        raw$Lineage = lineage[[org[o]]]
        table_raw = rbind(table_raw, raw)
      }
    }

    # Table Codons
    o_codons_files = list.files(paste0(dir,'/Codon/'))[grep(paste0(org[o],'.*'),list.files(paste0(dir,'/Codon/')))]
    list_genes = gsub(paste0(org[o],'_|_codon_preference.csv'),'',o_codons_files)
    for (gene in list_genes){
      codons=read.delim(file=paste0(dir,'/Codon/',org[o],'_',gene,'_codon_preference.csv'),sep=",",check.names = F,stringsAsFactors = F)
      codons=codons[,names(codons)!='']
      codons=na.omit(codons)
      if (nrow(codons) > 0){
        codons$Organism = org[o]
        codons$Gene = gene
        codons$Lineage = lineage[[org[o]]]
        table_codons = rbind(table_codons, codons)
      }
    }
    
    # Table Types
    for (pos in c('first', 'second', 'third')){
      types=read.delim(file=paste0(dir,'/Types/',org[o],'_master_editing_out_',pos,'_pos_editing_types.csv'),sep=",",check.names = F,stringsAsFactors = F)
      types=types[,names(types)!='']
      types=na.omit(types)
      types$position = pos
      types$Organism = org[o]
      types$Lineage = lineage[[org[o]]]
      types$Gene.family = gsub('[A-Z0-9]','',types$gene)
      types$Gene.type = ifelse(types$Gene.family %in% c('psa','psb','atp','pet'),'Photosynthetic','Housekeeping')
      table_types = rbind(table_types, types)
    }
    
    # Table GC Positions  
    GC_position=read.delim(file=paste0(dir,'/Masterfiles/',org[o],'_master_editing_out_GC_position.csv'),sep=",",check.names = F,stringsAsFactors = F) 
    GC_position=GC_position[,names(GC_position)!='']
    GC_position=na.omit(GC_position)
    GC_position$Organism = org[o]
    GC_position$Lineage = lineage[[org[o]]]
    GC_position$Gene.family = gsub('[A-Z0-9]','',GC_position$gene)
    GC_position$Gene.type = ifelse(GC_position$Gene.family %in% c('psa','psb','atp','pet'),'Photosynthetic','Housekeeping')
    table_GC_position = rbind(table_GC_position, GC_position)
    
    # Table Summary  
    summary=read.delim(file=paste0(dir,'/Masterfiles/',org[o],'_master_editing_out.csv'),sep=",",check.names = F,stringsAsFactors = F) 
    summary=summary[,names(summary)!='']
    summary=na.omit(summary)
    summary$Organism = org[o]
    summary$Lineage = lineage[[org[o]]]
    summary$Gene.family = gsub('[A-Z0-9]','',summary$gene)
    summary$Gene.type = ifelse(summary$Gene.family %in% c('psa','psb','atp','pet'),'Photosynthetic','Housekeeping')
    table_summary = rbind(table_summary, summary)
  }
  
  table_raw$Gene.family = gsub('[A-Z0-9]','',table_raw$Gene)
  table_raw$Gene.type = ifelse(table_raw$Gene.family %in% c('psa','psb','atp','pet'),'Photosynthetic','Housekeeping')
  table_codons$Gene.family = gsub('[A-Z0-9]','',table_codons$Gene)
  table_codons$Gene.type = ifelse(table_codons$Gene.family %in% c('psa','psb','atp','pet'),'Photosynthetic','Housekeeping')
  table_translate = read.delim('../DATA/nucleotide_table.csv', check.names = F, stringsAsFactors = F)
  for (c in unique(table_codons$codon)){
    AA.n = which(table_translate$Nucleotide == c)
    AA = unique(table_translate[AA.n, 'Amino Acid'])
    c.n = which(table_codons$codon == c)
    table_codons[c.n, 'amino acid'] = AA
  }
  
  table_raw[table_raw$genome.base == "TRUE", "genome.base"] = "T"
  write.table(table_raw,file="../DATA/pipeline_raw.csv",sep=",",row.names=F)
  write.table(table_codons,file="../DATA/pipeline_codon_usage.csv",sep=",",row.names=F)
  write.table(table_types,file="../DATA/pipeline_types.csv",sep=",",row.names=F)
  write.table(table_GC_position,file="../DATA/pipeline_GC_position.csv",sep=",",row.names=F)
  write.table(table_summary,file="../DATA/pipeline_summary.csv",sep=",",row.names=F)
}
format.sliding.window <- function(dir, org, ref, lineage){
  
  table_summary = NULL
  for (r in 1:5){
    
    table_cat = NULL
    for (o in 1:length(org)){
      table=read.delim(file=paste0(dir,'/',org[o],'_',ref[r],'_sw_trim.csv'),sep=",")
      table=na.omit(table) 
      table$Organism = org[o]
      table$Lineage = lineage[[org[o]]]
      table$Reference=ref[r]
      table$Gene = gsub(paste0(org[o],'_|_trimmed'),'',table$name)
      table_cat = rbind(table_cat, table)
    }
    table_summary = rbind(table_summary, table_cat)
  }
  
  table_summary$Gene.family = gsub('[A-Z0-9]','',table_summary$Gene)
  table_summary$Gene.type = ifelse(table_summary$Gene.family %in% c('psa','psb','atp','pet'),'Photosynthetic','Housekeeping')
  
  write.table(table_summary,"../DATA/SW_analysis.csv",sep=",",row.names=F)
}
format.amino.acid <- function(dir, org, ref, lineage){
  
  table_summary = NULL
  for (r in 1:length(ref)){
    table_cat = NULL
    for (o in 1:length(org)){
      o_files = list.files(dir)[grep(paste0(org[o],'_',ref[r],'.*'),list.files(dir))]
      genes = gsub(paste0(org[o],'_',ref[r],'_|_frame1_aminoacid_changes.csv'),'',o_files)
      for (g in genes){
        table = read.delim(file=paste0(dir,'/',org[o],'_',ref[r],'_',g,'_frame1_aminoacid_changes.csv'),sep=",",stringsAsFactors=F)
        table=na.omit(table)
        if (nrow(table) > 0){
          table$Gene = g
          table$Organism = org[o]
          table$Reference = ref[r]
          table$Lineage = lineage[[org[o]]]
          table_cat = rbind(table_cat, table)
        }
      }
    }
    table_summary = rbind(table_summary, table_cat)
  }
  table_summary$Gene.family = gsub('[A-Z0-9]','',table_summary$Gene)
  table_summary$Gene.type = ifelse(table_summary$Gene.family %in% c('psa','psb','atp','pet'),'Photosynthetic','Housekeeping')
  write.table(table_summary,"../DATA/AA_analysis.csv",sep=",",row.names=F)
}
format.gravy <- function(dir, org, lineage){
  
  table_summary = NULL
  for (o in 1:length(org)){
    table=read.delim(file=paste0(dir,'/',org[o],'_gravy.csv'),sep=",")
    table=na.omit(table) 
    table$Organism = org[o]
    table$Lineage = lineage[[org[o]]]
    table_summary = rbind(table_summary, table)
  }
  
  table_summary$Gene.family = gsub('[A-Z0-9]','',table_summary$gene)
  table_summary$Gene.type = ifelse(table_summary$Gene.family %in% c('psa','psb','atp','pet'),'Photosynthetic','Housekeeping')
  table_summary$weight = table_summary[,"rna.mol.weight"]-table_summary[,"gen.mol.weight"]
  table_summary$gravy = table_summary[,"rna.gravy"]-table_summary[,"gen.gravy"]
  
  write.table(table_summary,file="../DATA/gravy_analysis.csv",sep=",",row.names=F)
}
format.entropy <- function(dir, org, lineage){
  
  table_summary = NULL
  for (o in 1:length(org)){
    o_files = list.files(dir)[grep(paste0(org[o],'.*'),list.files(dir))]
    genes = gsub(paste0(org[o],'_|_entropy.csv'),'',o_files)
    for (g in genes){
      table=read.delim(file=paste0(dir,'/',org[o],'_',g,'_entropy.csv'),sep=",")
      table=na.omit(table) 
      table$Gene=g
      table$Organism=org[o]
      table$Lineage=lineage[[o]]
      table_summary = rbind(table_summary, table)
    }
  }
  
  table_summary$Gene.family = gsub('[A-Z0-9]','',table_summary$Gene)
  table_summary$Gene.type = ifelse(table_summary$Gene.family %in% c('psa','psb','atp','pet'),'Photosynthetic','Housekeeping')
  table_summary[,"edited.or.not"]=gsub(1,"yes",table_summary[,"edited.or.not"])
  table_summary[,"edited.or.not"]=gsub(0,"no",table_summary[,"edited.or.not"])
  
  write.table(table_summary,file="../DATA/entropy_analysis.csv",sep=",",row.names=F)
}
#############

####################
# Edit File structure (if needed)
####################
# format.amino.acid.directory(dir = '../DATA/Dino/AminoAcid/',
#                             organism = c("Kveneficum", "Kmikimotoi", "Plunula", "Sminutum", "Htriquetra"),
#                             reference = c("Ptricornutum", "Vbrassicaformis", "Ctobin", "Ehuxleyi", "Amphidinium"),
#                             org = c("Kv", "Km", "Pl", "Sm", "Ht"),
#                             ref = c("Pt", "Vb", "Ct", "Eh", "Ac"))
####################

####################
# General Pipeline 
####################

pi.dir = '../DATA/Dino/Pipeline'
format.pipeline(pi.dir, org.all, lineage)
####################

##################
# Sliding Window
##################

sw.dir = '../DATA/Dino/SlidingWindow/'
format.sliding.window(sw.dir, org.all, ref, lineage)
##################

################################
# Editing Score on Amino Acids
################################

aa.dir = '../DATA/Dino/AminoAcid/'
format.amino.acid(aa.dir, org.all, ref, lineage)
################################

################
# Gravy Scores
################

gr.dir = '../DATA/Dino/GRAVY/'
format.gravy(gr.dir, org.all, lineage)
################

###########
# Entropy 
###########

en.dir = '../DATA/Dino/Entropy/'
format.entropy(en.dir, org, lineage)
###########
