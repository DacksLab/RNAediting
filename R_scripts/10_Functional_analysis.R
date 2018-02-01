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
Short.names = c("Km","Kv","Pl","Sm")
Ext.names = c("K. mikimotoi", "K. veneficum", "P. lunula", "S. minutum")
input.dir = "/Users/Lucas/Dropbox/Chris/Very Final Stuff/New_Sheets_Lucas/"
genes = list.files(input.dir)

col2=c("#1f78b4","#fb9a99")
col2F=c("#1f78b4","#a6cee3")
col2P=c("#e31a1c","#fb9a99")
col4=c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3")
col5=c("#8dd3c7","#ffffb3","#bebada","#fb8072","#80b1d3")
###########

###########
# Libraries
library(tidyverse)
library(cowplot)
###########

#################################
# Function & Edits distribution #
#################################

################
# Distribution #
################

# Load and concatenate summary tables
summary_table = NULL
for (g in genes){
  temp=read.delim(paste0(input.dir,'/',g,'/',g,'_summary.csv'),row.names=NULL,sep=',',stringsAsFactors = F,header=F)
  summary_table = rbind.fill(summary_table, temp)
}

# Tidy table
summary_table = separate(data = summary_table, col = V1, into = c("Organism", "Gene"), sep = "_")
summary_table$Gene = gsub(".afa","", summary_table$Gene)
names(summary_table)=c('Organism', 'Genes', 'Number.Edits', 'Aligned.Positions', 'Edits.in.Functionals', 'Aligned.in.Functionals', 
                       'Lumen', 'Edits.Lumen', 'Aligned.Lumen', 
                       'Stroma', 'Edits.Stroma', 'Aligned.Stroma',
                       'TMD', 'Edits.TMD', 'Aligned.TMD')
summary_table$TMD = gsub('transmembrane domain.*','transmembrane domain',summary_table$TMD)
summary_table$Plastid = ifelse(summary_table$Organism == "Kv" | summary_table$Organism == "Km", 'Fucoxanthin', 'Peridinin')

# Subset for organisms of interest : Km, Kv, Pl, Sm

## Build sub df functional vs non functional
df = data.frame(value = c( (summary_table$Edits.in.Functionals / summary_table$Aligned.in.Functionals )*100,
                           (summary_table$Number.Edits - summary_table$Edits.in.Functionals) / (summary_table$Aligned.Positions - summary_table$Aligned.in.Functionals)*100),
                Type = c(rep('Functional',nrow(summary_table)),
                         rep('Non Functional',nrow(summary_table))),
                Plastid = summary_table$Plastid)
df$Type.Plastid = paste0(df$Type,':',df$Plastid)

# Check distribution
ggplot(df) + geom_density(aes(value)) + theme_grey()
summary(df$value ~ df$Type.Plastid)
pairwise.wilcox.test(df$value, df$Type.Plastid)

df$Plastid = gsub('Peridinin', 'Peridinin (**)', df$Plastid)

p1 = ggplot(df) + geom_boxplot(aes(x=Type, y=value, fill=Type.Plastid)) + theme_grey() + 
  theme(legend.position='none', axis.title.x=element_blank()) + facet_wrap(~Plastid) + scale_fill_manual(values=c("#1f78b4","#e31a1c","#a6cee3","#fb9a99"))
p1


## Build sub dataframe to investigate distrib in different domains
dom.df = data.frame(Organism = rep(summary_table$Organism,3),
                    Plastid = rep(summary_table$Plastid,3),
                    Gene = rep(summary_table$Genes,3),
                    Value = c(summary_table$Edits.Lumen/summary_table$Aligned.Lumen*100,
                              summary_table$Edits.Stroma/summary_table$Aligned.Stroma*100,
                              summary_table$Edits.TMD/summary_table$Aligned.TMD*100),
                    Env = c(summary_table$Lumen, summary_table$Stroma, summary_table$TMD))
dom.df$Env.Plastid = paste0(dom.df$Env,':',dom.df$Plastid)

ggplot(dom.df) + geom_density(aes(Value)) + theme_grey()
summary(dom.df$Value~dom.df$Env.Plastid)
pairwise.wilcox.test(dom.df$Value, dom.df$Env.Plastid)
pairwise.t.test(dom.df$Value, dom.df$Env.Plastid)

p2 = ggplot(dom.df) + geom_boxplot(aes(x = Plastid, y = Value, fill = Plastid)) + facet_wrap(~Env) + theme_grey()
p2 = p2 + scale_fill_manual(values=col2) + theme(legend.position = 'none', axis.title.x=element_blank())
p2

##################
# Correctiveness #
##################

# Load and cat tables
func_table = NULL
for (g in genes){
  for (org in Short.names){
    if (paste0(org,'_',g,'.csv') %in% list.files(paste0(input.dir,'/',g,'/'))){
      temp=read.delim(paste0(input.dir,'/',g,'/',org,'_',g,'.csv'),row.names=NULL,sep=',',stringsAsFactors = F)
      temp$gene = g
      temp$organism = org
      func_table = rbind(func_table, temp)
    }
  }
}

# Tidy table
func_table[func_table == ""]=NA
func_table$biochemical.environment=gsub('stroma','Stroma',func_table$biochemical.environment)
func_table$biochemical.environment=gsub(' Stroma','Stroma',func_table$biochemical.environment)
func_table$biochemical.environment=gsub('lumen','Lumen',func_table$biochemical.environment)
func_table$biochemical.environment[is.na(func_table$biochemical.environment)] = 'none'
table(func_table$biochemical.environment)


func_table$annotated = ifelse(!is.na(func_table$interface) | !is.na(func_table$known.function), 'Functional', 'Non Functional')

func_df = data.frame(Score = c(func_table$dino.score, func_table$ancestral.score, func_table$consensus.score),
                     Type = c(rep('Dino', nrow(func_table)), rep('Ancestral', nrow(func_table)), rep('Consensus', nrow(func_table))),
                     annotated = rep(func_table$annotated, 3),
                     Env = rep(func_table$biochemical.environment, 3),
                     Plastid = ifelse(func_table$organism == "Kv" | func_table$organism == "Km", 'Fucoxanthin', 'Peridinin'))
func_df=func_df[!(func_df$Env=='none'),]

ggplot(func_df) + geom_density(aes(Score)) + theme_grey()

summary(func_df$Score ~ func_df$Type)
summary(func_df$Score ~ paste0(func_df$Type,func_df$annotated))
summary(func_df$Score ~ paste0(func_df$Type,func_df$Env))
summary(func_df$Score ~ paste0(func_df$Type,func_df$Plastid))

score.lm = lm(data=func_df, Score ~ annotated*Type + Type*Env + Plastid*Type +
                Type*Plastid*Env + Type*annotated*Plastid + Type*annotated*Env)
summary(score.lm)
Anova(score.lm)
score.lm0 = lm(data=func_df, Score ~ annotated*Type + Type*Env + Plastid*Type)
summary(score.lm0)
Anova(score.lm0)

# No type effect -> Subselect Dino.

sub.func_df = subset(func_df, Type == 'Consensus')
summary(sub.func_df$Score ~ sub.func_df$Env)
summary(sub.func_df$Score ~ sub.func_df$Plastid)
summary(sub.func_df$Score ~ sub.func_df$annotated)
summary(sub.func_df$Score ~ paste0(sub.func_df$Plastid,sub.func_df$Env))
summary(sub.func_df$Score ~ paste0(sub.func_df$Plastid,sub.func_df$annotated))

dino.lm = lm(data=sub.func_df, Score ~ Env + Plastid + annotated +
               Plastid*Env + Plastid*annotated + annotated*Env +
               Plastid*Env*annotated)
summary(dino.lm)
Anova(dino.lm)
dino.lm0 = lm(data=sub.func_df, Score ~ Env + Plastid + annotated +
                Plastid*annotated)
summary(dino.lm0)
Anova(dino.lm0)

pairwise.t.test(sub.func_df$Score,paste0(sub.func_df$Plastid,sub.func_df$annotated))
pairwise.t.test(sub.func_df$Score,paste0(sub.func_df$Env))
pairwise.t.test(sub.func_df$Score,paste0(sub.func_df$Plastid,sub.func_df$Env))

p3.sub.func_df = sub.func_df
p3.sub.func_df$Plastid = ifelse(p3.sub.func_df$Plastid == 'Fucoxanthin', 'Fucoxanthin (*)', 'Peridinin (***)')

p3 = ggplot(p3.sub.func_df)+geom_boxplot(aes(x=annotated,y=Score, fill = Plastid))+theme_grey()+facet_wrap(~Plastid)
p3 = p3 + scale_fill_manual(values = col2) + theme(legend.position = 'none', axis.title.x=element_blank()) 
p3

p4 = ggplot(sub.func_df)+geom_boxplot(aes(x = Env, y = Score, fill = Plastid))+theme_grey()+facet_wrap(~Plastid)
p4 = p4 + scale_fill_manual(values = col2) + theme(legend.position = 'none', axis.title.x=element_blank()) 
p4
