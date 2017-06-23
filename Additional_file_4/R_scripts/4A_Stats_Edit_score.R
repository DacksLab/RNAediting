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
# EDITING SCORE ANALYSIS #
##########################

library(car)
library(lsmeans)
library(nlme)
library(lme4)
library(multcomp)
library(lmerTest)

setwd("/Users/Lucas/Documents/ENS/BIO_M1_STAGE/Paper_Dino/ANALYSIS")

table=read.delim(file="../DATA/table_analysis_global_lucas.csv",sep=",")

attach(table)

View(table) # Inspecting the Data

# "Reference" causes dependance in the data
# We ought to investigate this effect using a mixed model using (1|gene), gene being the random effect

# Complete model, with all potential effects
score.lmer1=lmer(average.edit.score.diff~num.AA.changes+gene_family+reference+plastid+
                   gene_family*plastid+plastid*reference+gene_family*reference+
                   (1|gene),REML=F)
summary(score.lmer1)
step(score.lmer1)
Anova(score.lmer1)
# Removing reference effects
score.lmer0=lmer(average.edit.score.diff~plastid+
                   (1|organism),REML=F)

# lmer0 is a more efficient model :
anova(score.lmer1,score.lmer0)
# Fixed effect significance
Anova(score.lmer0)
# Adding gene_family as it was signigicant in the first model
score.lmer=lmer(average.edit.score.diff~gene_family+plastid+(1|organism),REML=F)

# Validation : Plastid and gene_family are significant fixed effect
anova(score.lmer,score.lmer0)
summary(score.lmer)
Anova(score.lmer)
anova(score.lmer)

# Validation of the model
plot(score.lmer)
qqnorm(resid(score.lmer))
abline(0,1)

table(plastid,gene_family) # Contingence table

# Mean comparison, using lsmeans and glht 
# to perform the tests using Tukey post hoc tests.
lsmeans(score.lmer)
summary(glht(score.lmer,linfct=mcp(plastid="Tukey")))
summary(glht(score.lmer,linfct=mcp(gene_family="Tukey")))

r.squaredGLMM(score.lmer)

# Investigating potential organism effects.
###########################################
# As plastid is a linear combinaison of organism, this is a different model.

# Complete model, with all potential effects
org.lmer1=lmer(average.edit.score.diff~num.AA.changes+gene_family+reference+organism+
                 gene_family*organism+organism*reference+gene_family*reference+
                 (1|gene),REML=F)
summary(org.lmer1)
step(org.lmer1)
Anova(org.lmer1)

# Step result model
org.lmer0=lmer(average.edit.score.diff~organism+(1|gene),REML=F)

# comparing the two models : no significant differences
anova(org.lmer1,org.lmer0)

# Model based on the Anova() results
org.lmer=lmer(average.edit.score.diff~gene_family+organism+num.AA.changes+
                gene_family*organism+
                (1|gene),REML=F)

summary(org.lmer)
vcov(org.lmer)>0.5

# Validation : Anova() model slightly more efficient but has less meaning
anova(org.lmer,org.lmer0)
Anova(org.lmer)
anova(org.lmer)

# comparing with plastid organism
anova(score.lmer,org.lmer)

# Validation of the model
plot(org.lmer)
qqnorm(resid(org.lmer))
abline(0,1)

table(organism,gene_family) # Contingence table

# Mean comparison, using lsmeans
lsmeans(org.lmer)
lsmeans(org.lmer0)
summary(glht(org.lmer0,linfct=mcp(organism="Tukey")))

r.squaredGLMM(org.lmer0)
r.squaredGLMM(org.lmer)
