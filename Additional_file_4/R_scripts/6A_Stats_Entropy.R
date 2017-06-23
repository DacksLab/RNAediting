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

####################
# ENTROPY ANALYSIS #
####################

library(car)
library(lsmeans)

################
# LOADING DATA #
################

table = read.table(file="../DATA/entropy_global_table.csv",sep=",",header=T)

attach(table)

View(table) # Inspecting the Data

#########
# Stats #
#########

# Investing entropy relations with other variables.

hist(positional.entropy)

# Due to the distribution we will use statistical tests for each situation :

# Qualititative effects

# Test for a plastid effect : none
table(positional.entropy,plastid)
chisq.test(table(positional.entropy,plastid))
fisher.test(positional.entropy,plastid)

# Test for an organism effect : none
table(positional.entropy,organism)
chisq.test(table(positional.entropy,organism))

# Test for an edition effect : very significant
table(positional.entropy,edited.or.not)
chisq.test(table(positional.entropy,edited.or.not))

# Interactions

# Testing 2 by 2, plastid effects among edition categories : no effects
x=subset(table,edited.or.not=="yes")
chisq.test(table(x$positional.entropy,x$plastid))
y=subset(table,edited.or.not=="no")
chisq.test(table(y$positional.entropy,y$plastid))

# Testing 2 by 2, organism effects among edition categories : no effects
chisq.test(table(x$positional.entropy,x$organism))
chisq.test(table(y$positional.entropy,y$organism))

# Quantitative effects

# Correlation with score
cor.test(positional.entropy,score)

# Score:plastid
cor.test(subset(table,plastid=="fucoxanthin")$positional.entropy,
                subset(table,plastid=="fucoxanthin")$score)

cor.test(subset(table,plastid=="peridinin")$positional.entropy,
         subset(table,plastid=="peridinin")$score)

# Score:organism
cor.test(subset(table,organism=="kv")$positional.entropy,
         subset(table,organism=="kv")$score)
cor.test(subset(table,organism=="km")$positional.entropy,
         subset(table,organism=="km")$score)
cor.test(subset(table,organism=="sm")$positional.entropy,
         subset(table,organism=="sm")$score)
cor.test(subset(table,organism=="pl")$positional.entropy,
         subset(table,organism=="pl")$score)

# Test of the sign of the score :
x[,"score"]=x$score>0
chisq.test(table(x$positional.entropy,x$score))

########################
# EDIT SCORE ~ ENTROPY #
########################

# Selecting only edited position as they are the one of interest.
z=subset(table,edited.or.not=="yes")

# Creating a linear model
lm.score=lm(score~positional.entropy+plastid+positional.entropy*plastid,data=z)

summary(lm.score)

# Significance
Anova(lm.score)

# Validation
plot(lm.score)

# Distribution of score : normal so we can use t-tests
hist(z$score)

# Defining conserved entropy>0.95 regions based on the entropy analysis (see graphs)
z[,"positional.entropy"]=z$positional.entropy>=0.95

#Testing conserved (TRUE) versus non conserved (FALSE)
t.test(z$score~z$positional.entropy)
#Testing plastids
t.test(z$score~z$plastid)

#Testing for conservation*plastid
z=cbind(z,paste0(z$positional.entropy,z$plastid))

# Their is a significant interaction :
pairwise.t.test(z$score,z[,length(z)],p.adjust.method = "fdr")
mean(subset(z,z[,length(z)]=="TRUEfucoxanthin")$score)
mean(subset(z,z[,length(z)]=="TRUEperidinin")$score)
mean(subset(z,z[,length(z)]=="FALSEfucoxanthin")$score)
mean(subset(z,z[,length(z)]=="FALSEperidinin")$score)

##############################################
# FIGURE : EDIT SCORE ~ CONSERVATION*PLASTID #
##############################################

theme.set = theme_grey() + theme(text = element_text(size = 10), legend.position = "none", plot.title = element_text(face="bold", color = "black", size=12),
                                 legend.title = NULL, legend.text = element_text(size = 7), axis.title.x=element_blank(), plot.subtitle = element_text(size=8))

cols = c("#1f78b4","#fb9a99","#1f78b4","#fb9a99")

p = ggplot(z, aes(x=z[,length(z)],y=score)) + geom_boxplot(aes(fill=z[,length(z)])) + theme.set
p = p + ylim(-15,16) + ylab("Score") + ggtitle("Editing score depending on the conservation and the plastid type") +
  scale_fill_manual(values=cols)+scale_x_discrete(labels=c("Non conserved\nFucoxanthin",
                                                           "Non conserved\nPeredinin",
                                                           "Conserved\nFucoxanthin",
                                                           "Conserved\nPeredinin"))
p = p + annotate("text", x = c(1,2,3,4), y = rep(16,4), size = 3, fontface = "italic",
                 label=c(paste0("n = ",nrow(subset(z,z[,length(z)]=="FALSEfucoxanthin"))),
                         paste0("n = ",nrow(subset(z,z[,length(z)]=="FALSEperidinin"))),
                         paste0("n = ",nrow(subset(z,z[,length(z)]=="TRUEfucoxanthin"))),
                         paste0("n = ",nrow(subset(z,z[,length(z)]=="TRUEperidinin")))))
p = p + annotate("segment", x = 3, xend = 4, y = 14.4, yend = 14.4) + annotate("text", x = 3.5, y = 14.7,label="***") +
  annotate("segment", x = 2, xend = 4, y = 13.2, yend = 13.2) + annotate("text", x = 3, y = 13.5,label="***") +
  annotate("segment", x = 1, xend = 4, y = 12, yend = 12) + annotate("text", x = 2.5, y = 12.3,label="***")
p

pdf(file= "../RESULTS/6_Editing_Score_Entropy.pdf", width=7, height=7)
p
dev.off()

