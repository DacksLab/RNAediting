##################################################
# Statistical Analysis of Dinos ARN editing Data #
##################################################

library(lattice)
library(ggplot2)
library(cowplot)

################################
# Positional distrib in codons #
################################

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

#############
# Variables #
#############

setwd("/Users/Lucas/Documents/ENS/BIO_M1_STAGE/Paper_Dino/ANALYSIS")

Short.names = c("Ch","Ht","Km","Kv","Lp","Pl","Sm")
Ext.names = c("C. horridum", "H. triquetra", "K. mikimotoi",
              "K. veneficum", "L. polyedrum", "P. lunula", "S. minutum")

input.dir = "/Users/Lucas/Dropbox/Chris/Very Final Stuff/Pipeline_CSV_files_for_Lucas/"

summary.edit.table = get.synth.table(paste0(input.dir,"Gene_specific_edit_files/"))

###########
# Model : #
###########

results.model = NULL


# Number of codons for km, kv, pl, sm & Ch
num.codons = sum(c(996,5230,8117,4607,4465))
# Number of edits
number.edition = sum(c(197,1098,1086,748,389))
# Averaged position bias
position.rate = c(0.4846503695, 0.3931210915, 0.1222285389)
position.rate[2] = position.rate[1]+position.rate[2]

for (num in 1:100){
  print(num)
  # Simulate the distribution of edits
  simulation = data.frame(number = 1:num.codons, first = 0, second = 0, third = 0)
  
  for (i in 1:number.edition){
    x = round((runif(1,0,num.codons)))
    y = runif(1,0,1)
    if (y <= position.rate[1]){
      simulation[x,'first']=1
    }else if (y > position.rate[1]  & y <= position.rate[2]){
      simulation[x,'second']=1
    }else{
      simulation[x,'third']=1
    }
  }
  
  # Check the distribution accross positions -> expected
  sum(simulation$first) ; sum(simulation$second) ; sum(simulation$third)
  
  # is the next position edited ?
  simulation$first.next = (simulation$second == 1 & simulation$first == 1)
  simulation$second.next = (simulation$third == 1 & simulation$second == 1)
  simulation$third.next = c(simulation[2:nrow(simulation),'first'],NA)
  simulation$third.next = (simulation$third.next == 1 & simulation$third == 1)
  
  # is the previous position edited ?
  simulation$first.prev = c(NA,simulation[1:(nrow(simulation)-1),'third'])
  simulation$first.prev = (simulation$first.prev == 1 & simulation$first == 1)
  simulation$second.prev = (simulation$first == 1 & simulation$second == 1)               
  simulation$third.prev = (simulation$second == 1 & simulation$third == 1)
  
  # Number of edit in the codon
  simulation$number.in.codon = simulation$first + simulation$second + simulation$third
  
  t.num.codon = table(subset(simulation, number.in.codon >=1)$number.in.codon)
  
  # Neighbour :
  simulation$first.neighbour = (simulation$first.prev == TRUE | simulation$first.next == TRUE)
  simulation$second.neighbour = (simulation$second.prev == TRUE | simulation$second.next == TRUE)
  simulation$third.neighbour = (simulation$third.prev == TRUE | simulation$third.next == TRUE)
  
  t.neighbour.1 = table(subset(simulation, first == 1)$first.neighbour)
  t.neighbour.2 = table(subset(simulation, second == 1)$second.neighbour)
  t.neighbour.3 = table(subset(simulation, third == 1)$third.neighbour)
  
  # Close :
  simulation$first.close = (simulation$first.neighbour == TRUE | simulation$number.in.codon >= 2)
  simulation$second.close = (simulation$second.neighbour == TRUE | simulation$number.in.codon >= 2)
  simulation$third.close = (simulation$third.neighbour == TRUE | simulation$number.in.codon >= 2)
  
  t.close.1 = table(subset(simulation, first == 1)$first.close)
  t.close.2 = table(subset(simulation, second == 1)$second.close)
  t.close.3 = table(subset(simulation, third == 1)$third.close)
  
  temp.df = data.frame(codon.1.edit = t.num.codon[1]/sum(t.num.codon), 
                       codon.2.edit = t.num.codon[2]/sum(t.num.codon), 
                       codon.3.edit = t.num.codon[3]/sum(t.num.codon),
                       first.no.neighbour = t.neighbour.1[1]/sum(t.neighbour.1), 
                       second.no.neighbour =  t.neighbour.2[1]/sum(t.neighbour.2), 
                       third.no.neighbour = t.neighbour.3[1]/sum(t.neighbour.3),
                       first.no.close = t.close.1[1]/sum(t.close.1), 
                       second.no.close = t.close.2[1]/sum(t.close.2), 
                       third.no.close = t.close.3[1]/sum(t.close.3))
  
  results.model = rbind(results.model, temp.df)
}

# Extracting statistics :
model.stats = results.model[1:2,]
model.stats[] <- NA
row.names(model.stats)=c('mean','sd')
for (n in names(model.stats)){
  model.stats[1,n] = mean(results.model[,n],na.rm=T)
  model.stats[2,n] = sd(results.model[,n],na.rm=T)
}

#################
# Observation : #
#################

summary.edit.table$num.edit.codon = NA
summary.edit.table$is.next = NA
for (i in 1:nrow(summary.edit.table)){
  summary.edit.table[i,'num.edit.codon'] = adist(summary.edit.table[i,'genome.codon'],
                                                 summary.edit.table[i, 'mRNA.codon'])
  if(i != 1 & i != nrow(summary.edit.table)){
    is.next = (summary.edit.table[(i+1),'position']==summary.edit.table[i,'position']+1) | (summary.edit.table[(i-1),'position']==summary.edit.table[i,'position']-1)
    summary.edit.table[i,'is.next']=is.next
  }
}

# Number of edits per codon :
t.num.codon = table(summary.edit.table$num.edit.codon)/sum(table(summary.edit.table$num.edit.codon))
t.num.codon ; model.stats[,1:3]

chisq.test(t.num.codon, as.numeric(model.stats[1,1:3]))
pnorm(as.numeric(t.num.codon[1]), mean = model.stats['mean',1], sd = model.stats['sd',1])
pnorm(as.numeric(t.num.codon[2]), mean = model.stats['mean',2], sd = model.stats['sd',2],lower.tail = F)
pnorm(as.numeric(t.num.codon[3]), mean = model.stats['mean',3], sd = model.stats['sd',3],lower.tail = F)

obs.first = subset(summary.edit.table, codon.position == 1)
table(paste0(obs.first$is.next,'+',obs.first$num.edit.codon))
#Direct neighbour :
t.obs.neighbour.1 = table(obs.first$is.next)/sum(table(obs.first$is.next))
t.obs.neighbour.1 ; model.stats[,4] ; print('mean & sd')
pnorm(as.numeric(t.obs.neighbour.1[1]), mean = model.stats['mean',4], sd = model.stats['sd',4])
#Close :
t.obs.close.1 = table((obs.first$is.next | obs.first$num.edit.codon >= 2))/sum(table((obs.first$is.next | obs.first$num.edit.codon >= 2)))
t.obs.close.1 ; model.stats[,7] ; print('mean & sd')
pnorm(as.numeric(t.obs.close.1[1]), mean = model.stats['mean',7], sd = model.stats['sd',7])

obs.second = subset(summary.edit.table, codon.position == 2)
table(paste0(obs.second$is.next,'+',obs.second$num.edit.codon))
#Direct neighbour :
t.obs.neighbour.2 = table(obs.second$is.next)/sum(table(obs.second$is.next))
t.obs.neighbour.2 ; model.stats[,5] ; print('mean & sd')
pnorm(as.numeric(t.obs.neighbour.2[1]), mean = model.stats['mean',5], sd = model.stats['sd',5])
#Close :
t.obs.close.2 = table((obs.second$is.next | obs.second$num.edit.codon >= 2))/sum(table((obs.second$is.next | obs.second$num.edit.codon >= 2)))
t.obs.close.2 ; model.stats[,8] ; print('mean & sd')
pnorm(as.numeric(t.obs.close.2[1]), mean = model.stats['mean',8], sd = model.stats['sd',8])

obs.third = subset(summary.edit.table, codon.position == 3)
table(paste0(obs.third$is.next,'+',obs.third$num.edit.codon))
#Direct neighbour :
t.obs.neighbour.3 = table(obs.third$is.next)/sum(table(obs.third$is.next))
t.obs.neighbour.3 ; model.stats[,6] ; print('mean & sd')
pnorm(as.numeric(t.obs.neighbour.3[1]), mean = model.stats['mean',6], sd = model.stats['sd',6])
#Close :
t.obs.close.3 = table((obs.third$is.next | obs.third$num.edit.codon >= 2))/sum(table((obs.third$is.next | obs.third$num.edit.codon >= 2)))
t.obs.close.3 ; model.stats[,9] ; print('mean & sd')
pnorm(as.numeric(t.obs.close.3[1]), mean = model.stats['mean',9], sd = model.stats['sd',9])

##########
# Figure #
##########

model = as.numeric(model.stats[1,1:3])
model.sd = as.numeric(model.stats[2,1:3])
observation = as.numeric(t.num.codon)

plot.1 = data.frame(number.per.codon = rep(c('1 edit','2 edits','3 edits'),2),
                    value = c(model,observation),
                    interv = c(model.sd,rep(NA,3)),
                    type = c(rep("model",3),rep("observation",3))
)

p1 = ggplot(plot.1, aes(x = number.per.codon, y = value, fill = type)) + geom_bar(stat='identity',position=position_dodge()) + theme_gray() + theme(legend.position = 'none') +
  geom_errorbar(aes(ymin=value-interv, ymax=value+interv), width=.2,position=position_dodge(.9)) + xlab("Number of edits per codon") +
  ylab("Percentage of edited codons")+ggtitle("Number of edits per codon : model vs observations") + ylim(0,1.1)
p1 = p1 + annotate("text", x = 1, y = 1.05,
                   label= paste0("p-value = ",pnorm(as.numeric(t.num.codon[1]), mean = model.stats['mean',1], sd = model.stats['sd',1])))
p1 = p1 + annotate("text", x = 2, y = 0.3,
                   label= paste0("p-value = ",pnorm(as.numeric(t.num.codon[2]), mean = model.stats['mean',2], sd = model.stats['sd',2],lower.tail=F)))
p1 = p1 + annotate("text", x = 3, y = 0.15,
                   label= paste0("p-value = ",pnorm(as.numeric(t.num.codon[3]), mean = model.stats['mean',3], sd = model.stats['sd',3],lower.tail=F)))



model = 1-as.numeric(model.stats[1,7:9])
model.sd = as.numeric(model.stats[2,7:9])
observation = 1-as.numeric(c(t.obs.close.1[1],t.obs.close.2[1],t.obs.close.3[1]))

plot.2 = data.frame(number.per.codon = rep(c('First','Second','Third'),2),
                    value = c(model,observation),
                    interv = c(model.sd,rep(NA,3)),
                    type = c(rep("model",3),rep("observation",3))
)

p2 = ggplot(plot.2, aes(x = number.per.codon, y = value, fill = type)) + geom_bar(stat='identity',position=position_dodge()) + theme_gray() +
  geom_errorbar(aes(ymin=value-interv, ymax=value+interv), width=.2,position=position_dodge(.9)) + xlab("Position in codon") + ylim(0,0.5)+
  ylab("Percentage of edits located close to other edits\n(adjacent or in same codon)") + ggtitle("Edit clustering : model vs observations")
p2 = p2 + annotate("text", x = 1, y = 0.27,
                   label= paste0("p-value = ",prettyNum(pnorm(as.numeric(t.obs.close.1[1]), mean = model.stats['mean',7], sd = model.stats['sd',7]),digits=2)))
p2 = p2 + annotate("text", x = 2, y = 0.3,
                   label= paste0("p-value = ",prettyNum(pnorm(as.numeric(t.obs.close.2[1]), mean = model.stats['mean',8], sd = model.stats['sd',8]),digits=2)))
p2 = p2 + annotate("text", x = 3, y = 0.43,
                   label= paste0("p-value = ",prettyNum(pnorm(as.numeric(t.obs.close.3[1]), mean = model.stats['mean',9], sd = model.stats['sd',9]),digits=2)))


pdf(file = "../RESULTS/9_position.bias.pdf", width=14, height=6)
plot_grid(p1,p2,nrow=1,labels="AUTO")
dev.off()
