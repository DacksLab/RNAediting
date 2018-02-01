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
names.tab = data.frame(short=c("Kv","Km","Pl","Sm","Ch","Ht","Lp"),
                       long=c("K. veneficum",
                              "K. mikimotoi",
                              "P. lunula",
                              "S. minutum",
                              "C. horridum",
                              "H. triquetra",
                              "L. polyedrum"))
table.GC=read.delim(file="../DATA/pipeline_GC_position.csv",sep=',')
table.types=read.delim(file="../DATA/pipeline_types.csv",sep=',')
###########

###########
# Libraries
library(tidyverse)
library(cowplot)
###########

###########
# Functions
###########
get.increase = function(df){
  n = c("A.to.T", "A.to.G", "A.to.C", "T.to.A", "T.to.G", "T.to.C",
                   "G.to.A", "G.to.T", "G.to.C", "C.to.A", "C.to.T", "C.to.G")
  subset.names = names(df) %in% n
  return(mean(na.omit((df$A.to.C+df$A.to.G+df$T.to.C+df$T.to.G-
                         df$C.to.A-df$C.to.T-df$G.to.A-df$G.to.T)/
                        rowSums(df[,subset.names])*100)))
}

get.var = function(df){
  n = c("A.to.T", "A.to.G", "A.to.C", "T.to.A", "T.to.G", "T.to.C",
        "G.to.A", "G.to.T", "G.to.C", "C.to.A", "C.to.T", "C.to.G")
  subset.names = names(df) %in% n
  return(sqrt(var(na.omit((df$A.to.C+df$A.to.G+df$T.to.C+df$T.to.G-
                             df$C.to.A-df$C.to.T-df$G.to.A-df$G.to.T)/
                            rowSums(df[,subset.names])*100)))/
           sqrt(length(na.omit((df$A.to.C+df$A.to.G+df$T.to.C+df$T.to.G-
                                  df$C.to.A-df$C.to.T-df$G.to.A-df$G.to.T)/
                                 rowSums(df[,subset.names])*100))))
}

get.conf.int = function(df,sd.first,sd.second,sd.third){
  return(c(sd(df$X1st.GC.before)/sqrt(nrow(df)),
           sd(df$X1st.GC.after)/sqrt(nrow(df)),
           sd.first,
           sd(df$X2nd.GC.before)/sqrt(nrow(df)),
           sd(df$X2nd.GC.after)/sqrt(nrow(df)),
           sd.second,
           sd(df$X3rd.GC.before)/sqrt(nrow(df)),
           sd(df$X3rd.GC.after)/sqrt(nrow(df)),
           sd.third)*1.96)
}

get.pos.int = function(df,increase.first,increase.second,increase.third){
  return(c(mean(df$X1st.GC.before),mean(df$X1st.GC.after),increase.first,
           mean(df$X2nd.GC.before),mean(df$X2nd.GC.after),increase.second,
           mean(df$X3rd.GC.before),mean(df$X3rd.GC.after),increase.third))
}
###########

##############
# GC CONTENT #
##############
  
first=subset(table.types, position=="first")
increase.first=get.increase(first)
sd.first=get.var(first)

second=subset(table.types, position=="second")
increase.second=get.increase(second)
sd.second=get.var(second)

third=subset(table.types, position=="third")
increase.third=get.increase(third)
sd.third=get.var(third)

conf.int=get.conf.int(table.GC,sd.first,sd.second,sd.third)
pos.int=get.pos.int(table.GC,increase.first,increase.second,increase.third)

names = c("First, %GC before","First, %GC after","First, net GC increase",
          "Second, %GC before","Second, %GC after","Second, net GC increase",
          "Third, %GC before","Third, %GC after","Third, net GC increase")


plot.df = data.frame(Legend = names, mean = pos.int, interv = conf.int)

col = c("#4bb14b","#7fc97f","#b5e0b5",
        "#b14b4b","#c97f7f","#e0b5b5",
        "#4b4bb1","#7f7fc9","#b5b5e0")

p = ggplot(plot.df, aes(x=Legend, y = mean, fill = Legend)) + geom_bar(stat="identity") + 
  theme_grey() + theme(axis.text.x=element_text(angle = 45,hjust = 1),legend.position = 'none')+
  geom_errorbar(aes(ymin=mean-interv, ymax=mean+interv), width=.2,position=position_dodge(.9)) +
  ylab("%GC or\nNet % of edits increasing GC content")+xlab("")+ylim(-5,100)+scale_fill_manual(values=col)+
  ggtitle(paste0('Editing effect on the GC content per position'))

p = p + annotate("segment", x = 1, xend = 2, y = 51, yend = 51) + annotate("text", x = 1.5, y = 52,label="*")
p = p + annotate("segment", x = 4, xend = 5, y = 46, yend = 46) + annotate("text", x = 4.5, y = 47,label="*")
p = p + annotate("text", x = 2, y = 80,label="First Position",col="#4bb14b",fontface='bold') + 
  annotate("text", x = 5, y = 80,label="Second Position",col="#b14b4b",fontface='bold') + 
  annotate("text", x = 8, y = 80,label="Third Position",col="#4b4bb1",fontface='bold')

liste=c("Kv","Km","Pl","Sm","Ch","Ht","Lp")
for (i in liste){
  first=subset(table.types,organism==i & position=="first")
  increase.first=get.increase(first)
  sd.first=get.var(first)
  
  second=subset(table.types,organism==i & position=="second")
  increase.second=get.increase(second)
  sd.second=get.var(second)
  
  third=subset(table.types,organism==i & position=="third")
  increase.third=get.increase(third)
  sd.third=get.var(third)
  
  GC.sp = subset(table.GC, organism == i)
  
  conf.int=get.conf.int(GC.sp,sd.first,sd.second,sd.third)
  pos.int=get.pos.int(GC.sp,increase.first,increase.second,increase.third)
  
  plot.df = data.frame(Legend = names, mean = pos.int, interv = conf.int)
  
  sp=as.character(names.tab[which(names.tab$short==i),"long"])
  
  assign(paste0('p',i),ggplot(plot.df, aes(x=Legend, y = mean, fill = Legend)) + geom_bar(stat="identity") + 
           theme_grey() + theme(axis.text.x=element_text(angle = 45,hjust = 1),legend.position = 'none')+
           geom_errorbar(aes(ymin=mean-interv, ymax=mean+interv), width=.2,position=position_dodge(.9)) +
           ylab("%GC or\nNet % of edits increasing GC content")+xlab("")+scale_fill_manual(values=col)+
           coord_cartesian(ylim=c(ifelse(i=="Ch"|i=="Lp",-55,-5),100))+ggtitle(paste0('Editing effect on the GC content per position for ',sp)))
  
  if (i=="Kv"){
    pKv = pKv + annotate("segment", x = 1, xend = 2, y = 51, yend = 51) + annotate("text", x = 1.5, y = 52,label="*") +
      annotate("segment", x = 4, xend = 5, y = 46, yend = 46) + annotate("text", x = 4.5, y = 47,label="*")
  }
  if (i=="Pl"){
    pPl = pPl + annotate("segment", x = 1, xend = 2, y = 51, yend = 51) + annotate("text", x = 1.5, y = 52,label="*") +
      annotate("segment", x = 4, xend = 5, y = 51, yend = 51) + annotate("text", x = 4.5, y = 52,label="*")
  }
  if (i=="Ch"){
    pCh = pCh + annotate("segment", x = 4, xend = 5, y = 47, yend = 47) + annotate("text", x = 4.5, y = 48,label="*")
    }
}

p_plot = plot_grid(p,pKv,pKm,pPl,pSm,pHt,pCh,pLp,nrow=2,labels="AUTO",hjust=-3)
ggsave(filename = "../RESULTS/2_GC_barplot_review.pdf", p_plot, width=24, height=12)
