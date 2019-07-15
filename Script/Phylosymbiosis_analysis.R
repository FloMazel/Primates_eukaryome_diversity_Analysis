# Phlyosymbiosis for Primates 18S dataset 

rm(list=ls())
library(vegan)
library(dplyr)
library(ggplot2)
library(ape)
library(phytools)

source("/Users/fmazel/Desktop/Recherche/En_cours/MiscellanousRcode/Multiplot.R")
setwd("/Users/fmazel/Desktop/Recherche/En_cours/Analyses_en_cours/Primates/")
MyTheme=theme_bw()+theme(axis.text=element_text(size=9),axis.text.x =element_text(size=9),axis.title=element_text(size=14),legend.text=element_text(size=11),axis.ticks = element_line(size = 1.5)) #Define my ggplot theme


##################
#   Load data    #
##################

tree=read.tree("Data/MammalTrees_fix.tre")
Phydist=cophenetic(tree)

Metadata=read.table("Data/Metadata.txt",header = T,stringsAsFactors = F)
rownames(Metadata)=Metadata$SampleID
Beta=read.table("Data/weighted_unifrac_swarm_otus.wtax.final.txt")
Beta=as.matrix(Beta)
Beta=Beta[Metadata$SampleID,Metadata$SampleID]



####################################
#   Plot the general relationship  #
####################################

PhyDist_allsamples=Phydist[Metadata[row.names(Beta),"Taxa"],Metadata[row.names(Beta),"Taxa"]]
Data_Plot=data.frame(Phy_dist=c(PhyDist_allsamples),Beta=c(Beta))
Phylosymbiosis_plot=ggplot(data=Data_Plot,aes(y=Beta,x=Phy_dist))+
  geom_point()+MyTheme+
  xlab("Pairwise Host Phylogenetic distance (Millions years)")+ylab("Pairwise Beta-diversity (wUnifrac)")+
  ggtitle("A. Pairwise plots with all individuals")
Phylosymbiosis_plot



########################################################################
#          Measuting phylosymbiosns with the mantel test               #
########################################################################

x=100
results=list()
for (i in 1:x)
{
  #sub sample one individual per species
  sub_SP_Samples <- Metadata %>% group_by(Taxa) %>% sample_n(1) %>% as.data.frame()
  
  #subsample beta matrix 
  BetaSub=Beta[sub_SP_Samples$SampleID,sub_SP_Samples$SampleID]
  colnames(BetaSub)=rownames(BetaSub)=sub_SP_Samples$Taxa
  Phydist=Phydist[sub_SP_Samples$Taxa,sub_SP_Samples$Taxa]
  
  #Correlation 
  plot(Phydist,BetaSub)
  mantel_output=mantel(as.dist(BetaSub),as.dist(Phydist))
  mantel_outputS=mantel(as.dist(BetaSub),as.dist(Phydist),method = "spearman")
  results[[i]]=c(mantel_output$statistic,mantel_output$signif,mantel_outputS$statistic,mantel_outputS$signif)
  print(i)
}

Mantel_outputs_results=as.data.frame(do.call(rbind,results))
colnames(Mantel_outputs_results)=c("Statistic_Pearson","Significance_Pearson","Statistic_Spearman","Significance_Spearman")

########################################
#     Analyse and final plot results   #
########################################

#Not a single significnat correlation 
mean(Mantel_outputs_results$Statistic_Pearson)
table(Mantel_outputs_results$Significance_Pearson<.05) #no significant results
table(Mantel_outputs_results$Significance_Spearman<.05) #no significant results

#Plot the histograms of the statistics
Boxplot_Mantel=ggplot(Mantel_outputs_results,aes(y=Statistic_Pearson,x=NA))+
  geom_boxplot()+ylim(c(-.3,.3))+ylab("Mantel statistic (Pearson correlation)")+xlab("")+MyTheme+
  ggtitle("B. Mantel statistic distribution")
  
#General plot
pdf("Phylosymbiosis_plot.pdf",width = 13,height = 5)
multiplot(Phylosymbiosis_plot,Boxplot_Mantel,cols=2)
dev.off()






