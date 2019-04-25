################################################################################################################
##############  Correlation between 16S and 18S alpha diversity taking into account host phylogeny   ###########
################################################################################################################

library(MCMCglmm)
library(caper)
library(phytools)
library(dplyr)

setwd("/Users/fmazel/Desktop/Recherche/En_cours/Analyses_en_cours/Primates/PGLS/") #to be updtae with the publication of the data


###### Loa data
tree=read.tree("Data/MammalTrees_fix.tre")
tree=force.ultrametric(tree)
Alpha_div=read.table("Data/adiv_forFlo.txt",header = T)


####### Checking the distribution of the data
hist(Alpha_div$observed_otus_bac)
hist(Alpha_div$observed_otus_euks) # strongly left skewed, we should prpbably use a poisson distribution 


##### PLotting the correlation 
dev.off()
ggplot(aes(x=observed_otus_bac,y=observed_otus_euks,colour=Genus_species),data=Alpha_div)+geom_point()
ggplot(aes(x=observed_otus_bac,y=log(observed_otus_euks),colour=Genus_species),data=Alpha_div)+geom_point()
plot(tree)

##### MCMCglmm

#Build the model
inv.phylo<-inverseA(tree,nodes="TIPS",scale=TRUE) # Create a matrix of similarity between species ("pedigree")
Ginv=list(Genus_species=inv.phylo$Ainv) # Format it to use in MCMCglmm function

#Run the Poisson model without prior (with other uninformative prior, we got  similar results)
model1_Poisson<-MCMCglmm(fixed=observed_otus_euks~observed_otus_bac,random=~Genus_species,family="poisson",ginverse=Ginv,data=Alpha_div,nitt=1000000,burnin=10000,thin=1000) #run the mode

#Check that model converged 
par(mfrow=c(2,2))
plot(model1_Poisson$Sol, auto.layout=F)
slopeMCMCGlmm=ggplot(data.frame(model1_Poisson$Sol), aes(x=observed_otus_bac))+geom_histogram()+ labs(title="MCMCglmm Slope estimate")

# Outputs of the models 
summary_Poisson=summary.MCMCglmm(model1_Poisson) #paramters estimate and significance
lambdaPoisson <- model_Poisson$VCV[,'Genus_species']/(model_Poisson$VCV[,'Genus_species']+model_Poisson$VCV[,'units'])
median(lambdaPoisson) # estimation of phylo. signal (~0.71)
lambdaMCMCGlmm=ggplot(data.frame(lambdaPoisson), aes(x=lambdaPoisson))+geom_histogram()+ labs(title="MCMCglmm Lambda estimate")

#save the outputs: 
summary_Poisson=summary(model1_Poisson) 
write.csv(summary_model$solutions,file="/Results/MCMCglmm_PoissonFamily.csv")
