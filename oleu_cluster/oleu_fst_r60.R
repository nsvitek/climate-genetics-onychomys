resample<-1000
missing.data.levels<-seq(from=20,to=80,by=10)
mdat<-missing.data.levels[5]

#Phenotype as pseudolandmarks or PCs? Run one
phenotype.option <- "points"
source("oleu_cluster_base.R")