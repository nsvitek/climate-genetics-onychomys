resample<-1000
missing.data.levels<-seq(from=20,to=80,by=10)
mdat<-missing.data.levels[6]

#Phenotype as pseudolandmarks or PCs? Run one
phenotype.option <- "PCs"
source("oleu_cluster_base.R")