# !diagnostics off
#above added to try to stop getting warnings message about unknown factors
# source dependencies -----
#closer to the end, source this section so that you can switch between cluster & personal computer
#set file locations
datadir <- "/gpfs/scratch/nvitek/onychomys_leucogaster"
# datadir <- "D:/Dropbox/Documents/Dissertation/modern/onychomys_leucogaster"

gen_path_name<-paste('/radseq_analysis/populations.r',mdat,sep="") #for reading and writing

setwd(datadir)

library(dplyr,lib.loc="/gpfs/home/nvitek/R_packages")
library(ggplot2,lib.loc="/gpfs/home/nvitek/R_packages")
library(ggthemes,lib.loc="/gpfs/home/nvitek/R_packages")
library(scales,lib.loc="/gpfs/home/nvitek/R_packages") 
library(geosphere,lib.loc="/gpfs/home/nvitek/R_packages") #for distm: geographic distance
library(RRPP,lib.loc="/gpfs/home/nvitek/R_packages") 
library(geomorph,lib.loc="/gpfs/home/nvitek/R_packages") 
library(reshape2,lib.loc="/gpfs/home/nvitek/R_packages")
library(inlmisc,lib.loc="/gpfs/home/nvitek/R_packages") #for bright Paul Tol color scheme

source("oleu_pst_fst_functions.R") #mostly Pst functions
source("MMRR.R") #MMRR functions

#number of resamples
resample<-1000
# read in raw metadata ------
metadata<-read.csv("metadata.csv")
# create more readable labels for plotting -------
#sampling area labels
# two areas, one in OK and one in NE, cross multiple counties and are labelled differently
metadata$`Sampling Area`<-metadata$population
metadata$`Sampling Area`[which(metadata$population!="ip_ne"&
                                 metadata$population!="ip_ok")]<-metadata$county[which(metadata$population!="ip_ne"&
                                                                                         metadata$population!="ip_ok")]
metadata$`Sampling Area`[which(metadata$population=="ip_ne"|
                                 metadata$population=="ip_ok")]<-metadata$stateprovince[which(metadata$population=="ip_ne"|
                                                                                                metadata$population=="ip_ok")]

#genetic group labels
metadata$`Genetic Group`<-metadata$gengroup
metadata$`Genetic Group`[which(metadata$gengroup=="cp")]<-"CPBR"
metadata$`Genetic Group`[which(metadata$gengroup=="ip")]<-"IP"
metadata$`Genetic Group`[which(metadata$gengroup=="wb")]<-"WB"

# read in genotype ------
#if running output of oleu_gen.R, comment out the above lines and run lines below instead.
snps<-read.csv(paste(datadir,gen_path_name,"/snps.missing.interpolated.mean.csv",sep=""),row.names=1)
PCA.gen<-read.csv(paste(datadir,"/PCA_gen_mdat_",mdat,".csv", sep=""),row.names=1)
stop.gen<-ncol(PCA.gen)
#get global Fst
raw.mean.fst.wc<-read.table(paste(datadir,gen_path_name,"/fst.overall.tsv",sep=""),
                            sep=" ",header=TRUE)
# get pairwise Fst
pop.combo.wc<-read.table(paste(datadir,gen_path_name,"/pairwise.fst.tsv",sep=""),sep=" ",header=TRUE)

if("MSB-66111" %in% rownames(snps)){
  snps<-snps[-which(rownames(snps)=="MSB-66111"),]
}

# read environment -----
#read in raw environmental data per specimen
env.block<-read.csv("env.block.csv",row.names=1)
geo.block<-read.csv("geo.block.csv",row.names=1)

# MULTIPLE: decide on phenotype block -----
#phenotype has options
if (phenotype.option == "points") {
  outputdir.top<-paste(datadir,"/output/shape_point",sep="")
  mor.block<-read.csv(paste(outputdir.top,"/mor.block.csv",sep=""),row.names=1)
}

if (phenotype.option == "PCs")    {
  outputdir.top<-paste(datadir,"/output/shape_pc",sep="")
  mor.block<-read.csv(paste(outputdir.top,"/mor.block.csv",sep=""),row.names=1)
}

print("Data read into R.")
# genetic data: block ---------- 
#control for missing data level
outputdir.mid<-paste(outputdir.top,"/r.",mdat,sep="")
#this will change with multiple datasets
outputdir.choice<-outputdir.mid

# Pst-Fst --------------
outputdir.choice<-outputdir.mid

#Pst variables
# h2.levels<-seq(from=0.1, to=0.9, by=0.1)
# h.squared<-h2.levels[5] #set at 0.5 for code development
# add.gen.proportion<-0.25 #c, assumed additive genetic proportion of differences between populations
#the idea from Brommer 2011 is that any time c > h.squared your estimate is anti-conservative, especially if Pst>Fst
#so you might want to start with the null, c = h.squared
#but also see what happens when you alter both sides of the assumption -- plot as a surface? 
#try with the starting settings chosen above

print("started the Pst-Fst loop.")
source(paste(datadir,"/oleu_pst_fst.R",sep=""))

print("completed the Pst-Fst loops.")
# MMRR: pairwise distances between individuals-----
#see note below, OMNH-52913, MSB-66182 need removal: too much missing genetic data
sadness<-c("OMNH-52913","MSB-66182")
# mor.block<-mor.block[-which(row.names(mor.block) %in% sadness),]
# env.block<-env.block[-which(row.names(env.block) %in% sadness),]
# geo.block<-geo.block[-which(row.names(geo.block) %in% sadness),]
# metadata<-metadata[-which(metadata$filename %in% sadness),]
snps<-snps[-which(rownames(snps) %in% sadness),]
# genetic block of data
gen.block<-snps

#make fst.pairwise.wc
pops<-unique(c(levels(factor(pop.combo.wc$POP1)),levels(factor(pop.combo.wc$POP2))))
n.pops<-length(pops)
fst.pairwise.wc <- array(0, c(n.pops, n.pops), list(pops, pops))
i <- match(pop.combo.wc$POP1, pops)
j <- match(pop.combo.wc$POP2, pops)
fst.pairwise.wc[cbind(i,j)] <- fst.pairwise.wc[cbind(j,i)] <- pop.combo.wc$FST

print("Started the distance analyses.")
source(paste(datadir,"/oleu_distances.R",sep=""))
sink()

# 2B-PLS -------
#2-block partial least squares. Doesn't rely on different values of h.squared or c, only mdat
source(paste(datadir,"/oleu_pls.R",sep=""))
sink()

# try again without Wyoming Basin --------
outputdir.deep<-paste(outputdir.mid,"/no_WB",sep="")
if (file.exists(outputdir.deep)){
  print("no_wb directory exists")
  } else {
  dir.create(file.path(outputdir.deep))
}
#this will change with multiple datasets
outputdir.choice<-outputdir.deep
#save original data blocks, then make slimmed ones
mor.block.orig<-mor.block[,]
gen.block.orig<-gen.block[,]
env.block.orig<-env.block[,]
geo.block.orig<-geo.block[,]
metadata.orig<-metadata[,]
PCA.gen.orig<-PCA.gen
pop.combo.wc.orig<-pop.combo.wc

wyoming<-which(metadata$gengroup=="wb")
mor.block.2<-mor.block[-wyoming,]
gen.block.2<-gen.block[-wyoming,]
env.block.2<-env.block[-wyoming,]
geo.block.2<-geo.block[-wyoming,]
metadata.2<-metadata[-wyoming,]
PCA.gen.2<-PCA.gen[-wyoming,]

wyoming.pops<-c("wb_carbon","wb_sweetwater")
which.wyoming.pops<-which(pop.combo.wc$POP1 %in% wyoming.pops|pop.combo.wc$POP2 %in% wyoming.pops) %>% print()
pop.combo.wc.2<-pop.combo.wc[-which.wyoming.pops,]
  
mor.block<-mor.block.2
gen.block<-gen.block.2
env.block<-env.block.2
geo.block<-geo.block.2
metadata<-metadata.2
PCA.gen<-PCA.gen.2
pop.combo.wc<-pop.combo.wc.2

#make fst.pairwise.wc
pops<-unique(c(levels(factor(pop.combo.wc$POP1)),levels(factor(pop.combo.wc$POP2))))
n.pops<-length(pops)
fst.pairwise.wc <- array(0, c(n.pops, n.pops), list(pops, pops))
i <- match(pop.combo.wc$POP1, pops)
j <- match(pop.combo.wc$POP2, pops)
fst.pairwise.wc[cbind(i,j)] <- fst.pairwise.wc[cbind(j,i)] <- pop.combo.wc$FST

source(paste(datadir,"/oleu_distances.R",sep=""))
sink() 
print("If you're seeing this, then things are tripping up *after* No Wyoming distances")
source(paste(datadir,"/oleu_pls.R",sep=""))

sink()
print("Looks like we made it without Wyoming")

# just the 3 populations ------
outputdir.deep<-paste(outputdir.mid,"/key3",sep="")
if (file.exists(outputdir.deep)){
  print("key 3 directory exists")
} else {
  dir.create(file.path(outputdir.deep))
}
#this will change with multiple datasets
outputdir.choice<-outputdir.deep
focal.3<-which(metadata$population=="cp_otero"|metadata$population=="ip_roosevelt"|
                 metadata$population=="ip_winkler")
mor.block.2<-mor.block[focal.3,]
gen.block.2<-gen.block[focal.3,]
env.block.2<-env.block[focal.3,]
geo.block.2<-geo.block[focal.3,]
metadata.2<-metadata[focal.3,]
PCA.gen.2<-PCA.gen[focal.3,]

focal.3.pops<-c("cp_otero","ip_roosevelt","ip_winkler")
which.focal.3.pops<-which(pop.combo.wc$POP1 %in% focal.3.pops&pop.combo.wc$POP2 %in% focal.3.pops)
pop.combo.wc.2<-pop.combo.wc[which.focal.3.pops,]

mor.block<-mor.block.2
gen.block<-gen.block.2
env.block<-env.block.2
geo.block<-geo.block.2
metadata<-metadata.2
PCA.gen<-PCA.gen.2
pop.combo.wc<-pop.combo.wc.2

source(paste(datadir,"/oleu_pls.R",sep=""))
sink()