# !diagnostics off
#above added to try to stop getting warnings message about unknown factors
# source dependencies -----
#closer to the end, source this section so that you can switch between cluster & personal computer
#set file locations
# scriptsdir <- "C://cygwin/home/N.S/scripts/scripts"
scriptsdir <- "C://scripts/climate-genetics-onychomys"
datadir <- "D:/Dropbox/Documents/Dissertation/modern/onychomys_leucogaster"
# datadir <- "C:/Users/N.S/Dropbox/Documents/Dissertation/modern/onychomys_leucogaster"

# source(paste(scriptsdir,"/erinaceomorph_dependencies.R",sep=""))
library(readxl)
library(dplyr)
library(ggplot2)
library(ggthemes) #also ptol
library(reshape2) #used in error eval
library(geomorph) #reading, formatting, analyzing morphometric data, two.b.pls
library(Rvcg) #importing/exporting 3d files
library(molaR) #not just dental calcs, but also cSize
library(scales) #for rescale for colors, not cutting anymore
library(inlmisc) #for colors in heat map, error
library(adegenet) #for SNP data, but used in previously sourced code
library(hierfstat) #for SNP data. If not running oleu_gen again can possibly skip. 
library(geosphere) #for distm: geographic distance
# library(RColorBrewer)
# library(Morpho) #more two-block partial least squares

setwd(scriptsdir)
source("../observer-free-morphotype-characterization/rearrange3ddat.R")
source(paste(scriptsdir,"/oleu_pst_fst_functions.R",sep="")) #mostly Pst functions
source(paste(scriptsdir,"/MMRR.R",sep="")) #MMRR functions

# #if visualizing/calculating error or visualizing average, predicted, etc. shapes
# source("../../observer-free-morphotype-characterization/figpiece3d.R")
# source("../../observer-free-morphotype-characterization/find_repeatablePCs.R")
# source("../../scripts/calculate_error.R")

ncut<-256
raintable<-GetColors(ncut,scheme="sunset") #set gradient color scheme
ramp<-GetColors(, scheme="sunset")

# genetic variables ------
resample<-1000
missing.data.levels<-seq(from=40,to=80,by=10)
mdat<-missing.data.levels[4]
gen_path_name<-paste('/radseq_analysis/populations.r',mdat,sep="") #for reading and writing

# phenotype option -------
#Phenotype as pseudolandmarks or PCs?
phenotype.option <- "points"
# phenotype.option <- "PCs"

# read in raw metadata ------
setwd(paste(datadir,"/output",sep=""))
raw_metadata<-read_excel(paste(datadir,"onychomys_ct.xlsx",sep="/"))
raw_metadata$gengroup<-substr(raw_metadata$population,1,2)

# read in genotype ------
#read the output of oleu_gen.R from cluster, as opposed to executing oleu_gen script itself
#turns out that R can't read a csv that large without crashing, can go down to 40% missing data?

# #if running genetic data for first time, source oleu_gen.R:
# for (levelx in 1:5){
#   mdat<-missing.data.levels[levelx]
#   gen_path_name<-paste('/radseq_analysis/populations.r',mdat,sep="") #for reading and writing
#   source(paste(scriptsdir,"/oleu_gen.R",sep=""))
# }

#if running output of oleu_gen.R, comment out the above lines and run lines below instead.
snps<-read.csv(paste(datadir,gen_path_name,"/snps.missing.interpolated.mean.csv",sep=""),row.names=1)

if("MSB-66111" %in% rownames(snps)){
  snps<-snps[-which(rownames(snps)=="MSB-66111"),]
}

metadata<-raw_metadata[raw_metadata$filename %in% rownames(snps),]
#order to match genetic data
metadata<-metadata[match(rownames(snps),metadata$filename),]
#get global Fst
raw.mean.fst.wc<-read.table(paste(datadir,gen_path_name,"/fst.overall.tsv",sep=""),
                            sep=" ",header=TRUE)
# get pairwise Fst
pop.combo.wc<-read.table(paste(datadir,gen_path_name,"/pairwise.fst.tsv",sep=""),sep=" ",header=TRUE)

#make fst.pairwise.wc
pops<-unique(c(levels(factor(pop.combo.wc$POP1)),levels(factor(pop.combo.wc$POP2))))
n.pops<-length(pops)
fst.pairwise.wc <- array(0, c(n.pops, n.pops), list(pops, pops))
i <- match(pop.combo.wc$POP1, pops)
j <- match(pop.combo.wc$POP2, pops)
fst.pairwise.wc[cbind(i,j)] <- fst.pairwise.wc[cbind(j,i)] <- pop.combo.wc$FST

fst.pairwise.print<-fst.pairwise.wc
fst.pairwise.print[cbind(i,j)] <- round(pop.combo.wc$FST,3)
fst.pairwise.print[cbind(j,i)] <- paste(round(pop.combo.wc$CI_LOW,3),"-", round(pop.combo.wc$CI_HIGH,3),sep="")
matrix(fst.pairwise.print, nrow=n.pops, ncol=n.pops,dimnames=list(pops,pops)) %>% 
  write.csv(.,paste(datadir,gen_path_name,"/pairwise.fst.report.csv",sep=""))

#see note below, OMNH-52913, MSB-66182 need removal: too much missing genetic data
#OMNH-52913, MSB-66182 are essentially at (0,0), so all missing data?
sadness<-which(metadata$filename=="OMNH-52913"|metadata$filename=="MSB-66182")
metadata<-metadata[-sadness,]
snps<-snps[-which(! rownames(snps) %in% metadata$filename),]

# read environment and phenotype -----
#read and format phenotypic data, especially shapes
source(paste(scriptsdir,"/oleu_phenotype.R",sep=""))

#read in raw environmental data per specimen
enviro.raw<-read.csv(paste(datadir,"oleu_environment.csv",sep="/"))
# standardize among datasets -----
#standardize to the order of the aligned shapes
phenotype<-shapes[match(rownames(snps),rownames(shapes)),]
rm("shapes")

#reorganize rows of environmental data to match other tables
enviro.filtered<-enviro.raw[match(rownames(snps),enviro.raw$specimen),]
enviro.all<-enviro.filtered[match(metadata$filename,enviro.filtered$specimen),]
row.names(enviro.all)<-enviro.all$specimen

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

# # Sampling Reporting -------
# source(paste(scriptsdir,"/oleu_map.R",sep=""))

# genetic PCA -------

nrow(metadata)
nrow(snps)
#genetic Principal Components Analysis
PCA.genotype <- dudi.pca(snps,center=TRUE,scale=FALSE,nf=79,scannf=FALSE)

# #look for specimens with possibly too much missing data
# which(PCA.genotype$li[,1] < 0.01 & PCA.genotype$li[,1] > -0.01) %>% metadata$filename[.] %>% print()   

percents.gen<-(PCA.genotype$eig/sum(PCA.genotype$eig)*100) %>% round(1)
# sum(percents.gen[1:47])
for(i in 1:79){
  if(sum(percents.gen[1:i]) > 95){
    print(paste("PCs 1-",i," account for at least 95% of variation", sep=""))
    stop.gen<-i
    break()
  }
} 

master.variables<-cbind(metadata,PCA.genotype$li[,c(1:stop.gen)])
#Plot SNP PCA
# ggplot(data=master.variables,aes(x=Axis1,y=Axis2, color = `Sampling Area`,shape=`Genetic Group`))+
#   geom_point(size=2) +
#   theme_classic() + theme(legend.position="none") +
#   xlab(paste("PC 1 (",percents.gen[1],"%)",sep="")) +
#   ylab(paste("PC 2 (",percents.gen[2],"%)",sep=""))
# ggsave(paste("PCA_gen_mdat_",mdat,"_gen.pdf",sep=""), device = cairo_pdf, width = 8, height = 8,units="cm",dpi=600)
# 
# write.csv(PCA.genotype$li[,1:stop.gen],paste("PCA_gen_mdat_",mdat,".csv", sep=""))

# create phenotypic, environmental and geo block of data, explore structure ------------
source(paste(scriptsdir,"/oleu_environment.R",sep=""))

# MULTIPLE: decide on phenotype block -----
#env.block and geo.block already exist from oleu_environment
#write blocks to files
write.csv(env.block,"env.block.csv",quote=FALSE)
write.csv(geo.block,"geo.block.csv",quote=FALSE)
# write.csv(geo.dist, "geo.dist.csv", quote=FALSE)
write.csv(metadata,"metadata.csv")

#phenotype has options
if (phenotype.option == "points") {
  mor.block<-phenotype
  setwd(paste(datadir,"/output/shape_point",sep=""))
  #format: make sure row names are in there
  rownames(mor.block)<-metadata$filename
  #write these blocks to files
  write.csv(mor.block,"mor.block.csv",quote=FALSE)
}

if (phenotype.option == "PCs")    {
  mor.block<-PCA.phenotype$x[,c(1:12)]
  setwd(paste(datadir,"/output/shape_pc",sep=""))
  #format: make sure row names are in there
  rownames(mor.block)<-metadata$filename
  #write these blocks to files
  write.csv(mor.block,"mor.block.csv",quote=FALSE)
  
}
# setwd("../")
# #run disparity comparisons for each shape type
# source(paste(scriptsdir,"/oleu_disparity.R",sep=""))

# MULTIPLE: genetic data: block ---------- 
#control for missing data level
setwd(paste("r",mdat,sep="."))

# genetic block of data
gen.block<-snps

#At this point, your analyses should probably be done on the cluster. see oleu_cluster_base as a starting point
# troubleshooting: temporary: hierarchical ANOVA ------
dim(mor.block)
metadata %>% colnames
metadata$`Sampling Area`
metadata$`Genetic Group`
# Trying a hierarchical ANOVA as recommended by reviewer for submitted L305-307
#start with a geomorph data frame
gdf<-geomorph.data.frame(coords=arrayspecs(mor.block,1024,3),
                         SA=factor(metadata$`Sampling Area`),
                         GG=factor(metadata$`Genetic Group`))

gdf<-geomorph.data.frame(coords=as.matrix(PCA.genotype$li),
                         SA=factor(metadata$`Sampling Area`),
                         GG=factor(metadata$`Genetic Group`))

# dimnames(gdf$coords)<-list(
#   as.character(seq(1:1024)),
#   c("X","Y","Z"),
#   metadata$filename
# )


fit.gdf<-procD.lm(coords ~ GG/SA, data = gdf, print.progress = FALSE)

summary(fit.gdf)

anova(fit.gdf) #this doesn't seem right
anova(fit.gdf, error = c("GG:SA", "Residuals")) #example of correct adjustment

(76 * 4638) / (76 * 4639 + 10 * 18926)
(76 * 0.10968) / (76 * 0.10968 + 10 * 0.53920)

#may also be worth comparing models?
fit.gdf2<-procD.lm(coords ~ GG, data = gdf, print.progress = FALSE) #GG works, not SA. Doing something conceptually wrong

anova(fit.gdf2, fit.gdf) #adding info about sampling area gives significantly improved fit to model with physiogeographic area
#can't test the other way around...because one is completely redundant with the other?

fit.gdf3<-procD.lm(coords ~ GG + SA, data = gdf, print.progress = FALSE) #GG works, not SA. Doing something conceptually wrong
anova(fit.gdf3) #see example below of correct adjustment

#known working example:
# data("larvalMorph")
# Y.gpa <- gpagen(larvalMorph$tailcoords, 
#                 curves = larvalMorph$tail.sliders,
#                 ProcD = TRUE, print.progress = FALSE)
# gdf <- geomorph.data.frame(Y.gpa, treatment = larvalMorph$treatment, 
#                            family = larvalMorph$family)
# 
# fit <- procD.lm(coords ~ treatment/family, data = gdf, turbo = TRUE,
#                 print.progress = FALSE, iter = 999)
# anova(fit) # treatment effect not adjusted
# anova(fit, error = c("treatment:family", "Residuals"))


# # troubleshooting 2: hierarchical Fst (genetic) WORKING ------
# #also required lines 9-20 of oleu_gen.R
# library(hierfstat)
# x.matched<-x[indNames(x) %in% metadata$filename]
# x.df<-genind2hierfstat(x.matched,pop=x.matched@pop)
# 
# hierarchical.fst<-varcomp.glob(data.frame(metadata$`Genetic Group`,metadata$`Sampling Area`),x.df[,-1])



# troubleshooting 3: Pst --------------------------------
#choose how you want to measure "morphology"
sample.phenotype<-PCA.phenotype$x
sample.populations<-metadata$population

#get numbers to calculate the two variances used in Pst
pops.pd<-partition.distance(sample.phenotype,sample.populations,unique(sample.populations))
pops.pd<-partition.distance(as.matrix(PCA.genotype$li),sample.populations,unique(sample.populations))

#note: these sum of squares get pretty close to the geomorph procD.lm sum of squares which seems like a good sign.

#with a jackknife (leave out one specimen, recalculate)
# Numerator: (n-k-1)SSB*
# Denominator: (k-1)(SSW*+SSB*)
#SSB: sum of squared distance between group means and the grand mean
#SSW: sum of squared distances between individuals and their group mean
#k = # groups
#n = sample size

pst.polly<-function(SSB, SSW, k, n) {
  numerator<-(n - k - 1) * SSB
  denominator<- (k - 1) * (SSB + SSW)
  result<-numerator/denominator
  return(result)
}

pst.polly(pops.pd$SSB,pops.pd$SSW, pops.pd$k, pops.pd$n)

gg.pd<-partition.distance(sample.phenotype,metadata$`Genetic Group`,unique(metadata$`Genetic Group`))
pst.polly(gg.pd$SSB,gg.pd$SSW, gg.pd$k, gg.pd$n)

# individuals within sampling areas
pst.polly(0.10968,0.53920, 11, 77)
pst.polly(0.10968,0.53920, 11, 7)
#sampling areas within genetic groups
pst.polly(0.05005,0.10968, 3, 11)
pst.polly(0.05005,0.10968, 3, 3.7)
