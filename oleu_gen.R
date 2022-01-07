freeze<-ls() 
# read --------
#read metadata: should already be raw_metadata from oleu_base.R
#In case it is not:
# raw_metadata<-read_excel("onychomys_ct.xlsx")
#read the genepop file into the adegenet format, which by default stores alleles as integers [0, 1, 2]
#where 0 and 2 are homozygotes
x<-paste(datadir,gen_path_name,'/populations.snps.gen',sep="") %>%
  read.genepop(.)

metadata<-raw_metadata[raw_metadata$filename %in% indNames(x),]
#order to match genepop data
metadata<-metadata[match(indNames(x),metadata$filename),]


# create more readable labels for hierarchical groups -------
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
#add population data to genetic object
pop(x)<-metadata$population

#add both levels of hierarchy to genind object as strata
strata(x)<-data.frame(cbind(metadata$gengroup,metadata$population))
nameStrata(x)<- ~Genetic.Group/Sample.Area

#for each option below, the goal is to write the results

# calculate Fst for MMRR & Pst/Fst -------
source(paste(scriptsdir,"/oleu_gen_fst.R",sep=""))
# should write fst files.
# Interpolate missing data for PCA/PLS/individual genetic distances -------
source(paste(scriptsdir,"/oleu_gen_none_missing.R",sep=""))
#should write "snps.missing.interpolated.csv" file

rm(list = setdiff(ls(),freeze)) #clean up environment, also removes freeze