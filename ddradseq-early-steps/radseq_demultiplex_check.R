#check results of demultiplex step, see if any specimen has too few retained reads.

# Dependencies -----
library(dplyr)

# set file paths ----- 
# path_in<-"C://Users/N.S/Dropbox/Documents/Dissertation/modern/onychomys_leucogaster/radseq_analysis/02-demultiplex"
path_in<-"D://Dropbox/Documents/Dissertation/modern/onychomys_leucogaster/radseq_analysis/02-demultiplex"

name_samples<-"retained_reads.txt"

# read, clean ---------
specimen_df<-paste(path_in,name_samples,sep="/") %>% 
  read.table(., header = TRUE, sep = "\t")

#look at % of reads retained for each specimen
specimen_df<-specimen_df %>% mutate(.,percent = Retained/Total) 

specimen_df$Retained %>% hist()
specimen_df$percent %>% hist()

#look at which specimens have a very low % retained reads
specimen_df[which(specimen_df$percent<0.4),] #still at least 400k reads.

#order the specimens by # reads retained
specimen_df<-specimen_df[order(specimen_df$Retained),]

#look at specimens with fewest retained reads
specimen_df$Retained[1:20] %>% hist()
specimen_df$Filename[1:5]


#look at lowest number of reads
specimen_df[which(specimen_df$Filename=="MSB-140226"),]
specimen_df[which(specimen_df$Filename=="MSB-122896"),]


# reporting -----
#look at average number of reads per specimen
specimen_df$Total %>% mean
specimen_df$Total %>% range

specimen_df$Retained %>% mean
specimen_df$Retained %>% range

#remove the specimens that got cut:

got_cut<-c("MSB-140226","MSB-88629","OMNH-52913", "MSB-66182","MSB-66111")

specimen_report<-specimen_df[-which(specimen_df$Filename %in% got_cut),]
specimen_report$Retained %>% mean
specimen_report$Retained %>% range
