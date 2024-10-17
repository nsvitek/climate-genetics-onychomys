#match barcode sequences to barcode names used in library prep, write to a barcode file for stacks. 

# Dependencies -----
library(dplyr)
library(readxl)
library(writexl)

# set file paths ----- 
path_in<-"C://Users/N.S/Dropbox/Documents/Dissertation/modern/peromyscus_maniculatus"
# path_in<-"D://Dropbox/Documents/Dissertation/modern/onychomys_leucogaster"
path_out<-paste(path_in,"radseq_analysis",sep="/")
path_scripts<-"C://cygwin/home/N.S/scripts/cataloging"

name_samples<-"pman_dna.xlsx"
name_populations<-name_samples
name_barcodes<-"Possible_Barcodes.xlsx"
name_out<-"pman_barcodes.txt"

# read, clean ---------
specimen_df<-paste(path_in,name_samples,sep="/") %>% read_excel(.)
population_df<-paste(path_in,name_populations,sep="/") %>% read_excel(.)
barcode_df<-paste(path_in,name_barcodes,sep="/") %>% read_excel(.)

#remove samples not in library
specimen_df<-specimen_df[-which(is.na(specimen_df$EcoR1_barcode)),]

#set up columns for use later
specimen_df$EcoR1_code<-NA
specimen_df$Mse1_code<-NA
specimen_df$population<-NA
# specimen_df$ID<-paste(specimen_df$institution,specimen_df$catalog,sep="-")
specimen_df$ID<-specimen_df$full_specimenb

# match -----

for (i in 1:nrow(specimen_df)){
  #fill forward code
  barmatch1<-which(barcode_df$`Barcode name`==specimen_df$EcoR1_barcode[i])
  specimen_df$EcoR1_code[i]<-barcode_df$Sequence[barmatch1]
  #fill reverse code
  barmatch2<-which(barcode_df$`Barcode name`==specimen_df$Mse1_barcode[i])
  specimen_df$Mse1_code[i]<-barcode_df$Sequence[barmatch2]
  #fill population
  popmatch<-which(population_df$full_specimenb==specimen_df$full_specimenb[i])
  specimen_df$population[i]<-population_df$population[popmatch]
}

# format, write -----
# want to end with a .txt file, no headers, where each line looks like, for example:
# CTCCGAAGAA<tab><sample-name, ex: apal-204-15-an>

# Combinatorial barcodes are specified, one per column, separated by a tab:
# CGATA<tab>ACGTA<tab>sample_01

#population file popmap.tsv should be <sample name><tab><population>
out_pop_df<-specimen_df %>% select(ID,population)

# According to adaptor_schematic_mcdaniel_rad_seq.txt and a check of the raw fastq.gz files,
#EcoRI cut site mutated to add a "C" to the beginning, so a second version of EcoRI barcodes has the C
specimen_df$EcoR1_code_mutated<-paste(specimen_df$EcoR1_code,"C",sep="")
specimen_df$Mse1_code_mutated<-paste(specimen_df$Mse1_code,"G",sep="")

#make different combinations of variables for possible barcode files
# out_df<-specimen_df %>% select(EcoR1_code,ID) 
# out_df2<-specimen_df %>% select(EcoR1_code,Mse1_code,ID) 
out_df_mutated<-specimen_df %>% select(EcoR1_code_mutated,Mse1_code_mutated,ID) #use this one
# out_df_mutated2<-specimen_df %>% select(EcoR1_code_mutated,ID)

#write those tables to tab-separated text files
# write.table(out_df, file = paste(path_out,name_out,sep="/"), row.names=FALSE, col.names=FALSE,
#             quote=FALSE, sep = "\t")
write.table(out_df_mutated, file = paste(path_out,"/",name_out,"_mutated.txt",sep=""), row.names=FALSE, col.names=FALSE,
            quote=FALSE, sep = "\t")
# write.table(out_df2, file = paste(path_out,"/",name_out,"PE.txt",sep=""), row.names=FALSE, col.names=FALSE,
#             quote=FALSE, sep = "\t")
# write.table(out_df_mutated2, file = paste(path_out,"/",name_out,"_mutatedSE.txt",sep=""), row.names=FALSE, col.names=FALSE,
#             quote=FALSE, sep = "\t")
write.table(out_pop_df, file = paste(path_out,"/","popmap.tsv",sep=""), row.names=FALSE, col.names=FALSE,
            quote=FALSE, sep = "\t")