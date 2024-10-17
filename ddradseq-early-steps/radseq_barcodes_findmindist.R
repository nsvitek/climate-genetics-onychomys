####
#took out_df from radseq_barcodes_make.R

library(ape)
library(stringr)

#okay, min barcode is 8 bp and all need to be same length, so just work with first 8

x<-paste(out_df$EcoR1_code,"NNNNNNNN",sep="") %>% str_trunc(.,width = 14,side="right")
y <- t(sapply(strsplit(x,""), tolower))
z <- as.DNAbin(y)
distmat<-dist.dna(z,model = "N") %>% as.matrix
summary(distmat)
which(distmat==1, arr.ind = TRUE)
out_df$EcoR1_code[which(distmat==2, arr.ind = TRUE)[,1]]

seqlen<-str_length(out_df$EcoR1_code)
y <- t(sapply(strsplit(out_df$EcoR1_code[which(seqlen==8)],""), tolower))
z <- as.DNAbin(y)
distmat<-dist.dna(z,model = "N")
summary(distmat) #min = 4

y <- t(sapply(strsplit(out_df$EcoR1_code[which(seqlen==9)],""), tolower))
z <- as.DNAbin(y)
distmat<-dist.dna(z,model = "N")
summary(distmat) #min = 4

y <- t(sapply(strsplit(out_df$EcoR1_code[which(seqlen==10)],""), tolower))
z <- as.DNAbin(y)
distmat<-dist.dna(z,model = "N")
summary(distmat) #min = 4

y <- t(sapply(strsplit(out_df$EcoR1_code[which(seqlen==14)],""), tolower))
z <- as.DNAbin(y)
distmat<-dist.dna(z,model = "N")
summary(distmat) #min = 5

