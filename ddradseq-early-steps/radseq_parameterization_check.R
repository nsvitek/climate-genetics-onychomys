#From Rochette and Catchen 2017
#plus code by Yihan Wu on R-Bloggers: https://www.r-bloggers.com/parameter-testing-in-stacks-snp-extraction-and-visualization-in-r/
library(ggplot2)

path_in<-"D://RADseq_PM"
setwd(path_in)

# make a data frame combining all parameter distributions 
count <- 1
files <- list.files(pattern="*_snp_distribution.tsv", 
                    full.names = T)
for (i in files){
  table <- read.delim(i, skip=1, header=T)
  table$n_loci_percent<- table$n_loci/sum(table$n_loci)
  table$m<- count
  write.table(table, "distributions.tsv", append=T, row.names=F, col.names = F)
  snp_count <- data.frame("m"= count, "n_snps"=sum(table$n_loci))
  write.table(snp_count, "total_count.tsv", append=T, row.names=F, col.names = F)
  count <- count + 1
}
# total_count.tsv is used to display total SNP by parameter.

# visualize --------
#snp count plot
snp_count<-read.delim("total_count.tsv", sep=" ", header=F)
names(snp_count)<-c("m", "n_snps")
snp_count$m<-as.factor(snp_count$m)
pdf("paramaterize_m_snps.pdf",width=4,height=4)
ggplot(data=snp_count, aes(x=m, y=n_snps)) +
  geom_point() + theme_classic() #+ scale_y_continuous(limits = c(0, 3000))
dev.off()
#
snp_table<-read.delim("distributions.tsv", sep=" ", header=F)
names(snp_table)<- c("n_snps","n_loci", "n_loci_percent", "m") 
snp_table$n_loci_percent<-snp_table$n_loci_percent*100
snp_table$n_snps<-ifelse(snp_table$n_snps < 9, snp_table$n_snps, "9 +")
snp_table$n_snps<-as.factor(snp_table$n_snps)
snp_table$m<-as.factor(snp_table$m)

pdf("parameterize_snps_loci_percent.pdf",width=4,height=4)
ggplot(data = snp_table) + 
  geom_col(aes(x=n_snps, y=n_loci_percent, fill=m), position="dodge") + 
  theme_classic() + scale_fill_brewer(palette="RdBu")
dev.off()
