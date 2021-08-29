#calculate Fst with bootstrapped confidence interval
#relies on objects and libraries in oleu_gen.R

# # weir-cockerham estimator assigner -----
# # if (!require("devtools")) install.packages("devtools")
# # devtools::install_github("thierrygosselin/assigner")
# # library(assigner)
# #calculate pairwise Fst with the weir-cockerham estimator. use assigner package, which requires tidy data
# test.assigner<-radiator::tidy_genomic_data(
#   data = paste('radseq_analysis/populations.r',mdat,'/populations.snps.gen',sep=""), 
#   strata = "radseq_analysis/strata.tsv") %>% dplyr::glimpse()
# #run pairwise calculations
# fst <- assigner::fst_WC84(
#   data = test.assigner, 
#   pop.levels = NULL,
#   pairwise = TRUE,
#   ci = TRUE,
#   iteration.ci = 100,
#   filename = "testing_fst", #will trigger the function to write the results in directory as well
#   verbose = TRUE
# )

# compare to Hierfstat -------
#open previously saved results and compare to these:

#convert format into genetic data frame
x.df<-genind2hierfstat(x,pop=x@pop)

#repeat for genetic groups
x.df.gg<-genind2hierfstat(x,pop=x@strata$Genetic.Group)

# global Fst -------
global.stat<-wc(x.df)
global.stat.gg<-wc(x.df.gg)

#bootstrap for confidence interval
global.fst.boot<-boot.vc(levels=x.df[,1], loci=x.df[,-1], nboot = resample)
global.fst.boot.gg<-boot.vc(levels=x.df.gg[,1], loci=x.df.gg[,-1], nboot = resample)

#write using format FST	STRATA CI_LOW	CI_HIGH , matching assigner results as closely as possible
sample.area<-c(global.stat$FST,"Sample.Area",global.fst.boot$ci[1,2],global.fst.boot$ci[3,2])
genetic.group<-c(global.stat.gg$FST,"Genetic.Group",global.fst.boot.gg$ci[1,2],global.fst.boot.gg$ci[3,2])
results.global<-rbind(sample.area,genetic.group)
colnames(results.global)<-c("FST","STRATA","CI_LOW","CI_HIGH")

# mean pairwise Fst ----------------- 
#for sampling areas
pairwise.fst<-pairwise.WCfst(x.df)

#pairwise upper lower limits for sampling areas
pairwise.fst.boot<-boot.ppfst(x, nboot = resample)

# write using format POP1	POP2	FST	CI_LOW	CI_HIGH
dist.pop.pairs<-t(combn(row.names(pairwise.fst.boot$ll),2))
dist.pop<-data.frame(POP1=dist.pop.pairs[,1], POP2 = dist.pop.pairs[,2],
                      FST=pairwise.fst[dist.pop.pairs],
                      CI_LOW=pairwise.fst.boot$ll[dist.pop.pairs],
                      CI_HIGH=pairwise.fst.boot$ul[dist.pop.pairs])


#repeat for genetic groups
pairwise.fst.gg<-pairwise.WCfst(x.df.gg)

#pairwise upper lower limits for sampling areas
pairwise.fst.boot.gg<-boot.ppfst(x.df.gg, nboot = resample)

dist.gg.pairs<-t(combn(row.names(pairwise.fst.boot.gg$ll),2))
dist.gg<-data.frame(POP1=dist.gg.pairs[,1], POP2 = dist.gg.pairs[,2],
                     FST=pairwise.fst.gg[dist.gg.pairs],
                     CI_LOW=pairwise.fst.boot.gg$ll[dist.gg.pairs],
                     CI_HIGH=pairwise.fst.boot.gg$ul[dist.gg.pairs])

# write -----
#write to paste(datadir, gen_path_name,'/populations.snps.gen',sep="")
write.table(results.global,paste(datadir,gen_path_name,'/fst.overall.tsv',sep=""),
            quote=FALSE, sep="  ",row.names=FALSE)

write.table(dist.pop,paste(datadir,gen_path_name,'/pairwise.fst.tsv',sep=""),
            quote=FALSE, sep="  ",row.names=FALSE)

write.table(dist.gg,paste(datadir,gen_path_name,'/pairwise.fst.geneticgroups.tsv',sep=""),
            quote=FALSE, sep="  ",row.names=FALSE)

