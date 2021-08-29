
#genetic Principal Components Analysis
PCA.genotype <- dudi.pca(gen.block,center=TRUE,scale=FALSE,nf=79,scannf=FALSE)

#look for specimens with possibly too much missing data
which(PCA.genotype$li[,1] < 0.01 & PCA.genotype$li[,1] > -0.01) %>% metadata$filename[.] %>% print()   
#OMNH-52913, MSB-66182 are essentially at (0,0), so all missing data? Skip up a few liens and remove these from datasets

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
ggplot(data=master.variables,aes(x=Axis1,y=Axis2, color = `Sampling Area`,shape=`Genetic Group`))+
  geom_point(size=2) +
  theme_classic() + theme(legend.position="none") +
  xlab(paste("PC 1 (",percents.gen[1],"%)",sep="")) +
  ylab(paste("PC 2 (",percents.gen[2],"%)",sep=""))
ggsave("PCA_gen.pdf", device = cairo_pdf, width = 8, height = 8,units="cm",dpi=600)
