#Pst:Fst analyses: In two forms: Pst-Fst across populations, and a pairwise matrix-based approach

sink(paste(outputdir.choice,"/PstFst_output_md",mdat,"_h",h.squared,"_agp",add.gen.proportion,"_",phenotype.option,".txt",sep=""))
print(paste("r = ",mdat,", h.squared = ",h.squared,", c (proportion additive genetic variation) = ",
            add.gen.proportion,sep = ""))
# Fst Weir-Cockerham --------
#take the pair-column structure and reorganize into a distance matrix structure

print("Weir-Cockerham Pairwise Fst")
fst.pairwise.wc %>% print

#calculate observed mean Fst
print("Weir-Cockerham Mean Fst")
fst.mean.wc<-raw.mean.fst.wc$FST[1] %>% print
fst.modelled.variance<-(raw.mean.fst.wc$CI_HIGH - fst.mean.wc)/2 #larger interval, more conservative
bootstrap.fst.wc<-rnorm(resample,mean=fst.mean.wc,sd=sqrt(fst.modelled.variance))

# Pst: choose morphology --------
#choose how you want to measure "morphology"
sample.phenotype<-mor.block
sample.populations<-metadata$population

# Pst start ------------
#get numbers to calculate the two variances used in Pst
phenotype.pd<-partition.distance(sample.phenotype,sample.populations,unique(sample.populations))

#calculate those two variances
phenotype.v<-genetic.components.variance(add.gen.proportion,h.squared,phenotype.pd$n,phenotype.pd$k,
                                  phenotype.pd$SSB,phenotype.pd$SSW)

print("Observed Pst")
Pst.observed<-Pst(sample.phenotype,add.gen.proportion,h.squared,
                  sample.populations,unique(sample.populations)) %>% print

# # Pst:Fst comparison -------------
# sample V-a,within [additive genetic variance within populations].
#In the Pst world, V-a,within is (h^2)V-p,within. In this Qst CI world, you take your [(h^2)V-p] value
#divide by df, then multiply it by a random number from a chi-square distribution with appropriate d.f.
#[df = #individual-#groups], n-k
V.within.resample<-phenotype.v$V.within/(phenotype.pd$n-phenotype.pd$k) *
  rchisq(resample,phenotype.pd$n-phenotype.pd$k)

#estimate V-g, among, by algebraically reorganizing the Qst formula assuming Qst=Fst (null hypothesis)
#and therefore making V-g,among equal to 2*Fst*V-a,within/(1-Fst).
#DON'T need to multiply it by c to make it the Pst-world value of [c*V-p,among] b/c it's estimated, not observed
#divide by df, then multiply that whole value by a random number for a chi-square distribution with
#df = number of populations - 1
new.V.among<-(2*fst.mean.wc*phenotype.v$V.within/(1-fst.mean.wc))
V.among.resample<-new.V.among/(phenotype.pd$k-1) * rchisq(resample,phenotype.pd$k-1)
# 
# #After all of that, they'll repeatedly make simulated estimates of Qst-Fst. You resample Fst, you take 
# #this resample of V-a,within and V-g,among, and calculate the equation using the 3 simulated numbers.
Pst.resample<-V.among.resample/(V.among.resample + 2 * V.within.resample)
PstFst.distribution.WC<-Pst.resample-bootstrap.fst.wc
PstFst.observed.WC<-Pst.observed-fst.mean.wc 

p.wc<-which(PstFst.distribution.WC>=PstFst.observed.WC) %>% length/(resample+1) %>% print

# plot -----
pdf(file=paste(outputdir.choice,"/pst_wc_md",mdat,"_h",h.squared,"_agp",add.gen.proportion,".pdf",sep=""))
hist(PstFst.distribution.WC,xlim = range(c(PstFst.distribution.WC,PstFst.observed.WC)),
     main=paste("Pst = ", round(Pst.observed,3),", p = ",round(p.wc,3),sep=""))
abline(v=PstFst.observed.WC,col="red")
dev.off()