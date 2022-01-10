#Pst:Fst analyses: In two forms: Pst-Fst across populations, and a pairwise matrix-based approach
# For Lind et al. 2011 functions, need a global Fst value, Fst values for each locus, and number of loci
n.resamples<-resample #resample object from base or fst_rX script

# Observed Fst values -------
hierarchical.fst<-readRDS(paste(datadir,gen_path_name,'/hierarchical.fst.rds',sep=""))
colnames(hierarchical.fst$loc)<-c("Genetic.Group","Sampling.Area","Ind","Error") #not necessary but helpful reminder

n.loci<-nrow(hierarchical.fst$loc)

#calculate the global Fst values
global.Fst.GG<-hierarchical.fst$F[1,1]
global.Fst.SA<-hierarchical.fst$F[2,2]

# Observed Pst values --------
# Trying a hierarchical ANOVA as recommended by reviewer for submitted L305-307
#start with a geomorph data frame, assuming mor.block is points, not PCs

if (phenotype.option == "points") {
  gdf<-geomorph.data.frame(coords=arrayspecs(mor.block,1024,3), 
                           SA=factor(metadata$`Sampling Area`),
                           GG=factor(metadata$`Genetic Group`))
  #fit nested ANOVA model
  fit.gdf<-procD.lm(coords ~ GG/SA, data = gdf, print.progress = FALSE)
  #look at ANOVA table, with example of correct adjustment of F and p values
  anova.morph<-anova(fit.gdf, error = c("GG:SA", "Residuals"))
}

if (phenotype.option == "PCs")    {
  fit.gdf<-procD.lm(mor.block ~ metadata$`Genetic Group`/ metadata$`Sampling Area`, print.progress = FALSE)
  #look at ANOVA table, with example of correct adjustment of F and p values
  anova.morph<-anova(fit.gdf, error = c("metadata$`Genetic Group`:metadata$`Sampling Area`", "Residuals"))
}

#pull mean square values from the variance components
msGG<-anova.morph$table$MS[1] #MS are not the same as variance components.
msSA<-anova.morph$table$MS[2]
msIN<-anova.morph$table$MS[3]

#calculate variance components
vpGG<-(msGG-msIN)/ (3*11) #3 being number of genetic groups
vpSA<-(msSA-msIN)/11 #11 being number of sampling areas
vpIN<-msIN #within-individual variance component, or just the MS

#pull degrees of freedom
df.indv<-anova.morph$table$Df[3]
df.sa<-anova.morph$table$Df[2]
df.gg<-anova.morph$table$Df[1]
# 
# #calculate observed Pst
# Pst.observed.gg<-pst.equation(add.gen.proportion, h.squared, v.between = vpGG, v.within = vpSA+vpIN)
# Pst.observed.sa<-pst.equation(add.gen.proportion, h.squared, v.between = vpSA, v.within = vpIN)
# 
# # Test for difference between Pst and Fst --------
# #set up variables
# #for genetic groups
# var.comp.within<-(vpSA + vpIN) * h.squared
# df.within<-df.indv
# df.between<-df.gg
# observed.Fst<-global.Fst.GG
# 
# #for sampling areas
# # var.comp.within<-(vpIN) * h.squared
# # df.within<-df.indv
# # df.between<-df.sa
# # observed.Fst<-global.Fst.SA
# 
# #Step 5: Repeat Steps 1-4 many times to build a null distribution of Qst-Fst
# null.distribution<-NULL
# for (i in 1:resample){
#   #Step 1 of Whitlock and Guillame 2009: calculate a bootstrap confidence interval for Fst by randomly resampling loci
#   resampled.Fst<-resample.Fst(hierarchical.fst$loc,n.loci)[1]
#   #Steps 2-4
#   null.distribution[i]<-simulate.Pst.Fst(add.gen.proportion,var.comp.within,df.within,df.between,observed.Fst)
# }
# # hist(null.distribution)
# 
# #Step 6: Test for significance by calculating the number of simulated Qst-Fst values that are less than observed.
# #(can use as a p-value)
# #which observed value?
# Observed.Pst.Fst<-Pst.observed.gg - global.Fst.GG
# # Observed.Pst.Fst<-Pst.observed.sa - global.Fst.SA
# 
# p.val<-mean(null.distribution<Observed.Pst.Fst) #proportion of observations less than observed value
# 
# permute additive genetic variance for Sampling Area--------
h2.levels<-seq(from=0.1, to=1, by=0.1)
prop.sim.SA<-matrix(data=NA,nrow=10,ncol=10,
       dimnames=list(seq(from=0.1, to=1, by=0.1),
                     seq(from=0.1, to=1, by=0.1)))

print(paste("Observed Fst for Sampling Area:",global.Fst.SA))

for (j in 1:length(h2.levels)){
  for (k in 1: length(h2.levels)){
    h.squared<-h2.levels[j]
    add.gen.proportion<-h2.levels[k]  #c, assumed additive genetic proportion of differences between populations
    var.comp.within<-(vpIN) * h.squared
    df.within<-df.indv
    df.between<-df.sa
    observed.Fst<-global.Fst.SA
    null.distribution<-NULL
    for (i in 1:resample){
      #Step 1 of Whitlock and Guillame 2009: calculate a bootstrap confidence interval for Fst by randomly resampling loci
      resampled.Fst<-resample.Fst(hierarchical.fst$loc,n.loci)[1]
      #Steps 2-4
      null.distribution[i]<-simulate.Pst.Fst(add.gen.proportion,var.comp.within,df.within,df.between,observed.Fst)
    }
    # hist(null.distribution)
    
    #Step 6: Test for significance by calculating the number of simulated Qst-Fst values that are less than observed.
    #(can use as a p-value)
    #which observed value?
    Pst.observed.sa<-pst.equation(add.gen.proportion, h.squared, v.between = vpSA, v.within = vpIN)
    Observed.Pst.Fst<-Pst.observed.sa - global.Fst.SA
    prop.sim.SA[j,k]<-mean(null.distribution<Observed.Pst.Fst) #proportion of observations less than observed value
    print(paste("Where h.squared = ",h.squared, " and c = ", add.gen.proportion,":",sep=""))
    print(paste("Observed Pst for Sampling Area:",Pst.observed.sa))
    print(paste("Observed Pst-Fst for Sampling Area:",Observed.Pst.Fst))
    print(paste("Pst-Fst significance:", mean(null.distribution<Observed.Pst.Fst)))
  }
}

write.csv(prop.sim.SA,paste(outputdir.choice,"/significance_Pst_Fst_SA.csv",sep=""))

prop.sim.long<-prop.sim.SA %>% melt

plotA<-ggplot(data = prop.sim.long, 
       mapping = aes(x = Var1, y = Var2, fill = 1-value)) +
  geom_tile() +
  geom_text(aes(label = value), 
            vjust = .5, alpha = 1,size=2) + #fontface  = "bold", 
  scale_fill_gradientn(limits = c(0,1),
                       colours=GetColors(n=4,scheme="sunset"),
                       values=c(0,0.9,0.95,0.99,1)) +
  theme_minimal()  + xlab("Heritability") + ylab("C")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none", text = element_text(size = 8),
        axis.text.y=element_text(margin = margin(r=-10)),
        axis.text.x=element_text(margin = margin(t=-10))) +
  scale_x_continuous(breaks=seq(from=0.1, to=1, by=0.1)) +
  scale_y_continuous(breaks=seq(from=0.1, to=1, by=0.1))

ggsave(filename=paste(outputdir.choice,"/Pst-Fst_SA.pdf",sep=""), plotA, width = 8, height = 8,units="cm",dpi=600)

# permute additive genetic variance for Genetic Group--------
h2.levels<-seq(from=0.1, to=1, by=0.1)
prop.sim.GG<-matrix(data=NA,nrow=10,ncol=10,
                    dimnames=list(seq(from=0.1, to=1, by=0.1),
                                  seq(from=0.1, to=1, by=0.1)))

print(paste("Observed Fst for Genetic Group:",global.Fst.GG))

for (j in 1:length(h2.levels)){
  for (k in 1: length(h2.levels)){
    h.squared<-h2.levels[j]
    add.gen.proportion<-h2.levels[k]  #c, assumed additive genetic proportion of differences between populations
    var.comp.within<-(vpSA + vpIN) * h.squared
    df.within<-df.indv
    df.between<-df.gg
    observed.Fst<-global.Fst.GG
    null.distribution<-NULL
    for (i in 1:resample){
      #Step 1 of Whitlock and Guillame 2009: calculate a bootstrap confidence interval for Fst by randomly resampling loci
      resampled.Fst<-resample.Fst(hierarchical.fst$loc,n.loci)[1]
      #Steps 2-4
      null.distribution[i]<-simulate.Pst.Fst(add.gen.proportion,var.comp.within,df.within,df.between,observed.Fst)
    }
    # hist(null.distribution)
    
    #Step 6: Test for significance by calculating the number of simulated Qst-Fst values that are less than observed.
    #(can use as a p-value)
    #which observed value?
    Pst.observed.gg<-pst.equation(add.gen.proportion, h.squared, v.between = vpGG, v.within = (vpSA+vpIN))
    Observed.Pst.Fst<-Pst.observed.gg - global.Fst.GG
    prop.sim.GG[j,k]<-mean(null.distribution<Observed.Pst.Fst) #proportion of observations less than observed value
    print(paste("Where h.squared = ",h.squared, " and c = ", add.gen.proportion,":",sep=""))
    print(paste("Observed Pst for Genetic Group:",Pst.observed.gg))
    print(paste("Observed Pst-Fst for Genetic Group:",Observed.Pst.Fst))
    print(paste("Pst-Fst significance:", mean(null.distribution<Observed.Pst.Fst)))
  }
}

write.csv(prop.sim.GG,paste(outputdir.choice,"/significance_Pst_Fst_GG.csv",sep=""))

prop.sim.long<-prop.sim.GG %>% melt

plotB<-ggplot(data = prop.sim.long, 
       mapping = aes(x = Var1, y = Var2, fill = 1-value)) +
  geom_tile() +
  geom_text(aes(label = value), 
            vjust = .5, alpha = 1,size=2) + #fontface  = "bold", 
  scale_fill_gradientn(limits = c(0,1),
                       colours=GetColors(n=4,scheme="sunset"),
                       values=c(0,0.9,0.95,0.99,1)) +
  theme_minimal()  + xlab("Heritability") + ylab("C")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none", text = element_text(size = 8),
        axis.text.y=element_text(margin = margin(r=-10)),
        axis.text.x=element_text(margin = margin(t=-10))) +
  scale_x_continuous(breaks=seq(from=0.1, to=1, by=0.1)) +
  scale_y_continuous(breaks=seq(from=0.1, to=1, by=0.1))

ggsave(filename=paste(outputdir.choice,"/Pst-Fst_GG.pdf",sep=""), plotB, width = 8, height = 8,units="cm",dpi=600)
