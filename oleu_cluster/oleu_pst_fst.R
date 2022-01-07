#Pst:Fst analyses: In two forms: Pst-Fst across populations, and a pairwise matrix-based approach
# For Lind et al. 2011 functions, need a global Fst value, Fst values for each locus, and number of loci
n.resamples<-100 #or resample object from base

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
gdf<-geomorph.data.frame(coords=mor.block, #arrayspecs(mor.block,1024,3)
                         SA=factor(metadata$`Sampling Area`),
                         GG=factor(metadata$`Genetic Group`))
#fit nested ANOVA model
fit.gdf<-procD.lm(coords ~ GG/SA, data = gdf, print.progress = FALSE)
#look at ANOVA table, with example of correct adjustment of F and p values
anova.morph<-anova(fit.gdf, error = c("GG:SA", "Residuals")) 

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

#calculate observed Pst
Pst.observed.gg<-pst.equation(add.gen.proportion, h.squared, v.between = vpGG, v.within = vpSA+vpIN)
Pst.observed.sa<-pst.equation(add.gen.proportion, h.squared, v.between = vpSA, v.within = vpIN)

# Test for difference between Pst and Fst --------
#set up variables
#for genetic groups
var.comp.within<-(vpSA + vpIN) * h.squared
df.within<-df.indv
df.between<-df.gg
observed.Fst<-global.Fst.GG

#for sampling areas
# var.comp.within<-(vpIN) * h.squared
# df.within<-df.indv
# df.between<-df.sa
# observed.Fst<-global.Fst.SA

#Step 5: Repeat Steps 1-4 many times to build a null distribution of Qst-Fst
null.distribution<-NULL
for (i in 1:100){
  #Step 1 of Whitlock and Guillame 2009: calculate a bootstrap confidence interval for Fst by randomly resampling loci
  resampled.Fst<-resample.Fst(hierarchical.fst$loc,n.loci)[1]
  #Steps 2-4
  null.distribution[i]<-simulate.Pst.Fst(add.gen.proportion,var.comp.within,df.within,df.between,observed.Fst)
}
# hist(null.distribution)

#Step 6: Test for significance by calculating the number of simulated Qst-Fst values that are less than observed.
#(can use as a p-value)
#which observed value?
Observed.Pst.Fst<-Pst.observed.gg - global.Fst.GG
# Observed.Pst.Fst<-Pst.observed.sa - global.Fst.SA

p.val<-mean(null.distribution<Observed.Pst.Fst) #proportion of observations less than observed value

# permute additive genetic variance --------
h2.levels<-seq(from=0.1, to=1, by=0.1)
prop.sim.SA<-prop.sim.GG<-matrix(data=NA,nrow=10,ncol=10,
       dimnames=list(c("h0.1","h0.2,","h0.3","h0.4","h0.5","h0.6","h0.7","h0.8","h0.9","h1.0"),
                     c("c0.1","c0.2,","c0.3","c0.4","c0.5","c0.6","c0.7","c0.8","c0.9","c1.0")))

for (j in 1:length(h2.levels)){
  for (k in 1: length(h2.levels)){
    h.squared<-h2.levels[j]
    add.gen.proportion<-h2.levels[k]
    var.comp.within<-(vpIN) * h.squared
    df.within<-df.indv
    df.between<-df.sa
    observed.Fst<-global.Fst.SA
    null.distribution<-NULL
    for (i in 1:100){
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
  }
}

ggplot(data=prop.sim.SA)

prop.sim.long<-prop.sim.SA %>% melt
ggplot(data = prop.sim.long, 
       mapping = aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  geom_text(aes(label = value), 
            vjust = .5, fontface  = "bold", alpha = 1) +
  scale_fill_gradient2(low="red", high="blue",mid="white",midpoint=0.05) +
  theme_minimal()  + #xlab("Reference Genus") + ylab("Predicted Genus")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(face="italic"),legend.position = "none",
        axis.text.y = element_text(face="italic", angle=90), 
        axis.title.y=element_text(size=8),
        axis.title.x=element_text(size=8)) #+
  ylim(rev(levels(confusion.table$reference.class)))
?scale_fill_gradient2
  