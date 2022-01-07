#functions for calculating Pst and Fst
# # Post-Review Revised Functions -------
# Both Fst and Pst will be calculated using the same underlying equation, for comparability. 
pst.equation<-function(add.gen.proportion,h.squared,v.between,v.within){
  #Should look like the Fst hierarchial functions, but with additive genetic corrections
  # Pst.gg<-vpGG*c / (vpGG*c+ vpSA + vpIN) #an example for 2 hierarchical levels
  # Pst.sa<-vpSA*c / (vpSA*c + 2*h.squared*vpIN)
  # c = add.gen.proportion = proportion of total variance presumed to be because of additive genetic effects across populations (Brommer 2011)
  # h.squared = heritability
  # v.between = variance component from desired between-group level of hierarchy
  # v.within = variance component from desired within-group level of hierarchy (often sum of lower levels)
  #do the Pst calculation.combined equation:
  # Numerator: c * variance.between
  # Denominator: c * variance.between + h.squared  * 2 * variance.within
  V.among<-add.gen.proportion * v.between
  V.within<-h.squared * v.within
  Pst<-V.among/(V.among + 2 * V.within) #identical answer to Pstat, but can now do multidimensional distance!
  return(Pst)
}

#Alter Lind et al. approach to match the hierarchical checks below (taking mean of by-locus Fst does not give the correct value)
#per.locus.variance is a matrix structured like the $loc output of hierfstat package's varcomp.glob() function
#For now, function is tailored to a particular dataset and hierarchy, but could be modified to be more general.
resample.Fst<-function(per.locus.variance,n.loci = nrow(per.locus.variance)){
  rows2sample<-(sample(n.loci,replace=T))
  resampled.loci<-per.locus.variance[rows2sample,]
  vGG<-apply(resampled.loci,2,sum)[1]
  vSA<-apply(resampled.loci,2,sum)[2]
  vIN<-apply(resampled.loci,2,sum)[3]
  vER<-apply(resampled.loci,2,sum)[4]
  Fst.GG<-vGG/(vGG + vSA + vIN + vER)
  Fst.SA<-vSA/(vSA + vIN + vER)
  return(c(Fst.GG,Fst.SA))
}

simulate.Pst.Fst<-function(add.gen.proportion,var.comp.within,df.within,df.between,observed.Fst){
  #Step 2: construct a distribution for Va (additive genetic variance within populations) by 
  #multiplying it by a random number from a chi-squared distribution with the appropriate df
  #in the Pst world, Va = h.squared * V-p,within (within-population phenotypic variance component)
  #in an ANOVA-based framework, V-p,within = MS-residuals or individuals
  # Va.hat<-Va/df.within*rchisq(1,df=df.within) #Lind et al. version of this step
  v.within.hat<- var.comp.within / df.within * rchisq(1,df=df.within)
  #Step 3: estimate V-g,among by algebraically reorganizing the Qst formula assuming Qst=Fst (null hypothesis) 
  #and therefore making V-g,among [which is c*V-p,among in Pst world) equal to 2*Fst*V-a,within/(1-Fst). 
  #And then multiply that whole value by a random number for a chi-square distribution with df = number of populations - 1
  #In the Pst world, Vg = Vg*c, so might be good to include it here. 
  #Lind et al. version of this step
  # Vpop<-2*Global.Fst*Va/(1-Global.Fst) 
  # Vpop.hat<-Vpop/df.pop*rchisq(1,df.pop)
  v.between<- 2 * observed.Fst * var.comp.within / (1 - observed.Fst)
  v.between.hat<- v.between / df.between * rchisq(1,df.between)
  #Step 4: Use values from Steps 1-3 calculate a new, hypothetical Qst, and subtract it from the observed global Fst
  #Lind et al. version of this step
  # Qst<-Vpop.hat/(Vpop.hat+2*Va.hat) #estimated Qst value
  # return(Qst-Global.Fst)
  Qst.estimated<-v.between.hat / (v.between.hat + 2 * v.within.hat)
  return(Qst.estimated - observed.Fst)
}


# TO BE REMOVED: Lind et al. 2011 code, for reference and checking -------------
# The example below is written by, and uses data from, Lind MI, PK Ingvarsson, H Johansson, D Hall, and F Johansson. 2010. "Gene flow and selection on phenotypic plasticity in an island system of Rana temporaria." Evolution
# From Supplementary file 1. R code used for simulating the distribution of Qst-Fst for a neutral trait, and estimating the quantile of the simulated distribution that has more extreme values than the observed value of Qst-Fst. For more information on the theoretical foundations, see Whitlock & Guillaume (2009).
#  https://doi.org/10.1111/j.1558-5646.2010.01122.x

# #EXAMPLE
# #INPUT DATA
# 
# #From mixed model, partition of variation for quantitative trait
# Va<-1.03 #additive genetic variance within populations, includes the heritability tweaks
# Vpop<-0.9133 #additive variance between populations
# df.within<-70 #is this 81 - 10 - 1? Equivalent to n - k - 1
# df.pop<-9 #10 populations - 1? Equivalent to k-1
# QST<-0.3071671 #the observed Qst value, won't get used until near end, for p-value and plotting
# 
# #From microsatellite data
# Global.Fst<-0.051
# fst<-c(0.008,0.081,0.045,0.079,0.036) #vector of Fst values (for each loci)
# n.loci<-5 #number of loci

#SIMULATION
#Simulating distribution of Qst-Fst for a neutral trait

# Fst.Qst<-function(Va,df.within,Global.Fst,df.pop){
#   Va.hat<-Va/df.within*rchisq(1,df=df.within)
#   Vpop<-2*Global.Fst*Va/(1-Global.Fst)
#   Vpop.hat<-Vpop/df.pop*rchisq(1,df.pop)
#   Qst<-Vpop.hat/(Vpop.hat+2*Va.hat)
#   return(Qst-Global.Fst)
# }
# 
# sim.Qst.Fst<-NULL
# for(i in 1:10000){
#   sim.Qst.Fst[i]<-Fst.Qst(Va,df.within,mean(sample(fst,n.loci,replace=T)),df.pop)
# }

# #EXAMPLE
# #RESULTS
# hist(sim.Qst.Fst) # Distribution of Qst-FSt for simulated neutral trait
# 
# #P-TEST
# 1-mean(sim.Qst.Fst<QST-Global.Fst) 
# #The quantile of the simulated distribution that has more extreme values than the observed value of Qst-Fst
# 
# #FIGURE
# 
# par(mai=c(1,1,0.5,0.5))
# hist(sim.Qst.Fst,xlab="Qst-Fst",breaks=c(20), las=0,xlim=c(-0.1,QST+0.1),main="Trait name")
# arrows(QST-Global.Fst,800,QST-Global.Fst,100,length=0.1)
# #The arrow indicates the observed Qst-Fst
# Pst notes --------
# Polly:
#with a jackknife (leave out one specimen, recalculate)
# Numerator: (n-k-1)SSB*
# Denominator: (k-1)(SSW*+SSB*)
#SSB: sum of squared distance between group means and the grand mean
#SSW: sum of squared distances between individuals and their group mean
#k = # groups
#n = sample size
#distance for shapes is Procrustes distance but I'm using PCs so Euclidean makes sense

# Brommer: #Qst: 
# Numerator: sig^2B 
# Denominator: sig^2B + 2*sig^2W
#where  B is for additive genetic variance between groups
#and W is for additive genetic variance within groups

#to Pst:
#sig^2B becomes phenotypic variance between populations (symbol changes too)
#multiply the sig^2B every time by a scalar c [proportion of variance that is because of additive genetic effects]
#sig^2W becomes the phenotypic variance within populations
#multiply sig^2W times h^2 [narrow sense heritibility] 

#remember from basic stats: variance = average squared differences from the mean 
#(divided by N-1 instead of N for real world samples, which is why the extra numbers in Polly's equation)
#think of differences as distance and the core of the equations are similar

#last major difference is where the 2 comes from in denominator. See box 1 of Leinonen et al. 2013
#comes from how you convert Fst to Qst
# # Checks ----
# #Just so you can reassure yourself where the hierarchical Fst values are coming from 
# #This is the Fst for proportion of variance partioned at the level of genetic groups
# Fst.GG<-hierarchical.fst$F[1,1]
# 
# #variance copmonents
# vcGG<-hierarchical.fst$overall[1]
# 
# #how variance components come from by-locus values
# vGG<-apply(hierarchical.fst$loc,2,sum)[1]
# vcGG==vGG #should be TRUE
# 
# #now calculate the other four
# vSA<-apply(hierarchical.fst$loc,2,sum)[2]
# vIN<-apply(hierarchical.fst$loc,2,sum)[3]
# vER<-apply(hierarchical.fst$loc,2,sum)[4]
# 
# #use variance components to calculate Fst for Genetic Groups
# FGG<-vGG/(vGG + vSA + vIN + vER)
# FGG==Fst.GG #should be TRUE