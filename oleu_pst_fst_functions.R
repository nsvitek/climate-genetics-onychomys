#functions for calculating Pst and Fst

# Pst functions -------
group.means<-function(sample.morph, subpop.choice, subpop.vector){
  temp.match<-which(subpop.vector==subpop.choice) #find only the individuals in the subpopulation
  subpop.means<-sample.morph[temp.match,] %>% #take morphology from those individuals
    as.matrix(.,nrow=length(temp.match)) %>% #make sure it stays in matrix format (for unidimensional measures)
    colMeans
  return(subpop.means)
}

#calculates SSW, or sum of squared distances between individuals and their group mean (Polly 2007)
dist.individual.subpop<-function(sample.morph, subpop.options, subpop.choice, subpop.vector){
  temp.match<-which(subpop.vector==subpop.choice) #find only the individuals in the subpopulation
  temp.pop<-sample.morph[temp.match,] %>% #take morphology from those individuals
    as.matrix #make sure it stays in matrix format (for unidimensional measures)
  temp.pop<-rbind(temp.pop,colMeans(temp.pop)) #add subpop.mean to the end
  dist.individual2subpop.mean<-dist(temp.pop,method="euclidean") %>% #calculate distances
    as.matrix() %>% #make sure it stays in matrix format so you can pull rows
    .[nrow(temp.pop),1:(nrow(temp.pop)-1)] #and take just the last row, which should be distances from mean
  return(dist.individual2subpop.mean)
  # #comparing to Pst equation:
  # #if variance is (dist^2 %>% sum) / (n-1) and then you multiply by (n-1)
  # #you are left with (dist^2 %>% sum), where distance is distance of individuals from total mean
  # temp.match<-which(critters$population %in% pair)
  # temp.pop<-sample.morph[temp.match,] %>% as.matrix %>% rbind(.,grand.mean)
  # dist.ig<-dist(temp.pop) %>% as.matrix() %>% .[nrow(temp.pop),1:(nrow(temp.pop)-1)] #%>% .^2 *n.each[i]
  # SST<-(dist.ig^2) %>% sum
  # SSW2<-SST-SSB #is the same as SSW. This seems good. 
}

#calculates SSB, or sum of squared distances between group means and grand mean (Polly 2007), as well as
#pulling dist.individual.subpop to calculate SSW and the n, k variables you need for eq. 1 of Polly 2007
partition.distance<-function(sample.morph,subpop.vector,subpop.set){
  #shortening step: make object for the vector of possible subpopulations:
  pop.options<-unique(subpop.vector)
  k<-length(subpop.set) #total number of groups/subpopulations
  n<-which(subpop.vector %in% subpop.set) %>% length #total sample size
  n.each<-sapply(subpop.set, function(x) which(subpop.vector %in% x) %>% length) #sample size for each group
  grand.mean<-sample.morph[which(subpop.vector %in% subpop.set),] %>% as.matrix %>% colMeans #calculate grand mean
  #lapply not working, try a for loop instead
  subpop.means.pops<-matrix(0, nrow=k, ncol=ncol(sample.morph))
  for (i in 1:k){
    subpop.means.pops[i,]<-group.means(sample.morph, subpop.set[i], subpop.vector)
  }
  subpop.means<-rbind(subpop.means.pops,grand.mean) #add grand mean to the bottom of the matrix
  dist.group.grand<-dist(subpop.means,method="euclidean") %>% #calculate euclidean distances
    as.matrix() %>% #keep in matrix format so you can pull rows
    .[nrow(subpop.means),1:(nrow(subpop.means)-1)] #take last row, which should be distances from grand mean
  #SSB: sum of squared distances between group means and the grand mean
  SSB<-(dist.group.grand^2 * n.each) %>% sum
  #SSW: sum of squared distances between individuals and their group mean
  SSW<-sapply(subpop.set,function(x) dist.individual.subpop(sample.morph,pop.options,x,subpop.vector)) %>% #get distances
    unlist %>% #make into a single vector
    .^2 %>% sum #square and sum to get sum of squares. 
  result<-list(n,k,SSB,SSW)
  names(result)<-c("n","k","SSB","SSW")
  return(result)
}

genetic.components.variance<-function(add.gen.proportion,h.squared,n,k,SSB,SSW){
  V.among<-add.gen.proportion * SSB/(k-1)
  V.within<-h.squared * SSW/(n-k)
  result<-list(V.among,V.within)
  names(result)<-c("V.among","V.within")
  return(result)
}


pst.equation<-function(add.gen.proportion,h.squared,n,k,SSB,SSW){
  #do the Pst calculation. 
  #combined equation
  # Numerator: c * SSB/(k-1)
  # Denominator: c * SSB/(k-1) + h.squared  * 2 * SSW/(n-k)
  #SSW divided by n-k according to Pstat package function. Reason related to df?
  #a form closer to Psts package (but algebraically equivalent) is:
  # Numerator: c/h.squared * SSB/(k-1), or V.among
  # Denominator: c/h.squared * SSB/(k-1) + 2* SSW/(n-k), or V.among + V.within
  V<-genetic.components.variance(add.gen.proportion,h.squared,n,k,SSB,SSW)
  pst.est<-V$V.among/(V$V.among + 2 * V$V.within) #identical answer to Pstat, but can now do multidimensional distance!
  return(pst.est)
}

Pst<-function(sample.morph,add.gen.proportion,h.squared,subpop.vector,subpop.set){
  pd<-partition.distance(sample.morph,subpop.vector,subpop.set)
  pst<-pst.equation(add.gen.proportion,h.squared,pd[[1]],pd[[2]],pd[[3]],pd[[4]])
  return(pst)
}

Pst.pairwise<-function(sample.morph,add.gen.proportion,h.squared,subpop.vector){
  pop.options<-unique(subpop.vector)
  pairwise.pops<-combn(pop.options,2) %>% t
  vector.pst<-apply(pairwise.pops,1,function(x) Pst(sample.morph,add.gen.proportion,h.squared,subpop.vector,x))
  pw.pst<-matrix(0, nrow = length(pop.options), ncol = length(pop.options), 
                 dimnames = list(pop.options,pop.options)) #make empty matrix for pairwise Pst
  pw.pst[lower.tri(pw.pst,diag=FALSE)]<-vector.pst #fill in lower triangle
  pw.pst<-t(pw.pst) #fill in upper triangle by transposing. Thank you Dave Tang
  pw.pst[lower.tri(pw.pst,diag=FALSE)]<-vector.pst #re-fill in lower triangle
  return(pw.pst)
}
# #example
# Pst(sample.morph,csh,critters$population,pair)

# better Pst notes --------
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

