# disparity functions -------
############################################################################
#
#   Creates a function called individual.disparity() that calculates the 
#   morphological disparity among several specimens 
#   from continuous data, including PC scores.  Disparity is calculated as 
#   the mean squared distance among the specimens.  There is no standardization
#   because individual specimens do not have variances like groups do.
#
#   The format is: individual.disparity( data )
#
#   where data are the continuous data with individual specimens on
#   rows and variables in columns. Written by P. David Polly, 2008
#
############################################################################

individual.disparity <- function(d) {
  dists <-( dist(d))^2
  return(mean(dists))
}

############################################################################
#
#   Creates a function called group.disparity() that calculates the 
#   morphological disparity among the means of several groups of specimen 
#   from continuous data, including PC scores.  Disparity is calculated as 
#   the mean squared distance among the means of the taxa, either raw or 
#   standardized by the  variance within the taxa.
#
#   The format is: group.disparity( data , labels , standardize=T or F)
#
#   where data are the continuous data for the several taxa, with individual specimens on
#   rows and variables in columns, labels are row labels with the taxon name, and the 
#   option standardize determines whether the output is a raw mean distance squared or 
#   whether it is standardized by dividing by the mean variance within the taxa.  Written by
#   P. David Polly, 2008
#
############################################################################

group.disparity <- function(d, l,standardize=F) {
  mns <- array(data=NA,dim=c(length(d[1,]), length(levels(factor(l))) ))
  for (i in 1:length(mns[,1])) { mns[i,] <-  tapply(d[,i],INDEX=factor(l),FUN=mean) }
  dists <-( dist(aperm(mns,c(2,1))))^2
  
  vrs <- array(data=NA,dim=c(length(d[1,]), length(levels(factor(l))) ))
  for (i in 1:length(vrs[,1])) { vrs[i,] <-  tapply(d[,i],INDEX=factor(l),FUN=var) }
  vrs <- sum(vrs)/length(vrs[,1])
  
  if(standardize) { return(c(mean(dists)/vrs)) }
  else {return(mean(dists))}
}

# calculate disparity -------
#choose how size is being measured
mor.block<-mor.block %>% as.data.frame
var.factor<-factor(metadata$population)

#this value is the avg. distance of group means from grand mean
group.disp<-group.disparity(d = mor.block,l = var.factor, standardize = FALSE)

#individual disparity for roosevelt county
roosevelt.disp<-individual.disparity(mor.block[which(var.factor=="ip_roosevelt"),])

#individual disparity for each population
pop.indv.disp<-rep(NA,length(levels(var.factor)))
for (ind in 1:length(levels(var.factor))){ 
  pop.indv.disp[ind]<-individual.disparity(mor.block[which(var.factor==levels(var.factor)[ind]),])
}

# bootstrap -------
#okay, now how to build a reasonable bootstrap resample of each sample? 
#what's the hypothesis?
#that disparity between individuals of a population is higher than disparity between group means
#so a confidence interval around group means would randomly resample the total sample into 3 groups
#take the mean shape of each, and calculate individual disparity
nperm<-resample-1

#create a starter for random group disparity distribution
random.group.disparity<-group.disp

#created an empty matrix for random population-level individual disparity distribution
random.pop.disparity<-matrix(data=NA,nrow=length(pop.indv.disp),ncol=nperm+1,dimnames=list(c(levels(var.factor))))
#make the first column the observd values
random.pop.disparity[,1]<-pop.indv.disp

#randomly resample total mouse sample to populations with same original population sample size
#then calculate group disparity for each random sample
#then calculate individual disparity of each population for each random sample
for(j in 1:nperm){
  mor.block.random<-mor.block 
  for(i in c(1:length(levels(var.factor)))){
    which.var<-which(var.factor == levels(var.factor)[i]) 
    rando.sample<-sample(c(1:nrow(mor.block)), size = length(which.var))
    mor.block.random[which.var,]<-mor.block[rando.sample,]
  }
  random.group.disparity[j+1]<-group.disparity(d = mor.block.random,l = var.factor, standardize = F)
  for (ind in 1:length(levels(var.factor))){ 
    random.pop.disparity[ind,j+1]<-individual.disparity(mor.block.random[which(var.factor==levels(var.factor)[ind]),])
  }
}

#center the distributions around the observed value
random.distribution<-scale(random.group.disparity,center = TRUE, scale = FALSE) + group.disp

#center the distributions around the observed values
random.pop.disparity.scaled<-random.pop.disparity - apply(random.pop.disparity,1,mean) + pop.indv.disp

#put all resamples in a dataframe together
disparity.all<-data.frame(t(random.pop.disparity.scaled), group = random.distribution)

#report significance values
which(disparity.all$cp_apache < group.disp) %>% length/resample

apply(disparity.all, 2, function(x) which(x < group.disp) %>% length/resample) %>% 
  cbind(.,c(pop.indv.disp,group.disp)) %>%
  write.csv(.,"disparity_p_vals.csv")


#do some renaming
cbind(colnames(disparity.all),levels(factor(metadata$`Sampling Area`))[c(1,4,7,9,3,5,6,8,11,2,10)])
colnames(disparity.all)<-c(levels(factor(metadata$`Sampling Area`))[c(1,4,7,9,3,5,6,8,11,2,10)],"Z. Group Disparity")

#reformat all of those resamples to be more conducive to a ggplot boxplot
disparity2plot<-melt(disparity.all)

#disparity shouldn't be negative -- convert to zeroes
disparity2plot$value[which(disparity2plot$value < 0)]<- 0

#this is silly, but I have to unfactor in order to alphebatize
disparity2plot$variable<-disparity2plot$variable %>% as.character
disparity2plot<-arrange(disparity2plot,variable)
disparity2plot$variable<-disparity2plot$variable %>% factor

#and then rename the group disparity that I artificially set to be last
levels(disparity2plot$variable)[12]<-"Group Disparity"
disparity2plot$variable[which(disparity2plot$variable=="Z. Group Disparity")]<-"Group Disparity"
#set color palette to match previous plots
match.col<-c(hue_pal()(11),"black")
# plot ------
ggplot(data=disparity2plot,aes(x=variable,y=value)) + 
  geom_boxplot(aes(fill=variable))  +   
  scale_fill_manual(values = match.col) +theme_classic() + ylab("Disparity") +
  theme(axis.title.x=element_blank(), axis.text.x = element_text(angle = 45,hjust=1), legend.position="none")
  
ggsave("population_shape_disparity.pdf", device = cairo_pdf, width = 8, height = 8,units="cm",dpi=600)

ggplot(data=disparity2plot,aes(x=variable,y=value)) + 
  geom_boxplot(aes(fill=variable))  +   
  scale_fill_manual(values = match.col) +theme_classic() + ylab("Disparity") +
ggsave("population_shape_disparity_legend.pdf", device = cairo_pdf, width = 8, height = 16,units="cm",dpi=300)