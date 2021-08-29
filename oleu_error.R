# #Evaluate cropping error
# (based on erinaceomorphs, pseudolandmark placement error seems to have been solved)
# general approach: calculate repeatability values from those subset, averaging replicates across dataset

# create labels ---------
#label which specimens are repeate crops vs. everything else
repeated<-critters$filename[which(critters$cropped==3)] 
label_repeats<-NULL
for(i in 1:length(label_error)){
  if (specimen[i] %in% repeated){label_repeats[i]<-specimen[i]}
  else (label_repeats[i]<-"other")
}
label_repeats<-factor(label_repeats)

# shape error ----------
#replicates in context of total dataset
PCA_er<-prcomp(molars$m2d[,],scale.=FALSE) #PCA of all surfaces
testgdf2<-geomorph.data.frame(coords=PCA_er$x[,],specimen=specimen)
errorANOVA2<-procD.lm(coords~specimen,data=testgdf2,iter=999,RRPP=TRUE) %>% .$aov.table
mean_rep<-mean(shapemetadata$cropped)
err_cr2<-error3d(errorANOVA2,mean_rep,f1=1,f2=2)

# plot shape error -------
jpeg(filename=paste("error_crop_pca_all.jpg",sep=""))
plot(PCA_er$x[,1:2],bg=GetColors(n=length(levels(label_repeats)))[label_repeats],pch=21,cex=1.5)
title(paste("Cropping Repeatability =", round(err_cr2$repeatability,3)))
dev.off()

jpeg(filename=paste("error_crop_by_pc.jpg",sep=""))
repeatPCs_result<-find_repeatablePCs(PCA_er$x,specimen,rep=mean_rep) #due to cropping error
plot(repeatPCs_result,xlab="principal components")
lines(repeatPCs_result)
abline(h=0.90,col="red",lty=2)
abline(h=0.80,col="blue",lty=3)
title(paste("# repeatable PCs =", (min(which(round(repeatPCs_result,2)<0.8))-1)))
dev.off()

# example shape error --------
#example shapes
newshapes1<-rotateMorphologika(molars$scaled[,,which(specimen==repeated[1])],degreesX=0,degreesY=0,degreesZ=0) #compare recrops
newshapes2<-rotateMorphologika(molars$scaled[,,which(specimen==repeated[2])],degreesX=0,degreesY=0,degreesZ=0)
newshapes3<-rotateMorphologika(molars$scaled[,,which(specimen==repeated[3])],degreesX=0,degreesY=0,degreesZ=0)

#plot first
meanshape<-apply(newshapes1,1:2,mean) #calculate mean shape
variances<-apply(newshapes1,1:2,var) %>% apply(.,1,mean) #find mean variance for all points, unrotated
colbydiff<-rescale(variances,to=c(1,ncut)) %>% round
result<-raintable[colbydiff]

open3d()
plot3d(meanshape,axes=F,col=result,size=10,xlab="",ylab="",zlab="") #some bright colors near cropping margin, but also at local minima
writePLY(paste("error_example_crop_1.ply",sep=""),format="ascii",pointRadius=0.005)
rgl.close()

#plot second
meanshape<-apply(newshapes2,1:2,mean) #calculate mean shape
variances<-apply(newshapes2,1:2,var) %>% apply(.,1,mean) #find mean variance for all points, unrotated
colbydiff<-rescale(variances,to=c(1,ncut)) %>% round
result<-raintable[colbydiff]

open3d()
plot3d(meanshape,axes=F,col=result,size=10,xlab="",ylab="",zlab="") #some bright colors near cropping margin, but also at local minima
writePLY(paste("error_example_crop_2.ply",sep=""),format="ascii",pointRadius=0.005)
rgl.close()

#plot third
meanshape<-apply(newshapes3,1:2,mean) #calculate mean shape
variances<-apply(newshapes3,1:2,var) %>% apply(.,1,mean) #find mean variance for all points, unrotated
colbydiff<-rescale(variances,to=c(1,ncut)) %>% round
result<-raintable[colbydiff]

open3d()
plot3d(meanshape,axes=F,col=result,size=10,xlab="",ylab="",zlab="") #some bright colors near cropping margin, but also at local minima
writePLY(paste("error_example_crop_3.ply",sep=""),format="ascii",pointRadius=0.005)
rgl.close()
# linear measurement error --------
error.data<-shapemetadata %>% 
  filter(!is.na(m1_length_1)) %>%
  select(filename,m1_length_1,m1_length_2,m1_length_3,m1_width_1,m1_width_2,m1_width_3)

in_df_long<-melt(error.data,id=c("filename")) 
in_df_long$orientation<-gsub("m[1-3]_(length|width)_[1-3]","\\1",in_df_long$variable)
in_df_long$replicate<-gsub("m[1-3]_(length|width)_([1-3])","\\2",in_df_long$variable)

ANOVA.length<-anova(lm(value~filename, data = filter(in_df_long,orientation=="length")))
ANOVA.width<-anova(lm(value~filename, data = filter(in_df_long,orientation=="width")))

#number of replicates
r<-max(in_df_long$replicate) %>% as.numeric
error2d(ANOVA.length,r)
error2d(ANOVA.width,r)