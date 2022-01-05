sink(paste(outputdir.choice,"/distance_output_md",mdat,"_",phenotype.option,".txt",sep=""))
print(paste("r = ",mdat,sep = ""))

# calculate individual  distances -------
#genetic distance as PCA distance
# pw.gen<-dist(PCA.gen,diag=TRUE,upper=TRUE) %>% as.matrix
pw.morph<-dist(mor.block,method="euclidean",diag=TRUE,upper=TRUE) %>% as.matrix
pw.env<-dist(env.block,method="euclidean",diag=TRUE,upper=TRUE) %>% as.matrix
pw.geo<-distm(cbind(geo.block$decimallongitude, geo.block$decimallatitude),
              fun=distVincentyEllipsoid)
row.names(pw.geo)<-colnames(pw.geo)<-row.names(geo.block)

#POST-REVIEW: suggested change
#population genetic distance should be complete, not PCA distances
pw.gen<-dist(gen.block,diag=TRUE,upper=TRUE) %>% as.matrix

# calculate population distances -----

# #fiddly renaming, reorganizing
# cbind(colnames(fst.pairwise.wc),levels(factor(metadata$population))[c(1,4,7,9,3,5,6,8,11,2,10)])
# colnames(fst.pairwise.wc)<-rownames(fst.pairwise.wc)<-levels(factor(metadata$population))[c(1,4,7,9,3,5,6,8,11,2,10)]
# fst.pairwise.wc<-fst.pairwise.wc[order(colnames(fst.pairwise.wc)),order(colnames(fst.pairwise.wc))]

#take column means for each sampling area and neaten up resulting data frames
pop.ids<-metadata$population 

pop.mean.morph<-data.frame(pop.ids, mor.block) %>% 
  group_by(., pop.ids) %>% summarise_at(vars(-group_cols()),mean) %>% data.frame
pop.mean.env<-data.frame(pop.ids, env.block) %>% 
  group_by(., pop.ids) %>% summarise_at(vars(-group_cols()),mean) %>% data.frame
pop.mean.geo<-data.frame(pop.ids, geo.block) %>% 
  group_by(., pop.ids) %>% summarise_at(vars(-group_cols()),mean) %>% data.frame

sample.names<-pop.mean.morph$pop.ids 
pop.mean.morph<-subset(pop.mean.morph, select = -pop.ids)
pop.mean.env<-subset(pop.mean.env, select = -pop.ids)
pop.mean.geo<-subset(pop.mean.geo, select = -pop.ids)
row.names(pop.mean.morph)<-row.names(pop.mean.env)<-row.names(pop.mean.geo)<-sample.names

#calculate distance matrices for populations
pw.gen.pop<-fst.pairwise.wc
pw.morph.pop<-dist(pop.mean.morph,method="euclidean",diag=TRUE,upper=TRUE) %>% as.matrix
pw.env.pop<-dist(pop.mean.env,method="euclidean",diag=TRUE,upper=TRUE) %>% as.matrix
pw.geo.pop<-distm(cbind(pop.mean.geo$decimallongitude, pop.mean.geo$decimallatitude),
              fun=distVincentyEllipsoid)
row.names(pw.geo.pop)<-colnames(pw.geo.pop)<-row.names(pop.mean.geo)

# MMRR --------
#print all of these analyses
print("Individual morphology MMRR")
(results.morph<-MMRR(pw.morph,list(pw.gen,pw.env,pw.geo),nperm=resample))
print(results.morph)

print("Individual genetic MMRR")
(results.gen<-MMRR(pw.gen,list(pw.env,pw.geo),nperm=resample))
print(results.gen)

print("Population morphology MMRR")
(results.morph.pop<-MMRR(pw.morph.pop,list(pw.gen.pop,pw.env.pop,pw.geo.pop),nperm=resample))
print(results.morph.pop)

print("Population genetic MMRR")
(results.gen.pop<-MMRR(pw.gen.pop,list(pw.env.pop,pw.geo.pop),nperm=resample))
print(results.gen.pop)
# Reformat individuals for plotting -------
# dist.indv.gen<-distgenEUCL
dist.indv.gen<-t(combn(colnames(pw.gen),2))
# dist.indv.gen<-data.frame(dist.indv.gen,genetic=dist.indv.gen.mat[dist.indv.gen])
dist.indv<-data.frame(dist.indv.gen,genetic=pw.gen[dist.indv.gen],environment=pw.env[dist.indv.gen],
                          shape=pw.morph[dist.indv.gen],geography=pw.geo[dist.indv.gen])

# linear methods? -----
print("Phenotype-Genotype distance model")
dist.lm<-lm(shape~genetic, data = dist.indv)
summary(dist.lm) %>% print
#does the slope include 0?
confint(dist.lm,'genetic',level=0.95) %>% print

print("Genotype-Geography distance model")
dist.lm<-lm(genetic~geography, data = dist.indv)
summary(dist.lm) %>% print
#does the slope include 0?
confint(dist.lm,'geography',level=0.95) %>% print

print("Phenotype-Geography distance model")
dist.lm<-lm(shape~geography, data = dist.indv)
summary(dist.lm) %>% print
#does the slope include 0?
confint(dist.lm,'geography',level=0.95) %>% print

print("Environment-Phenotype distance model")
dist.lm<-lm(shape~environment, data = dist.indv)
summary(dist.lm) %>% print
#does the slope include 0?
confint(dist.lm,'environment',level=0.95) %>% print

# Plot -----
#MSB 75789, 75790, 75791, 75792, 75799, 75808, 75809
# dist.indv[which(dist.indv$genetic<0.2),]
# dist.indv[which(dist.indv$genetic>17.5),]

r.critters<-metadata$filename[which(metadata$population=="ip_roosevelt")]
pick.r<-which(dist.indv$X1 %in% r.critters & dist.indv$X2 %in% r.critters)

gr.critters<-metadata$filename[which(metadata$population=="cp_grand")]
pick.g<-which(dist.indv$X1 %in% gr.critters & dist.indv$X2 %in% gr.critters)

wb.critters<-metadata$filename[which(metadata$gengroup=="wb")]
pick.w<-which((! dist.indv$X1 %in% wb.critters & dist.indv$X2 %in% wb.critters)|
                (dist.indv$X1 %in% wb.critters & ! dist.indv$X2 %in% wb.critters))

ex.cols<-GetColors(5,scheme="bright")

plot1<-ggplot(data=dist.indv,aes(x=genetic,y=shape)) + 
  geom_point(alpha=0.3,shape=16,size=0.5) +  
  geom_point(data=dist.indv[pick.w,],color=ex.cols[5],shape=17,size=1) +
  geom_point(data=dist.indv[pick.r,],color=ex.cols[2],shape=17,size=1) +
  geom_point(data=dist.indv[pick.g,],color=ex.cols[4],shape=17,size=1) +
  # geom_smooth(method=lm, se=FALSE) + 
  theme_classic()
ggsave(filename=paste(outputdir.choice,"/dist_gen_mor.pdf", sep=""),plot1,  width = 8, height = 8,units="cm",dpi=600)

plot2<-ggplot(data=dist.indv,aes(x=environment,y=shape)) + 
  geom_point(alpha=0.3,shape=16,size=0.5) +  
  geom_point(data=dist.indv[pick.w,],color=ex.cols[5],shape=17,size=1) +
  geom_point(data=dist.indv[pick.r,],color=ex.cols[2],shape=17,size=1) +
  geom_point(data=dist.indv[pick.g,],color=ex.cols[4],shape=17,size=1) +
  # geom_smooth(method=lm, se=FALSE) + 
  theme_classic()
ggsave(paste(outputdir.choice,"/dist_env_mor.pdf", sep=""), plot2, width = 8, height = 8,units="cm",dpi=600)

plot3<-ggplot(data=dist.indv,aes(x=geography,y=shape)) + 
  geom_point(alpha=0.3,shape=16,size=0.5) +  
  geom_point(data=dist.indv[pick.w,],color=ex.cols[5],shape=17,size=1) +
  geom_point(data=dist.indv[pick.r,],color=ex.cols[2],shape=17,size=1) +
  geom_point(data=dist.indv[pick.g,],color=ex.cols[4],shape=17,size=1) +
  # geom_smooth(method=lm, se=FALSE) + 
  theme_classic()
ggsave(paste(outputdir.choice,"/dist_geo_mor.pdf", sep=""),plot3, width = 8, height = 8,units="cm",dpi=600)

plot4<-ggplot(data=dist.indv,aes(x=geography,y=genetic)) + 
  geom_point(alpha=0.3,shape=16,size=0.5) +  
  geom_point(data=dist.indv[pick.w,],color=ex.cols[5],shape=17,size=1) +
  geom_point(data=dist.indv[pick.r,],color=ex.cols[2],shape=17,size=1) +
  geom_point(data=dist.indv[pick.g,],color=ex.cols[4],shape=17,size=1) +
  # geom_smooth(method=lm, se=FALSE) + 
  theme_classic()
ggsave(paste(outputdir.choice,"/dist_geo_gen.pdf", sep=""), plot4, width = 8, height = 8,units="cm",dpi=600)

# # compare mean shapes of sample areas, genetic groups -----
# source(paste(scriptsdir,"/oleu_visualize_pops.R",sep=""))