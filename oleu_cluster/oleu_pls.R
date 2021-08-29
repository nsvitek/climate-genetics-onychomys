
sink(paste(outputdir.choice,"/PLS_output_md",mdat,"_",phenotype.option,".txt",sep=""))
print(paste("r = ",mdat,sep = ""))

print(paste("PCs 1-",stop.gen," account for at least 95% of variation", sep=""))

# 2-b pls --------
gen.env.b.pls<-two.b.pls(gen.block, env.block, iter = 999, seed = NULL, print.progress = FALSE) #genotype x envrionment
gen.geo.b.pls<-two.b.pls(gen.block, geo.block, iter = 999, seed = NULL, print.progress = FALSE) #genotype x geography
env.geo.b.pls<-two.b.pls(env.block, geo.block, iter = 999, seed = NULL, print.progress = FALSE) #geography x environment
summary(gen.env.b.pls) #genetics and environment strongly spatially related
summary(gen.geo.b.pls) #population structure strongly spatially related
summary(env.geo.b.pls) #environmental structure strongly spatially related

# #example of comparison of block-comparison
com.2b<-compare.pls(gen.env.b.pls,gen.geo.b.pls,env.geo.b.pls)
summary(com.2b) #result: genxenv explains signficantly less than geoxenv explains

# men.b.pls<-two.b.pls(arrayspecs(morgenblock,p=ncol(morgenblock)/3,k=3), envmorblock, iter = 999, seed = NULL, print.progress = TRUE) #morphology x envrionment
men.b.pls<-two.b.pls(mor.block, env.block, iter = 999, seed = NULL, print.progress = FALSE) #morphology x envrionment
mge.b.pls<-two.b.pls(mor.block, gen.block, iter = 999, seed = NULL, print.progress = FALSE) #morphology x genotype
mgo.b.pls<-two.b.pls(mor.block, geo.block, iter = 999, seed = NULL, print.progress = FALSE) #morphology x geography
summary(men.b.pls)
summary(mge.b.pls)
summary(mgo.b.pls)

com.2b<-compare.pls(men.b.pls,mge.b.pls,mgo.b.pls)
com.2b<-compare.pls(men.b.pls,mge.b.pls)
summary(com.2b)

# explore results -------
# str(mge.b.pls)
default.cols<-hue_pal()(11)
# ptol_pal()(3)

pdf(file=paste(outputdir.choice,"/PLS_MorGen_md",mdat,".pdf",sep=""))
plot(mge.b.pls,col=hue_pal()(11)[factor(metadata$`Sampling Area`)],
     pch=c(16,17,15)[factor(metadata$`Genetic Group`)],cex=2)
dev.off()

pdf(file=paste(outputdir.choice,"/PLS_MorEnv_md",mdat,".pdf",sep=""))
plot(men.b.pls,col=hue_pal()(11)[factor(metadata$`Sampling Area`)],
     pch=c(16,17,15)[factor(metadata$`Genetic Group`)],cex=2)
dev.off()

# pdf(file=paste(outputdir.choice,"/legend_PLS_MorEnv_md",mdat,".pdf",sep=""))
# plot(c(1,1),col=hue_pal()(11)[factor(metadata$`Sampling Area`)],
#      pch=c(16,17,15)[factor(metadata$`Genetic Group`)])
# legend("center",levels(factor(metadata$gengroup)),pch=16,col=ptol_pal()(3),pt.cex=2)
# dev.off()

#what is driving environmental connection with morphology?
print("")
print("environmental loadings on morphology")
men.b.pls$right.pls.vectors[,1] %>% print()
print("")
print("spatial loadings on morphology")
mgo.b.pls$right.pls.vectors[,1] %>% print()

#only run if phenotype.option == "points"
#   # 2B-PLS shape prediction type 1: environment ------
#   #below only works if you use coordinates, not PCs
#   pred.limits<-men.b.pls$XScores %>% range
#   preds <- shape.predictor(arrayspecs(mor.block,p=ncol(mor.block)/3,k=3),
#                            env.block, Intercept = FALSE,
#                            method = "PLS",
#                            pred1 = pred.limits[1], pred2 = pred.limits[2]) # using PLS plot as a guide
#   colors1<-shpdif(preds$pred1,preds$pred2,ramp,alter="none")
#   
#   open3d()
#   plot3d(preds$pred2,axes=F,col=colors1,size=10,xlab="",ylab="",zlab="")
#   writePLY(paste(outputdir.choice,"/PLS_MorEnv_md",mdat,"_max.ply",sep=""),format="ascii",pointRadius=0.005)
#   rgl.close()
#   
#   open3d()
#   plot3d(preds$pred1,axes=F,col=colors1,size=10,xlab="",ylab="",zlab="")
#   writePLY(paste(outputdir.choice,"/PLS_MorEnv_md",mdat,"_min.ply",sep=""),format="ascii",pointRadius=0.005)
#   rgl.close()
#   
#   # 2B-PLS shape prediction type 1: genotype ------
#   pred.limits<-mge.b.pls$XScores %>% range
#   preds <- shape.predictor(arrayspecs(mor.block,p=ncol(mor.block)/3,k=3),
#                            gen.block, Intercept = FALSE,
#                            method = "PLS",
#                            pred1 = pred.limits[1], pred2 = pred.limits[2]) # using PLS plot as a guide
#   colors1<-shpdif(preds$pred1,preds$pred2,ramp,alter="none")
#   
#   open3d()
#   plot3d(preds$pred2,axes=F,col=colors1,size=10,xlab="",ylab="",zlab="")
#   writePLY(paste(outputdir.choice,"/PLS_MorGen_md",mdat,"_max.ply",sep=""),format="ascii",pointRadius=0.005)
#   rgl.close()
#   
#   open3d()
#   plot3d(preds$pred1,axes=F,col=colors1,size=10,xlab="",ylab="",zlab="")
#   writePLY(paste(outputdir.choice,"/PLS_MorGen_md",mdat,"_min.ply",sep=""),format="ascii",pointRadius=0.005)
#   rgl.close()
#   
#   # shape prediction type 2 -------
#   comp.cols<-ptol_pal()(2)
#   pls2.env<-pls2B(arrayspecs(mor.block,p=ncol(mor.block)/3,k=3), env.block)
#   plsEffects <- plsCoVar(pls2.env,i=1)
#   
#   open3d() #colors changed to match PLS a priori group
#   deformGrid3d(plsEffects$x[,,2],plsEffects$x[,,1],type="p",col1=comp.cols[1],col2=comp.cols[2])
#   writePLY(paste(outputdir.choice,"/PLS_MorEnv_md",mdat,"_comparison.ply",sep=""),format="ascii",pointRadius=0.03)
#   rgl.close()
#   
#   pls2.gen<-pls2B(arrayspecs(mor.block,p=ncol(mor.block)/3,k=3), gen.block)
#   deformGrid3d(plsEffects$x[,,1],plsEffects$x[,,2],type="p",col1=comp.cols[1],col2=comp.cols[2])
#   writePLY(paste(outputdir.choice,"/PLS_MorGen_md",mdat,"_comparison.ply",sep=""),format="ascii",pointRadius=0.03)
#   rgl.close()