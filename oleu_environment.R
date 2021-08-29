# standardize, summarize, explore the environmental data for sample of O. leucogaster
# geographic data -------
#create scaled block of geographic data
geo.block<-enviro.all %>% 
  dplyr::select(c(decimallatitude,decimallongitude))

# environmental data -------
env.block<-enviro.all %>% 
  dplyr::select(-specimen, -decimallatitude, -decimallongitude,
                -geodeticdatum, -coordinateuncertaintyinmeters)  %>%
  scale()

#run correlation PCA w/ pre-scaled data 
PCA.environment<- prcomp(env.block,center=FALSE,scale=FALSE)

# phenotype PCA ------------------
PCA.phenotype<-prcomp(phenotype,scale.=FALSE)


# # Run once, then comment out: make PC plots ------
# master.variables<-cbind(metadata,enviro.all,PCA.phenotype$x[,1:10])
# source(paste(scriptsdir,"/oleu_PCA_non_gen.R",sep=""))

# clean --------
rm("enviro.filtered","enviro.raw", "PCA.e.rot","PCA.e.perc","PCA.e.2plot","PCA.p.2plot","PCA.p.perc","hull.wb.1","hull.not.1",
   "hull.wb.2","hull.not.2","hull.plot1","hull.plot2")