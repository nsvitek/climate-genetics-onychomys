# Script to average 10 years of eMODIS phenology data into single raster. Hopefully, only done once.

# #Phenology from eMODIS. Projection, etc. info included on load.
#make mean phenology map for each layer, then save.
phen_abbr<-c("AMP","DUR","EOSN","EOST","MAXN","MAXT","SOSN",
             "SOST","TIN")
phen_names<-c("NDVI_amplitude","growing_season_duration","season_end_NDVI","season_end_time","NDVI_max",
              "NDVI_max_time","season_start_NDVI","season_start_time","NDVI_cumulative")
# 1. Start of Season Time (SOST): starting time of the onset of the growing season (in day of the year). 
# 2. Start of Season NDVI (SOSN): NDVI value at the starting time of the onset of the growing season (unitless- based on NDVI units). 
# 3. End of Season Time (EOST): ending time of the growing season (in day of the year). 
# 4. End of Season NDVI (EOSN): NDVI value at the ending time of the growing season (unitless-based on NDVI units).
# 5. Maximum Time (MAXT): the day of the year when the NDVI reaches its maximum during the growing season (in day of the year). 
# 6. Maximum NDVI (MAXN): the highest (or peak) value in NDVI observed in a growing season (unitless-based on NDVI units). 
# 7. Duration (DUR): the length of the growing season-the time between the start of season and end of season (in number of days). 
# 8. Amplitude (AMP): the difference between the Maximum NDVI and NDVI at the day of start of season (unitless-based on NDVI units). 
# 9. Time Integrated NDVI (TIN): the cumulative value of NDVI from the start to the end of the growing season (unitless-based on accumulated NDVI units).


phenology_files<-files2<-list.files(paste(getwd(),"/phenology/",sep=""),full.names=TRUE)
for (var in 1:length(phen_abbr)){
  for(i in 1:length(phenology_files)){
    files2[i]<-list.files(phenology_files[i],full.names=TRUE,pattern = paste("^",phen_abbr[var],".*.tif$",sep=""))
  }
  stackE<-stack(files2[1:10])
  stackW<-stack(files2[11:20])
  print(paste("calculating mean for East stack", var))
  stackE<-mean(stackE)
  print(paste("calculating mean for West stack",var))
  stackW<-mean(stackW)
  print(paste("merging stacks",var))
  stackEW<-merge(stackE,stackW)
  writeRaster(stackEW, filename = paste("phen_2001-2010_mean_",phen_abbr[var],".grd",sep=""))
}



