# Script to load GIS data, especially climate rasters, into R.
# # From the PRISM 30-year Normals Dataset (1981-2010):
# # precipitation, mean temperature, maximum temperature, minimum temperature,
# # mean dewpoint temperature (measure of humidity), minimum vapor pressure deficit,
# # maximum vapor pressure deficit (both also related to humidity), elevation
# # From here: http://sedac.ciesin.columbia.edu/data/set/hanpp-net-primary-productivity/data-download#close
# # Net Primary Productivity
# # From USGS: National Land Cover and eMODIS Phenology (annual, avg'd 2001-2010) 

#create list of all objects currently there
# freeze<-ls()

# install.packages("rgdal")
setwd(paste(datadir,"/../environment",sep=""))

# #USA National Land Cover from USGS. Projection, extent data read in automatically
landcover<-raster("nlcd_2001_landcover_2011_edition_2014_10_10/nlcd_2001_landcover_2011_edition_2014_10_10.img")
# CRS.landcover<-"+proj=aea +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +6378137 +ellps=GRS80 +datum=NAD83 +no_defs"

# #Net Primary Productivity, from CIESIN. Extent, projection, etc. read in automatically.
productivity<-raster("npp-geotiff/npp_geotiff.tif") 
# CRS.productivity<-"+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

# #PRISM Data. Projection, extent, coordinate reference, etc read in automatically. 
# prism_files<-list.files(paste(getwd(),"/prism/",sep=""),full.names=TRUE)
# for(i in 1:length(prism_files)){
#   prism_files[i]<-list.files(prism_files[i],full.names=TRUE,pattern = "\\.bil$")
# }
# prism<-stack(prism_files)
# names(prism)<-c("precipitation_mean","dew_point_temperature_mean",
#                 "temperature_max","temperature_mean",
#                 "temperature_min","elevation",
#                 "vapor_pressure_max","vapor_pressure_min")
# # CRS.prism<-"+proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs"
#read in projected files
prism<-list.files(getwd(),pattern = "^prism-proj*.*\\.tif$") %>% stack

# #Phenology from eMODIS. Projection, etc. info included on load.
# phenology<- list.files(getwd(),pattern = "\\.grd$") %>% stack
# phen_names<-c("NDVI_amplitude","growing_season_duration","season_end_NDVI","season_end_time","NDVI_max",
#               "NDVI_max_time","season_start_NDVI","season_start_time","NDVI_cumulative")
# names(phenology)<-phen_names
#read in projected files 
phenology<-list.files(getwd(),pattern = "^phenology_proj*.*\\.tif$") %>% stack

productivity<-projectRaster(from=productivity,crs=proj4string(landcover))
#reprojecting phenology and prism rasters is slow, don't wish to re-do, work from saved versions above
# prism<-projectRaster(from=prism,crs=proj4string(landcover),useNA=FALSE) 
# phenology.proj<-projectRaster(from=phenology,crs=proj4string(landcover))
# writeRaster(phenology.proj,filename="phenology_proj.tif",progress='text',
#             bylayer=TRUE,suffix=names(phenology.proj),format="GTiff")
# writeRaster(prism,filename="prism-proj.tif",progress='text',
#             bylayer=TRUE,suffix=names(prism),format="GTiff")

#set extent on all?
# make a matrix out of extents, each column represents a raster, rows the values
# productivity.c<-crop(productivity, extent(prism))
# extent_objects<-c(prism,productivity.c,landcover,phenology)
# extent_list<-lapply(extent_objects,extent) %>% 
#   lapply(.,as.matrix) %>% unlist %>% matrix(.,ncol=length(extent_objects))
# 
# # create an extent with the extrem values of your extent
# best_extent<-extent(min(extent_list[1,]), max(extent_list[3,]),
#                     min(extent_list[2,]), max(extent_list[4,]))
# 
# # the range of your extent in degrees
# ranges<-apply(as.matrix(best_extent), 1, diff)
# # the resolution of your raster (pick one) or add a desired resolution
# dummy<-list()
# for (i in 1:length(extent_objects)){
#   nrow_ncol<-res(extent_objects[[i]][[1]]) #from here on needs to be done as a loop
#   # deviding the range by your desired resolution gives you the number of rows and columns
#   nrow_ncol<-ranges/reso
#   dummy[[i]]<-raster(best_extent, nrows=nrow_ncol[2], ncols=nrow_ncol[1], crs=extent_objects[[i]][[1]]@crs)
# }
# # nrow_ncol<-res(extent_objects[[4]][[1]]) #from here on needs to be done as a loop
# # # deviding the range by your desired resolution gives you the number of rows and columns
# # nrow_ncol<-ranges/reso
# # 
# # # create your raster with the following
# # dummy<-raster(best_extent, nrows=nrow_ncol[2], ncols=nrow_ncol[1], crs=extent_objects[1][[1]]@crs)
# productivity_r<-resample(productivity.c, dummy[[2]], method="ngb")
# phenology_r<-resample(phenology, dummy[[4]], method="ngb")
# prism_r<-resample(prism, dummy[[1]], method="ngb")
# landcover_r<-resample(landcover, dummy[[3]], method="bilinear")
# 
# phenology_m<-mask(phenology_r,prism_r[[1]])
# productivity_m<-mask(productivity_r,prism_r[[1]])
# landcover_m<-mask(landcover,prism_r[[1]])

# plot(productivity_m, add=TRUE,legend=FALSE)  
# plot(phenology_m[[1]])
# plot(phenology_r[[1]])
# plot(prism_r[[1]], add=TRUE) 
# plot(landcover_m,add=TRUE,legend=F)
# plot(landcover)
# plot(landcover_r)

# environment<-stack(prism_r,productivity_m,phenology_m,landcover_m)
#for visualization only. seems to have lost a lot of information.
# writeRaster(environment,filename="environment.tif",progress='text',
#   bylayer=TRUE,suffix=names(environment),format="GTiff")
environment<-list.files(getwd(),pattern = "environment*.*\\.tif$") %>% stack

# freeze<-c(freeze,"environment")
# rm(list = setdiff(ls(),freeze))