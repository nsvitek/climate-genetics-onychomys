# # Base script of family to process environmental data from a series of localities.

# # Dependencies ------
# library(dismo)
library(dplyr)
# library(geosphere)
# library(gtools)
# library(maptools)
# library(raster)
# library(rasterVis)
# library(rgdal)
# library(rgeos)
library(sp)
# library(vegan)
# library(ggplot2)
# library(reshape2)
library(readxl) #used in sourced scripts: GIS_specimens_[critter name]

scriptsdir<-"C://scripts/scripts/GIS-sample-PRISM"
datadir<-"D:/Dropbox/Documents/Dissertation/modern/onychomys_leucogaster"

# # Load locality coordinates -------
setwd(scriptsdir)
source("GIS_specimens_oleu.R") #makes specimen_coords, or cleaned coordinates for specimens

# # no longer need environmental data, can skip all below ------
# # also need to get landcover object -----
# source("GIS_map-load.R") #load maps of environmental data
# 
# # # deal with uncertainty (coordinateuncertaintyinmeters) -----
# # # Put everything in one datum & projection: NAD83, Albers Equal Area
# specimen_circles<-list()
# for (i in 1:nrow(specimen_coords)){
#   specimen_circles[[i]]<-c(specimen_coords$decimallongitude[i],specimen_coords$decimallatitude[i]) %>%
#     matrix(.,ncol=2) %>% 
#     circles(.,d=specimen_coords$coordinateuncertaintyinmeters[i],lonlat=T) %>%
#     .@polygons %>% spTransform(.,crs(landcover))
# }
# specimen_circles<-do.call(bind,specimen_circles)
# 
# # # Extract raster data for coordinates: pull, average, and turn into a table
# specimen_prism<-extract(prism, specimen_circles) %>% 
#   lapply(.,function(x) apply(x,2,mean)) %>%
#   unlist %>% matrix(.,byrow=TRUE,ncol=dim(prism)[3])
# colnames(specimen_prism)<-c("dew_point_temperature_mean","elevation",
#                  "precipitation_mean","temperature_max",
#                  "temperature_mean","temperature_min",
#                  "vapor_pressure_max","vapor_pressure_min")
# specimen_phenology<-extract(phenology, specimen_circles) %>% 
#   lapply(.,function(x) apply(x,2,mean)) %>%
#   unlist %>% matrix(.,byrow=TRUE,ncol=dim(phenology)[3])
# colnames(specimen_phenology)<-c("growing_season_duration","NDVI_amplitude",
#                                 "NDVI_cumulative","NDVI_max","NDVI_max_time",
#                                 "season_end_NDVI","season_end_time",
#                                 "season_start_NDVI","season_start_time")
# specimen_productivity<-extract(productivity, specimen_circles) %>% 
#   lapply(.,mean) %>%
#   unlist %>% matrix(.,byrow=TRUE,ncol=dim(productivity)[3])
# colnames(specimen_productivity)<-"net_primary_productivity"
# specimen_environment<-cbind(specimen_prism,specimen_phenology,specimen_productivity)
# head(specimen_environment)
# 
# #put specimen names as row names
# row.names(specimen_environment)<-specimen_coords$full_specimenb
# specimen_environment<-cbind(specimen_coords[,2:ncol(specimen_coords)],specimen_environment)
# # # Save environmental data as a table. 
# setwd(datadir)
# write.csv(specimen_environment,"oleu_environment.csv",quote=FALSE,row.names=TRUE)
