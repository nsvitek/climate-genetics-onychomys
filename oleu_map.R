# Script to report specimen locality sampling
# Sourced early in oleu_base.R, uses objectes enviro.all, metadata


# extra mapping dependencies -----

library(gtools)
# library(maptools) #map2SpatialPolygons()
# library(raster) #stack function?
# library(rasterVis)
# library(rgdal)
# library(rgeos)
# library(sf) #all the spatial transforms, crs
# library(maps) #map() function

# setwd(datadir)
# setwd("../environment")
# environment<-list.files(getwd(),pattern = "environment*.*\\.tif$")[-1] %>% stack
# Combine in state and county data with environmental data. -------
specimen_environment<-cbind(enviro.all,metadata)

# Report sample size by sampling area ------
group_by(specimen_environment,population) %>% 
  summarise(.,n(),MAT=mean(temperature_mean),MAP=mean(precipitation_mean),MNDVI=mean(NDVI_cumulative)) %>% 
  write.csv(.,"sampling_population.csv")


# Deal with distances. -----------------
#slim down to only unique points
# modelmice<-distinct(specimen_environment,decimallatitude,decimallongitude,.keep_all=TRUE) %>% droplevels
modelmice<-distinct(specimen_environment,county,.keep_all=TRUE) %>% droplevels
locality_coords<-modelmice
coordinates(locality_coords)<-c("decimallongitude", "decimallatitude")
proj4string(locality_coords)<-CRS("+init=epsg:4326") #WGS84

# Note from GIS processing: climate layers were in NAD83/Albers North American projectsion:
raster.CRS<-"+proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs "
locality_coords<-spTransform(locality_coords,CRS(raster.CRS))
locality_coords$population
#envrionment

# # make regular map of US states & localities ------------
US <- map("state",fill=TRUE, plot=FALSE)
US.names <- US$names
US.IDs <- sapply(strsplit(US.names,":"),function(x) x[1])
US_poly_sp <- map2SpatialPolygons(US,IDs=US.IDs,proj4string=CRS("+proj=longlat + datum=WGS84"))
US_poly_sp.albers <- spTransform(US_poly_sp, CRS=CRS(raster.CRS))

# read in physiogeographic provinces -----
#old, cut this block once working
phys.prov.raw<-readOGR(dsn = "physio_shp")

phys.prov.raw$

# set_RGEOS_CheckValidity(phys.prov.raw)
# #now only want province level
# phys.prov.raw$PROVCODE
# phys.prov.raw@data

groups = aggregate(phys.prov.raw, by = "PROVINCE")
plot(phys.prov.raw)
plot(groups)

# see localities colored -----
ylims2<-range(specimen_environment$decimallatitude)+c(-0.8,+0.8)
xlims2<-range(specimen_environment$decimallongitude)+c(-0.8,+0.8)

#create a map
ggplot(specimen_environment,aes(x=decimallongitude,y=decimallatitude))+
  geom_polygon(data=US,aes(x=long,y=lat,group=group),fill=NA,color="gray",size=1.0) +
  geom_polygon(data=groups, aes(x=long,y=lat,group=group),fill=NA,color="black",size=0.5) +
  coord_cartesian(xlim=xlims2,ylim=ylims2)+
  theme_classic() + xlab("Longitude") + ylab("Latitude") +
  geom_point(aes(color=`Sampling Area`),size=2)
ggsave("specimen_map_v2.pdf", device = cairo_pdf, width = 4, height = 4,units="in",dpi=600)

# quick compare environments -----
# #range of temperature in sample?
# specimen_environment$temperature_mean %>% range %>% round(1)
# 18.5-4.4
# 
# #range of precipitation?
# specimen_environment$precipitation_mean %>% range %>% round(1)
# (720.7-223.2)/10/2.2

specimen_environment %>% 
  group_by(.,population) %>% write.csv(.,"average_population_environment.csv")
