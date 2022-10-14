# Onychomys leucogaster-specific script for reading, cleaning, and processing locality information from VertNet.

setwd(datadir)
#create list of all objects currently there
freeze<-ls()

# create output from vertnet only for specimens scanned or sequenced -----
sampled_metadata<-read_excel("onychomys_ct.xlsx") %>% filter(vial!="WORN FLAT"&vial!="TOO WORN"&vial!="NOT LOANED")
raw_metadata<-read.csv("planning_notes/Onychomys_potential.csv",header=T)
raw_metadata$filename<-paste(raw_metadata$institutioncode,raw_metadata$catalognumber,sep="-")
raw_specimens<-raw_metadata[raw_metadata$filename %in% sampled_metadata$filename,]
#enter datum for specimens georeferenced using Google Earth, GEOlocate, or Wikipedia
#Wikipedia Idaho specimens appear to be in WGS84, based on following the Wikipedia article to GeoHack
raw_specimens$geodeticdatum[which(raw_specimens$georeferencesources=="Google Earth"|raw_specimens$georeferencesources=="Wikipedia")]<-"WGS84"
#According to HerpNet (http://www.herpnet.org/Gazetteer/GeorefResources.htm), GEOlocate is in WGS84
raw_specimens$geodeticdatum[which(raw_specimens$georeferenceprotocol=="GEOLocate")]<-"WGS84"
#Remaining specimens missing datum: likely either NAD83 or WGS84, which are within 1 m of each other
#Solution: write datum as WGS1984
raw_specimens$geodeticdatum[which(raw_specimens$geodeticdatum==""|raw_specimens$geodeticdatum=="not recorded"|raw_specimens$geodeticdatum=="unknown")]<-"WGS84"


# Explore data, see what needs doing: -------
easy_specimens<-raw_specimens %>%
  filter(.,geodeticdatum!=""&geodeticdatum!="not recorded"&geodeticdatum!="unknown") %>%
  dplyr::select(.,filename,decimallatitude,decimallongitude,geodeticdatum,coordinateuncertaintyinmeters)

# #get the specimens that have no datum :( But no specimens have no datum!
# problem_specimens<-raw_specimens %>%
#   filter(.,geodeticdatum==""|geodeticdatum=="not recorded"|geodeticdatum=="unknown") %>%
#   select(.,c(2,74:75,77:78,89:107))
# # check to make sure all specimens in one object or the other
# nrow(problem_specimens)+nrow(easy_specimens)==nrow(raw_specimens)

# Fill in coordinate uncertainty -----
easy_specimens$coordinateuncertaintyinmeters[is.na(easy_specimens$coordinateuncertaintyinmeters)]<-0

# Standardize datums ------
easy_specimens<-droplevels(easy_specimens)
levels(easy_specimens$geodeticdatum)
easy_specimens$geodeticdatum[easy_specimens$geodeticdatum=="World Geodetic System 1984"]<-"WGS84"
#add small uncertainty for NAD83 to WGS84 conversion
easy_specimens$coordinateuncertaintyinmeters[easy_specimens$geodeticdatum=="North American Datum 1983"]<-1
easy_specimens$geodeticdatum[easy_specimens$geodeticdatum=="North American Datum 1983"]<-"WGS84"
#convert coordinates from NAD27 to WGS84
convert<-easy_specimens %>% filter(.,geodeticdatum=="North American Datum 1927") %>% 
  dplyr::select(.,decimallongitude,decimallatitude) %>% data.frame
coordinates(convert)<-c("decimallongitude","decimallatitude")
proj4string(convert)<-CRS("+init=epsg:4267") # NAD27
CRS.WGS84<-CRS("+init=epsg:4326") #WGS 1984
converted<-spTransform(convert,CRS.WGS84) %>% as.data.frame(.@coords)
easy_specimens$decimallatitude[easy_specimens$geodeticdatum=="North American Datum 1927"]<-converted$decimallatitude
easy_specimens$decimallongitude[easy_specimens$geodeticdatum=="North American Datum 1927"]<-converted$decimallongitude
easy_specimens$geodeticdatum[easy_specimens$geodeticdatum=="North American Datum 1927"]<-"WGS84"

#Give tiny extent so that circles can exist -----
easy_specimens$coordinateuncertaintyinmeters[which(easy_specimens$coordinateuncertaintyinmeters==0.0)]<-1

#final output ------
specimen_coords<-easy_specimens

#remove all temporary objects, i.e., everything except the final produce and objects that were here before. 
freeze<-c(freeze,"specimen_coords")
rm(list = setdiff(ls(),freeze))
