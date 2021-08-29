freeze<-ls() 

# read metadata/linear measurements ----
#takes raw_metadata from oleu_base. If not running those scrips, uncomment line below
# source(paste(scriptsdir,"/oleu_linear_error.R",sep=""))
# linear_metadata<-metadata %>% filter(is.na(m1_length_1)==FALSE) %>% 
#   mutate(m1.length=(m1_length_1+m1_length_2+m1_length_3)/3,
#                  m1.width=(m1_width_1 + m1_width_2 + m1_width_3)/3) 
# linear_metadata<-mutate(linear_metadata,m1.area=m1.length*m1.width,m1.proportion=m1.width/m1.length) %>% 
#   select(.,filename,stateprovince,county,population,m1.length,m1.width,m1.area,m1.proportion,cropped)
#  
# # #exploration figures for linear data
# n<-nrow(raw_metadata) #number of specimens
# ggplot(data=linear_metadata,aes(x=m1.length,y=m1.width, color = population, fill=gengroup))+
#   geom_point(alpha=1,size=2)
# 
# ggplot(data=linear_metadata,aes(x=population,y=m1.area, fill=gengroup))+
#   geom_boxplot()
# 
# ggplot(data=linear_metadata,aes(x=population,y=m1.proportion, fill=gengroup))+
#   geom_violin()

# #proportions look slightly different between Colorado Plateao and Interior Plains. Really?
# aov(data=linear_metadata[which(linear_metadata$gengroup!="wb"),], m1.proportion~gengroup) %>% summary

# read shape coordinate data -----
shapedir<-"morphology/oleu_crop4x_191231_1024" #no cs

#read in pseudolandmark data
molars<-read.morphologika(paste(datadir,shapedir,
                                "morphologika_unscaled_high.txt",sep="/")) %>%
  preprocess(.)

# match metadata and coordinate data ---------
label_error<-dimnames(molars$scaled)[[3]] %>% strsplit(.,split="[_-]")
#create variables based on file names
specimen<-crop<-NULL
for (i in 1:length(label_error)){
  specimen[i]<-paste(label_error[[i]][1],label_error[[i]][2],sep="-")
  crop[i]<-label_error[[i]][6]
}

shapes<-molars$m2d[which(crop==1),]
rownames(shapes)<-specimen[which(crop==1)]

#clean space
freeze<-c(freeze, "shapes")
rm(list = setdiff(ls(),freeze))
# # error -----
# setwd(paste(datadir,"/output/results_error",sep="")) #move to folder for keeping error-related results
# freeze<-ls() #take a snapshop of objects in environment
# source(paste(scriptsdir,"/oleu_error.R",sep=""))
# rm(list = setdiff(ls(),freeze)) #clean up environment, also removes freeze
# setwd("../")