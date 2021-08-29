# phenotype PCA ------------------
PCA.phenotype<-prcomp(phenotype,scale.=FALSE)

#for reference:
plot(PCA.phenotype$x[,c(1:2)],xlim=c(-0.20,0.10),ylim=c(-0.1,0.1))

# choose which points on PC1-2 to sample -----
coords<-matrix(c(-0.2,0.1,
                 -0.1,0.1,
                 0,0.1,
                 0.1,0.1,
                 -0.2,0,
                 -0.1,0,
                 0,0,
                 0.1,0,
                 -0.2,-0.1,
                 -0.1,-0.1,
                 0,-0.1,
                 0.1,-0.1),byrow=TRUE, ncol=2) %>% as.data.frame

# example for calculating single shape -----

# ### calculate the shapes at the maxima and minima of principal component axes
# #input: object output from prcomp() function, mshp output, desired PCs (ex: 1, 2:4)
# #output: k x m shapes for the maximum and minimum of each desired PC
# pcs<-c(1,2)
# coords<-c(-0.1,-0.06) #example in empty space
# #to transform to PC1 space, take the location along pc1 axis * the $rotation values from prcomp object
# #then make it as a matrix with 3 rows and #LM of columns (which will then be transformed)
# pc.trans.1<-t(matrix(coords[1]*PCA.phenotype$rotation[,pcs[1]],3,dim(phenotype)[2]/3))
# pc.trans.2<-t(matrix(coords[2]*PCA.phenotype$rotation[,pcs[2]],3,dim(phenotype)[2]/3))
# #calculate mean shape
# meanshape<-t(matrix(colMeans(phenotype),3,dim(phenotype)[2]/3))
# hypothetical.shape<-meanshape+pc.trans.1+pc.trans.2
# 
# library(plot3D)
# plot3d(x=hypothetical.shape,size=10)
# writePLY(paste("hypothetical_shape_01.ply",sep=""),format="ascii",pointRadius=0.005)
# rgl.close()

# loop through all hypotethetical shapes ----
### calculate the shapes at the maxima and minima of principal component axes
#input: object output from prcomp() function, mshp output, desired PCs (ex: 1, 2:4)
#output: k x m shapes for the maximum and minimum of each desired PC
pcs<-c(1,2)

for (i in c(1:nrow(coords))){
  #to transform to PC1 space, take the location along pc1 axis * the $rotation values from prcomp object
  #then make it as a matrix with 3 rows and #LM of columns (which will then be transformed)
  pc.trans.1<-t(matrix(coords[i,1]*PCA.phenotype$rotation[,pcs[1]],3,dim(phenotype)[2]/3))
  pc.trans.2<-t(matrix(coords[i,2]*PCA.phenotype$rotation[,pcs[2]],3,dim(phenotype)[2]/3))
  #calculate mean shape
  meanshape<-t(matrix(colMeans(phenotype),3,dim(phenotype)[2]/3))
  hypothetical.shape<-meanshape+pc.trans.1+pc.trans.2
  
  library(plot3D)
  plot3d(x=hypothetical.shape,size=10)
  writePLY(paste("hypothetical_shape_",i,".ply",sep=""),format="ascii",pointRadius=0.005)
  rgl.close()
}


# plot PCA -----
#for theoretical shapes?
PCA.p.perc<-round(summary(PCA.phenotype)$importance[2,]*100,1)
PCA.p.2plot<-cbind(metadata,PCA.phenotype$x)


ggplot(data=PCA.p.2plot,aes(x=PC1,y=PC2, color = `Sampling Area`,shape=`Genetic Group`))+
  geom_point(size=2) +
  theme_classic() + theme(legend.position="none") +
  xlab(paste("PC 1 (",PCA.p.perc[1],"%)",sep="")) +
  ylab(paste("PC 2 (",PCA.p.perc[2],"%)",sep="")) +
  geom_point(data=coords, aes(x=V1, y=V2), size=3,color = "black", shape = 16)
ggsave("PCA_phen_theory.pdf", device = cairo_pdf, width = 8, height = 8,units="cm",dpi=600)