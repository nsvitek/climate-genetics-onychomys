# describe mean shapes -------
#use repeatable PCs only to describe differences between groups

ip<-which(metadata$gengroup=="ip")
cp<-which(metadata$gengroup=="cp")
otero<-which(metadata$population=="cp_otero")
winkler<-which(metadata$population=="ip_winkler")
roosevelt<-which(metadata$population=="ip_roosevelt")

ip_PC<-PC2shape(PCA.phenotype,ip,c(1:11),colMeans(phenotype)) %>% 
    matrix(.,ncol=molars$m,nrow=molars$k,byrow=TRUE) 

cp_PC<-PC2shape(PCA.phenotype,cp,c(1:11),colMeans(phenotype)) %>% 
  matrix(.,ncol=molars$m,nrow=molars$k,byrow=TRUE) 

#calculate heat map to paint shape
color_ip_cp<-shpdif(ip_PC,cp_PC,ramp,alter="square",outlier=TRUE)
  
open3d()
plot3d(ip_PC,col=color_ip_cp,size=10,axes=F,xlab="",ylab="",zlab="")
writePLY(paste("ip_vs_cp.ply",sep=""),format="ascii",pointRadius=0.005)
rgl.close()

open3d()
plot3d(cp_PC,col=color_ip_cp,size=10,axes=F,xlab="",ylab="",zlab="")
writePLY(paste("cp_vs_ip.ply",sep=""),format="ascii",pointRadius=0.005)
rgl.close()

otero_PC<-PC2shape(PCA.phenotype,otero,c(1:11),colMeans(phenotype)) %>% 
  matrix(.,ncol=molars$m,nrow=molars$k,byrow=TRUE) 
roosevelt_PC<-PC2shape(PCA.phenotype,roosevelt,c(1:11),colMeans(phenotype)) %>% 
  matrix(.,ncol=molars$m,nrow=molars$k,byrow=TRUE) 
winkler_PC<-PC2shape(PCA.phenotype,winkler,c(1:11),colMeans(phenotype)) %>% 
  matrix(.,ncol=molars$m,nrow=molars$k,byrow=TRUE) 

color_ot_ro<-shpdif(otero_PC,roosevelt_PC,ramp,alter="square",outlier=TRUE)
color_ot_wi<-shpdif(otero_PC,winkler_PC,ramp,alter="square",outlier=TRUE)
color_ro_wi<-shpdif(roosevelt_PC,winkler_PC,ramp,alter="square",outlier=TRUE)

open3d()
plot3d(otero_PC,col=color_ot_ro,size=10,axes=F,xlab="",ylab="",zlab="")
writePLY(paste("otero_vs_roosevelt.ply",sep=""),format="ascii",pointRadius=0.005)
rgl.close()

open3d()
plot3d(roosevelt_PC,col=color_ot_ro,size=10,axes=F,xlab="",ylab="",zlab="")
writePLY(paste("roosevelt_vs_otero.ply",sep=""),format="ascii",pointRadius=0.005)
rgl.close()

open3d()
plot3d(otero_PC,col=color_ot_wi,size=10,axes=F,xlab="",ylab="",zlab="")
writePLY(paste("otero_vs_winkler.ply",sep=""),format="ascii",pointRadius=0.005)
rgl.close()

open3d()
plot3d(winkler_PC,col=color_ot_wi,size=10,axes=F,xlab="",ylab="",zlab="")
writePLY(paste("winkler_vs_otero.ply",sep=""),format="ascii",pointRadius=0.005)
rgl.close()

open3d()
plot3d(winkler_PC,col=color_ro_wi,size=10,axes=F,xlab="",ylab="",zlab="")
writePLY(paste("winkler_vs_roosevelt.ply",sep=""),format="ascii",pointRadius=0.005)
rgl.close()

open3d()
plot3d(roosevelt_PC,col=color_ro_wi,size=10,axes=F,xlab="",ylab="",zlab="")
writePLY(paste("roosevelt_vs_winkler.ply",sep=""),format="ascii",pointRadius=0.005)
rgl.close()